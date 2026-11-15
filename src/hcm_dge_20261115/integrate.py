"""Join post-lock evidence without modifying or re-ranking the discovery list."""
import argparse
import hashlib
import json
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import pandas as pd


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def tab(path):
    return pd.read_csv(path, sep='\t', keep_default_na=False, na_values=['NA'])


def write(df, path):
    df.to_csv(path, sep='\t', index=False, na_rep='NA', float_format='%.15g')


def relation(effect, reference):
    values = np.full(len(effect), 'not_evaluable', dtype=object)
    ok = np.isfinite(effect)
    values[ok & effect.eq(0)] = 'zero'
    values[ok & (effect * reference > 0)] = 'same'
    values[ok & (effect * reference < 0)] = 'opposite'
    return values


def main():
    p = argparse.ArgumentParser()
    p.add_argument('work', type=Path)
    w = p.parse_args().work
    for tag in ['fit_d', 'lock', 'pb_s', 'coverage_s', 'fit_s']:
        assert (w/f'logs/{tag}.exit').read_text().strip() == '0', tag
    lock_path = w/'lock/CANDIDATES.tsv'
    digest = sha(lock_path)
    assert (w/'lock/CANDIDATES.tsv.sha256').read_text().split()[0] == digest
    assert (w/'coverage/support.lock.sha256').read_text().split()[0] == digest
    assert tab(w/'stats/validation_design.tsv').lock_sha256.iloc[0] == digest
    locked = tab(lock_path)
    ids = locked.gene_id.tolist()
    assert locked.gene_id.is_unique and len(ids) > 0
    common = ['rank', 'gene_id', 'gene_symbol', 'association_direction', 'token_id']
    d = locked[common].copy()
    primary = [
        'logFC', 'FDR', 'PValue', 'cb_logFC', 'cb_FDR',
        'loo_min_abs_logFC', 'loo_same_direction_fraction', 'loo_min_logFC',
        'loo_max_logFC', 'loo_max_abs_delta', 'loo_max_influence_omitted_donor',
        'hhatl_logFC', 'hhatl_FDR', 'Control_donors_CPM_ge1', 'HCM_donors_CPM_ge1',
        'Control_median_CPM', 'HCM_median_CPM', 'hcm_donors_token_ge10pct',
        'hcm_token_fraction', 'hcm_expr_fraction',
    ]
    for c in primary:
        d['scp_' + c] = locked[c]
    d['scp_hhatl_direction'] = relation(d.scp_hhatl_logFC, d.scp_logFC)
    for study in ['liu', 'chin', 'whole']:
        s = tab(w/f'stats/{study}.tsv')
        assert s.gene_id.is_unique and s.gene_id.tolist() == ids
        d = d.merge(s.rename(columns={c: study+'_'+c for c in s if c != 'gene_id'}),
                    on='gene_id', how='left', validate='1:1', sort=False)
    d['liu_direction'] = relation(d.liu_logFC, d.scp_logFC)
    omit = [c for c in d if c.startswith('liu_logFC_omit_')]
    assert sorted(omit) == ['liu_logFC_omit_Liu2023_CL0412', 'liu_logFC_omit_Liu2023_CL2']
    for c in omit:
        d[c.replace('logFC_', 'direction_')] = relation(d[c], d.scp_logFC)
    main_same = d.liu_direction.eq('same')
    both_same = ((d[omit].to_numpy() * d.scp_logFC.to_numpy()[:, None]) > 0).all(axis=1)
    sensitive_estimable = np.isfinite(d[omit]).all(axis=1)
    category = d.liu_direction.map({
        'not_evaluable': 'not_evaluable', 'zero': 'zero', 'opposite': 'opposite',
        'same': 'same_control_sensitive',
    })
    category.loc[main_same & ~sensitive_estimable] = 'same_sensitivity_not_evaluable'
    category.loc[main_same & both_same] = 'same_both_single_controls'
    d['liu_evidence_category'] = category
    d['liu_not_evaluable_reason'] = 'evaluated'
    missing = ~d.liu_source_feature_present
    filtered = d.liu_source_feature_present & ~d.liu_passes_expression_filter
    d.loc[missing, 'liu_not_evaluable_reason'] = 'feature_absent_or_unmapped'
    d.loc[filtered, 'liu_not_evaluable_reason'] = 'below_expression_filter'
    d['liu_FDR_lt05_same_direction'] = main_same & d.liu_FDR.lt(.05)
    for study in ['chin', 'whole']:
        d[study+'_descriptive_direction'] = relation(
            d[study+'_descriptive_log2_cpm_difference'], d.scp_logFC)
    cov = tab(w/'coverage/support.tsv')
    assert not cov.duplicated(['dataset', 'disease', 'gene_id']).any()
    coverage_cols = ['donors', 'cells', 'native_cells', 'expr_fraction',
                     'token_fraction', 'donors_native', 'donors_token_ge10pct',
                     'median_donor_token_fraction']
    coverage_summary = {}
    for ds, prefix in [('liu2023', 'liu'), ('GSE255296', 'gse255296'), ('GSE292754', 'gse292754')]:
        s = cov[(cov.dataset == ds) & (cov.disease == 'HCM')].copy()
        assert len(s) == len(ids) and set(s.gene_id) == set(ids)
        s = s.set_index('gene_id').loc[ids].reset_index()
        assert (s.token_id.to_numpy() == d.token_id.to_numpy()).all()
        assert s.native_cells.gt(0).equals(s.expr_fraction.notna())
        d = d.merge(s[['gene_id'] + coverage_cols].rename(
            columns={c: prefix+'_HCM_'+c for c in coverage_cols}), on='gene_id', validate='1:1', sort=False)
        # Aggregate evidence uses donor prevalence, not a new inclusion threshold.
        coverage_summary[ds] = {
            'HCM_donors': int(s.donors.iloc[0]), 'HCM_cells': int(s.cells.iloc[0]),
            'genes_measured_any_cells': int(s.native_cells.gt(0).sum()),
            'genes_token_in_any_cells': int(s.token_fraction.gt(0).sum()),
            'genes_token_ge10pct_in_at_least_half_donors': int(
                s.donors_token_ge10pct.ge(np.ceil(s.donors/2)).sum()),
            'median_gene_token_fraction': float(s.token_fraction.median()),
        }
    d['interpretation'] = 'HCM_associated_candidate;beneficial_perturbation_direction_undetermined'
    assert d.gene_id.tolist() == ids
    assert d['rank'].tolist() == list(range(1, len(d)+1))
    assert sha(lock_path) == digest
    summary = {
        'generated_utc': datetime.now(timezone.utc).isoformat(),
        'candidate_sha256': digest,
        'locked_candidates': len(d),
        'higher_in_HCM': int(d.scp_logFC.gt(0).sum()),
        'lower_in_HCM': int(d.scp_logFC.lt(0).sum()),
        'liu_evaluable': int(d.liu_logFC.notna().sum()),
        'liu_direction_counts': d.liu_direction.value_counts().to_dict(),
        'liu_evidence_category_counts': d.liu_evidence_category.value_counts().to_dict(),
        'liu_FDR_lt05_same_direction': int(d.liu_FDR_lt05_same_direction.sum()),
        'liu_FDR_lt05_opposite': int((d.liu_FDR.lt(.05) & d.liu_direction.eq('opposite')).sum()),
        'liu_omit_direction_counts': {c: d[c].value_counts().to_dict() for c in d if c.startswith('liu_direction_omit_')},
        'scp_hhatl_direction_counts': d.scp_hhatl_direction.value_counts().to_dict(),
        'descriptive_direction_counts': {study: d[study+'_descriptive_direction'].value_counts().to_dict() for study in ['chin', 'whole']},
        'coverage': coverage_summary,
        'interpretation': 'Association candidates, not proven therapeutic knockout targets. Validation never changes rank.',
    }
    out = w/'final'
    out.mkdir(exist_ok=True)
    write(d, out/'CANDIDATE_EVIDENCE.tsv')
    (out/'SUMMARY.json').write_text(json.dumps(summary, indent=2, allow_nan=False)+'\n', encoding='utf8')
    print(json.dumps(summary, ensure_ascii=False, allow_nan=False))


if __name__ == '__main__':
    main()
