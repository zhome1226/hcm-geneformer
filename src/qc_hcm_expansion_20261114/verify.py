"""Independently reopen delivered objects and compare every count with original GEO matrices."""
from datetime import datetime, timezone
from pathlib import Path
import gzip
import hashlib
import json

import anndata as ad
import numpy as np
import pandas as pd
from scipy.io import mmread
from threadpoolctl import threadpool_limits

ROOT = Path('/home/zhome/hcm-geneformer-data')
WORK = Path('/home/zhome/hcm-geneformer-analysis/hcm_expansion_20261114/nuclear_qc')
samples = pd.read_csv(WORK / 'samples.tsv', sep='\t')
official=json.loads((Path(__file__).parent/"SAMPLES.json").read_text())
geo={}
for source in Path(__file__).parent.glob('*_selected_samples.soft'):
    for line in source.read_text().splitlines():
        if line.startswith('^SAMPLE = '):
            gsm=line.split(' = ',1)[1];geo[gsm]={}
        elif line.startswith('!Sample_') and ' = ' in line:
            key,value=line.removeprefix('!Sample_').split(' = ',1)
            geo[gsm].setdefault(key,[]).append(value)
for row in official:
    source=geo[row['sample_id']]
    ch=dict(value.split(': ',1) for value in source['characteristics_ch1'])
    assert row['donor_original']==(ch['patient id'] if row['dataset']=='GSE292754' else source['title'][0].split(',')[0])
    if row['dataset']=='GSE292754':
        assert row['age_years']==int(ch['age']) and row['sex']==ch['Sex']


def check(checks, name, result):
    checks[name] = bool(result)


def exact(left, right):
    return left.shape == right.shape and (left != right).nnz == 0


summary = []
for dataset, group in samples.groupby('dataset', sort=True):
    output = ROOT / f'processed/{dataset}_qc_20261114'
    checks = {}
    exports = json.loads((output / 'EXPORT_SUCCESS.json').read_text())
    policy = json.loads((output / 'qc_policy.json').read_text())
    names = ['all_cells_qc.h5ad', 'cm_qc_full.h5ad', 'cm_qc_light.h5ad',
             'cm_review_pending.h5ad', 'suspected_doublets.h5ad']
    objects = {name: ad.read_h5ad(output / name) for name in names}
    allcell = objects['all_cells_qc.h5ad']
    dbl = objects['suspected_doublets.h5ad']
    core = objects['cm_qc_full.h5ad']
    pending = objects['cm_review_pending.h5ad']
    light = objects['cm_qc_light.h5ad']
    identity = pd.read_csv(output / 'cell_identity_audit.tsv.gz', sep='\t', index_col='barcode', low_memory=False)
    check(checks, 'policy_strict_limits', policy['min_genes'] == 500 and policy['max_genes'] == 5000 and policy['mt_percent_exclusive_max'] == 5)
    check(checks, 'core_all_identity_retained', core.obs.review_action.eq('retain_CM').all())
    check(checks, 'core_best_and_pruned_reference_CM', core.obs.reference_label.eq('Cardiomyocyte').all() and core.obs.reference_pruned.eq('Cardiomyocyte').all())
    check(checks, 'core_multigene_and_program_support', core.obs.Cardiomyocyte_markers.ge(2).all() and core.obs.cm_program_dominates.all())
    check(checks, 'core_no_strong_mixed_program', not core.obs.coherent_mixed_panel_flag.any())
    check(checks, 'core_reference_tail_or_context_support', (core.obs.reference_tail_pass | (core.obs.Cardiomyocyte_markers.ge(3) & core.obs.cm_anchor_neighbor_fraction.ge(.5))).all())
    check(checks, 'pending_explicit_and_disjoint_core', pending.obs.review_action.eq('pending_CM').all() and not pending.obs_names.isin(core.obs_names).any())
    check(checks, 'doublet_exclusions_have_two_evidence_sources', dbl.obs.scDblFinder_class.eq('doublet').all() and dbl.obs.coherent_mixed_panel_flag.all())
    check(checks, 'allcell_no_excluded_doublets', not allcell.obs.exclude_suspected_doublet.any())
    check(checks, 'core_matches_audit_membership', set(core.obs_names) == set(identity.index[identity.review_action.eq('retain_CM')]))
    check(checks, 'pending_matches_audit_membership', set(pending.obs_names) == set(identity.index[identity.review_action.eq('pending_CM')]))
    check(checks, 'all_strict_audit_members_exported', set(identity.index) == set(allcell.obs_names) | set(dbl.obs_names))
    check(checks, 'doublets_disjoint_allcell', not dbl.obs_names.isin(allcell.obs_names).any())
    combined = ad.concat([allcell, dbl], join='inner', merge='same')[identity.index].copy()
    refdir = Path('/home/zhome/xiaoe-internship-data/processed/SCP1303_review_20261111/identity')
    rmeta = pd.read_csv(refdir / 'reference_metadata.tsv', sep='\t')
    cmeta = pd.read_csv(refdir / 'calibration_metadata.tsv', sep='\t')
    check(checks, 'reference_training_and_holdout_donors_disjoint', not rmeta.donor_id.isin(cmeta.donor_id).any())
    rg = pd.read_csv(refdir / 'reference_genes.tsv', sep='\t')
    with threadpool_limits(limits=2):
        reference_counts = mmread(refdir / 'calibration_counts.mtx.gz').T.tocsr()
    reference_totals = np.asarray(reference_counts.sum(axis=1)).ravel()
    reference_truth = cmeta.broad_type.replace({'Pericyte': 'Mural', 'Smooth_muscle': 'Mural'})
    count_total = np.asarray(combined.X.sum(axis=1)).ravel()
    programs, detected, repeated, bounds = {}, {}, {}, {}
    for lineage, geneset in policy['panels'].items():
        vals = combined.X[:, combined.var.gene_symbol.isin(geneset).to_numpy()].toarray()
        detected[lineage] = (vals > 0).sum(axis=1)
        repeated[lineage] = (vals >= 2).sum(axis=1)
        programs[lineage] = np.log1p(vals * (10000 / np.maximum(count_total, 1))[:, None]).mean(axis=1)
        check(checks, lineage + '_markers_recomputed', np.array_equal(detected[lineage], combined.obs[lineage + '_markers']) and np.array_equal(repeated[lineage], combined.obs[lineage + '_markers_ge2_UMI']))
        check(checks, lineage + '_program_recomputed', np.allclose(programs[lineage], combined.obs[lineage + '_program']))
        if lineage in policy['mixed_program_bounds']:
            rv = reference_counts[:, rg.gene_symbol.isin(geneset).to_numpy()].toarray()
            rs = np.log1p(rv * (10000 / np.maximum(reference_totals, 1))[:, None]).mean(axis=1)
            threshold = max(np.quantile(rs[reference_truth.eq(lineage)], .1), np.quantile(rs[reference_truth.ne(lineage)], .99))
            bounds[lineage] = threshold
            check(checks, lineage + '_strong_threshold_from_reference', np.isclose(threshold, policy['mixed_program_bounds'][lineage]['threshold']))
    weak = np.zeros(combined.n_obs, dtype=bool)
    strong = weak.copy()
    for left, right in policy['doublet_pairs']:
        ratio = np.minimum(programs[left], programs[right]) / np.maximum(np.maximum(programs[left], programs[right]), 1e-12)
        pair = (detected[left] >= 3) & (detected[right] >= 3) & (ratio >= .25)
        weak |= pair
        strong |= pair & (repeated[left] >= 3) & (repeated[right] >= 3) & (programs[left] >= bounds[left]) & (programs[right] >= bounds[right])
    check(checks, 'weak_flags_diagnostic_only', np.array_equal(weak, combined.obs.weak_coherent_mixed_panel_flag))
    check(checks, 'strong_mixed_flags_recomputed', np.array_equal(strong, combined.obs.coherent_mixed_panel_flag))
    check(checks, 'doublet_exclusion_rule_recomputed', np.array_equal(strong & combined.obs.scDblFinder_class.eq('doublet'), combined.obs.exclude_suspected_doublet))
    dominates = programs['Cardiomyocyte'] >= np.max(np.vstack([v for k, v in programs.items() if k != 'Cardiomyocyte']), axis=0)
    check(checks, 'CM_program_dominance_recomputed', np.array_equal(dominates, combined.obs.cm_program_dominates))
    ctable = pd.read_csv(WORK / 'reference_calibration.tsv', sep='\t', keep_default_na=False)
    check(checks, 'reference_truth_traceable_to_heldout', np.array_equal(ctable.truth, cmeta.set_index('reference_barcode').loc[ctable.barcode, 'broad_type']))
    tail = np.zeros(combined.n_obs, dtype=bool)
    for lineage, frame in ctable.groupby('truth'):
        correct = frame[frame.reference_label.eq(lineage)]
        if len(correct) >= 30:
            tail |= (combined.obs.reference_label.eq(lineage) & combined.obs.reference_pruned.eq(lineage) & combined.obs.reference_delta.ge(correct.reference_delta.quantile(.05)) & combined.obs.reference_delta_next.ge(correct.reference_delta_next.quantile(.05))).to_numpy()
    check(checks, 'reference_tail_flags_recomputed', np.array_equal(tail, combined.obs.reference_tail_pass))
    core_rule = combined.obs.reference_label.eq('Cardiomyocyte') & combined.obs.reference_pruned.eq('Cardiomyocyte') & (detected['Cardiomyocyte'] >= 2) & dominates & ~strong
    core_rule &= tail | ((detected['Cardiomyocyte'] >= 3) & combined.obs.cm_anchor_neighbor_fraction.ge(.5))
    check(checks, 'core_rule_complete_membership', set(combined.obs_names[core_rule]) == set(core.obs_names))
    candidate = combined.obs.reference_label.eq('Cardiomyocyte') | ((detected['Cardiomyocyte'] >= 3) & dominates)
    pending_rule = candidate & ~core_rule & ~combined.obs.exclude_suspected_doublet
    check(checks, 'pending_rule_complete_membership', set(combined.obs_names[pending_rule]) == set(pending.obs_names))
    diagnostic = ad.read_h5ad(output / 'diagnostic_normalized.h5ad')
    check(checks, 'diagnostic_explicitly_analysis_only', diagnostic.uns.get('analysis_only', False) and 'Diagnostic log1p CP10k' in diagnostic.uns.get('count_source', ''))
    check(checks, 'diagnostic_full_strict_axes', set(diagnostic.obs_names) == set(identity.index) and diagnostic.var_names.equals(allcell.var_names))
    original_diagnostic = ad.concat([allcell, dbl], join='inner', merge='same')[diagnostic.obs_names].X.astype(float).tocsr()
    size = np.asarray(original_diagnostic.sum(axis=1)).ravel()
    expected_diagnostic = original_diagnostic.multiply((10000 / np.maximum(size, 1))[:, None]).tocsr()
    expected_diagnostic.data = np.log1p(expected_diagnostic.data)
    diagnostic_difference = expected_diagnostic - diagnostic.X
    check(checks, 'diagnostic_values_match_declared_log_CP10k', np.max(np.abs(diagnostic_difference.data), initial=0) < 1e-5)
    check(checks, 'light_unique_no_padding', light.obs_names.is_unique and light.n_obs == min(20000, core.n_obs))
    check(checks, 'light_subset_core', light.obs_names.isin(core.obs_names).all())
    for name in ['cm_qc_full.h5ad', 'cm_qc_light.h5ad', 'cm_review_pending.h5ad']:
        obj = objects[name]
        check(checks, name + '_matches_allcell_X', exact(obj.X, allcell[obj.obs_names].X))
    for name, obj in objects.items():
        for score_col in ['reference_delta', 'reference_delta_next', 'reference_best_score', 'scDblFinder_score']:
            check(checks, name + '_' + score_col + '_numeric_dtype', pd.api.types.is_numeric_dtype(obj.obs[score_col]))
        check(checks, name + '_identity_scores_finite', np.isfinite(obj.obs[['reference_delta', 'reference_delta_next', 'reference_best_score']].to_numpy(dtype=float)).all())
        check(checks, name + '_unique_axes', obj.obs_names.is_unique and obj.var_names.is_unique)
        check(checks, name + '_strict_bounds', obj.obs.n_genes_by_counts.between(500,5000).all() and obj.obs.pct_counts_mt.lt(5).all())
        check(checks, name + '_cell_calling_uncertainty_explicit', obj.obs.cell_calling_status.eq('unverified_released_candidate_no_ambient_pool').all() and 'no independent EmptyDrops' in obj.uns.get('cell_calling_status',''))
        check(checks, name + '_finite_integer_nonnegative', np.isfinite(obj.X.data).all() and (obj.X.data >= 0).all() and (obj.X.data == np.round(obj.X.data)).all())
        check(checks, name + '_X_equals_counts', exact(obj.X, obj.layers['counts']))
        check(checks, name + '_no_training_split', 'split' not in obj.obs)
    strict_expected = set()
    all_source_features = {}
    release_total = 0
    for row in group.itertuples():
        sample = row.sample_id
        stage = json.loads((WORK / sample / 'STAGED.json').read_text())
        for feature, record in stage['sources'].items():
            source = Path(record['path'])
            with source.open('rb') as handle:
                sha = hashlib.file_digest(handle, 'sha256').hexdigest()
            check(checks, sample + '_' + feature + '_source_unchanged', source.stat().st_size == record['bytes'] and sha == record['sha256'])
        prefix = row.source_prefix
        var = pd.read_csv(prefix + '_features.tsv.gz', sep='\t', header=None, dtype=str, keep_default_na=False)
        barcodes = pd.read_csv(prefix + '_barcodes.tsv.gz', sep='\t', header=None, dtype=str)[0]
        with gzip.open(prefix + '_matrix.mtx.gz', 'rb') as handle, threadpool_limits(limits=2):
            raw = mmread(handle).T.tocsr()
        raw.sum_duplicates(); raw.eliminate_zeros()
        counts = np.asarray(raw.sum(axis=1)).ravel()
        genes = np.diff(raw.indptr)
        mt_counts = np.asarray(raw[:, var[1].str.startswith('MT-').to_numpy()].sum(axis=1)).ravel()
        mt = 100 * mt_counts / np.maximum(counts, 1)
        release_total += len(barcodes)
        source_space = hashlib.sha256(var.to_csv(sep='\t', header=False, index=False).encode()).hexdigest()[:12]
        all_source_features[source_space] = var
        check(checks, sample + '_feature_space_hash_from_original', source_space == row.feature_space_id)
        prefixed = pd.Index(sample + ':' + barcodes)
        check(checks, sample + '_no_standard_ambient_pool', not ((counts>0)&(counts<=100)).any())
        strict = (genes >= 500) & (genes <= 5000) & (mt < 5)
        strict_expected.update(prefixed[strict])
        pred = pd.read_csv(WORK / sample / 'identity.tsv', sep='\t', index_col='barcode', keep_default_na=False)
        part_identity = identity.loc[identity.sample_id.eq(sample)]
        check(checks, sample + '_reference_predictions_cover_strict', set(pred.index) == set(prefixed[strict]))
        for col in ['reference_label', 'reference_pruned']:
            check(checks, sample + '_' + col + '_prediction_trace', np.array_equal(part_identity[col].fillna('').astype(str), pred.loc[part_identity.index, col]))
        model = pd.read_csv(WORK / sample / 'doublets_v2.tsv', sep='\t')
        check(checks, sample + '_doublet_model_complete_ge500_capture', set(model.barcode) == set(prefixed[genes >= 500]) and model.barcode.is_unique)
        if (genes>=500).sum()>=100:
            check(checks, sample + '_doublet_scores_finite', np.isfinite(model.scDblFinder_score).all() and model.scDblFinder_class.isin(['singlet', 'doublet']).all())
        else:
            check(checks, sample + '_small_input_explicitly_unscored', model.scDblFinder_score.isna().all() and model.scDblFinder_class.eq('not_scored_insufficient_nuclei').all())
        audit = pd.read_csv(output / f'{sample}_released_barcode_audit.tsv.gz', sep='\t', index_col='barcode')
        check(checks, sample + '_complete_droplet_audit', audit.index.equals(prefixed) and np.array_equal(audit.total_counts, counts) and np.array_equal(audit.n_genes_by_counts, genes) and np.allclose(audit.pct_counts_mt, mt))
        check(checks, sample + '_numeric_exclusions_accounted', np.array_equal(audit.pass_strict_numeric, strict) and audit.loc[~strict, 'final_action'].eq('outside_strict_numeric').all())
        original_record=[r for r in official if r['sample_id']==sample]
        check(checks,sample+'_unique_frozen_source_record',len(original_record)==1)
        record=original_record[0]
        original_donor=record['donor_original'];expected_subtype=record['disease_subtype'];expected_disease=record['disease'];expected_label=record['label']
        feature_axis = pd.read_csv(output / 'feature_axes' / (row.feature_space_id + '.tsv'),sep='\t',header=None,dtype=str,keep_default_na=False)
        check(checks, sample + '_original_feature_axis_preserved', feature_axis.equals(var))
        for name, obj in objects.items():
            part = obj[obj.obs.sample_id.eq(sample)]
            positions = barcodes.to_numpy()
            positions = pd.Index(positions).get_indexer(part.obs.barcode_original.astype(str))
            stem = sample + '_' + name
            union_positions = obj.var_names.get_indexer(var[0])
            native_present = obj.var_names.isin(var[0])
            capture_positions = list(obj.uns['feature_availability_captures']).index(sample)
            check(checks, stem + '_all_original_genes_retained', (union_positions>=0).all() and len(union_positions)==len(var))
            check(checks, stem + '_missingness_mask_exact', np.array_equal(obj.varm['feature_available_by_capture'][:,capture_positions],native_present))
            check(checks, stem + '_unavailable_zero_padding_only', part.X[:,~native_present].nnz==0)
            check(checks, stem + '_all_barcodes_original', (positions >= 0).all() and np.array_equal(part.obs_names.astype(str), sample + ':' + part.obs.barcode_original.astype(str)))
            check(checks, stem + '_all_counts_original', exact(part.X[:,union_positions], raw[positions]))
            check(checks, stem + '_recomputed_QC', np.array_equal(part.obs.total_counts, counts[positions]) and np.array_equal(part.obs.n_genes_by_counts, genes[positions]) and np.allclose(part.obs.pct_counts_mt, mt[positions]))
            check(checks, stem + '_source_metadata', part.obs.donor_id.eq(record['donor_id']).all() and part.obs.disease.eq(expected_disease).all() and part.obs.label.eq(expected_label).all() and part.obs.capture_id.eq(sample).all() and part.obs.disease_subtype.eq(expected_subtype).all() and part.obs.feature_space_id.eq(row.feature_space_id).all())
            check(checks,stem+'_clinical_and_age_stratum',all(part.obs[c].astype(str).eq(str(record[c])).all() for c in ['clinical_stage','age_stratum','sex','age_years','technology','tissue']))
        print('Verified original count values: ' + sample, flush=True)
    check(checks, 'every_original_strict_candidate_accounted_for', strict_expected == set(allcell.obs_names) | set(dbl.obs_names))
    check(checks, 'study_not_merged', all(obj.obs.dataset.eq(dataset).all() for obj in objects.values()))
    check(checks, 'expected_original_sample_count', len(group) == (4 if dataset=='GSE255296' else 2) and group.donor_id.nunique()==len(group))
    union_ids = sorted(set.union(*[set(v[0]) for v in all_source_features.values()]))
    common_ids = set.intersection(*[set(v[0]) for v in all_source_features.values()])
    check(checks, 'all_union_genes_no_inner_loss', len(union_ids) == 36601 and all(np.array_equal(x.var_names, union_ids) for x in objects.values()))
    check(checks, 'shared_gene_mask_from_original', len(common_ids) == 36601 and np.array_equal(diagnostic.var.present_in_all_captures, diagnostic.var_names.isin(common_ids)))
    check(checks, 'graph_uses_shared_nonMT_features_only', not (diagnostic.var.highly_variable & (~diagnostic.var.present_in_all_captures | diagnostic.var.mt)).any())
    symbols = {}
    for key in reversed(sorted(all_source_features, key=lambda k: (-len(all_source_features[k]), k))):
        symbols.update(dict(zip(all_source_features[key][0], all_source_features[key][1])))
    check(checks, 'canonical_symbols_recomputed', np.array_equal(allcell.var.gene_symbol.astype(str), [symbols[g] for g in allcell.var_names]))
    anchor = combined.obs.reference_pruned.eq('Cardiomyocyte') & (detected['Cardiomyocyte'] >= 3) & dominates & ~strong
    anchor_series = pd.Series(anchor.to_numpy(dtype=float), index=combined.obs_names).loc[diagnostic.obs_names]
    graph = diagnostic.obsp['connectivities'].tocsr()
    support = np.asarray(graph @ anchor_series.to_numpy()).ravel() / np.maximum(np.asarray(graph.sum(axis=1)).ravel(), 1e-12)
    check(checks, 'CM_neighborhood_support_recomputed', np.allclose(pd.Series(support, index=diagnostic.obs_names).loc[combined.obs_names], combined.obs.cm_anchor_neighbor_fraction))
    hashes = {}
    for name, record in exports['hashes'].items():
        path = output / name
        md5 = hashlib.md5(); sha = hashlib.sha256()
        with path.open('rb') as handle:
            while chunk := handle.read(32 * 1024**2):
                md5.update(chunk); sha.update(chunk)
        hashes[name] = {'bytes': path.stat().st_size, 'md5': md5.hexdigest(), 'sha256': sha.hexdigest()}
        check(checks, name + '_size_and_export_hash', hashes[name]['bytes'] == record['bytes'] and hashes[name]['sha256'] == record['sha256'])
    status = 'passed' if all(checks.values()) else 'failed'
    result = {'status': status, 'dataset': dataset, 'verified_utc': datetime.now(timezone.utc).isoformat(),
              'checks': checks, 'hashes': hashes, 'all_cells_qc': allcell.n_obs, 'core_CM': core.n_obs,
              'light_CM': light.n_obs, 'pending_CM': pending.n_obs, 'suspected_doublets': dbl.n_obs,
              'input_donors': len(group), 'core_donors': int(core.obs.donor_id.nunique()), 'genes': core.n_vars,
              'input_released_barcodes': release_total, 'strict_candidates': len(strict_expected),
              'core_disease_counts': {d: int(core.obs.disease.eq(d).sum()) for d in ['Control','HCM']},
              'allcell_disease_counts': {d: int(allcell.obs.disease.eq(d).sum()) for d in ['Control','HCM']},
              'clinical_stage_core': core.obs.groupby('clinical_stage', observed=True).size().to_dict(),
              'feature_spaces': {key: len(value) for key, value in all_source_features.items()},
              'suitability': 'HCM-only supplementary resource; ambient cell calling unavailable; adult-reference domain shift and cross-study independence unresolved',
              'verifier_sha256': hashlib.sha256(Path(__file__).read_bytes()).hexdigest()}
    (output / 'INDEPENDENT_VERIFICATION.json').write_text(json.dumps(result, indent=2) + '\n')
    if status != 'passed':
        raise AssertionError([k for k, v in checks.items() if not v])
    for algorithm in ['md5', 'sha256']:
        (output / ('MD5SUMS' if algorithm == 'md5' else 'SHA256SUMS')).write_text(''.join(f"{r[algorithm]}  {n}\n" for n, r in hashes.items()))
    (output/'REPORT.md').write_text(f"# {dataset} QC verification\n\n{len(checks)} checks passed; {core.n_obs} core cardiomyocyte candidates, {pending.n_obs} pending; {core.obs.donor_id.nunique()}/{len(group)} source donors retained.\n\nProject nuclear baseline 500<=genes<=5000 and MT<5%; released counts preserved. Adult reference domain shift and missing ambient background remain explicit. HCM-only, split unassigned.\n")
    summary.append({k: v for k, v in result.items() if k not in ['checks', 'hashes']})
    print(f'VERIFIED {dataset}: {len(checks)} checks, core={core.n_obs}', flush=True)
(WORK / 'VERIFICATION_SUMMARY.json').write_text(json.dumps(summary, indent=2) + '\n')
