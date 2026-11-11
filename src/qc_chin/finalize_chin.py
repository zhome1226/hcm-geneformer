"""Join independent QC evidence, freeze explicit rules, and export unchanged counts."""
from datetime import datetime, timezone
from importlib.metadata import version
from pathlib import Path
import gzip
import hashlib
import json

import anndata as ad
import numpy as np
import pandas as pd
import scanpy as sc

ROOT = Path.home() / 'xiaoe-internship-data'
WORK = ROOT / '_processing/all_data_cleaning_20261111/chin_pair'
SEED = 20261111
PANELS = {
    'Cardiomyocyte': ['TNNT2', 'TTN', 'TNNI3', 'ACTC1', 'MYH7', 'MYBPC3', 'MYL2', 'TNNC1', 'RYR2', 'PLN', 'ACTN2', 'MYH6'],
    'Fibroblast': ['COL1A1', 'COL1A2', 'DCN', 'LUM', 'PDGFRA'],
    'Endothelial': ['PECAM1', 'VWF', 'KDR', 'EMCN', 'CLDN5'],
    'Immune': ['PTPRC', 'LST1', 'TYROBP', 'CD3D', 'CD79A'],
    'Mural': ['RGS5', 'PDGFRB', 'CSPG4', 'MYH11', 'ACTA2', 'TAGLN'],
    'Neural': ['NRXN1', 'SNAP25', 'SYT1', 'SLC1A3', 'PLP1'],
    'Adipocyte': ['ADIPOQ', 'FABP4', 'PLIN1', 'LPL', 'CIDEC'],
    'Mesothelial': ['WT1', 'MSLN', 'UPK3B', 'KRT8', 'KRT18'],
}
# Mural/fibroblast overlap is deliberately not used as a doublet combination.
PAIRS = [('Cardiomyocyte', 'Fibroblast'), ('Cardiomyocyte', 'Endothelial'),
         ('Cardiomyocyte', 'Immune'), ('Cardiomyocyte', 'Mural'),
         ('Fibroblast', 'Endothelial'), ('Fibroblast', 'Immune'),
         ('Endothelial', 'Immune'), ('Immune', 'Mural')]
samples = pd.read_csv(WORK / 'samples.tsv', sep='\t')
calibration = pd.read_csv(WORK / 'reference_calibration.tsv', sep='\t').fillna({'reference_pruned': ''})
tails = {}
for label, group in calibration.groupby('truth'):
    correct = group.loc[group['reference_label'].eq(label)]
    if len(correct) >= 30:
        tails[label] = {'n_correct': len(correct), 'delta_q05': float(correct.reference_delta.quantile(.05)),
                        'delta_next_q05': float(correct.reference_delta_next.quantile(.05))}
calibration_report = {'source': 'Litvinukova2020 adult left ventricular nuclei; same fixed training/held-out donor partitions as SCP1303 reference',
                      'accuracy': float(calibration.truth.eq(calibration.reference_label).mean()),
                      'common_genes': len((WORK / 'reference_common_genes.txt').read_text().splitlines()),
                      'thresholds': tails, 'per_label': [],
                      'limitation': 'Empirical reference tails are not posterior probabilities or independent disease validation; threshold selection and apparent tail metrics share the same held-out reference set.'}
for label in sorted(calibration.truth.unique()):
    truth = calibration.truth.eq(label)
    predicted = calibration.reference_label.eq(label)
    calibration_report['per_label'].append({'label': label, 'n_truth': int(truth.sum()),
                                            'precision': float((truth & predicted).sum() / max(1, predicted.sum())),
                                            'recall': float((truth & predicted).sum() / max(1, truth.sum()))})
(WORK / 'REFERENCE_CALIBRATION.json').write_text(json.dumps(calibration_report, indent=2) + '\n')

POLICY = {
    'version': 'Chin_pair_strict_v2', 'seed': SEED, 'min_genes': 500, 'max_genes': 5000,
    'mt_percent_exclusive_max': 5, 'count_source': 'Original GEO Cell Ranger gene-expression counts; no background correction, normalization, imputation or integration applied to delivered counts',
    'cell_calling': 'EmptyDrops lower=100, FDR<=.001, automatic knee retention; 10k iterations increased to 100k if strict candidate Limited and not significant',
    'doublet_model': 'Per GEO capture; all EmptyDrops-called nuclei with >=500 genes, including >5000 genes and high MT; scDblFinder dbr=NULL, dbr.per1k=.008 default loading-rate assumption',
    'doublet_input_revision': 'Initial >=200 input retained too many low-information called droplets and inflated the count-dependent loading rate (up to about .53). Revised before export to the project gene lower bound >=500; old calls archived and never used for final filtering.',
    'doublet_exclusion': 'scDblFinder doublet AND at least one listed distinct lineage pair with >=3 detected panel genes in both and smaller/larger mean log1p(CP10k) >=.25. These are project heuristics; calls are suspected, not experimentally confirmed.',
    'doublet_pairs': PAIRS, 'panels': PANELS,
    'identity': 'SingleR trained on adapted common genes from adult LV nuclear reference; fixed reference donors held out. Cell type assigned only with pruned reference label and >=2 matched panel genes; otherwise Unknown.',
    'CM_core': 'Pruned and best reference label both Cardiomyocyte, >=2 CM panel genes, CM mean log1p(CP10k) >= every other panel, and (empirical reference-tail pass OR >=3 CM genes and >=.5 weighted CM-anchor neighbor support). Strong multi-evidence doublet exclusions override.',
    'CM_pending': 'Strict numeric/cell-called, no strong doublet exclusion; best reference CM OR >=3 CM panel genes and CM program dominates, but core identity rule fails.',
    'graph': 'Each cohort separately; log-CP10k diagnostic copy, 3000 HVGs, up to 40 PCs, 15 neighbors; no disease labels or test labels used. CM anchors require pruned CM plus >=3 CM panel genes and panel dominance. Neighborhood evidence is correlated contextual support, not independent confidence.',
    'light': 'At most 20000 unique core nuclei, balanced allocations across available donors with fixed seed, never duplicated to fill a target',
    'disease_map': {'Control': 0, 'HCM': 1}, 'split_status': 'No train/test split assigned; external cohort status pending donor-overlap audit',
    'scope': 'Each study saved separately. No disease differential expression, tokenization or training.',
}


def marker_evidence(a):
    total = np.asarray(a.X.sum(axis=1)).ravel()
    for label, genes in PANELS.items():
        positions = np.flatnonzero(a.var.gene_symbol.isin(genes).to_numpy())
        assert len(positions) >= 3, label
        values = a.X[:, positions].toarray()
        a.obs[label + '_markers'] = (values > 0).sum(axis=1)
        a.obs[label + '_program'] = np.log1p(values * (1e4 / np.maximum(total, 1))[:, None]).mean(axis=1)
    cm = a.obs['Cardiomyocyte_program']
    a.obs['cm_program_dominates'] = cm.ge(a.obs[[k + '_program' for k in PANELS if k != 'Cardiomyocyte']].max(axis=1))
    mixed = np.zeros(a.n_obs, dtype=bool)
    reasons = np.full(a.n_obs, '', dtype=object)
    for left, right in PAIRS:
        first = a.obs[left + '_program'].to_numpy()
        second = a.obs[right + '_program'].to_numpy()
        ratio = np.minimum(first, second) / np.maximum(np.maximum(first, second), 1e-12)
        valid = a.obs[left + '_markers'].ge(3).to_numpy() & a.obs[right + '_markers'].ge(3).to_numpy() & (ratio >= .25)
        reasons[valid & ~mixed] = left + '+' + right
        mixed |= valid
    a.obs['coherent_mixed_panel_flag'] = mixed
    a.obs['mixed_panel_pair'] = reasons
    a.obs['exclude_suspected_doublet'] = a.obs.scDblFinder_class.eq('doublet') & mixed


def prepare_sample(row):
    directory = WORK / row.sample_id
    assert (directory / 'R_V2_SUCCESS.txt').exists()
    a = ad.read_h5ad(directory / 'prefilter.h5ad')
    ed = pd.read_csv(directory / 'emptydrops.tsv.gz', sep='\t', index_col='barcode',
                     usecols=['barcode', 'Total', 'FDR', 'Limited', 'emptydrops_called'])
    assert ed.index.is_unique and a.obs_names.is_unique
    ev = ed.reindex(a.obs_names)
    assert ev.Total.notna().all() and np.array_equal(a.obs.total_counts, ev.Total)
    a.obs['emptydrops_FDR'] = ev.FDR.to_numpy()
    a.obs['emptydrops_Limited'] = ev.Limited.fillna(False).to_numpy(dtype=bool)
    a.obs['emptydrops_called'] = ev.emptydrops_called.to_numpy(dtype=bool)
    strict = a.obs['pass_strict_numeric'] & a.obs.emptydrops_called
    a = a[strict].copy()
    identity = pd.read_csv(directory / 'identity.tsv', sep='\t', index_col='barcode').fillna({'reference_pruned': ''})
    doublets = pd.read_csv(directory / 'doublets_v2.tsv', sep='\t', index_col='barcode')
    assert set(identity.index) == set(a.obs_names) and a.obs_names.isin(doublets.index).all()
    a.obs = a.obs.join(identity).join(doublets)
    a.obs['reference_tail_pass'] = False
    for label, bounds in tails.items():
        passed = a.obs.reference_label.eq(label) & a.obs.reference_pruned.eq(label)
        passed &= a.obs.reference_delta.ge(bounds['delta_q05']) & a.obs.reference_delta_next.ge(bounds['delta_next_q05'])
        a.obs.loc[passed, 'reference_tail_pass'] = True
    marker_evidence(a)
    return a


def light_ids(obs):
    target = min(20000, len(obs))
    groups = {str(k): list(v.index) for k, v in obs.groupby('donor_id', observed=True)}
    allocation = {k: 0 for k in groups}
    for _ in range(target):
        available = [k for k in groups if allocation[k] < len(groups[k])]
        chosen = min(available, key=lambda k: (allocation[k], k))
        allocation[chosen] += 1
    rng = np.random.default_rng(SEED)
    selected = []
    for key in sorted(groups):
        selected.extend(rng.choice(groups[key], size=allocation[key], replace=False))
    return sorted(selected)


exported = []
for dataset, group in samples.groupby('dataset', sort=True):
    output = ROOT / f'processed/{dataset}_qc_20261111'
    output.mkdir(parents=True, exist_ok=True)
    if (output / 'EXPORT_SUCCESS.json').exists():
        existing = json.loads((output / 'EXPORT_SUCCESS.json').read_text())
        assert json.loads((output / 'qc_policy.json').read_text()) == json.loads(json.dumps(POLICY))
        for filename, record in existing['hashes'].items():
            path = output / filename
            with path.open('rb') as handle:
                actual = hashlib.file_digest(handle, 'sha256').hexdigest()
            assert path.stat().st_size == record['bytes'] and actual == record['sha256']
        exported.append(existing)
        print('Reusing existing count-preserving export with matching policy and digests: ' + dataset, flush=True)
        continue
    arrays = [prepare_sample(row) for row in group.itertuples()]
    assert all(a.var_names.equals(arrays[0].var_names) for a in arrays)
    a = ad.concat(arrays, join='inner', merge='same', uns_merge='same')
    assert a.n_vars == arrays[0].n_vars and a.obs_names.is_unique
    # NUMBA_NUM_THREADS=2 is the deliberate job limit in the background driver.
    # The approximate-neighbor backend enforces that ceiling on larger cohorts.
    sc.settings.n_jobs = 2
    diagnostic = a.copy()
    sc.pp.normalize_total(diagnostic, target_sum=1e4)
    sc.pp.log1p(diagnostic)
    sc.pp.highly_variable_genes(diagnostic, n_top_genes=3000, flavor='seurat')
    diagnostic.var.loc[diagnostic.var.mt, 'highly_variable'] = False
    pcs = min(40, diagnostic.n_obs - 1, int(diagnostic.var.highly_variable.sum()) - 1)
    sc.pp.pca(diagnostic, n_comps=pcs, mask_var='highly_variable', random_state=SEED)
    # The approximate backend produces non-contiguous float32 distances rejected
    # by this installed UMAP/Numba combination. Exact sklearn neighbors avoid that
    # backend and are tractable for these separately processed cohorts.
    sc.pp.neighbors(diagnostic, n_neighbors=min(15, a.n_obs - 1), n_pcs=pcs,
                    transformer='sklearn', random_state=SEED)
    sc.tl.leiden(diagnostic, resolution=.6, flavor='igraph', n_iterations=2, random_state=SEED)
    anchor = (a.obs.reference_pruned.eq('Cardiomyocyte') & a.obs.Cardiomyocyte_markers.ge(3)
              & a.obs.cm_program_dominates & ~a.obs.exclude_suspected_doublet)
    graph = diagnostic.obsp['connectivities'].tocsr()
    support = np.asarray(graph @ anchor.to_numpy(dtype=float)).ravel() / np.maximum(np.asarray(graph.sum(axis=1)).ravel(), 1e-12)
    a.obs['cm_anchor_neighbor_fraction'] = support
    a.obs['diagnostic_leiden'] = diagnostic.obs['leiden'].astype(str)
    a.obsm['X_pca_qc'] = diagnostic.obsm['X_pca']
    cm = a.obs.reference_label.eq('Cardiomyocyte') & a.obs.reference_pruned.eq('Cardiomyocyte')
    cm &= a.obs.Cardiomyocyte_markers.ge(2) & a.obs.cm_program_dominates
    cm &= a.obs.reference_tail_pass | (a.obs.Cardiomyocyte_markers.ge(3) & a.obs.cm_anchor_neighbor_fraction.ge(.5))
    cm &= ~a.obs.exclude_suspected_doublet
    candidate = a.obs.reference_label.eq('Cardiomyocyte') | (a.obs.Cardiomyocyte_markers.ge(3) & a.obs.cm_program_dominates)
    pending = candidate & ~cm & ~a.obs.exclude_suspected_doublet
    a.obs['cell_type_qc'] = 'Unknown'
    for label in sorted(a.obs.reference_label.unique()):
        panel = 'Mural' if label in ['Pericyte', 'Smooth_muscle'] else label
        if panel in PANELS:
            passed = a.obs.reference_pruned.eq(label) & a.obs[panel + '_markers'].ge(2)
            a.obs.loc[passed, 'cell_type_qc'] = label
    a.obs.loc[a.obs.cell_type_qc.eq('Cardiomyocyte'), 'cell_type_qc'] = 'Unknown'
    a.obs.loc[cm, 'cell_type_qc'] = 'Cardiomyocyte'
    a.obs['review_action'] = np.select([a.obs.exclude_suspected_doublet, cm, pending],
                                      ['exclude_suspected_doublet', 'retain_CM', 'pending_CM'],
                                      default='retain_other_or_unknown')
    a.obs.to_csv(output / 'cell_identity_audit.tsv.gz', sep='\t', index_label='barcode', compression='gzip')
    counts_qc = []
    for row in group.itertuples():
        directory = WORK / row.sample_id
        audit = pd.read_csv(directory / 'droplet_qc.tsv.gz', sep='\t', index_col='barcode')
        ed = pd.read_csv(directory / 'emptydrops.tsv.gz', sep='\t', index_col='barcode',
                         usecols=['barcode', 'FDR', 'Limited', 'emptydrops_called'])
        assert audit.index.equals(ed.index)
        audit = audit.join(ed)
        status = a.obs.review_action.reindex(audit.index)
        audit['final_action'] = np.where(~audit.emptydrops_called, 'not_cell_called',
                                         np.where(~audit.pass_strict_numeric, 'outside_strict_numeric', status.fillna('UNRESOLVED')))
        assert not audit.final_action.eq('UNRESOLVED').any()
        audit.to_csv(output / f'{row.sample_id}_droplet_audit.tsv.gz', sep='\t', compression='gzip')
        ro = pd.read_csv(directory / 'R_V2_SUMMARY.tsv', sep='\t').iloc[0].to_dict()
        current = a.obs.loc[a.obs.sample_id.eq(row.sample_id)]
        counts_qc.append({**ro, 'donor_id': row.donor_id, 'disease': row.disease,
                          'numeric_before_calling': int(audit.pass_strict_numeric.sum()),
                          'cell_called_ge500': int((audit.emptydrops_called & audit.n_genes_by_counts.ge(500)).sum()),
                          'called_500to5000': int((audit.emptydrops_called & audit.n_genes_by_counts.between(500, 5000)).sum()),
                          'core_CM': int(current.review_action.eq('retain_CM').sum()),
                          'pending_CM': int(current.review_action.eq('pending_CM').sum()),
                          'suspected_doublets': int(current.exclude_suspected_doublet.sum()),
                          'algorithmic_doublet_calls': int(current.scDblFinder_class.eq('doublet').sum()),
                          'median_genes_strict': float(current.n_genes_by_counts.median()),
                          'median_MT_strict': float(current.pct_counts_mt.median())})
    coverage = pd.DataFrame(counts_qc)
    coverage.to_csv(output / 'sample_qc_summary.tsv', sep='\t', index=False)
    group.to_csv(output / 'sample_metadata.tsv', sep='\t', index=False)
    a.obs.groupby(['sample_id', 'diagnostic_leiden', 'reference_label', 'review_action'], observed=True).size().to_csv(
        output / 'identity_group_summary.tsv', sep='\t', header=['nuclei'])
    a.layers['counts'] = a.X.copy()
    a.uns['qc_policy_json'] = json.dumps(POLICY)
    a.uns['reference_calibration_json'] = json.dumps(calibration_report)
    a.uns['source_file_digests_json'] = json.dumps([json.loads((WORK / row.sample_id / 'STAGED.json').read_text())['sources'] for row in group.itertuples()])
    outputs = {'all_cells_qc.h5ad': ~a.obs.exclude_suspected_doublet,
               'cm_qc_full.h5ad': cm, 'cm_review_pending.h5ad': pending,
               'suspected_doublets.h5ad': a.obs.exclude_suspected_doublet}
    for filename, mask in outputs.items():
        obj = a[mask].copy()
        obj.write_h5ad(output / filename, compression='gzip')
    obj = a[light_ids(a.obs.loc[cm])].copy()
    obj.write_h5ad(output / 'cm_qc_light.h5ad', compression='gzip')
    diagnostic.obs['review_action'] = a.obs.review_action
    diagnostic.uns['count_source'] = 'Diagnostic log1p CP10k expression; original counts are in the separately delivered count-preserving objects, not this X'
    diagnostic.uns['analysis_only'] = True
    diagnostic.write_h5ad(output / 'diagnostic_normalized.h5ad', compression='gzip')
    (output / 'qc_policy.json').write_text(json.dumps(POLICY, indent=2) + '\n')
    (output / 'REFERENCE_CALIBRATION.json').write_text(json.dumps(calibration_report, indent=2) + '\n')
    hashes = {}
    for path in sorted(output.glob('*.h5ad')):
        with path.open('rb') as handle:
            sha = hashlib.file_digest(handle, 'sha256').hexdigest()
        hashes[path.name] = {'bytes': path.stat().st_size, 'sha256': sha}
    result = {'status': 'exported_pending_independent_verification', 'dataset': dataset,
              'completed_utc': datetime.now(timezone.utc).isoformat(), 'strict_nuclei': a.n_obs,
              'all_cells_qc': int((~a.obs.exclude_suspected_doublet).sum()), 'core_CM': int(cm.sum()),
              'light_CM': obj.n_obs, 'pending_CM': int(pending.sum()),
              'suspected_doublets': int(a.obs.exclude_suspected_doublet.sum()),
              'genes': a.n_vars, 'input_donors': len(group),
              'core_donors': int(a.obs.loc[cm, 'donor_id'].nunique()), 'hashes': hashes,
              'neighbors_transformer': 'sklearn (exact)',
              'python_versions': {name: version(name) for name in ['scanpy', 'anndata', 'numpy', 'scipy', 'pandas']}}
    (output / 'EXPORT_SUCCESS.json').write_text(json.dumps(result, indent=2) + '\n')
    exported.append(result)
    print(json.dumps({k: v for k, v in result.items() if k != 'hashes'}, indent=2), flush=True)
(WORK / 'EXPORT_SUMMARY.json').write_text(json.dumps(exported, indent=2) + '\n')
print('CHIN_PAIR_EXPORTED', flush=True)
