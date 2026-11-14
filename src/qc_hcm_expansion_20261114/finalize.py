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

ROOT = Path('/home/zhome/hcm-geneformer-data')
WORK = Path('/home/zhome/hcm-geneformer-analysis/hcm_expansion_20261114/nuclear_qc')
SEED = 20261114
MIXED_REVIEW = json.loads((WORK / 'MIXED_PROGRAM_REVIEW.json').read_text())
MIXED_BOUNDS = MIXED_REVIEW['bounds']
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
union = pd.read_csv(WORK / 'UNION_FEATURES.tsv', sep='\t').set_index('gene_id')
common_ids = set((WORK / 'COMMON_GENE_IDS.txt').read_text().splitlines())
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

POLICY = json.loads('{"version": "HCM_expansion_snRNA_20261114_v1", "seed": 20261114, "min_genes": 500, "max_genes": 5000, "mt_percent_exclusive_max": 5, "count_source": "Released GEO nonnegative integer matrices; upstream SoupX or other correction not established by integer appearance. Deliverable X/counts remain unchanged; all native features retained.", "cell_calling": "No independent EmptyDrops: no positive-UMI<=100 background barcodes in any released matrix. Provenance-screened released candidate barcodes; ambient correction not rerun without unfiltered background.", "author_QC": "GSE255296 GEO Cell Ranger6.0.1; paper PMC12054467 uses SoupX1.6.1 and DoubletFinder2.0.3 but application to deposited matrices is unverified. GSE292754 only new P186/P194 retained. This batch applies agreed project nuclear numeric baseline and audits sample loss.", "doublet_model": "scDblFinder per capture on all >=500-gene released barcodes, including high MT and >5000 genes; dbr=NULL,dbr.per1k=.008. Screening heuristic: true loading and upstream doublet history unknown.", "doublet_exclusion": "Algorithmic doublet AND strong mixed pair (3 genes each >=2 UMI in each lineage, both program means >=max(reference positive q10,negative q99), smaller/larger>=.25). Algorithm-only calls retained if other CM evidence passes; same-lineage doublets not excluded reliably.", "identity": "Fixed adult-LV reference and donor-held-out calibration; SingleR raw/pruned labels, marker panels and contextual support. Pediatric/region domain shift remains explicit.", "CM_core": "Best/pruned reference both CM, >=2 detected CM markers, CM program dominance, reference-tail pass OR >=3 CM markers and >=.5 CM-anchor neighborhood support; no strong mixed flag.", "CM_pending": "Best-reference CM or >=3 CM markers and CM program dominance, failing core and not jointly excluded; not tokenized.", "graph": "Separate log1p CP10k diagnostic copy; <=3000 shared non-MT HVGs, <=40 PCs, 15 exact neighbors, Leiden .6; no disease-based integration.", "feature_space": "Native Ensembl union, full feature axes and per-capture availability mask; structural missingness is not observed zero.", "light": "<=20000 unique CM cores, balanced across source donor identifiers; no duplication.", "disease_map": {"Control": 0, "HCM": 1}, "split_status": "unassigned;model_history_and_cross_study_donor_review_pending", "replicate_policy": "GSE255296 has four deposited human samples only; reused Chin controls/HCM excluded from new download. GSE292754 P86 and reused controls excluded. Cross-study donor independence not proven.", "clinical_stage": "GSE255296 septal myectomy; age/sex unresolved. GSE292754 pediatric end-stage LV explants P186/P194. Both released subsets HCM-only; not standalone HCM-v-Control cohorts.", "scope": "QC then separate Geneformer V1 serialization only; no DGE, training or in-silico perturbation."}')
POLICY.update(mixed_program_bounds=MIXED_BOUNDS, doublet_pairs=PAIRS, panels=PANELS)


def marker_evidence(a):
    total = np.asarray(a.X.sum(axis=1)).ravel()
    for label, genes in PANELS.items():
        positions = np.flatnonzero(a.var.gene_symbol.isin(genes).to_numpy())
        assert len(positions) >= 3, label
        values = a.X[:, positions].toarray()
        a.obs[label + '_markers'] = (values > 0).sum(axis=1)
        a.obs[label + '_markers_ge2_UMI'] = (values >= 2).sum(axis=1)
        a.obs[label + '_program'] = np.log1p(values * (1e4 / np.maximum(total, 1))[:, None]).mean(axis=1)
    cm = a.obs['Cardiomyocyte_program']
    a.obs['cm_program_dominates'] = cm.ge(a.obs[[k + '_program' for k in PANELS if k != 'Cardiomyocyte']].max(axis=1))
    mixed = np.zeros(a.n_obs, dtype=bool)
    weak = mixed.copy()
    reasons = np.full(a.n_obs, '', dtype=object)
    weak_reasons = reasons.copy()
    for left, right in PAIRS:
        first = a.obs[left + '_program'].to_numpy()
        second = a.obs[right + '_program'].to_numpy()
        ratio = np.minimum(first, second) / np.maximum(np.maximum(first, second), 1e-12)
        valid = a.obs[left + '_markers'].ge(3).to_numpy() & a.obs[right + '_markers'].ge(3).to_numpy() & (ratio >= .25)
        weak_reasons[valid & ~weak] = left + '+' + right
        weak |= valid
        valid &= a.obs[left + '_markers_ge2_UMI'].ge(3).to_numpy() & a.obs[right + '_markers_ge2_UMI'].ge(3).to_numpy()
        valid &= first >= MIXED_BOUNDS[left]['threshold']
        valid &= second >= MIXED_BOUNDS[right]['threshold']
        reasons[valid & ~mixed] = left + '+' + right
        mixed |= valid
    a.obs['coherent_mixed_panel_flag'] = mixed
    a.obs['weak_coherent_mixed_panel_flag'] = weak
    a.obs['weak_mixed_panel_pair'] = weak_reasons
    a.obs['mixed_panel_pair'] = reasons
    a.obs['exclude_suspected_doublet'] = a.obs.scDblFinder_class.eq('doublet') & mixed


def prepare_sample(row):
    directory = WORK / row.sample_id
    assert (directory / 'R_SUCCESS.txt').exists()
    a = ad.read_h5ad(directory / 'prefilter.h5ad')
    a.obs['cell_calling_status'] = 'unverified_released_candidate_no_ambient_pool'
    strict = a.obs['pass_strict_numeric']
    a = a[strict].copy()
    # A sample with zero strict candidates still has a header-only identity table.
    # Explicit numeric dtypes prevent that empty table from turning real scores
    # from the other captures into object columns during outer concatenation.
    score_types = {key: float for key in ['reference_delta', 'reference_delta_next', 'reference_best_score']}
    identity = pd.read_csv(directory / 'identity.tsv', sep='\t', index_col='barcode',
                           dtype=score_types).fillna({'reference_pruned': ''})
    doublets = pd.read_csv(directory / 'doublets_v2.tsv', sep='\t', index_col='barcode',
                          dtype={'scDblFinder_score': float})
    assert set(identity.index) == set(a.obs_names) and a.obs_names.isin(doublets.index).all()
    a.obs = a.obs.join(identity).join(doublets)
    a.obs['reference_tail_pass'] = False
    for label, bounds in tails.items():
        passed = a.obs.reference_label.eq(label) & a.obs.reference_pruned.eq(label)
        passed &= a.obs.reference_delta.ge(bounds['delta_q05']) & a.obs.reference_delta_next.ge(bounds['delta_next_q05'])
        a.obs.loc[passed, 'reference_tail_pass'] = True
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
    output = ROOT / f'processed/{dataset}_qc_20261114'
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
    a = ad.concat(arrays, join='outer', fill_value=0, merge='first')
    assert set(a.var_names) == set(union.index) and a.obs_names.is_unique
    a = a[:, union.index].copy()
    a.var = union[['gene_symbol']].copy()
    a.var['mt'] = a.var.gene_symbol.str.startswith('MT-')
    a.var['present_in_all_captures'] = a.var_names.isin(common_ids)
    availability = np.column_stack([union['available_' + row.feature_space_id].to_numpy(dtype=bool) for row in group.itertuples()])
    a.varm['feature_available_by_capture'] = availability
    a.uns['feature_availability_captures'] = group.sample_id.to_numpy(dtype=str)
    a.uns['feature_availability_note'] = 'Unavailable source features are structural missingness, represented by sparse zero padding plus this explicit mask; do not use padded entries as observed zero targets.'
    a.uns['cell_calling_status'] = 'Released candidates only; no independent EmptyDrops/ambient-based nucleus calls possible from provided truncated inputs'
    marker_evidence(a)
    # NUMBA_NUM_THREADS=2 is the deliberate job limit in the background driver.
    # The approximate-neighbor backend enforces that ceiling on larger cohorts.
    sc.settings.n_jobs = 2
    diagnostic = a.copy()
    sc.pp.normalize_total(diagnostic, target_sum=1e4)
    sc.pp.log1p(diagnostic)
    sc.pp.highly_variable_genes(diagnostic, n_top_genes=3000, flavor='seurat')
    diagnostic.var.loc[diagnostic.var.mt | ~diagnostic.var.present_in_all_captures, 'highly_variable'] = False
    pcs = min(40, diagnostic.n_obs - 1, int(diagnostic.var.highly_variable.sum()) - 1)
    sc.pp.pca(diagnostic, n_comps=pcs, mask_var='highly_variable', random_state=SEED)
    sc.pp.neighbors(diagnostic, n_neighbors=min(15, a.n_obs - 1), n_pcs=pcs,
                    transformer='sklearn', random_state=SEED)
    sc.tl.leiden(diagnostic, resolution=.6, flavor='igraph', n_iterations=2, random_state=SEED)
    anchor = (a.obs.reference_pruned.eq('Cardiomyocyte') & a.obs.Cardiomyocyte_markers.ge(3)
              & a.obs.cm_program_dominates & ~a.obs.coherent_mixed_panel_flag)
    graph = diagnostic.obsp['connectivities'].tocsr()
    support = np.asarray(graph @ anchor.to_numpy(dtype=float)).ravel() / np.maximum(np.asarray(graph.sum(axis=1)).ravel(), 1e-12)
    a.obs['cm_anchor_neighbor_fraction'] = support
    a.obs['diagnostic_leiden'] = (diagnostic.obs['leiden'] if 'leiden' in diagnostic.obs else diagnostic.obs['diagnostic_leiden']).astype(str)
    a.obsm['X_pca_qc'] = diagnostic.obsm['X_pca']
    cm = a.obs.reference_label.eq('Cardiomyocyte') & a.obs.reference_pruned.eq('Cardiomyocyte')
    cm &= a.obs.Cardiomyocyte_markers.ge(2) & a.obs.cm_program_dominates
    cm &= a.obs.reference_tail_pass | (a.obs.Cardiomyocyte_markers.ge(3) & a.obs.cm_anchor_neighbor_fraction.ge(.5))
    cm &= ~a.obs.exclude_suspected_doublet
    cm &= ~a.obs.coherent_mixed_panel_flag
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
        status = a.obs.review_action.reindex(audit.index)
        audit['final_action'] = np.where(~audit.pass_strict_numeric, 'outside_strict_numeric', status.fillna('UNRESOLVED'))
        assert not audit.final_action.eq('UNRESOLVED').any()
        audit['cell_calling_status'] = 'unverified_released_candidate_no_ambient_pool'
        audit.to_csv(output / f'{row.sample_id}_released_barcode_audit.tsv.gz', sep='\t', compression='gzip')
        ro = pd.read_csv(directory / 'R_SUMMARY.tsv', sep='\t').iloc[0].to_dict()
        current = a.obs.loc[a.obs.sample_id.eq(row.sample_id)]
        counts_qc.append({**ro, 'donor_id': row.donor_id, 'disease': row.disease,
                          'input_strict_numeric': int(audit.pass_strict_numeric.sum()),
                          'input_ge500': int(audit.n_genes_by_counts.ge(500).sum()),
                          'input_500to5000': int(audit.n_genes_by_counts.between(500, 5000).sum()),
                          'disease_subtype': row.disease_subtype, 'clinical_stage': row.clinical_stage,
                          'feature_space_id': row.feature_space_id, 'source_genes': row.source_gene_count,
                          'released_below_500_genes': int(audit.n_genes_by_counts.lt(500).sum()),
                          'released_above_5000_genes': int(audit.n_genes_by_counts.gt(5000).sum()),
                          'released_MT_ge5': int(audit.pct_counts_mt.ge(5).sum()),
                          'core_CM': int(current.review_action.eq('retain_CM').sum()),
                          'pending_CM': int(current.review_action.eq('pending_CM').sum()),
                          'suspected_doublets': int(current.exclude_suspected_doublet.sum()),
                          'algorithmic_doublet_calls': int(current.scDblFinder_class.eq('doublet').sum()),
                          'median_genes_strict': float(current.n_genes_by_counts.median()),
                          'median_MT_strict': float(current.pct_counts_mt.median())})
    coverage = pd.DataFrame(counts_qc)
    coverage.to_csv(output / 'sample_qc_summary.tsv', sep='\t', index=False)
    group.to_csv(output / 'sample_metadata.tsv', sep='\t', index=False)
    a.obs.groupby(['donor_id', 'disease_subtype', 'clinical_stage', 'review_action'], observed=True).size().to_csv(
        output / 'clinical_stage_retention.tsv', sep='\t', header=['nuclei'])
    axes = output / 'feature_axes'
    axes.mkdir(exist_ok=True)
    for path in (WORK / 'feature_axes').glob('*.tsv'):
        (axes / path.name).write_bytes(path.read_bytes())
    union.to_csv(output / 'union_feature_metadata.tsv', sep='\t')
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
    (output / 'MIXED_PROGRAM_REVIEW.json').write_text(json.dumps(MIXED_REVIEW, indent=2) + '\n')
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
              'genes': a.n_vars, 'input_donors': int(group.donor_id.nunique()), 'input_samples': len(group),
              'core_donors': int(a.obs.loc[cm, 'donor_id'].nunique()), 'hashes': hashes,
              'neighbors_transformer': 'sklearn (exact)',
              'python_versions': {name: version(name) for name in ['scanpy', 'anndata', 'numpy', 'scipy', 'pandas']}}
    (output / 'EXPORT_SUCCESS.json').write_text(json.dumps(result, indent=2) + '\n')
    exported.append(result)
    print(json.dumps({k: v for k, v in result.items() if k != 'hashes'}, indent=2), flush=True)
(WORK / 'EXPORT_SUMMARY.json').write_text(json.dumps(exported, indent=2) + '\n')
print('EXPANSION_EXPORTED', flush=True)
