"""Auditable HCM/NF cardiomyocyte filtering of an author-QC'd SCP1303 H5AD.

Custom workflow is required for paired count layers, reused author QC,
percent/fraction auditing, full-gene preservation and donor-aware outputs.
Does not perform differential expression, integration, tokenization or training.
"""
from __future__ import annotations

import argparse
import gc
import hashlib
import importlib.metadata
import json
import os
from pathlib import Path
import sys
from datetime import datetime, timezone

import anndata as ad
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc
from scipy import sparse

from qc_core import baseline_mask, balanced_donor_sample, make_donor_splits


def log(message):
    print(datetime.now(timezone.utc).isoformat(), message, flush=True)


def write_json(path, value):
    Path(path).write_text(json.dumps(value, indent=2, ensure_ascii=False, default=str) + '\n')


def audit_matrix(matrix, name, require_integer=False):
    if not sparse.issparse(matrix):
        raise ValueError(f'{name}: expected sparse source counts')
    matrix = matrix.tocsr()
    matrix.eliminate_zeros()
    fractional = 0
    for start in range(0, matrix.data.size, 5_000_000):
        block = matrix.data[start:start+5_000_000]
        if not np.isfinite(block).all() or (block < 0).any():
            raise ValueError(f'{name}: negative/nonfinite count values')
        fractional += int(np.count_nonzero(block != np.rint(block)))
    if require_integer and fractional:
        raise ValueError(f'{name}: noninteger raw counts')
    return matrix, {'nnz': int(matrix.nnz), 'dtype': str(matrix.dtype),
                    'fractional_nonzero_values': fractional}


def grouped_stages(obs, stages, columns):
    tables = []
    for name, keep in stages.items():
        table = obs.loc[keep].groupby(columns, observed=True, dropna=False).size().rename('nuclei').reset_index()
        table['stage'] = name
        tables.append(table)
    return pd.concat(tables, ignore_index=True)


def marker_audit(a, obs, reports):
    marker_sets = {'cardiomyocyte': ['TNNT2', 'TTN', 'ACTC1', 'TNNI3', 'MYH7'],
                   'fibroblast': ['COL1A1', 'DCN'], 'endothelial': ['PECAM1', 'VWF'],
                   'immune': ['PTPRC']}
    markers = [g for genes in marker_sets.values() for g in genes]
    present = [g for g in markers if g in a.var_names]
    if not {'TNNT2', 'TTN', 'MYH7'}.issubset(present):
        raise ValueError('Required cardiomyocyte marker symbols missing')
    values = a[:, present].X.toarray()
    counts = obs['total_counts'].to_numpy()
    logged = np.log1p(values * (1e4 / np.maximum(counts, 1))[:, None])
    groups = obs['cell_type_original'].astype(str)
    means = pd.DataFrame(logged, columns=present).groupby(groups.to_numpy(), observed=True).mean()
    detection = pd.DataFrame(values > 0, columns=present).groupby(groups.to_numpy(), observed=True).mean()
    means.to_csv(reports/'marker_mean_log1p_cp10k_by_author_type.tsv', sep='\t')
    detection.to_csv(reports/'marker_detection_fraction_by_author_type.tsv', sep='\t')
    mean_scores = {}
    for name, genes in marker_sets.items():
        positions = [present.index(g) for g in genes if g in present]
        mean_scores[name] = logged[:, positions].mean(axis=1)
    obs['cm_marker_mean_log1p_cp10k'] = mean_scores['cardiomyocyte']
    obs['noncm_marker_max_mean_log1p_cp10k'] = np.max(
        np.vstack([v for k, v in mean_scores.items() if k != 'cardiomyocyte']), axis=0)
    obs['marker_mixed_warning'] = obs['noncm_marker_max_mean_log1p_cp10k'] > obs['cm_marker_mean_log1p_cp10k']
    cm_positions = [present.index(g) for g in marker_sets['cardiomyocyte'] if g in present]
    obs['n_cm_markers_detected'] = (values[:, cm_positions] > 0).sum(axis=1)
    fig, ax = plt.subplots(figsize=(10, max(5, .28*len(means))))
    im = ax.imshow(means.to_numpy(), aspect='auto', cmap='viridis')
    ax.set_xticks(range(len(present)), present, rotation=60, ha='right')
    ax.set_yticks(range(len(means)), means.index)
    ax.set_title('Author cell types: marker identity audit')
    fig.colorbar(im, ax=ax, label='Mean log1p(counts / total * 10,000)')
    fig.tight_layout(); fig.savefig(reports/'marker_identity.png', dpi=160); plt.close(fig)
    return {'present': present, 'missing': sorted(set(markers)-set(present)),
            'mixed_warning_policy': 'diagnostic only; ambient RNA and marker dropout preclude a hard single-gene rule'}


def qc_plots(obs, final, reports):
    target = obs['disease'].isin(['NF', 'HCM']) & obs['pass_cm_annotation']
    fig, axes = plt.subplots(2, 3, figsize=(14, 8))
    metrics = ['n_genes_by_counts', 'total_counts', 'pct_counts_mt']
    for row, mask in enumerate([target, final]):
        for col, metric in enumerate(metrics):
            for disease, color in [('NF', '#247BA0'), ('HCM', '#D95D39')]:
                vals = obs.loc[mask & (obs['disease'] == disease), metric].to_numpy()
                if metric == 'total_counts':
                    vals = np.log10(1 + vals)
                axes[row, col].hist(vals, bins=60, histtype='step', density=True,
                                    color=color, label=disease, linewidth=1.4)
            axes[row, col].set_xlabel('log10(1 + total_counts)' if metric == 'total_counts' else metric)
            axes[row, col].set_ylabel('Density')
            axes[row, col].set_title(('Before project filters' if row == 0 else 'After project filters'))
            if metric == 'n_genes_by_counts':
                axes[row, col].axvline(500, color='grey', ls='--')
                axes[row, col].axvline(5000, color='grey', ls='--')
            if metric == 'pct_counts_mt':
                axes[row, col].axvline(5, color='grey', ls='--')
    axes[0, 0].legend(); fig.tight_layout()
    fig.savefig(reports/'qc_before_after.png', dpi=160); plt.close(fig)
    before = obs.loc[target].groupby('donor_id', observed=True).size()
    after = obs.loc[final].groupby('donor_id', observed=True).size().reindex(before.index, fill_value=0)
    fig, ax = plt.subplots(figsize=(12, 5))
    ax.bar(before.index.astype(str), 100*after/before, color='#247BA0')
    ax.set_ylabel('Retained cardiomyocyte nuclei (%)'); ax.set_ylim(0, 100)
    ax.tick_params(axis='x', rotation=70); fig.tight_layout()
    fig.savefig(reports/'donor_retention.png', dpi=160); plt.close(fig)


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--input', required=True, type=Path)
    p.add_argument('--output', required=True, type=Path)
    p.add_argument('--seed', type=int, default=20261111)
    p.add_argument('--light-size', type=int, default=20000)
    args = p.parse_args()
    out = args.output.resolve()
    out.mkdir(parents=True, exist_ok=True)
    if (out/'SUCCESS.json').exists() or any(out.glob('cm_qc_*.h5ad')):
        raise RuntimeError('Refusing to overwrite existing processed objects; choose a new run directory')
    reports = out/'reports'; reports.mkdir(exist_ok=True)
    sc.settings.n_jobs = 4
    params = {'min_genes_inclusive': 500, 'max_genes_inclusive': 5000,
        'mt_percent_strict_upper': 5.0, 'gene_count_source': 'CellBender-adjusted X',
        'mt_count_source': "layers['cellranger_raw']; recomputed from MT- symbols",
        'seed': args.seed, 'light_target': args.light_size,
        'doublet_qc': 'reuse author per-sample Scrublet QC; no new global threshold or new call',
        'annotation': 'author cell_type_leiden0.6 containing cardiomyocyte; require nonmissing SubCluster',
        'split_status': 'provisional_70_15_15_by_donor; repository ratio pending confirmation',
        'cross_study_donor_overlap': 'not yet resolved; external validation not merged',
        'feature_policy': 'all original genes retained; no imputation, HVG restriction or normalization of saved X'}
    write_json(reports/'PARAMETERS.json', params)
    log('Loading source working H5AD')
    a = ad.read_h5ad(args.input)
    required = ['donor_id', 'biosample_id', 'disease', 'cell_type_leiden0.6',
                'SubCluster', 'cellranger_doublet_scores']
    if not set(required).issubset(a.obs):
        raise ValueError(f'Missing expected author metadata: {set(required)-set(a.obs)}')
    if not a.obs_names.is_unique or a.obs['donor_id'].isna().any():
        raise ValueError('Duplicate cell IDs or missing donors')
    if (a.obs.groupby('donor_id', observed=True)['disease'].nunique() != 1).any():
        raise ValueError('Conflicting donor diagnosis')
    if 'cellranger_raw' not in a.layers:
        raise ValueError('Missing raw counts for mitochondrial QC')
    if not a.var_names.is_unique:
        raise ValueError('Gene symbols duplicated: resolve mapping explicitly before QC')
    a.X, cb_stats = audit_matrix(a.X, 'CellBender X')
    a.layers['cellranger_raw'], raw_stats = audit_matrix(a.layers['cellranger_raw'], 'CellRanger raw', True)
    a.var['mt'] = a.var_names.str.startswith('MT-')
    if int(a.var['mt'].sum()) < 10:
        raise ValueError('Too few mitochondrial symbols; cannot use a zero/invalid MT fraction')
    log('Calculating QC on CellBender and Cell Ranger count layers')
    sc.pp.calculate_qc_metrics(a, qc_vars=['mt'], percent_top=None, log1p=False, inplace=True)
    raw_qc, _ = sc.pp.calculate_qc_metrics(a, qc_vars=['mt'], layer='cellranger_raw',
                                         percent_top=None, log1p=False, inplace=False)
    o = a.obs.copy()
    o['cell_type_original'] = o['cell_type_leiden0.6'].astype(str)
    o['pct_counts_mt_cellbender'] = o['pct_counts_mt']
    o['pct_counts_mt'] = raw_qc['pct_counts_mt']
    o['total_counts_cellranger'] = raw_qc['total_counts']
    o['n_genes_cellranger'] = raw_qc['n_genes_by_counts']
    o['dataset'] = 'SCP1303'
    o['sample_id'] = o['biosample_id'].astype(str)
    o['capture_id'] = o['sample_id']
    o['barcode_original'] = o.index.astype(str)
    o['disease_original'] = o['disease'].astype(str)
    o['doublet_score'] = o['cellranger_doublet_scores'].astype(float)
    o['doublet_qc_status'] = 'author_qc_reused_no_new_call'
    o['pass_author_qc'] = True
    o['pass_min_genes'] = np.isfinite(o['n_genes_by_counts']) & (o['n_genes_by_counts'] >= 500)
    o['pass_max_genes'] = np.isfinite(o['n_genes_by_counts']) & (o['n_genes_by_counts'] <= 5000)
    o['pass_mito'] = np.isfinite(o['pct_counts_mt']) & (o['pct_counts_mt'] < 5)
    o['pass_baseline'] = baseline_mask(o['n_genes_by_counts'], o['pct_counts_mt'])
    o['pass_target_disease'] = o['disease'].isin(['NF', 'HCM'])
    o['pass_cm_annotation'] = o['cell_type_original'].str.contains('cardiomyocyte', case=False, regex=False)
    sub = o['SubCluster'].astype(str).str.strip().str.lower()
    o['pass_author_subcluster'] = o['SubCluster'].notna() & ~sub.isin(['', 'nan', 'none', 'na'])
    log('Auditing cardiomyocyte and non-cardiomyocyte markers')
    marker_info = marker_audit(a, o, reports)
    keep = np.ones(a.n_obs, dtype=bool)
    stages = {'00_author_qc_input': keep.copy()}
    exclusion = np.full(a.n_obs, 'kept', dtype=object)
    for name, col in [('01_target_disease', 'pass_target_disease'), ('02_min_genes', 'pass_min_genes'),
                      ('03_max_genes', 'pass_max_genes'), ('04_raw_mito', 'pass_mito'),
                      ('05_author_cm', 'pass_cm_annotation'), ('06_author_subcluster', 'pass_author_subcluster')]:
        passed = o[col].to_numpy(dtype=bool)
        exclusion[keep & ~passed] = name
        keep &= passed
        stages[name] = keep.copy()
    o['exclusion_primary'] = pd.Categorical(exclusion)
    o['kept_full'] = keep
    o['label'] = o['disease'].astype(str).map({'NF': 0, 'HCM': 1}).fillna(-1).astype('int8')
    if not keep.any() or set(o.loc[keep, 'label']) != {0, 1}:
        raise ValueError('QC did not retain both target diseases')
    log(f'Project filters retain {int(keep.sum())} nuclei')
    grouped_stages(o, stages, ['disease']).to_csv(reports/'filter_counts_by_disease.tsv', sep='\t', index=False)
    grouped_stages(o, stages, ['donor_id', 'disease']).to_csv(reports/'filter_counts_by_donor.tsv', sep='\t', index=False)
    grouped_stages(o, stages, ['sample_id', 'donor_id', 'disease']).to_csv(reports/'filter_counts_by_sample.tsv', sep='\t', index=False)
    grouped_stages(o, stages, ['cell_type_original']).to_csv(reports/'filter_counts_by_celltype.tsv', sep='\t', index=False)
    sensitivity = []
    cm_target = o['pass_target_disease'] & o['pass_cm_annotation'] & o['pass_author_subcluster']
    for low in [150, 300, 500]:
        for high in [5000, 6000, np.inf]:
            candidate = cm_target & baseline_mask(o['n_genes_by_counts'], o['pct_counts_mt'], low, high)
            for disease in ['NF', 'HCM']:
                group = o['disease'] == disease
                denominator = int((cm_target & group).sum())
                sensitivity.append({'min_genes': low, 'max_genes': high, 'mt_lt': 5, 'disease': disease,
                    'retained': int((candidate & group).sum()), 'eligible_author_cm': denominator,
                    'retention_pct': 100*int((candidate & group).sum())/denominator})
    pd.DataFrame(sensitivity).to_csv(reports/'threshold_sensitivity.tsv', sep='\t', index=False)
    author_mito = pd.to_numeric(o.get('cellranger_percent_mito'), errors='coerce')
    differences = {str(factor): float(np.nanmedian(np.abs(o['pct_counts_mt'].to_numpy()-author_mito.to_numpy()*factor)))
                   for factor in [1, 100]}
    comparison = {'author_mito_to_percent_scale_median_abs_errors': differences,
                  'best_supported_scale': min(differences, key=differences.get),
                  'comparison_note': 'Recomputed raw count fraction is used; author field units are audited, never assumed.'}
    if 'cellbender_ngenes' in o:
        comparison['author_cb_ngenes_mismatches'] = int((o['cellbender_ngenes'] != o['n_genes_by_counts']).sum())
    write_json(reports/'SOURCE_COMPARISON.json', comparison)
    target_donors = set(o.loc[o['pass_target_disease'], 'donor_id'].astype(str))
    kept_donors = set(o.loc[keep, 'donor_id'].astype(str))
    donor_splits = make_donor_splits(o.loc[keep], args.seed)
    o['split'] = o['donor_id'].astype(str).map(donor_splits).fillna('excluded')
    o['split_status'] = 'provisional_pending_project_ratio'
    selected = balanced_donor_sample(o.loc[keep], args.light_size, args.seed)
    o['selected_light'] = o.index.isin(selected)
    o.to_csv(reports/'cell_qc_audit.tsv.gz', sep='\t', index_label='barcode', compression='gzip')
    donor_table = o.loc[keep].groupby(['donor_id', 'disease', 'label', 'split'], observed=True).agg(
        nuclei_full=('dataset', 'size'), nuclei_light=('selected_light', 'sum'),
        samples=('sample_id', 'nunique')).reset_index()
    donor_table.to_csv(reports/'donor_split_and_sampling.tsv', sep='\t', index=False)
    qc_plots(o, keep, reports)
    validation = {'count_matrices': {'cellbender': cb_stats, 'cellranger': raw_stats},
                  'mitochondrial_genes': int(a.var['mt'].sum()), 'markers': marker_info,
                  'dropped_target_donors': sorted(target_donors-kept_donors),
                  'kept_donors': len(kept_donors), 'split_ratio_is_provisional': True,
                  'new_doublet_calls': False,
                  'light_marker_mixed_warnings': int(o.loc[selected, 'marker_mixed_warning'].sum())}
    write_json(reports/'PRE_EXPORT_AUDIT.json', validation)
    if target_donors != kept_donors:
        raise ValueError('At least one target donor lost all nuclei; review written QC tables before export')
    log('Creating clean AnnData with full original gene space and count provenance')
    var = a.var.copy()
    var['gene_symbol'] = var.index.astype(str)
    if 'gene_ids' in var:
        var['ensembl_id'] = var['gene_ids'].astype(str)
    obs = o.loc[keep].copy()
    obs.index = 'SCP1303:' + obs.index.astype(str)
    full = ad.AnnData(X=a.X[keep].copy(), obs=obs, var=var)
    full.layers['counts'] = full.X.copy()
    full.layers['cellranger_raw'] = a.layers['cellranger_raw'][keep].copy()
    source_path = str(Path.home()/'xiaoe-internship-data/training/SCP1303/human_dcm_hcm_scportal_03.17.2022.h5ad.gz')
    full.uns['qc_parameters_json'] = json.dumps(params, ensure_ascii=False)
    full.uns['source_files_json'] = json.dumps({'archive': source_path,
        'archive_sha256': '088da9ed54b3d45e039a00b7f3db58e9a032496fb1e45b1ac32019c2280c9cb1',
        'author_readme': 'human_dcm_hcm_scportal_03.17.2022_README.txt'})
    full.uns['count_layer_provenance'] = {'X': 'Unnormalized CellBender-adjusted counts',
        'counts': 'Same CellBender-adjusted counts as X', 'cellranger_raw': 'Original Cell Ranger counts'}
    full.uns['label_map'] = {'NF': 0, 'HCM': 1}
    full.uns['author_qc_reused'] = 'CellBender ambient removal and author sample/cluster Scrublet QC; no new global doublet call.'
    full.uns['split_status'] = params['split_status']
    full.uns['external_validation_status'] = params['cross_study_donor_overlap']
    full.uns['created_utc'] = datetime.now(timezone.utc).isoformat()
    full.uns['software_versions_json'] = json.dumps({k: importlib.metadata.version(k)
        for k in ['scanpy', 'anndata', 'numpy', 'scipy', 'pandas', 'h5py', 'matplotlib']})
    light_ids = ['SCP1303:'+str(i) for i in selected]
    light = full[light_ids].copy()
    light.uns['sampling'] = {'seed': args.seed, 'target': args.light_size,
        'actual': light.n_obs, 'method': 'Balanced disease labels; near-equal distinct-nucleus allocation per donor'}
    del a
    gc.collect()
    for name, data in [('cm_qc_full.h5ad', full), ('cm_qc_light.h5ad', light)]:
        log(f'Writing {name}: {data.shape}')
        target = out/name
        temp = out/(name+'.tmp')
        data.write_h5ad(temp, compression='gzip', compression_opts=4)
        temp.replace(target)
    summary = {'status': 'complete', 'dataset': 'SCP1303', 'input_nuclei': len(o),
        'full_nuclei': full.n_obs, 'light_nuclei': light.n_obs, 'genes': full.n_vars,
        'donors': len(kept_donors), 'full_by_disease': full.obs['disease'].value_counts().to_dict(),
        'light_by_disease': light.obs['disease'].value_counts().to_dict(),
        'parameters': params, 'source_comparison': comparison, 'validation': validation,
        'completed_utc': datetime.now(timezone.utc).isoformat(), 'outputs': []}
    for path in sorted(out.glob('cm_qc_*.h5ad')):
        log('Hashing output '+path.name)
        with path.open('rb') as f:
            md5 = hashlib.file_digest(f, 'md5').hexdigest()
        with path.open('rb') as f:
            sha = hashlib.file_digest(f, 'sha256').hexdigest()
        summary['outputs'].append({'file': path.name, 'bytes': path.stat().st_size, 'md5': md5, 'sha256': sha})
    write_json(reports/'RUN_SUMMARY.json', summary)
    write_json(out/'SUCCESS.json', summary)
    (out/'MD5SUMS').write_text(''.join(r['md5']+'  '+r['file']+'\n' for r in summary['outputs']))
    (out/'SHA256SUMS').write_text(''.join(r['sha256']+'  '+r['file']+'\n' for r in summary['outputs']))
    log(f'QC EXPORT COMPLETE: full={full.n_obs}, light={light.n_obs}, donors={len(kept_donors)}')


if __name__ == '__main__':
    main()
