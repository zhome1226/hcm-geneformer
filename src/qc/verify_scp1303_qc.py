"""Independently reopen SCP1303 exports and check membership, counts and hashes."""
from __future__ import annotations

import argparse
from datetime import datetime, timezone
import hashlib
import json
from pathlib import Path

import anndata as ad
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


def require(condition, message):
    if not condition:
        raise AssertionError(message)


def same_sparse(left, right, message):
    require(left.shape == right.shape, message + ': shape')
    require((left != right).nnz == 0, message + ': values')


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--input', required=True, type=Path)
    parser.add_argument('--output', required=True, type=Path)
    args = parser.parse_args()
    out = args.output
    summary = json.loads((out / 'SUCCESS.json').read_text())
    print('Reopening source and both exported H5AD files', flush=True)
    source = ad.read_h5ad(args.input)
    full = ad.read_h5ad(out / 'cm_qc_full.h5ad')
    light = ad.read_h5ad(out / 'cm_qc_light.h5ad')
    require(full.n_obs == summary['full_nuclei'], 'Full shape differs from report')
    require(light.n_obs == summary['light_nuclei'], 'Light shape differs from report')
    require(full.var_names.equals(source.var_names), 'Original gene order/space differs')
    require(light.var_names.equals(full.var_names), 'Light gene order/space differs')
    require(source.obs_names.is_unique, 'Duplicate source cell IDs')
    require(full.obs_names.is_unique and light.obs_names.is_unique, 'Duplicate output cell IDs')

    # Independently reconstruct eligible source membership, including all target CM.
    o = source.obs
    cm = o['cell_type_leiden0.6'].astype(str).str.contains('cardiomyocyte', case=False, regex=False)
    sub = o['SubCluster'].astype(str).str.strip().str.lower()
    candidates = np.flatnonzero((o['disease'].isin(['HCM', 'NF']) & cm &
        o['SubCluster'].notna() & ~sub.isin(['', 'nan', 'none', 'na'])).to_numpy())
    mt = source.var_names.str.startswith('MT-')
    retained = []
    candidate_metrics = []
    for start in range(0, len(candidates), 4096):
        rows = candidates[start:start + 4096]
        cb = source.X[rows]
        raw = source.layers['cellranger_raw'][rows]
        genes = np.asarray((cb > 0).sum(axis=1)).ravel()
        totals = np.asarray(raw.sum(axis=1)).ravel()
        pct = 100 * np.asarray(raw[:, mt].sum(axis=1)).ravel() / np.maximum(totals, 1)
        selected = (genes >= 500) & (genes <= 5000) & (pct < 5) & (totals > 0)
        retained.extend(rows[selected].tolist())
        candidate_metrics.append(pd.DataFrame({'donor_id': o.iloc[rows]['donor_id'].astype(str).to_numpy(),
            'disease': o.iloc[rows]['disease'].astype(str).to_numpy(), 'nuclei_before': 1,
            'below_500': genes < 500, 'above_5000': genes > 5000,
            'mt_at_least_5': pct >= 5, 'nuclei_after': selected}))
    expected = 'SCP1303:' + source.obs_names[np.sort(retained)].astype(str)
    require(full.obs_names.equals(expected), 'Export is not exactly the independently reconstructed QC set')
    retention = pd.concat(candidate_metrics, ignore_index=True).groupby(['donor_id', 'disease'], observed=True).sum()
    retention['retention_pct'] = 100 * retention['nuclei_after'] / retention['nuclei_before']
    retention.to_csv(out / 'reports' / 'cm_retention_by_donor.tsv', sep='\t')
    # Draw positions and values explicitly from the independently verified table.
    positions = np.arange(len(retention))
    fig, ax = plt.subplots(figsize=(12, 5))
    ax.bar(positions, retention['retention_pct'].to_numpy(), color='#247BA0')
    ax.set_xticks(positions, retention.index.get_level_values('donor_id'), rotation=70)
    ax.set_ylabel('Retained cardiomyocyte nuclei (%)')
    ax.set_ylim(0, 100)
    ax.set_title('Independently verified donor retention')
    fig.tight_layout()
    fig.savefig(out / 'reports' / 'donor_retention.png', dpi=160)
    plt.close(fig)
    warnings = full.obs.groupby(['donor_id', 'disease'], observed=True).agg(
        nuclei=('label', 'size'), marker_mixed_warnings=('marker_mixed_warning', 'sum'),
        no_panel_cm_marker=('n_cm_markers_detected', lambda x: (x == 0).sum()))
    warnings.to_csv(out / 'reports' / 'marker_warnings_by_donor.tsv', sep='\t')

    for name, obj in [('full', full), ('light', light)]:
        require(set(obj.obs['disease'].astype(str)) == {'HCM', 'NF'}, name + ': diseases')
        expected_labels = obj.obs['disease'].astype(str).map({'NF': 0, 'HCM': 1}).to_numpy()
        require(np.array_equal(obj.obs['label'].to_numpy(), expected_labels), name + ': labels')
        require((obj.obs.groupby('donor_id', observed=True)['split'].nunique() == 1).all(), name + ': donor leakage')
        require(set(obj.obs['split'].astype(str)) == {'train', 'validation', 'test'}, name + ': split coverage')
        require(obj.obs['barcode_original'].is_unique, name + ': duplicate original barcodes')
        require(obj.raw is None and len(obj.obsm) == 0 and len(obj.obsp) == 0, name + ': stale analysis slots')
        genes = np.asarray((obj.X > 0).sum(axis=1)).ravel()
        require(((genes >= 500) & (genes <= 5000)).all(), name + ': gene thresholds')
        require(np.array_equal(genes, obj.obs['n_genes_by_counts']), name + ': gene metric mismatch')
        raw = obj.layers['cellranger_raw']
        totals = np.asarray(raw.sum(axis=1)).ravel()
        pct = 100 * np.asarray(raw[:, mt].sum(axis=1)).ravel() / np.maximum(totals, 1)
        require(((totals > 0) & (pct < 5)).all(), name + ': mitochondrial threshold')
        require(np.allclose(pct, obj.obs['pct_counts_mt'], rtol=1e-6, atol=1e-6), name + ': mt metric mismatch')
        same_sparse(obj.X, obj.layers['counts'], name + ': X/counts equality')
        for matrix in (obj.X, obj.layers['counts'], raw):
            require(np.isfinite(matrix.data).all() and (matrix.data >= 0).all(), name + ': illegal count values')

    require(light.n_obs == 20000, 'Light object must contain 20,000 distinct nuclei')
    require(light.obs['label'].value_counts().to_dict() == {0: 10000, 1: 10000}, 'Light disease balance')
    require(set(light.obs['donor_id']) == set(full.obs['donor_id']), 'Light donor coverage')
    light_rows = full.obs_names.get_indexer(light.obs_names)
    require((light_rows >= 0).all(), 'Light contains cells absent from full object')
    require(full.obs.iloc[light_rows]['split'].astype(str).tolist() == light.obs['split'].astype(str).tolist(), 'Light/full split mismatch')
    for layer in ['X', 'counts', 'cellranger_raw']:
        left = full.X if layer == 'X' else full.layers[layer]
        right = light.X if layer == 'X' else light.layers[layer]
        same_sparse(left[light_rows], right, 'All light/full counts: ' + layer)

    rng = np.random.default_rng(20261111)
    sampled = np.sort(rng.choice(full.n_obs, size=min(128, full.n_obs), replace=False))
    source_rows = source.obs_names.get_indexer(full.obs.iloc[sampled]['barcode_original'].astype(str))
    require((source_rows >= 0).all(), 'Output barcodes missing in source')
    same_sparse(full.X[sampled], source.X[source_rows], '128 output/source rows, CellBender all genes')
    same_sparse(full.layers['cellranger_raw'][sampled], source.layers['cellranger_raw'][source_rows],
                '128 output/source rows, Cell Ranger all genes')

    print('Membership, thresholds, donor splits and count consistency passed; checking files', flush=True)
    checked_files = []
    for item in summary['outputs']:
        path = out / item['file']
        require(path.stat().st_size == item['bytes'], path.name + ': byte size')
        for algorithm in ['md5', 'sha256']:
            with path.open('rb') as handle:
                require(hashlib.file_digest(handle, algorithm).hexdigest() == item[algorithm], path.name + ': ' + algorithm)
        checked_files.append(item)
    scripts = {}
    for path in sorted(Path(__file__).parent.glob('*')):
        if path.suffix in {'.py', '.sh'} or path.name == 'environment.freeze.txt':
            scripts[path.name] = hashlib.sha256(path.read_bytes()).hexdigest()
    result = {'status': 'passed', 'verified_utc': datetime.now(timezone.utc).isoformat(),
        'independent_source_membership': 'exact match across all target cardiomyocytes',
        'full_shape': list(full.shape), 'light_shape': list(light.shape),
        'full_by_disease': full.obs['disease'].value_counts().to_dict(),
        'light_by_disease': light.obs['disease'].value_counts().to_dict(),
        'donors': int(full.obs['donor_id'].nunique()),
        'full_marker_mixed_warnings': int(full.obs['marker_mixed_warning'].sum()),
        'light_marker_mixed_warnings': int(light.obs['marker_mixed_warning'].sum()),
        'full_no_panel_cm_marker': int((full.obs['n_cm_markers_detected'] == 0).sum()),
        'donor_retention_pct_range': [float(retention['retention_pct'].min()), float(retention['retention_pct'].max())],
        'source_count_comparison': '128 seeded rows, both matrices, every gene; exact',
        'light_count_comparison': 'all 20,000 rows, all three matrices; exact',
        'thresholds_recomputed_from_exported_counts': True, 'donor_leakage': False,
        'checked_files': checked_files, 'processing_script_sha256': scripts}
    (out / 'reports' / 'INDEPENDENT_VERIFICATION.json').write_text(json.dumps(result, indent=2) + '\n')
    print(json.dumps(result, indent=2), flush=True)


if __name__ == '__main__':
    main()
