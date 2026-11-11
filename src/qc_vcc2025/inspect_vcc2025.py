"""Inspect official split metadata without materializing the large count matrices."""
from datetime import datetime, timezone
from pathlib import Path
import json
import anndata as ad
import h5py
import numpy as np
import pandas as pd

ROOT = Path.home() / 'xiaoe-internship-data'
SOURCE = ROOT / 'benchmarks/competitions/vcc2025'
WORK = ROOT / '_processing/all_data_cleaning_20261111/VCC2025'
FILES = {'train': ('train', 'Training'), 'validation': ('validation', 'Validation'), 'test': ('test', 'Test')}


def main():
    WORK.mkdir(parents=True, exist_ok=True)
    gene_names = pd.read_csv(SOURCE / 'gene_names.csv', header=None)[0].astype(str).tolist()
    report = {'checked_utc': datetime.now(timezone.utc).isoformat(), 'gene_names': len(gene_names), 'splits': {}, 'overlap': {}}
    observations = {}
    for split, (directory, suffix) in FILES.items():
        path = SOURCE / directory / ('adata_' + suffix + '.h5ad')
        with h5py.File(path) as handle:
            obs = ad.io.read_elem(handle['obs'])
            var = ad.io.read_elem(handle['var'])
            observations[split] = obs
            values = handle['X']['data'][:100000]
            counts = pd.read_csv(SOURCE / directory / ('pert_counts_' + suffix + '.csv'))
            actual = obs.groupby('target_gene', observed=True).size()
            report['splits'][split] = {'cells': len(obs), 'genes': len(var), 'obs_fields': list(obs.columns),
                'var_fields': list(var.columns), 'X_dtype': str(handle['X']['data'].dtype),
                'indices_dtype': str(handle['X']['indices'].dtype), 'nnz': len(handle['X']['data']),
                'encoding': str(handle['X'].attrs['encoding-type']), 'has_raw': 'raw' in handle,
                'layers': list(handle['layers']), 'gene_axis_matches_official': var.index.tolist() == gene_names,
                'cell_id_examples': obs.index[:3].tolist(), 'target_genes': int(obs.target_gene.nunique()),
                'guides': int(obs.guide_id.nunique()), 'batches': int(obs.batch.nunique()),
                'control_cells': int(obs.target_gene.eq('non-targeting').sum()),
                'obs_missing': obs.isna().sum().to_dict(),
                'guide_target_conflicts': int((obs.groupby('guide_id', observed=True).target_gene.nunique() > 1).sum()),
                'pert_counts_columns': list(counts.columns), 'pert_counts_rows': len(counts),
                'count_file_labels_match_noncontrol': set(counts.target_gene) == set(actual.index) - {'non-targeting'},
                'count_file_sizes_match': all(actual.loc[row.target_gene] == row.n_cells for row in counts.itertuples()),
                'sample_only_numeric': {'min': float(values.min()), 'max': float(values.max()),
                    'finite': bool(np.isfinite(values).all()), 'integer': bool((values == np.rint(values)).all())}}
    for left, right in [('train', 'validation'), ('train', 'test'), ('validation', 'test')]:
        a, b = observations[left], observations[right]
        shared = a.index.intersection(b.index)
        report['overlap'][left + '__' + right] = {'shared_cells': len(shared),
            'all_shared_non_targeting': bool(a.loc[shared, 'target_gene'].eq('non-targeting').all() and b.loc[shared, 'target_gene'].eq('non-targeting').all()),
            'shared_metadata_exact': bool(a.loc[shared].astype(str).equals(b.loc[shared].astype(str))),
            'shared_noncontrol_targets': sorted((set(a.target_gene) & set(b.target_gene)) - {'non-targeting'})}
    report['unique_cell_identifiers'] = len(set().union(*(set(obs.index) for obs in observations.values())))
    (WORK / 'INPUT_INSPECTION.json').write_text(json.dumps(report, indent=2) + '\n')
    print(json.dumps(report, indent=2), flush=True)


if __name__ == '__main__':
    main()
