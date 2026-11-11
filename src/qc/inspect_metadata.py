"""Inspect author metadata and count layout without changing the downloaded data."""
import json
from pathlib import Path

import anndata as ad
import h5py
import numpy as np

work = Path(__file__).resolve().parent
path = work / 'author_input.h5ad'
a = ad.read_h5ad(path, backed='r')
o = a.obs
summary = {
    'shape': list(a.shape), 'obs_columns': list(o.columns), 'var_columns': list(a.var.columns),
    'layers': list(a.layers), 'cell_types': o['cell_type_leiden0.6'].value_counts(dropna=False).to_dict(),
    'disease': o['disease'].value_counts(dropna=False).to_dict(),
    'donors': int(o['donor_id'].nunique()), 'samples': int(o['biosample_id'].nunique()),
    'duplicate_barcodes': int(o.index.duplicated().sum()),
    'duplicate_gene_symbols': int(a.var.index.duplicated().sum()),
    'donor_disease_conflicts': int((o.groupby('donor_id', observed=True)['disease'].nunique() != 1).sum()),
    'numeric_ranges': {}, 'var_preview': a.var.head(6).astype(str).to_dict(),
}
for c in o.select_dtypes(include='number'):
    v = o[c]
    summary['numeric_ranges'][c] = {'missing': int(v.isna().sum()),
        'quantiles': {str(k): float(vv) for k, vv in v.quantile([0, .01, .25, .5, .75, .99, 1]).items()}}
for group in ['disease', 'donor_id', 'biosample_id', 'cell_type_leiden0.6']:
    o.groupby(group, observed=True).agg(nuclei=('disease', 'size'),
        donors=('donor_id', 'nunique')).to_csv(work / ('input_by_' + group.replace('.', '_') + '.tsv'), sep='\t')
o.groupby(['donor_id', 'disease', 'cell_type_leiden0.6'], observed=True).size().rename('nuclei').to_csv(work / 'input_donor_disease_celltype.tsv', sep='\t')
o.groupby(['biosample_id', 'donor_id', 'disease'], observed=True).size().rename('nuclei').to_csv(work / 'input_samples.tsv', sep='\t')
with h5py.File(path, 'r') as f:
    summary['matrix_layout'] = {k: {'attrs': {ak: str(av) for ak,av in f[k].attrs.items()},
        'children': {ck: list(cv.shape) for ck,cv in f[k].items() if isinstance(cv,h5py.Dataset)}}
        for k in ['X', 'layers/cellranger_raw']}
(work / 'INPUT_INSPECTION.json').write_text(json.dumps(summary, indent=2, default=str) + '\n')
print(json.dumps(summary, indent=2, default=str))
a.file.close()
