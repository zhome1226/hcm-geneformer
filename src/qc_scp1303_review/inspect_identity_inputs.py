"""Inspect reference layers and metadata before selecting an annotation reference."""
from pathlib import Path
import json
import anndata as ad
import h5py
import numpy as np
import pandas as pd

root = Path.home() / 'xiaoe-internship-data'
ref = root / 'references/Litvinukova2020/heart_atlas_49312f8c.h5ad'
out = root / 'processed/SCP1303_review_20261111'
report = {}
with h5py.File(ref) as f:
    obs = ad.io.read_elem(f['obs'])
    report['expression'] = {}
    for name in ['X', 'raw/X']:
        node = f[name]
        values = node['data'][:100000] if isinstance(node, h5py.Group) else node[:50, :1000].ravel()
        report['expression'][name] = {'encoding': str(node.attrs.get('encoding-type')),
            'shape': [int(x) for x in node.attrs.get('shape', getattr(node, 'shape', []))],
            'sample_min': float(values.min()), 'sample_max': float(values.max()),
            'sample_fractional_values': int((np.abs(values - np.round(values)) > 1e-5).sum())}
    for name in ['var', 'raw/var']:
        var = ad.io.read_elem(f[name])
        report[name] = {'shape': list(var.shape), 'index_example': list(var.index[:4]),
            'columns': list(var.columns)}
    report['uns_keys'] = list(f['uns'].keys())
    report['extra_categories'] = {k: obs[k].astype(str).value_counts().to_dict()
                                for k in ['source', 'cell_source', 'type', 'NRP']}
    lv = obs.loc[obs['suspension_type'].eq('nucleus') & obs['tissue'].eq('heart left ventricle')]
    report['lv_by_donor_type'] = lv.groupby(['donor_id', 'cell_type'], observed=True).size().reset_index(name='n').to_dict('records')
query = pd.read_csv(out / 'query_metadata.tsv', sep='\t', nrows=2)
report['query_fields'] = list(query.columns)
report['query_examples'] = query.fillna('').to_dict('records')
(out / 'IDENTITY_INPUT_INSPECTION.json').write_text(json.dumps(report, indent=2) + '\n')
print(json.dumps({k: report[k] for k in ['expression', 'var', 'raw/var', 'extra_categories', 'query_fields']}, indent=2))
