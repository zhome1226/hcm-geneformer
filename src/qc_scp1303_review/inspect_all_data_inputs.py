"""Read-only schema inventory for subsequent modality-specific cleaning."""
from datetime import datetime, timezone
from pathlib import Path
import gzip
import json
import tarfile

import anndata as ad
import h5py
import numpy as np
import pandas as pd

root = Path.home() / 'xiaoe-internship-data'
out = root / 'processed/all_data_cleaning_20261111'
out.mkdir(exist_ok=True)
manifest = pd.read_csv(out / 'ALL_FILES.tsv', sep='\t')
results = {'created_utc': datetime.now(timezone.utc).isoformat(),
    'scope': 'Schema/sample inspection only; not a completion audit or cleaning result.', 'datasets': {}}
for dataset, group in manifest.groupby('dataset', sort=True):
    if '/' in dataset:
        continue
    records = []
    for _, entry in group.iterrows():
        path = Path(entry['server_path'])
        name = path.name
        if not path.is_file():
            records.append({'path': str(path), 'error': 'missing'})
            continue
        record = {'path': str(path), 'bytes': path.stat().st_size}
        try:
            if name.endswith(('.h5ad', '.h5mu', '.h5')):
                with h5py.File(path, 'r') as handle:
                    record['hdf5_keys'] = list(handle.keys())
                    groups = {'root': handle} if 'obs' in handle and 'var' in handle else {}
                    if 'mod' in handle:
                        groups.update({key: handle['mod'][key] for key in handle['mod']})
                    record['modalities'] = {}
                    for key, node in groups.items():
                        obs = ad.io.read_elem(node['obs'])
                        var = ad.io.read_elem(node['var'])
                        item = {'shape': [len(obs), len(var)], 'obs_columns': list(obs.columns),
                            'var_columns': list(var.columns), 'index_unique': bool(obs.index.is_unique),
                            'gene_index_unique': bool(var.index.is_unique), 'gene_examples': list(var.index[:3]),
                            'layers': list(node['layers'].keys()) if 'layers' in node else [],
                            'has_raw': 'raw' in node,
                            'categories': {column: obs[column].astype(str).value_counts().head(12).to_dict()
                                for column in ['cell_type', 'celltype', 'disease', 'sample', 'donor', 'perturbation', 'batch', 'split', 'feature_types']
                                if column in obs}}
                        if 'X' in node:
                            x = node['X']
                            values = x['data'][:100000] if isinstance(x, h5py.Group) else x[:min(25, x.shape[0]), :min(1000, x.shape[1])].ravel()
                            item['X_sample_only'] = {'encoding': str(x.attrs.get('encoding-type')),
                                'min': float(values.min()) if len(values) else None,
                                'max': float(values.max()) if len(values) else None,
                                'nonfinite': int((~np.isfinite(values)).sum()),
                                'fractional': int((np.abs(values-np.round(values))>1e-5).sum())}
                        record['modalities'][key] = item
                records.append(record)
            elif name.endswith(('.tar', '.tar.gz')):
                with tarfile.open(path, 'r:*') as archive:
                    members = archive.getmembers()
                    record['archive_member_count'] = len(members)
                    record['member_examples'] = [m.name for m in members[:30]]
                records.append(record)
            elif '/raw/' in str(path) and name.endswith(('.mtx.gz', '.txt.gz', '.csv.gz', '.tsv.gz')):
                with gzip.open(path, 'rt') as handle:
                    record['head'] = [handle.readline().rstrip()[:400] for _ in range(4)]
                records.append(record)
            elif name.endswith(('.h5ad.gz', '.Rdata', '.soft')):
                record['note'] = 'Needs decompression/R/clinical-metadata parsing in cleaning stage'
                records.append(record)
        except Exception as exc:
            record['error'] = type(exc).__name__ + ': ' + str(exc)
            records.append(record)
    results['datasets'][dataset] = records
    print(f'Inspected schema: {dataset} ({len(records)} records)', flush=True)
(out / 'INPUT_SCHEMAS.json').write_text(json.dumps(results, indent=2, ensure_ascii=False) + '\n')
print('ALL_INPUT_SCHEMAS_READY', flush=True)
