"""Inspect all downloaded releases and modalities without modifying the inputs."""
from datetime import datetime, timezone
from pathlib import Path
import json
import anndata as ad
import h5py
import numpy as np
import pandas as pd

ROOT = Path.home() / 'xiaoe-internship-data'
SOURCE = ROOT / 'benchmarks/competitions/neurips2023'
WORK = ROOT / '_processing/all_data_cleaning_20261111/NeurIPS2023'


def native(value):
    if isinstance(value, np.generic):
        return value.item()
    if isinstance(value, np.ndarray):
        return value.tolist()
    return str(value)


def matrix_info(matrix):
    sparse = isinstance(matrix, h5py.Group)
    data = matrix['data'][:100000] if sparse else matrix[:min(25, matrix.shape[0])].ravel()
    return {'shape': list(matrix.attrs['shape'] if sparse else matrix.shape),
        'encoding': matrix.attrs.get('encoding-type', 'array'), 'dtype': str(data.dtype),
        'nnz_or_entries': len(matrix['data']) if sparse else int(np.prod(matrix.shape)),
        'sample_min': float(data.min()) if len(data) else None,
        'sample_max': float(data.max()) if len(data) else None,
        'sample_nonfinite': int((~np.isfinite(data)).sum()),
        'sample_fractional': int(np.count_nonzero(data != np.rint(data))),
        'sample_examples': data[:12].tolist()}


def inspect_group(group):
    obs = ad.io.read_elem(group['obs'])
    var = ad.io.read_elem(group['var'])
    result = {'keys': list(group), 'obs_shape': list(obs.shape), 'var_shape': list(var.shape),
        'obs_index_unique': obs.index.is_unique, 'var_index_unique': var.index.is_unique,
        'obs_examples': obs.head(2).reset_index().to_dict('records'),
        'var_examples': var.head(2).reset_index().to_dict('records'),
        'obs_missing': obs.isna().sum().to_dict(), 'var_missing': var.isna().sum().to_dict(),
        'categories': {}, 'matrices': {}, 'uns': {}, 'obsm_shapes': {}, 'obsm_tables': {}}
    for name in ['cell_type', 'cell_type_orig', 'split', 'donor_id', 'control', 'sm_name', 'group', 'library', 'feature_types']:
        if name in obs:
            result['categories'][name] = obs[name].value_counts(dropna=False).to_dict()
        if name in var:
            result['categories']['var_' + name] = var[name].value_counts(dropna=False).to_dict()
    if 'X' in group:
        result['matrices']['X'] = matrix_info(group['X'])
    if 'layers' in group:
        for name in group['layers']:
            result['matrices']['layers/' + name] = matrix_info(group['layers'][name])
    if 'raw' in group and 'X' in group['raw']:
        rawvar = ad.io.read_elem(group['raw']['var'])
        result['raw_var_shape'] = list(rawvar.shape)
        result['raw_var_examples'] = rawvar.head(2).reset_index().to_dict('records')
        result['raw_var_equal_main'] = rawvar.index.equals(var.index)
        result['matrices']['raw/X'] = matrix_info(group['raw']['X'])
    if 'uns' in group:
        for name in group['uns']:
            try:
                value = ad.io.read_elem(group['uns'][name])
                text = json.dumps(value, default=native)
                result['uns'][name] = text[:1800]
            except Exception as error:
                result['uns'][name] = repr(error)
    if 'obsm' in group:
        for name, value in group['obsm'].items():
            result['obsm_shapes'][name] = list(value.shape) if isinstance(value, h5py.Dataset) else str(value.attrs.get('shape'))
            if isinstance(value, h5py.Group) and value.attrs.get('encoding-type') == 'dataframe':
                table = ad.io.read_elem(value)
                array = table.to_numpy()
                result['obsm_tables'][name] = {'shape': list(table.shape), 'columns': table.columns.tolist(),
                    'index_aligns_obs': table.index.equals(obs.index), 'missing': int(table.isna().sum().sum()),
                    'min': float(np.nanmin(array)), 'max': float(np.nanmax(array)),
                    'fractional': int(np.count_nonzero(array != np.rint(array))),
                    'nonfinite': int(np.count_nonzero(~np.isfinite(array)))}
    return result, obs, var


def main():
    WORK.mkdir(parents=True, exist_ok=True)
    report = {'checked_utc': datetime.now(timezone.utc).isoformat(), 'files': {}, 'relationships': {}}
    observations, features = {}, {}
    baseline = json.loads((WORK / 'ORIGINAL_DOWNLOAD_BASELINE.json').read_text())
    for item in baseline:
        path = Path(item['server_path'])
        if path.suffix not in {'.h5ad', '.h5mu'}:
            continue
        with h5py.File(path, 'r') as handle:
            content = {}
            for name, group in [('root', handle)] + ([('mod/' + key, value) for key, value in handle['mod'].items()] if 'mod' in handle else []):
                info, obs, var = inspect_group(group)
                content[name] = info
                key = path.name + ':' + name
                observations[key], features[key] = obs, var
            report['files'][path.name] = content
    original = observations['sc_counts.h5ad:root']
    updated = observations['sc_counts_reannotated_with_counts.h5ad:root']
    shared = original.index.intersection(updated.index)
    report['relationships']['sc_release_versions'] = {
        'shared_cells': len(shared), 'original_only': len(original.index.difference(updated.index)),
        'updated_only': len(updated.index.difference(original.index)),
        'gene_axis_equal': features['sc_counts.h5ad:root'].index.equals(features['sc_counts_reannotated_with_counts.h5ad:root'].index),
        'shared_metadata_differences': {name: int((original.loc[shared, name].astype(str) != updated.loc[shared, name].astype(str)).sum()) for name in set(original.columns) & set(updated.columns)},
    }
    for name in ['root', 'mod/rna', 'mod/atac']:
        obs = observations['multiome_counts.h5mu:' + name]
        report['relationships']['multiome_' + name + '_aligns_root'] = obs.index.equals(observations['multiome_counts.h5mu:root'].index)
    report['auxiliary'] = {}
    for relative in ['moa_annotations.csv', 'workflow-resources/neurips-2023-kaggle/id_map.csv']:
        table = pd.read_csv(SOURCE / relative)
        report['auxiliary'][relative] = {'shape': list(table.shape), 'columns': list(table.columns),
            'missing': table.isna().sum().to_dict(), 'head': table.head(3).to_dict('records'),
            'exact_duplicate_rows': int(table.duplicated().sum())}
    report['dataset_info_yaml'] = (SOURCE / 'workflow-resources/neurips-2023-kaggle/dataset_info.yaml').read_text()
    (WORK / 'INPUT_INSPECTION.json').write_text(json.dumps(report, indent=2, default=native) + '\n')
    print('Inspection completed; full results: ' + str(WORK / 'INPUT_INSPECTION.json'), flush=True)


if __name__ == '__main__':
    main()
