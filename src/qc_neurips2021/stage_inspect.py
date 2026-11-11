"""Verify the two originals, expand gzip streams, and inspect their native assays."""
from datetime import datetime, timezone
from pathlib import Path
import gzip
import hashlib
import json
import os
import urllib.request
import anndata as ad
import h5py
import numpy as np
import pandas as pd

ROOT = Path.home() / 'xiaoe-internship-data'
WORK = ROOT / '_processing/all_data_cleaning_20261111/NeurIPS2021'
EXPANDED = WORK / 'expanded_inputs'


def digest(path, algorithm):
    with path.open('rb') as handle:
        return hashlib.file_digest(handle, algorithm).hexdigest()


def json_default(value):
    if isinstance(value, np.generic):
        return value.item()
    if isinstance(value, np.ndarray):
        return value.tolist()
    return str(value)


def matrix_info(node):
    encoding = node.attrs.get('encoding-type', 'array')
    if isinstance(node, h5py.Group) and encoding == 'dataframe':
        table = ad.io.read_elem(node)
        values = table.to_numpy()
        return {'shape': list(table.shape), 'encoding': encoding, 'columns': table.columns.tolist(),
            'dtype': str(values.dtype), 'sample_min': float(np.nanmin(values)), 'sample_max': float(np.nanmax(values)),
            'nonfinite': int(np.count_nonzero(~np.isfinite(values))),
            'fractional': int(np.count_nonzero(values != np.rint(values)))}
    if isinstance(node, h5py.Group):
        values = node['data'][:100000]
        shape = list(node.attrs['shape'])
        entries = len(node['data'])
    else:
        shape = list(node.shape)
        values = node[:min(len(node), 16)].ravel()
        entries = int(np.prod(shape))
    return {'shape': shape, 'encoding': encoding, 'dtype': str(values.dtype), 'stored_entries': entries,
        'sample_min': float(values.min()) if len(values) else None, 'sample_max': float(values.max()) if len(values) else None,
        'sample_nonfinite': int(np.count_nonzero(~np.isfinite(values))),
        'sample_fractional': int(np.count_nonzero(values != np.rint(values))), 'sample_values': values[:8].tolist()}


def inspect(path):
    with h5py.File(path, 'r') as handle:
        obs, var = ad.io.read_elem(handle['obs']), ad.io.read_elem(handle['var'])
        result = {'hdf5_keys': list(handle), 'obs_shape': list(obs.shape), 'var_shape': list(var.shape),
            'obs_columns': list(obs), 'var_columns': list(var),
            'obs_index_unique': bool(obs.index.is_unique), 'var_index_unique': bool(var.index.is_unique),
            'obs_missing': obs.isna().sum().to_dict(), 'var_missing': var.isna().sum().to_dict(),
            'obs_examples': obs.head(2).reset_index().to_dict('records'),
            'var_examples': var.head(3).reset_index().to_dict('records'), 'categories': {}, 'matrices': {}, 'uns': {}}
        for name in obs:
            if obs[name].nunique(dropna=False) <= 50:
                result['categories']['obs:' + name] = obs[name].astype(str).value_counts(dropna=False).to_dict()
        for name in var:
            if var[name].nunique(dropna=False) <= 50:
                result['categories']['var:' + name] = var[name].astype(str).value_counts(dropna=False).to_dict()
        result['matrices']['X'] = matrix_info(handle['X'])
        for name in handle.get('layers', {}):
            result['matrices']['layers/' + name] = matrix_info(handle['layers'][name])
        if 'raw' in handle:
            result['matrices']['raw/X'] = matrix_info(handle['raw/X'])
            rawvar = ad.io.read_elem(handle['raw/var'])
            result['raw_var_equal_native'] = rawvar.index.equals(var.index)
            result['raw_var_shape'] = list(rawvar.shape)
        for name in handle.get('obsm', {}):
            result['matrices']['obsm/' + name] = matrix_info(handle['obsm'][name])
        for name in handle.get('uns', {}):
            try:
                value = ad.io.read_elem(handle['uns'][name])
                result['uns'][name] = json.dumps(value, default=json_default)[:6000]
            except Exception as error:
                result['uns'][name] = {'read_error': repr(error)}
        result['varm_shapes'] = {name: list(node.shape) if isinstance(node, h5py.Dataset) else str(node.attrs.get('shape')) for name, node in handle.get('varm', {}).items()}
        result['obsp_shapes'] = {name: list(node.attrs.get('shape', [])) if isinstance(node, h5py.Group) else list(node.shape) for name, node in handle.get('obsp', {}).items()}
        return result


def main():
    EXPANDED.mkdir(parents=True, exist_ok=True)
    assert not (WORK / 'STAGE_SUCCESS.json').exists(), 'Successful stage exists; do not repeat expansion'
    baseline = json.loads((WORK / 'ORIGINAL_DOWNLOAD_BASELINE.json').read_text())
    assert len(baseline) == 2
    staged, inspections = [], {}
    for item in baseline:
        source = Path(item['server_path'])
        md5, sha256 = digest(source, 'md5'), digest(source, 'sha256')
        assert source.stat().st_size == int(item['bytes']) and md5 == item['md5'] and sha256 == item['sha256']
        print('Original checksum verified:', source.name, flush=True)
        output = EXPANDED / source.stem
        partial = output.with_suffix(output.suffix + '.partial')
        assert not output.exists() and not partial.exists(), 'Inspect existing expansion before resuming'
        expanded_sha, expanded_md5, size = hashlib.sha256(), hashlib.md5(), 0
        with gzip.open(source, 'rb') as compressed, partial.open('xb') as expanded:
            while block := compressed.read(16 * 1024 * 1024):
                expanded.write(block)
                expanded_sha.update(block)
                expanded_md5.update(block)
                size += len(block)
            expanded.flush()
            os.fsync(expanded.fileno())
        partial.rename(output)
        record = {'original_file': source.name, 'original_path': str(source),
            'original_bytes': source.stat().st_size, 'original_MD5': md5, 'original_SHA256': sha256,
            'expanded_path': str(output), 'expanded_bytes': size, 'expanded_MD5': expanded_md5.hexdigest(),
            'expanded_SHA256': expanded_sha.hexdigest(), 'gzip_CRC_verified_at_EOF': True, 'source_url': item['source_url']}
        staged.append(record)
        (WORK / 'STAGED_INPUTS.json').write_text(json.dumps(staged, indent=2) + '\n')
        inspections[source.name] = inspect(output)
        print('Expanded and inspected:', source.name, size, flush=True)
    evidence = {'checked_utc': datetime.now(timezone.utc).isoformat(), 'files': inspections,
        'scope': 'Input inspection only, not a cleaning completion gate.'}
    (WORK / 'INPUT_INSPECTION.json').write_text(json.dumps(evidence, indent=2, default=json_default) + '\n')
    url = 'https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE194122&targ=self&form=text&view=full'
    try:
        with urllib.request.urlopen(url, timeout=30) as response:
            payload = response.read()
        (WORK / 'GEO_SOURCE.txt').write_bytes(payload)
        lookup = {'url': url, 'status': 'retrieved', 'bytes': len(payload)}
    except Exception as error:
        lookup = {'url': url, 'status': 'unavailable', 'error': repr(error)}
    (WORK / 'SOURCE_LOOKUP.json').write_text(json.dumps(lookup, indent=2) + '\n')
    (WORK / 'STAGE_SUCCESS.json').write_text(json.dumps({'status': 'staged_and_inspected_only', 'staged_files': 2,
        'completed_utc': datetime.now(timezone.utc).isoformat(), 'not_cleaning_completion': True}, indent=2) + '\n')
    print('All originals expanded and inspected; cleaning not yet complete.', flush=True)


if __name__ == '__main__':
    main()
