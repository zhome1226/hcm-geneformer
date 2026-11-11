"""Correct inherited count-source wording on analysis-only copies; never edit delivered counts."""
from datetime import datetime, timezone
from pathlib import Path
import hashlib
import json

import anndata as ad
import h5py

ROOT = Path.home() / 'xiaoe-internship-data/processed'


def matrix_digest(group):
    digest = hashlib.sha256()
    for key in ['data', 'indices', 'indptr']:
        array = group[key]
        for start in range(0, len(array), 1024**2):
            digest.update(array[start:start + 1024**2].tobytes())
    return digest.hexdigest()


for dataset in ['chin_GSE161921', 'chin_GSE174691']:
    output = ROOT / f'{dataset}_qc_20261111'
    path = output / 'diagnostic_normalized.h5ad'
    if (output / 'DIAGNOSTIC_METADATA_FIX.json').exists():
        continue
    exported = json.loads((output / 'EXPORT_SUCCESS.json').read_text())
    original_record = dict(exported['hashes'][path.name])
    with path.open('rb') as handle:
        assert hashlib.file_digest(handle, 'sha256').hexdigest() == original_record['sha256']
    with h5py.File(path, 'r+') as handle:
        before = matrix_digest(handle['X'])
        ad.io.write_elem(handle['uns'], 'count_source',
                         'Diagnostic log1p CP10k expression; original counts are in separately delivered count-preserving objects, not this X')
        ad.io.write_elem(handle['uns'], 'analysis_only', True)
        handle.flush()
        after = matrix_digest(handle['X'])
        assert before == after
    with path.open('rb') as handle:
        new_sha = hashlib.file_digest(handle, 'sha256').hexdigest()
    (output / 'EXPORT_RECORD_BEFORE_DIAGNOSTIC_LABEL.json').write_text(json.dumps(exported, indent=2) + '\n')
    exported['hashes'][path.name] = {'bytes': path.stat().st_size, 'sha256': new_sha}
    exported['diagnostic_metadata_note'] = 'Inherited raw count-source wording corrected; X, indices and indptr digests unchanged. No delivered original-count object edited.'
    (output / 'EXPORT_SUCCESS.json').write_text(json.dumps(exported, indent=2) + '\n')
    (output / 'DIAGNOSTIC_METADATA_FIX.json').write_text(json.dumps({
        'utc': datetime.now(timezone.utc).isoformat(), 'old_file_record': original_record,
        'new_file_record': exported['hashes'][path.name], 'matrix_sha256_before': before,
        'matrix_sha256_after': after, 'matrix_unchanged': True,
        'change': 'Only uns/count_source and uns/analysis_only on the explicitly normalized diagnostic copy'}, indent=2) + '\n')
    print('Diagnostic metadata labeled, matrix unchanged: ' + dataset, flush=True)
