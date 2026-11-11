"""Verify the downloaded archive and unpack a separate working copy, atomically."""
import csv
import gzip
import hashlib
import os
from pathlib import Path
import shutil

work = Path(__file__).resolve().parent
source_dir = Path.home() / 'xiaoe-internship-data/training/SCP1303'
source = source_dir / 'human_dcm_hcm_scportal_03.17.2022.h5ad.gz'
record = next(r for r in csv.DictReader((source_dir / 'MANIFEST.tsv').open(), delimiter='\t')
              if r['file'] == source.name)
if source.stat().st_size != int(record['bytes']):
    raise ValueError('Input archive size differs from verified manifest')
sha = hashlib.file_digest(source.open('rb'), 'sha256').hexdigest()
if sha != record['sha256']:
    raise ValueError('Input archive SHA256 differs from verified official manifest')
print('Input archive SHA256 verified:', sha, flush=True)
target = work / 'author_input.h5ad'
if target.exists():
    print('Existing working H5AD retained:', target, flush=True)
else:
    if shutil.disk_usage(work).free < 30 * 1024**3:
        raise RuntimeError('Need at least 30 GiB free before unpacking')
    temp = target.with_suffix('.h5ad.tmp')
    with gzip.open(source, 'rb') as src, temp.open('wb') as dst:
        shutil.copyfileobj(src, dst, length=16 * 1024**2)
        dst.flush()
        os.fsync(dst.fileno())
    temp.replace(target)
    print('Unpacked working H5AD:', target.stat().st_size, 'bytes', flush=True)
