"""Download a catalogued reference with resumable verified HTTP ranges."""
from concurrent.futures import ThreadPoolExecutor, as_completed
from datetime import datetime, timezone
import argparse
import hashlib
import json
from pathlib import Path
import re
import shutil
import subprocess

p = argparse.ArgumentParser(description=__doc__)
p.add_argument('--proxy', required=True)
p.add_argument('--workers', type=int, default=8)
args = p.parse_args()
root = Path.home()/'xiaoe-internship-data/references/Litvinukova2020'
catalog = json.loads((root/'collection.json').read_text())
dataset = next(x for x in catalog['datasets'] if x['dataset_id'] == 'd4e69e01-3ba2-4d6b-a15d-e7048f78f22e')
asset = next(x for x in dataset['assets'] if x['filetype'] == 'H5AD')
expected = int(asset['filesize'])
target = root/'heart_atlas_49312f8c.h5ad'
parts = root/'range_parts'
parts.mkdir(exist_ok=True)
chunk = 32*1024**2
if shutil.disk_usage(root).free < expected*2 + 10*1024**3:
    raise RuntimeError('Insufficient space for range assembly')
prefix_size = target.stat().st_size if target.exists() else 0
if prefix_size > expected:
    raise ValueError('Oversized prefix')
# Retain bytes from the previously interrupted sequential transfer.
if prefix_size:
    with target.open('rb') as source:
        for number, start in enumerate(range(0, expected, chunk)):
            available = min(chunk, max(0, prefix_size-start))
            if not available:
                break
            part = parts/f'{number:04}.part'
            if not part.exists():
                source.seek(start)
                with part.open('wb') as handle:
                    handle.write(source.read(available))

def fetch(number, start):
    length = min(chunk, expected-start)
    part = parts/f'{number:04}.part'
    for attempt in range(8):
        current = part.stat().st_size if part.exists() else 0
        if current == length:
            return number
        if current > length:
            raise ValueError(f'Oversized range {number}')
        first, last = start+current, start+length-1
        temporary = parts/f'{number:04}.incoming'
        header = parts/f'{number:04}.headers'
        with (parts/f'{number:04}.log').open('ab') as log:
            result = subprocess.run(['curl', '--fail', '--location', '--silent', '--show-error',
                '--proxy', args.proxy, '--connect-timeout', '20', '--speed-limit', '1024',
                '--speed-time', '90', '--range', f'{first}-{last}',
                '--dump-header', str(header), '--output', str(temporary), asset['url']],
                stdout=log, stderr=log)
        text = header.read_text(errors='replace') if header.exists() else ''
        matches = re.findall(r'content-range:\s*bytes\s+(\d+)-(\d+)/(\d+)', text, re.I)
        if matches and tuple(map(int, matches[-1])) == (first, last, expected) and temporary.exists():
            received = temporary.stat().st_size
            if received > length-current:
                raise ValueError('Range payload exceeds requested size')
            if received:
                with part.open('ab') as dst, temporary.open('rb') as src:
                    shutil.copyfileobj(src, dst)
                temporary.unlink()
        if result.returncode == 0 and part.exists() and part.stat().st_size == length:
            return number
    raise RuntimeError(f'Range {number} failed after bounded retries')

total = (expected+chunk-1)//chunk
with ThreadPoolExecutor(max_workers=args.workers) as pool:
    jobs = [pool.submit(fetch, number, start) for number, start in enumerate(range(0, expected, chunk))]
    for done, job in enumerate(as_completed(jobs), 1):
        job.result()
        print(f'Completed verified HTTP ranges: {done}/{total}', flush=True)
assembled = target.with_suffix('.assembled')
with assembled.open('wb') as dst:
    for number in range(total):
        with (parts/f'{number:04}.part').open('rb') as src:
            shutil.copyfileobj(src, dst)
if assembled.stat().st_size != expected:
    raise AssertionError('Assembled reference size mismatch')
assembled.replace(target)
digests = {}
for algorithm in ['md5', 'sha256']:
    with target.open('rb') as handle:
        digests[algorithm] = hashlib.file_digest(handle, algorithm).hexdigest()
record = {'status': 'downloaded_size_verified', 'completed_utc': datetime.now(timezone.utc).isoformat(),
    'file': target.name, 'bytes': expected, **digests, 'source_url': asset['url'],
    'collection_id': catalog['collection_id'], 'dataset_id': dataset['dataset_id'],
    'dataset_version_id': dataset['dataset_version_id'], 'doi': catalog['doi'],
    'purpose': 'cell-type annotation reference only; not disease training or external validation',
    'integrity_note': 'Official catalog size and HTTP Content-Range checked; digests are measured, not official references.',
    'donor_independence': 'Different published study and identifiers; biological donor overlap not established.'}
(root/'DOWNLOAD_VERIFIED.json').write_text(json.dumps(record, indent=2)+'\n')
for algorithm in ['md5', 'sha256']:
    (root/(algorithm.upper()+'SUMS')).write_text(f"{digests[algorithm]}  {target.name}\n")
print('REFERENCE_DOWNLOAD_COMPLETE', flush=True)
