"""Resume the six selected human GEO samples, then verify gzip and 10x axes."""
import argparse
import concurrent.futures
import csv
import gzip
import hashlib
import json
import os
from pathlib import Path
import shutil
import subprocess
from datetime import datetime, timezone
from urllib.parse import urlparse
from urllib.request import Request, urlopen

os.environ.setdefault('OPENBLAS_NUM_THREADS', '2')
os.environ.setdefault('OMP_NUM_THREADS', '2')
import numpy as np
from scipy.io import mmread


def save_json(path, value):
    tmp = path.with_suffix(path.suffix + '.tmp')
    tmp.write_text(json.dumps(value, ensure_ascii=False, indent=2)+'\n', encoding='utf-8')
    tmp.replace(path)


def download_one(task):
    url, directory = task
    assert urlparse(url).hostname == 'ftp.ncbi.nlm.nih.gov'
    name = Path(urlparse(url).path).name
    target = directory/name
    headers_dir = directory/'source_headers'
    headers_dir.mkdir(exist_ok=True)
    with urlopen(Request(url, method='HEAD'), timeout=60) as response:
        headers = dict(response.headers)
        length = response.headers.get('Content-Length')
    assert length and int(length) > 0, 'Missing upstream length: '+url
    expected = int(length)
    save_json(headers_dir/(name+'.json'), {'url': url, 'headers': headers, 'expected_bytes': expected})
    if not target.exists():
        partial = directory/(name+'.part')
        with (directory/(name+'.wget.log')).open('a') as log:
            subprocess.run(['wget', '-c', '--tries=10', '--timeout=60', '--waitretry=5',
                            '--no-verbose', '-O', str(partial), url], check=True, stdout=log, stderr=log)
        assert partial.stat().st_size == expected, name+': upstream length mismatch'
        with gzip.open(partial, 'rb') as stream:
            while stream.read(4*1024*1024):
                pass  # Full decompression verifies gzip CRC/truncation.
        partial.replace(target)
    assert target.stat().st_size == expected
    with target.open('rb') as f:
        md5 = hashlib.file_digest(f, 'md5').hexdigest()
    with target.open('rb') as f:
        sha = hashlib.file_digest(f, 'sha256').hexdigest()
    with gzip.open(target, 'rb') as stream:
        while stream.read(4*1024*1024):
            pass
    return {'file': name, 'server_path': str(target), 'bytes': expected, 'MD5': md5, 'SHA256': sha,
            'source_url': url, 'upstream_size_match': True, 'gzip_CRC_pass': True,
            'upstream_MD5_status': 'not_provided_in_checked_source_headers; local_MD5_recorded'}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--samples', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    args.output.mkdir(parents=True, exist_ok=True)
    assert shutil.disk_usage(args.output).free >= 10*1024**3, 'Less than 10 GiB free'
    samples = list(csv.DictReader(args.samples.open(encoding='utf-8-sig'), delimiter='\t'))
    assert len(samples) == 6 and all(r['organism'] == 'Homo sapiens' and r['disease'] == 'HCM' for r in samples)
    tasks, lookup = [], {}
    for row in samples:
        directory = args.output/row['accession']/row['sample']
        directory.mkdir(parents=True, exist_ok=True)
        for url in row['matrix_urls'].split(';'):
            tasks.append((url, directory))
            lookup[url] = row
    status = {'started_utc': datetime.now(timezone.utc).isoformat(), 'complete': False,
              'expected_files': len(tasks), 'completed_files': [], 'errors': [], 'sample_checks': [],
              'new_QC_run': False, 'count_matrices_modified': False,
              'upstream_MD5_claim': 'No upstream digest provided; report size/CRC/axes separately from local MD5'}
    save_json(args.output/'DOWNLOAD_STATUS.json', status)
    with concurrent.futures.ThreadPoolExecutor(max_workers=2) as executor:
        futures = {executor.submit(download_one, task): task[0] for task in tasks}
        for future in concurrent.futures.as_completed(futures):
            url = futures[future]
            try:
                record = future.result()
                record.update(dataset=lookup[url]['accession'], sample=lookup[url]['sample'],
                              source_donor_id=lookup[url]['source_donor_id'], disease='HCM',
                              age_stratum=lookup[url]['age_stratum'])
                status['completed_files'].append(record)
                print(json.dumps({'downloaded': record['file'], 'bytes': record['bytes']}), flush=True)
            except Exception as exc:
                status['errors'].append({'url': url, 'error': repr(exc)})
            save_json(args.output/'DOWNLOAD_STATUS.json', status)
    if status['errors']:
        raise RuntimeError('Download failed; inspect DOWNLOAD_STATUS.json; partial files kept for resume')
    for row in samples:
        directory = args.output/row['accession']/row['sample']
        features = next(directory.glob('*_features.tsv.gz'))
        barcodes = next(directory.glob('*_barcodes.tsv.gz'))
        matrix = next(directory.glob('*_matrix.mtx.gz'))
        with gzip.open(features, 'rt') as f:
            genes = [line.rstrip('\n').split('\t')[0] for line in f]
        with gzip.open(barcodes, 'rt') as f:
            cells = [line.rstrip('\n') for line in f]
        with gzip.open(matrix, 'rb') as f:
            x = mmread(f).tocsr()
        assert x.shape == (len(genes), len(cells)), row['sample']+': axes mismatch'
        assert len(genes) == len(set(genes)) and len(cells) == len(set(cells))
        assert np.isfinite(x.data).all() and (x.data >= 0).all() and np.equal(x.data, np.rint(x.data)).all()
        check = {'dataset': row['accession'], 'sample': row['sample'], 'source_donor_id': row['source_donor_id'],
                 'published_barcodes': len(cells), 'features': len(genes), 'stored_entries': int(x.nnz),
                 'all_values_finite_nonnegative_integer': True, 'axes_unique_and_match': True,
                 'status': 'download_verified_not_QC_cleaned', 'disease': 'HCM', 'age_stratum': row['age_stratum']}
        status['sample_checks'].append(check)
        save_json(args.output/'DOWNLOAD_STATUS.json', status)
    records = sorted(status['completed_files'], key=lambda x: (x['dataset'], x['file']))
    with (args.output/'DOWNLOAD_MANIFEST.tsv').open('w', encoding='utf-8', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=list(records[0]), delimiter='\t', lineterminator='\n')
        writer.writeheader(); writer.writerows(records)
    with (args.output/'SAMPLE_MATRIX_CHECKS.tsv').open('w', encoding='utf-8', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=list(status['sample_checks'][0]), delimiter='\t', lineterminator='\n')
        writer.writeheader(); writer.writerows(status['sample_checks'])
    for accession in sorted({r['accession'] for r in samples}):
        selected = [r for r in records if r['dataset'] == accession]
        (args.output/accession/'MD5SUMS').write_text(''.join(
            r['MD5']+'  '+str(Path(r['server_path']).relative_to(args.output/accession))+'\n' for r in selected))
    status['complete'] = len(records) == 18 and len(status['sample_checks']) == 6
    status['completed_utc'] = datetime.now(timezone.utc).isoformat()
    save_json(args.output/'DOWNLOAD_STATUS.json', status)
    assert status['complete']
    print(json.dumps({'status': 'complete', 'files': len(records), 'samples': 6,
                      'bytes': sum(r['bytes'] for r in records), 'QC_run': False}), flush=True)


if __name__ == '__main__':
    main()
