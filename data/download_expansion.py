"""Resume public data downloads and verify bytes/hashes; never load expression data."""
import argparse
import base64
import concurrent.futures
import csv
import fcntl
import gzip
import hashlib
import json
from pathlib import Path
import shutil
import subprocess
import time
from datetime import datetime, timezone


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--catalog', required=True, type=Path)
    parser.add_argument('--root', type=Path, default=Path.home()/'xiaoe-internship-data')
    parser.add_argument('--proxy')
    parser.add_argument('--workers', type=int, default=3)
    args = parser.parse_args()
    root = args.root.resolve()
    work = root/'_download_expansion'
    work.mkdir(parents=True, exist_ok=True)
    assets = json.loads(args.catalog.read_text(encoding='utf-8-sig'))
    needed = 0
    for a in assets:
        p = (root/a['path']).resolve()
        if not p.is_relative_to(root) or not a.get('bytes'):
            raise ValueError('Require bounded target and known nonzero size')
        needed += max(0, a['bytes']-(p.stat().st_size if p.exists() else 0))
    if shutil.disk_usage(root).free < needed + 10*1024**3:
        raise RuntimeError('Insufficient space: require payload plus 10 GiB reserve')

    def fetch(a):
        path = root/a['path']
        path.parent.mkdir(parents=True, exist_ok=True)
        key = hashlib.sha256(a['path'].encode()).hexdigest()[:20]
        verified = work/(key+'.verified.json')
        print('START '+a['path'], flush=True)
        if path.exists() and path.stat().st_size > a['bytes']:
            raise ValueError('Oversized existing file '+str(path))
        with (work/(key+'.curl.log')).open('ab') as log:
            for attempt in range(12):
                if path.exists() and path.stat().st_size == a['bytes']:
                    break
                cmd = ['curl','--fail','--location','--continue-at','-',
                       '--retry','0','--connect-timeout','20','--speed-limit','1024',
                       '--speed-time','90','--output',str(path)]
                if a.get('proxy',False) and args.proxy:
                    cmd += ['--proxy',args.proxy]
                cmd.append(a['url'])
                result = subprocess.run(cmd, stdout=log, stderr=log)
                if result.returncode == 0:
                    break
                print(f'RESUME {a["path"]} attempt {attempt+1}', flush=True)
                if attempt == 11:
                    result.check_returncode()
                time.sleep(min(2**attempt,30))
        if path.stat().st_size != a['bytes']:
            raise ValueError('Size mismatch '+str(path))
        md5, sha = hashlib.md5(), hashlib.sha256()
        part_hashes = []
        crc = None
        if a.get('crc32c'):
            import google_crc32c
            crc = google_crc32c.Checksum()
        with path.open('rb') as f:
            for block in iter(lambda:f.read(8*1024**2),b''):
                md5.update(block)
                sha.update(block)
                if a.get('etag'):
                    part_hashes.append(hashlib.md5(block).digest())
                if crc is not None:
                    crc.update(block)
        actual_crc = base64.b64encode(crc.digest()).decode() if crc is not None else ''
        etag = ''
        if a.get('etag'):
            etag = (hashlib.md5(b''.join(part_hashes)).hexdigest()+'-'+str(len(part_hashes))
                    if '-' in a['etag'] else md5.hexdigest())
            if etag != a['etag']:
                raise ValueError('S3 ETag mismatch (8 MiB multipart assumption): '+str(path))
        for expected,actual,label in [(a.get('md5'),md5.hexdigest(),'MD5'),
                                      (a.get('sha256'),sha.hexdigest(),'SHA256'),
                                      (a.get('crc32c'),actual_crc,'CRC32C')]:
            if expected and expected != actual:
                raise ValueError(label+' mismatch '+str(path))
        if path.name.endswith('.gz'):
            with gzip.open(path,'rb') as f:
                while f.read(8*1024**2):
                    pass
        row = dict(dataset=a['dataset'],relative_path=a['path'],bytes=path.stat().st_size,
                   md5=md5.hexdigest(),sha256=sha.hexdigest(),crc32c=actual_crc,
                   source_url=a['url'],official_hash=('MD5' if a.get('md5') else
                   'SHA256' if a.get('sha256') else 'CRC32C' if a.get('crc32c') else
                   'S3 ETag (8 MiB multipart)' if a.get('etag') else 'none'),
                   verified_utc=datetime.now(timezone.utc).isoformat())
        verified.write_text(json.dumps(row,indent=2)+'\n')
        with (work/'inventory.lock').open('a') as lock:
            fcntl.flock(lock,fcntl.LOCK_EX)
            rows = [json.loads(p.read_text()) for p in work.glob('*.verified.json')]
            rows.sort(key=lambda r:r['relative_path'])
            with (work/'MANIFEST.tsv').open('w',newline='') as out:
                writer=csv.DictWriter(out,fieldnames=list(row),delimiter='\t')
                writer.writeheader()
                writer.writerows(rows)
            (work/'MD5SUMS').write_text(''.join(r['md5']+'  '+r['relative_path']+'\n' for r in rows))
        print('VERIFIED '+a['path'],flush=True)
        return row

    rows, errors = [], []
    with concurrent.futures.ThreadPoolExecutor(max_workers=args.workers) as pool:
        futures = {pool.submit(fetch,a):a for a in assets}
        for future in concurrent.futures.as_completed(futures):
            try:
                rows.append(future.result())
            except Exception as e:
                error=dict(path=futures[future]['path'],error=str(e))
                errors.append(error)
                print('FAILED '+json.dumps(error),flush=True)
    summary = dict(status='complete' if not errors else 'partial_failure',
                   files_verified=len(rows),bytes_verified=sum(r['bytes'] for r in rows),
                   errors=errors,updated_utc=datetime.now(timezone.utc).isoformat())
    (work/(args.catalog.stem+'.summary.json')).write_text(json.dumps(summary,indent=2)+'\n')
    print(json.dumps(summary),flush=True)
    if errors:
        raise SystemExit(1)


if __name__ == '__main__':
    main()
