#!/usr/bin/env python3
"""Archive competition code and small public data catalogs; never run third-party code."""
import argparse
import base64
import concurrent.futures
import csv
import gzip
import hashlib
import json
from pathlib import Path
import subprocess
import urllib.request


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--base', type=Path, default=Path.home()/'xiaoe-internship-data/benchmarks/competitions')
    parser.add_argument('--repository-inventory', type=Path, required=True)
    parser.add_argument('--vcc-catalog', type=Path, required=True)
    parser.add_argument('--verify-existing-code', action='store_true', help='Verify code archives already transferred from the PC')
    args = parser.parse_args()
    base = args.base.expanduser()
    repos = json.loads(args.repository_inventory.read_text(encoding='utf-8-sig'))
    catalog = json.loads(args.vcc_catalog.read_text())
    assets = []
    for repo in repos:
        if repo['license'] not in ('MIT', 'Apache-2.0'):
            continue
        folder = repo['repository'].replace('/', '__')
        assets.append(dict(kind='code_archive', name=repo['repository'],
                           relative_path=f'code/{folder}/{repo["sha"]}.tar.gz',
                           url=repo['source_archive'], reference=repo['sha'],
                           license=repo['license'], expected_size=None, expected_md5=None))
    for obj in catalog:
        if not obj['name'].endswith('.csv'):
            continue
        assets.append(dict(kind='dataset_metadata', name='VCC2025',
                           relative_path='vcc2025/'+obj['name'].split('2025/',1)[1],
                           url='https://storage.googleapis.com/arc-institute-virtual-cell-atlas/'+obj['name'],
                           reference='GCS object MD5', license='dataset terms are separate from code',
                           expected_size=int(obj['size']),
                           expected_md5=base64.b64decode(obj['md5Hash']).hex()))

    def fetch(asset):
        path = base/asset['relative_path']
        path.parent.mkdir(parents=True, exist_ok=True)
        if asset['kind']=='code_archive' and args.verify_existing_code:
            if not path.is_file() or not path.stat().st_size:
                raise FileNotFoundError(path)
        elif (asset['expected_md5'] is not None and path.is_file()
              and path.stat().st_size==asset['expected_size']
              and hashlib.md5(path.read_bytes()).hexdigest()==asset['expected_md5']):
            pass
        else:
            subprocess.run(['curl','-fL','-sS','-C','-','--connect-timeout','15','--max-time','300',
                            '--retry','2','--retry-delay','2','-o',str(path),asset['url']],check=True)
        h = hashlib.md5()
        with path.open('rb') as f:
            for block in iter(lambda:f.read(8*1024*1024), b''):
                h.update(block)
        md5, size = h.hexdigest(), path.stat().st_size
        if asset['expected_size'] is not None and size != asset['expected_size']:
            raise ValueError(f'Size mismatch: {path}')
        if asset['expected_md5'] is not None and md5 != asset['expected_md5']:
            raise ValueError(f'Upstream MD5 mismatch: {path}')
        if path.name.endswith('.gz'):
            with gzip.open(path,'rb') as f:
                for block in iter(lambda:f.read(8*1024*1024),b''):
                    pass
        return [asset['kind'],asset['name'],asset['relative_path'],size,md5,asset['url'],
                asset['reference'],asset['license'],
                'upstream_MD5_match' if asset['expected_md5'] else 'local_MD5_and_gzip_CRC_only']

    with concurrent.futures.ThreadPoolExecutor(max_workers=3) as pool:
        rows=list(pool.map(fetch,assets))
    base.mkdir(parents=True,exist_ok=True)
    with (base/'MANIFEST.tsv').open('w',encoding='utf-8',newline='') as f:
        writer=csv.writer(f,delimiter='\t',lineterminator='\n')
        writer.writerow(['type','resource','relative_path','bytes','MD5','source_url','upstream_reference','code_license_or_note','verification'])
        writer.writerows(rows)
    (base/'MD5SUMS').write_text(''.join(f'{row[4]}  {row[2]}\n' for row in rows))
    subprocess.run(['md5sum','-c','MD5SUMS'],cwd=base,check=True)
    summary={'files':len(rows),'code_archives':sum(row[0]=='code_archive' for row in rows),
             'metadata_files':sum(row[0]=='dataset_metadata' for row in rows),
             'bytes':sum(row[3] for row in rows),'expression_matrices_downloaded':0,'third_party_code_executed':False}
    (base/'_catalog/download_summary.json').write_text(json.dumps(summary,indent=2)+'\n')
    print(json.dumps(summary),flush=True)


if __name__ == '__main__':
    main()
