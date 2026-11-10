#!/usr/bin/env python3
"""Download and inventory selected public HCM supplements; no expression analysis."""
import argparse
import csv
import gzip
import hashlib
import json
from pathlib import Path
import shutil
import subprocess
import tarfile
import urllib.request

SPECS = [
    dict(dataset='chin_GSE181764', accession='GSE181764', filename='GSE181764_RAW.tar', size=301056000,
         group='1 梗阻性 HCM + 6 非梗阻性 HCM + 2 供体对照', tissue='成人室间隔', purpose='补充验证（Chin 同实验室扩展）'),
    dict(dataset='cox7b_GSE227734', accession='GSE227734', filename='GSE227734_hcm_raw.h5ad.gz', size=2398031,
         group='5 名 HCM 患者；338 个心肌细胞；无配套健康单细胞对照', tissue='成人室间隔', purpose='HCM 心肌细胞特征复核（不单独评估二分类）'),
    dict(dataset='bulk_GSE249925', accession='GSE249925', filename='GSE249925_20231208-For-GEO-Raw-Counts.tsv.gz', size=6968101,
         group='97 例 HCM + 23 例对照（bulk RNA-seq）', tissue='心肌组织（GEO: heart muscle）', purpose='靶点方向验证（bulk；不用于单细胞分类器）'),
]


def digest(path):
    h = hashlib.md5()
    with path.open('rb') as f:
        for block in iter(lambda: f.read(8 * 1024 * 1024), b''):
            h.update(block)
    return h.hexdigest()


def download(url, path, expected=None):
    with urllib.request.urlopen(urllib.request.Request(url, method='HEAD'), timeout=45) as response:
        source_size = response.headers.get('Content-Length')
    if source_size is not None:
        source_size = int(source_size)
        if expected is not None and source_size != expected:
            raise ValueError(f'Source size changed: {url}: {source_size} != {expected}')
    subprocess.run(['wget', '-c', '-nv', '--timeout=45', '--tries=8', '--waitretry=3', '-O', str(path), url], check=True)
    if expected is not None and path.stat().st_size != expected:
        raise ValueError(f'Size mismatch: {path}')
    if source_size is not None and path.stat().st_size != source_size:
        raise ValueError(f'Content-Length mismatch: {path}')


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--base', type=Path, default=Path.home() / 'xiaoe-internship-data/validation')
    parser.add_argument('--metadata-cache', type=Path, required=True, help='Directory of previously fetched unmodified GEO metadata')
    args = parser.parse_args()
    base, cache = args.base.expanduser(), args.metadata_cache.expanduser()
    if shutil.disk_usage(base).free < 10 * 1024**3:
        raise RuntimeError('Less than 10 GiB free')
    manifest_path = base / 'MANIFEST.tsv'
    with manifest_path.open(encoding='utf-8', newline='') as f:
        reader = csv.reader(f, delimiter='\t')
        header = next(reader)
        previous = list(reader)
    dataset_ids = {spec['dataset'] for spec in SPECS}
    retained = [row for row in previous if row[0] not in dataset_ids]
    added, catalog, report = [], [], []
    for spec in SPECS:
        folder = base / spec['dataset']
        folder.mkdir(exist_ok=True)
        accession = spec['accession']
        root_url = f'https://ftp.ncbi.nlm.nih.gov/geo/series/{accession[:-3]}nnn/{accession}/suppl/'
        payload = folder / spec['filename']
        source_urls = {spec['filename']: root_url + spec['filename']}
        download(source_urls[spec['filename']], payload, spec['size'])
        # Preserve published text metadata without rewriting or filtering it.
        for local_name, cached_suffix, target in [('series_metadata.soft', 'series.txt', 'self'), ('sample_metadata.soft', 'samples.soft', 'gsm')]:
            shutil.copyfile(cache / f'{accession}_{cached_suffix}', folder / local_name)
            view = 'full' if target == 'self' else 'brief'
            source_urls[local_name] = f'https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc={accession}&targ={target}&form=text&view={view}'
        sample_text = (folder / 'sample_metadata.soft').read_text()
        for record in sample_text.split('^SAMPLE = ')[1:]:
            lines = record.splitlines()
            gsm = lines[0]
            title = next(line.split(' = ', 1)[1] for line in lines if line.startswith('!Sample_title = '))
            source = next(line.split(' = ', 1)[1] for line in lines if line.startswith('!Sample_source_name_ch1 = '))
            characteristics = '; '.join(line.split(' = ', 1)[1] for line in lines if line.startswith('!Sample_characteristics_ch1 = '))
            catalog.append([accession, gsm, title, source, characteristics, f'https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc={gsm}'])
        if spec['filename'].endswith('.tar'):
            filelist = folder / 'filelist.txt'
            download(root_url + 'filelist.txt', filelist)
            source_urls['filelist.txt'] = root_url + 'filelist.txt'
            with filelist.open() as f:
                rows = list(csv.reader(f, delimiter='\t'))
            archive = next(row for row in rows if row[0] == 'Archive')
            if int(archive[3]) != payload.stat().st_size:
                raise ValueError('Archive size differs from GEO filelist')
            expected_members = {row[1]: int(row[3]) for row in rows if row[0] == 'File'}
            with tarfile.open(payload) as tar:
                members = tar.getmembers()
                if any(not member.isfile() or member.name not in expected_members for member in members):
                    raise ValueError('Unexpected archive member')
                if len(members) != len(expected_members) or {member.name for member in members} != set(expected_members):
                    raise ValueError('Archive member list differs from GEO filelist')
                for member in members:
                    if member.size != expected_members[member.name]:
                        raise ValueError(f'Archive member size mismatch: {member.name}')
            raw = folder / 'raw'
            raw.mkdir(exist_ok=True)
            subprocess.run(['tar', '-xf', str(payload), '-C', str(raw), '--no-same-owner'], check=True)
            for name, expected_size in expected_members.items():
                if (raw / name).stat().st_size != expected_size:
                    raise ValueError(f'Extracted file size mismatch: {name}')
                source_urls['raw/' + name] = root_url + spec['filename']
            actual_raw = {p.relative_to(raw).as_posix() for p in raw.rglob('*') if p.is_file()}
            if actual_raw != set(expected_members):
                raise ValueError('Unexpected extracted file')
        elif spec['filename'].endswith('.gz'):
            # Validate the gzip CRC without saving decompressed data or loading expression matrices.
            with gzip.open(payload, 'rb') as f:
                for block in iter(lambda: f.read(8 * 1024 * 1024), b''):
                    pass
        sums = ['# Actual MD5 and byte counts; verify with: md5sum -c MD5SUMS',
                '# No publisher reference MD5 was available; these are locally measured checksums.',
                '# Payload byte counts match NCBI Content-Length; GEO tar members match official filelist.']
        for path in sorted(folder.rglob('*')):
            if not path.is_file() or path.name == 'MD5SUMS':
                continue
            within = path.relative_to(folder).as_posix()
            checksum, size = digest(path), path.stat().st_size
            sums.extend([f'# bytes={size}  {within}', f'{checksum}  {within}'])
            added.append([spec['dataset'], path.name, path.relative_to(base).as_posix(), size, checksum,
                          source_urls[within], spec['group'], spec['tissue'], spec['purpose']])
        (folder / 'MD5SUMS').write_text('\n'.join(sums) + '\n', encoding='utf-8')
        subprocess.run(['md5sum', '-c', 'MD5SUMS'], cwd=folder, check=True)
        report.append(dict(dataset=spec['dataset'], file=spec['filename'], bytes=payload.stat().st_size,
                           md5=digest(payload), reference_md5_available=False, local_md5_selfcheck=True,
                           http_content_length_match=True))
    with manifest_path.open('w', encoding='utf-8', newline='') as f:
        writer = csv.writer(f, delimiter='\t', lineterminator='\n')
        writer.writerow(header)
        writer.writerows(retained + added)
    with (base / 'SUPPLEMENTAL_SAMPLES.tsv').open('w', encoding='utf-8', newline='') as f:
        writer = csv.writer(f, delimiter='\t', lineterminator='\n')
        writer.writerow(['accession', 'GSM', 'sample_title', 'source', 'published_characteristics', 'source_url'])
        writer.writerows(catalog)
    summary = dict(payloads=report, new_manifest_rows=len(added), total_manifest_rows=len(retained) + len(added), catalog_rows=len(catalog))
    (base / '_download/supplemental_verification.json').write_text(json.dumps(summary, ensure_ascii=False, indent=2) + '\n', encoding='utf-8')
    print(json.dumps(summary, ensure_ascii=False), flush=True)


if __name__ == '__main__':
    main()
