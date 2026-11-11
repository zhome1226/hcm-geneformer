"""Download author-corrected Wehrens counts; verify Zenodo bytes/MD5 before extraction."""
from datetime import datetime, timezone
from pathlib import Path
import csv
import hashlib
import json
import shutil
import subprocess
import zipfile

ROOT = Path.home() / 'xiaoe-internship-data'
WORK = ROOT / '_processing/all_data_cleaning_20261111/gse138262'
SOURCE = ROOT / 'validation/hcm_GSE138262/author_corrected_zenodo18206554'
SOURCE.mkdir(parents=True, exist_ok=True)
record = json.loads((WORK / 'ZENODO_RECORD_18206553.json').read_text())
assert record['id'] == 18206554
assert record['metadata']['doi'] == '10.5281/zenodo.18206554'
assert shutil.disk_usage(SOURCE).free >= 10 * 1024**3
files = []
for item in record['files']:
    name = item['key']
    assert Path(name).name == name
    expected = item['checksum'].removeprefix('md5:')
    path = SOURCE / name
    def valid():
        if not path.exists() or path.stat().st_size != item['size']:
            return False
        with path.open('rb') as handle:
            return hashlib.file_digest(handle, 'md5').hexdigest() == expected
    if not valid():
        if path.exists() and path.stat().st_size >= item['size']:
            raise RuntimeError(f'Existing complete-sized file failed checksum: {path}')
        print(f'Downloading {name}: {item["size"]} bytes', flush=True)
        subprocess.run(['curl', '--fail', '--location', '--retry', '4', '--retry-delay', '5',
                        '--connect-timeout', '30', '--continue-at', '-', '--output', str(path),
                        item['links']['self']], check=True)
    assert valid(), f'Bytes/MD5 mismatch: {name}'
    with path.open('rb') as handle:
        sha256 = hashlib.file_digest(handle, 'sha256').hexdigest()
    files.append({'file': name, 'bytes': path.stat().st_size, 'md5': expected,
                  'sha256': sha256, 'source_url': item['links']['self'], 'verified': True})
    print(f'Verified {name}: MD5 {expected}', flush=True)
    if name.endswith('.zip'):
        target = SOURCE / name.removesuffix('.zip')
        target.mkdir(exist_ok=True)
        with zipfile.ZipFile(path) as archive:
            assert archive.testzip() is None, f'ZIP CRC failed: {name}'
            for info in archive.infolist():
                dest = (target / info.filename).resolve()
                assert dest.is_relative_to(target.resolve()), info.filename
            archive.extractall(target)
shutil.copy2(WORK / 'ZENODO_RECORD_18206553.json', SOURCE / 'ZENODO_RECORD.json')
(SOURCE / 'MD5SUMS').write_text(''.join(f'{f["md5"]}  {f["file"]}\n' for f in files))
with (SOURCE / 'MANIFEST.tsv').open('w') as handle:
    writer = csv.DictWriter(handle, fieldnames=list(files[0]), delimiter='\t')
    writer.writeheader()
    writer.writerows(files)
inspection = []
for path in sorted(SOURCE.rglob('*counts.table.tsv')):
    if path.name.startswith('._') or '__MACOSX' in path.parts:
        continue
    with path.open() as handle:
        header = next(handle).split()
        first = next(handle).split()
        rows = 1 + sum(1 for _ in handle)
    inspection.append({'path': str(path), 'bytes': path.stat().st_size,
                       'columns': len(header), 'rows_excluding_header': rows,
                       'header_first_12': header[:12], 'first_row_first_12': first[:12],
                       'first_row_columns': len(first)})
assert len(inspection) == 32
result = {'checked_utc': datetime.now(timezone.utc).isoformat(), 'record_id': record['id'],
          'concept_doi_requested': '10.5281/zenodo.18206553',
          'version_doi_resolved': record['metadata']['doi'],
          'source_directory': str(SOURCE), 'files': files, 'count_table_inspection': inspection}
(WORK / 'DOWNLOAD_VERIFICATION.json').write_text(json.dumps(result, indent=2) + '\n')
print(json.dumps({'verified_files': len(files), 'count_tables': len(inspection),
                  'total_download_bytes': sum(f['bytes'] for f in files)}), flush=True)
