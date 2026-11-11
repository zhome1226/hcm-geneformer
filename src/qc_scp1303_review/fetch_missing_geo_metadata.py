"""Fetch official clinical/sample metadata needed to clean existing GEO datasets."""
from concurrent.futures import ThreadPoolExecutor
from datetime import datetime, timezone
from pathlib import Path
import hashlib
import json
import subprocess

root = Path.home() / 'xiaoe-internship-data'
out = root / 'processed/all_data_cleaning_20261111/geo_metadata'
out.mkdir(exist_ok=True)
accessions = ['GSE174691', 'GSE161921', 'GSE138262', 'GSE203274']

def fetch(accession):
    url = f'https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc={accession}&targ=gsm&form=text&view=full'
    destination = out / f'{accession}_samples.soft'
    if not destination.exists() or not destination.read_text(errors='replace').startswith('^SAMPLE'):
        temporary = destination.with_suffix('.incoming')
        subprocess.run(['curl', '--fail', '--location', '--retry', '3', '--connect-timeout', '20',
            '--max-time', '180', '--output', str(temporary), url], check=True, capture_output=True)
        text = temporary.read_text(errors='replace')
        if '^SAMPLE = GSM' not in text or '!Sample_geo_accession' not in text:
            raise ValueError('Unexpected GEO metadata response: ' + accession)
        temporary.replace(destination)
    raw = destination.read_bytes()
    return {'accession': accession, 'path': str(destination), 'source_url': url,
        'bytes': len(raw), 'sha256': hashlib.sha256(raw).hexdigest(),
        'sample_records': raw.count(b'^SAMPLE = GSM')}

with ThreadPoolExecutor(max_workers=4) as pool:
    records = list(pool.map(fetch, accessions))
(out / 'SOURCES.json').write_text(json.dumps({'downloaded_utc': datetime.now(timezone.utc).isoformat(), 'files': records}, indent=2) + '\n')
print(json.dumps(records, indent=2), flush=True)
