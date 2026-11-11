"""Read the complete bulk release and retain authoritative source evidence."""
from collections import Counter
from datetime import datetime, timezone
from pathlib import Path
import csv
import gzip
import hashlib
import json
import re
import urllib.request
import urllib.parse
import xml.etree.ElementTree as ET

import numpy as np
import pandas as pd

ROOT = Path.home() / 'xiaoe-internship-data'
SOURCE = ROOT / 'validation/bulk_GSE249925'
WORK = ROOT / '_processing/all_data_cleaning_20261111/GSE249925'
WORK.mkdir(parents=True, exist_ok=True)
LOOKUP = WORK / 'source_lookup'
LOOKUP.mkdir(exist_ok=True)


def parse_soft(path, section='SAMPLE'):
    records = []
    item = None
    for line in path.read_text().splitlines():
        if line.startswith('^' + section + ' = '):
            if item is not None:
                records.append(item)
            item = {'accession': line.split(' = ', 1)[1]}
        elif item is not None and line.startswith('!') and ' = ' in line:
            key, value = line[1:].split(' = ', 1)
            item.setdefault(key, []).append(value)
    if item is not None:
        records.append(item)
    return records


def fetch(url, destination):
    request = urllib.request.Request(url, headers={'User-Agent': 'HCM-source-audit/1.0'})
    with urllib.request.urlopen(request, timeout=25) as response:
        data = response.read()
    destination.write_bytes(data)
    return data


def main():
    baseline = json.loads((WORK / 'ORIGINAL_DOWNLOAD_BASELINE.json').read_text())
    hashes = {}
    for row in baseline:
        path = Path(row['server_path'])
        with path.open('rb') as handle:
            md5 = hashlib.file_digest(handle, 'md5').hexdigest()
        assert path.stat().st_size == int(row['bytes']) and md5 == row['md5'], path
        hashes[path.name] = {'bytes': path.stat().st_size, 'MD5': md5}
    original = SOURCE / 'GSE249925_20231208-For-GEO-Raw-Counts.tsv.gz'
    with gzip.open(original, 'rt') as handle:
        header = next(csv.reader(handle, delimiter='\t'))
    frame = pd.read_csv(original, sep='\t', dtype=str, keep_default_na=False)
    samples = parse_soft(SOURCE / 'sample_metadata.soft')
    series = parse_soft(SOURCE / 'series_metadata.soft', 'SERIES')
    values = frame.iloc[:, 2:].apply(pd.to_numeric, errors='raise').to_numpy()
    titles = [row['Sample_title'][0] for row in samples]
    columns = list(frame.columns[2:])
    characteristics = []
    for row in samples:
        record = {'accession': row['accession'], 'title': row['Sample_title'][0]}
        record.update(dict(item.split(': ', 1) for item in row.get('Sample_characteristics_ch1', [])))
        characteristics.append(record)
    pd.DataFrame(characteristics).to_csv(WORK / 'SOURCE_SAMPLE_CHARACTERISTICS.tsv', sep='\t', index=False)
    processing = sorted({value for row in samples for value in row.get('Sample_data_processing', [])})
    (LOOKUP / 'GEO_PROCESSING_STATEMENTS.txt').write_text('\n\n'.join(processing) + '\n')
    lookups = []
    try:
        query = urllib.parse.urlencode({'query': 'EXT_ID:39255898 AND SRC:MED', 'format': 'json'})
        metadata = json.loads(fetch('https://www.ebi.ac.uk/europepmc/webservices/rest/search?' + query,
                                    LOOKUP / 'PAPER_METADATA.json'))
        hits = metadata['resultList']['result']
        assert len(hits) == 1
        hit = hits[0]
        lookups.append({'source': 'EuropePMC', 'status': 'metadata_retrieved', 'article': hit})
        if hit.get('pmcid'):
            article = fetch('https://www.ebi.ac.uk/europepmc/webservices/rest/' + hit['pmcid'] + '/fullTextXML',
                            LOOKUP / 'PAPER_EFETCH.xml')
            doc = ET.fromstring(article)
            paragraphs = [' '.join(element.itertext()) for element in doc.iter('p')]
            relevant = [text for text in paragraphs if re.search(r'TMM|trimmed.mean|normaliz|HTSeq|HT.Seq|RNA.seq|249925|STAR align', text, re.I)]
            (LOOKUP / 'METHODS_EXCERPTS.txt').write_text('\n\n'.join(relevant) + '\n')
            lookups.append({'source': 'EuropePMC_fullTextXML', 'status': 'retrieved', 'paragraphs': len(relevant)})
    except Exception as error:
        lookups.append({'source': 'EuropePMC', 'status': 'lookup_incomplete', 'error': repr(error)})
    feature_parts = frame['Gene_Index'].str.split('@')
    result = {'dataset': 'GSE249925', 'checked_utc': datetime.now(timezone.utc).isoformat(),
              'shape_genes_by_samples': list(values.shape), 'header_unique': len(header) == len(set(header)),
              'header_matches_parsed': header == list(frame.columns),
              'finite': bool(np.isfinite(values).all()), 'nonnegative': bool((values >= 0).all()),
              'integer_valued': bool(np.equal(values, np.rint(values)).all()),
              'all_zero_features': int((values.sum(axis=1) == 0).sum()),
              'all_zero_samples': int((values.sum(axis=0) == 0).sum()),
              'duplicate_full_feature_identifiers': int(frame.Gene_Index.duplicated(keep=False).sum()),
              'feature_component_counts': dict(Counter(map(len, feature_parts))),
              'count_min': float(values.min()), 'count_max': float(values.max()),
              'library_sum_min': float(values.sum(axis=0).min()), 'library_sum_max': float(values.sum(axis=0).max()),
              'source_GEO_records': len(samples), 'GEO_titles_unique': len(titles) == len(set(titles)),
              'matrix_without_GEO': sorted(set(columns) - set(titles)), 'GEO_without_matrix': sorted(set(titles) - set(columns)),
              'source_groups': dict(Counter(record.get('group', '<missing>') for record in characteristics)),
              'characteristic_fields': sorted(set().union(*(set(record) for record in characteristics))),
              'processing_statements': processing, 'series': series, 'original_hashes': hashes,
              'literature_lookup': lookups}
    (WORK / 'INPUT_INSPECTION.json').write_text(json.dumps(result, indent=2) + '\n')
    print(json.dumps({key: value for key, value in result.items()
                      if key not in ['series', 'original_hashes', 'processing_statements', 'literature_lookup']}, indent=2))
    print('Literature source evidence saved:', json.dumps(lookups), flush=True)


if __name__ == '__main__':
    main()
