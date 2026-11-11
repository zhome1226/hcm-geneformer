"""Safely stage existing data archives and official clinical metadata; originals remain read-only."""
from datetime import datetime, timezone
from pathlib import Path
import gzip
import hashlib
import json
import shutil
import tarfile

import h5py
import pandas as pd

root = Path.home() / 'xiaoe-internship-data'
work = root / '_processing/all_data_cleaning_20261111'
out = root / 'processed/all_data_cleaning_20261111'
records = []
for dataset, source in [
    ('liu2023', root/'validation/liu2023/snRNA-seq.tar.gz'),
    ('GSE203274', root/'validation/pediatric_GSE203274/GSE203274_AllNuclei_snRNA_rawCount.tar.gz')]:
    target = work / dataset
    target.mkdir(exist_ok=True)
    with tarfile.open(source, 'r:*') as archive:
        members = archive.getmembers()
        for member in members:
            resolved = (target / member.name).resolve()
            if not resolved.is_relative_to(target.resolve()):
                raise ValueError('Unsafe archive member: '+member.name)
        existing = all((target/m.name).exists() and (not m.isfile() or (target/m.name).stat().st_size==m.size) for m in members)
        if not existing:
            archive.extractall(target, filter='data')
    records.append({'dataset':dataset, 'source':str(source), 'target':str(target), 'member_count':len(members)})
    print('Staged archive: '+dataset, flush=True)

sources = [
    ('GSE227734', root/'validation/cox7b_GSE227734/GSE227734_hcm_raw.h5ad.gz'),
    ('NeurIPS2021_CITE', root/'benchmarks/competitions/neurips2021/GSE194122_openproblems_neurips2021_cite_BMMC_processed.h5ad.gz'),
    ('NeurIPS2021_multiome', root/'benchmarks/competitions/neurips2021/GSE194122_openproblems_neurips2021_multiome_BMMC_processed.h5ad.gz')]
for dataset, source in sources:
    target_dir = work / dataset
    target_dir.mkdir(exist_ok=True)
    destination = target_dir / source.name.removesuffix('.gz')
    if not destination.exists():
        temporary = destination.with_suffix('.incoming')
        digest = hashlib.sha256()
        total = 0
        with gzip.open(source, 'rb') as src, temporary.open('wb') as dst:
            while chunk := src.read(32*1024**2):
                if shutil.disk_usage(target_dir).free < 10*1024**3:
                    raise OSError('Insufficient staging space; minimum 10GiB reserve')
                total += len(chunk)
                if total > 70*1024**3:
                    raise OSError('Unexpected decompressed file size')
                dst.write(chunk); digest.update(chunk)
        with h5py.File(temporary) as handle:
            if not {'obs','var','X'}.issubset(handle):
                raise ValueError('Decompressed file is not expected H5AD')
        temporary.replace(destination)
        (target_dir/'DECOMPRESSION.json').write_text(json.dumps({'source':str(source),'bytes':total,'sha256':digest.hexdigest()},indent=2)+'\n')
    records.append({'dataset':dataset, 'source':str(source), 'target':str(destination), 'bytes':destination.stat().st_size})
    print('Staged H5AD: '+dataset, flush=True)

softs = list((out/'geo_metadata').glob('*_samples.soft')) + [
    root/'validation/chin_GSE181764/sample_metadata.soft',
    root/'validation/cox7b_GSE227734/sample_metadata.soft',
    root/'validation/bulk_GSE249925/sample_metadata.soft']
samples = []
for path in softs:
    current = None
    for line in path.read_text(errors='replace').splitlines():
        if line.startswith('^SAMPLE = '):
            if current:
                samples.append(current)
            current = {'accession':line.split(' = ',1)[1], 'source_path':str(path)}
        elif current is not None and line.startswith('!Sample_') and ' = ' in line:
            key, value = line.split(' = ',1)
            if key in ['!Sample_title','!Sample_source_name_ch1','!Sample_organism_ch1','!Sample_characteristics_ch1','!Sample_description','!Sample_data_processing','!Sample_library_strategy']:
                current.setdefault(key.removeprefix('!Sample_'),[]).append(value)
    if current:
        samples.append(current)
(out/'GEO_SAMPLES.json').write_text(json.dumps(samples,indent=2,ensure_ascii=False)+'\n')
pd.DataFrame([{key: (' | '.join(value) if isinstance(value,list) else value) for key,value in row.items()} for row in samples]).to_csv(out/'GEO_SAMPLES.tsv',sep='\t',index=False)
(out/'STAGED_INPUTS.json').write_text(json.dumps({'prepared_utc':datetime.now(timezone.utc).isoformat(),'files':records},indent=2)+'\n')
print('EXTERNAL_INPUTS_STAGED',flush=True)
