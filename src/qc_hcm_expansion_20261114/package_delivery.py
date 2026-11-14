"""Package the three accepted QC cores and complete HF folders without touching older deliveries."""
from pathlib import Path
from datetime import datetime,timezone
import hashlib,json,shutil,tarfile
import pandas as pd
ROOT=Path('/home/zhome/hcm-geneformer-data')
TOKEN=ROOT/'tokenized/geneformer_v1_20261114'
DEST=ROOT/'deliveries/hcm_expansion_20261114'
assert not DEST.exists(),DEST
assert json.loads((TOKEN/'FULL_VERIFICATION.json').read_text())['status']=='passed'
assert json.loads((TOKEN/'SMOKE_VERIFICATION.json').read_text())['status']=='passed'
DEST.mkdir(parents=True)
records=json.loads((TOKEN/'prepared_inputs/PREPARATION_VERIFICATION.json').read_text())['datasets']
report=json.loads((TOKEN/'FULL_VERIFICATION.json').read_text())
for rec in records:
    ds=rec['dataset'];src=Path(rec['source']);out=DEST/'processed'/src.parent.name
    out.mkdir(parents=True)
    shutil.copy2(src,out/src.name)
    for name in ['INDEPENDENT_VERIFICATION.json','FINAL_COMPLETION_GATE.json','POLICY.json','qc_policy.json','DONOR_RETENTION.tsv','SAMPLE_RETENTION.tsv','sample_qc_summary.tsv','sample_metadata.tsv','REPORT.md','REFERENCE_CALIBRATION.json','MIXED_PROGRAM_REVIEW.json','PENDING_REASON_SUMMARY.tsv']:
        if (src.parent/name).is_file():shutil.copy2(src.parent/name,out/name)
    source=TOKEN/'datasets'/(ds+'.dataset');destination=DEST/'tokenized/geneformer_v1_20261114/datasets'/(ds+'.dataset')
    assert (source/'state.json').is_file() and (source/'dataset_info.json').is_file() and list(source.glob('*.arrow'))
    shutil.copytree(source,destination)
small=DEST/'tokenized/geneformer_v1_20261114'
for name in ['FULL_VERIFICATION.json','SMOKE_VERIFICATION.json','FULL_MANIFEST.tsv','ENCODING_ENVIRONMENT.json']:
    shutil.copy2(TOKEN/name,small/name)
for p in TOKEN.glob('full_*_DONOR_COVERAGE.tsv'):shutil.copy2(p,small/p.name)
shutil.copy2(TOKEN/'prepared_inputs/PREPARATION_VERIFICATION.json',small/'PREPARATION_VERIFICATION.json')
readme=['# HCM Geneformer V1 delivery','',f'Packaged UTC: {datetime.now(timezone.utc).isoformat()}','',
    'Destination root: D:/hcm-geneformer-data; each dataset is separate; no training/validation roles assigned.',
    'QC cores are unchanged source-count AnnData files, with all native genes and any native missingness masks.',
    'Geneformer V1: 2048-token maximum, no CLS/EOS, official median-scaled ranking. All tokens and metadata independently verified before packaging.',
    'New nuclear cohorts use 500<=detected genes<=5000 and MT<5%; whole-cell GSE138262 uses author-specific non-MT UMI>1000 and identity consensus.',
    'GSE255296 HCM-only (age unresolved); GSE292754 HCM-only pediatric females 14/16; whole-cell disease/technology/study confounded.',
    'No ambient-based independent droplet calls, no verified ambient correction history for new matrices, no independent validation/training approval.',
    '', '| Dataset | CM cells/nuclei | Donor identifiers | Disease counts |', '|---|---:|---:|---|']
for r in records:readme.append(f"| {r['dataset']} | {r['cells']} | {r['donors']} | {r['disease_counts']} |")
(DEST/'DELIVERY_README.md').write_text('\n'.join(readme)+'\n')
manifest=[]
for p in sorted(DEST.rglob('*')):
    if p.is_file():
        raw=p.read_bytes();manifest.append({'relative_path':p.relative_to(DEST).as_posix(),'bytes':len(raw),'md5':hashlib.md5(raw).hexdigest(),'sha256':hashlib.sha256(raw).hexdigest()})
pd.DataFrame(manifest).to_csv(DEST/'DELIVERY_MANIFEST.tsv',sep='\t',index=False)
archive=DEST.with_suffix('.tar.gz');assert not archive.exists()
with tarfile.open(archive,'w:gz',compresslevel=1) as tar:tar.add(DEST,arcname=DEST.name)
with archive.open('rb') as stream:sha=hashlib.file_digest(stream,'sha256').hexdigest()
result={'status':'packaged','utc':datetime.now(timezone.utc).isoformat(),'archive':str(archive),'archive_bytes':archive.stat().st_size,'archive_sha256':sha,'payload_files':len(manifest),'payload_bytes':sum(r['bytes'] for r in manifest),'dataset_cells':{r['dataset']:r['cells'] for r in records}}
(DEST.parent/'PACKAGE_STATUS.json').write_text(json.dumps(result,indent=2)+'\n')
print(json.dumps(result),flush=True)
