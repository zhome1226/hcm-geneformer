"""Read-only final integrity audit; does not rerun cleaning or alter data."""
import csv
import hashlib
import json
import sys
import time
from datetime import datetime, timezone
from pathlib import Path
import h5py

HERE=Path(__file__).resolve().parent
baseline=json.loads((HERE/'BASELINE.json').read_text())
results=[]
errors=[]
seen={}
start=time.time()

def record(kind,ok,**details):
    r={'kind':kind,'passed':bool(ok),**details}
    results.append(r)
    with (HERE/'PROGRESS.jsonl').open('a') as out:
        out.write(json.dumps(r)+'\n')
    if not ok: errors.append(r)
    return r

def digest_file(path):
    p=Path(path)
    stat=p.stat()
    md5=hashlib.md5();sha=hashlib.sha256()
    with p.open('rb') as f:
        for block in iter(lambda:f.read(8*1024*1024),b''):
            md5.update(block);sha.update(block)
    assert stat.st_size==p.stat().st_size and stat.st_mtime_ns==p.stat().st_mtime_ns,'file_changed_during_read'
    return {'bytes':stat.st_size,'md5':md5.hexdigest(),'sha256':sha.hexdigest()}

def check_artifact(a):
    p=a['path']
    try:
        if p not in seen:
            seen[p]=digest_file(p)
            if p.lower().endswith(('.h5ad','.h5','.hdf5')):
                with h5py.File(p,'r') as f:
                    record('HDF5_structure_readable',True,path=p,top_level_keys=list(f.keys()))
        actual=seen[p]
        for field in ('bytes','md5','sha256'):
            if a.get(field) is not None:
                assert str(actual[field]).lower()==str(a[field]).lower(),field+'_mismatch'
        record('file_integrity',True,dataset=a['dataset'],path=p,actual=actual,baseline=a.get('baseline','anchored_source_manifest'))
    except Exception as e:
        record('file_integrity',False,dataset=a['dataset'],path=p,error=repr(e))

def read_table(path):
    with path.open(encoding='utf-8-sig') as f:
        return list(csv.DictReader(f,delimiter='\t'))

for ds,d in baseline['datasets'].items():
    p=Path(d['output'])
    try:
        vfile=p/'INDEPENDENT_VERIFICATION.json'
        sha=hashlib.sha256(vfile.read_bytes()).hexdigest()
        assert sha==d['verification_sha256'],'verification_record_changed'
        v=json.loads(vfile.read_text())
        c=v.get('checks',v.get('all_checks'))
        flags=list(c.values()) if isinstance(c,dict) else [x['passed'] for x in c]
        assert len(flags)==d['independent_checks'] and all(x is True for x in flags)
        assert v['status'] in ('passed','complete_verified')
        record('independent_verification_record',True,dataset=ds,checks=len(flags),sha256=sha)
    except Exception as e:
        record('independent_verification_record',False,dataset=ds,error=repr(e))

for a in baseline['artifacts']:
    check_artifact(a)

# ORIGINAL_FILES tables are already hashed by the baseline above. Their records
# anchor the unchanged large inputs used by index-based competition deliveries.
for ds in ('NeurIPS2021','CausalBench2023','ACDC2017'):
    p=Path(baseline['datasets'][ds]['output'])/'ORIGINAL_FILES.tsv'
    try:
        assert str(p) in seen,'source_manifest_not_in_hashed_baseline'
        for r in read_table(p):
            r={k.lower():v for k,v in r.items()}
            check_artifact({'dataset':ds,'path':r.get('path',r.get('source_path')),
                            'bytes':int(r['bytes']),'md5':r['md5'],'sha256':r['sha256'],
                            'baseline':str(p)})
    except Exception as e:
        record('anchored_source_manifest',False,dataset=ds,error=repr(e))

# Check path references in small, delivered index tables. Full row/feature
# membership was independently verified earlier and is protected by these hashes.
for ds,d in baseline['datasets'].items():
    p=Path(d['output'])
    for name in ('ASSAY_INDEX.tsv','DATASET_INDEX.tsv','USABLE_IMAGE_MASK_PAIRS.tsv',
                 'USABLE_VOXEL_GRID_IMAGE_MASK_PAIRS.tsv'):
        f=p/name
        if not f.exists():continue
        try:
            refs=set()
            for row in read_table(f):
                for key,value in row.items():
                    if not value:continue
                    if value.startswith('/home/zhome/xiaoe-internship-data/'):
                        refs.add(value)
                    elif key in ('row_index','column_index'):
                        refs.add(str(p/value))
            missing=[v for v in sorted(refs) if not Path(v).is_file()]
            assert not missing,missing
            record('index_references_exist',True,dataset=ds,path=str(f),references=len(refs))
        except Exception as e:
            record('index_references_exist',False,dataset=ds,path=str(f),error=repr(e))

# Preserve historical failed attempts; inspect current terminal exit records only.
execution=[]
for ds,d in baseline['datasets'].items():
    q=Path(d['execution'])
    codes=[]
    for f in sorted(q.glob('*exit_code.txt')):
        text=f.read_text().strip()
        codes.append({'file':str(f),'code':text})
    execution.append({'dataset':ds,'current_exit_records':codes})
(HERE/'EXECUTION_RECORDS.json').write_text(json.dumps(execution,indent=2)+'\n')

summary={'status':'passed' if not errors else 'failed','checked_utc':datetime.now(timezone.utc).isoformat(),
         'datasets_in_scope':len(baseline['datasets']),'excluded':baseline['excluded'],
         'unique_files_rehashed':len(seen),'total_bytes_rehashed':sum(x['bytes'] for x in seen.values()),
         'historical_independent_checks_confirmed':sum(d['independent_checks'] for d in baseline['datasets'].values()),
         'fresh_checks_passed':sum(r['passed'] for r in results),'fresh_checks_failed':len(errors),
         'elapsed_seconds':round(time.time()-start,2),'errors':errors,
         'scope':'Integrity, exact independent-verification records, HDF5 structure and index targets; no recleaning or training-eligibility certification.'}
(HERE/'FINAL_INTEGRITY_AUDIT.json').write_text(json.dumps(summary,indent=2)+'\n')
print(json.dumps(summary),flush=True)
sys.exit(0 if not errors else 1)
