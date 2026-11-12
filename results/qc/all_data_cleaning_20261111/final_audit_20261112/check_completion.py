"""Small terminal-record and ACDC geometry supplement to the integrity audit."""
import csv, hashlib, json, sys
from datetime import datetime, timezone
from pathlib import Path

HERE=Path(__file__).resolve().parent
b=json.loads((HERE/'BASELINE.json').read_text())
s=json.loads((HERE/'SUPPLEMENT_BASELINE.json').read_text())
results=[]
def check(label,ok,**details):
    results.append({'check':label,'passed':bool(ok),**details})
for ds,names in s['exit_files'].items():
    for name in names:
        p=Path(b['datasets'][ds]['execution'])/name
        value=p.read_text().strip() if p.exists() else None
        check('terminal_exit_zero',value=='0',dataset=ds,path=str(p),value=value)
for a in s['geometry_artifacts']:
    p=Path(a['path'])
    try:
        data=p.read_bytes()
        ok=len(data)==a['bytes'] and hashlib.md5(data).hexdigest()==a['md5'] and hashlib.sha256(data).hexdigest()==a['sha256']
        check('geometry_artifact_integrity',ok,path=str(p))
    except Exception as e:check('geometry_artifact_integrity',False,path=str(p),error=repr(e))
for r in s['records']:
    d=b['datasets'][r['dataset']]
    paths=[Path(d['output'])/r['name'],Path(d['execution'])/r['name']]
    matches=[p for p in paths if p.is_file() and hashlib.sha256(p.read_bytes()).hexdigest()==r['sha256']]
    check('completion_record_unchanged',bool(matches),dataset=r['dataset'],name=r['name'],matching_paths=list(map(str,matches)))
    if matches and r['name']=='GEOMETRY_VIEWS_VERIFICATION.json':
        j=json.loads(matches[0].read_text())
        check('geometry_2100_checks_passed',j['status']=='passed' and len(j['checks'])==2100 and all(c['passed'] is True for c in j['checks']))
result={'status':'passed' if all(r['passed'] for r in results) else 'failed',
        'checked_utc':datetime.now(timezone.utc).isoformat(),'checks':results,
        'passed':sum(r['passed'] for r in results),'failed':sum(not r['passed'] for r in results)}
(HERE/'COMPLETION_AUDIT.json').write_text(json.dumps(result,indent=2)+'\n')
print(json.dumps({k:v for k,v in result.items() if k!='checks'}))
sys.exit(0 if result['status']=='passed' else 1)
