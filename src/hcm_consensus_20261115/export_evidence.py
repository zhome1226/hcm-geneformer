"""Archive checked analysis evidence; keep source cell and donor count matrices on Linux."""
import hashlib
import json
import sys
import tarfile
from datetime import datetime, timezone
from pathlib import Path

w=Path(sys.argv[1]).resolve()
def sha(p):
    with p.open('rb') as f:return hashlib.file_digest(f,'sha256').hexdigest()
assert json.loads((w/'audit/FINAL_AUDIT.json').read_text())['status']=='passed'
audit_tag=sys.argv[2] if len(sys.argv)>2 else 'audit'
assert audit_tag.replace('_','').isalnum()
for name in ['prepare','fit','coverage','integrate','review',audit_tag]:
    assert (w/f'logs/{name}.exit').read_text().strip()=='0',name
excluded={f'input/pb/{n}.tsv.gz' for n in ['scp_raw','scp_cb','scp_hhatl','liu','chin','whole']}
files=[]
for directory in ['code','input','stats','coverage','final','review','audit','logs']:
    for p in sorted((w/directory).rglob('*')):
        rel=p.relative_to(w).as_posix()
        if not p.is_file() or p.is_symlink() or '__pycache__' in p.parts or rel in excluded:continue
        if directory=='logs' and p.stem.startswith('export'):continue
        assert p.suffix in {'.py','.R','.c','.md','.tsv','.csv','.gz','.json','.txt','.log','.exit','.sha256'},rel
        files.append(p)
out=w/'export';out.mkdir(exist_ok=True)
assert not (out/'evidence.tar.gz').exists()
manifest={'created_utc':datetime.now(timezone.utc).isoformat(),'server_workspace':str(w),
          'excluded_donor_matrices':sorted(excluded),
          'files':[{'path':p.relative_to(w).as_posix(),'bytes':p.stat().st_size,'sha256':sha(p)} for p in files]}
mp=out/'EVIDENCE_MANIFEST.json';mp.write_text(json.dumps(manifest,indent=2)+'\n')
archive=out/'evidence.tar.gz'
with tarfile.open(archive,'x:gz') as tf:
    for p in files:tf.add(p,arcname=p.relative_to(w).as_posix(),recursive=False)
    tf.add(mp,arcname='EVIDENCE_MANIFEST.json',recursive=False)
receipt={'file':archive.name,'bytes':archive.stat().st_size,'sha256':sha(archive),'files':len(files),'manifest_sha256':sha(mp)}
(out/'EXPORT_RECEIPT.json').write_text(json.dumps(receipt,indent=2)+'\n')
print(json.dumps(receipt))
