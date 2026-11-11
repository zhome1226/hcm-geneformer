"""Require terminal jobs, exact source SHA256, and all verified-output hashes."""
from datetime import datetime, timezone
from pathlib import Path
import hashlib
import json
import pandas as pd

ROOT=Path.home()/'xiaoe-internship-data'
WORK=ROOT/'_processing/all_data_cleaning_20261111/ACDC2017'
OUT=ROOT/'processed/ACDC2017_qc_20261111'


def sha(path,algorithm='sha256'):
    with path.open('rb') as handle:
        return hashlib.file_digest(handle,algorithm).hexdigest()


def main():
    gate=json.loads((OUT/'FINAL_COMPLETION_GATE.json').read_text())
    verify=json.loads((OUT/'INDEPENDENT_VERIFICATION.json').read_text())
    assert gate['status']==verify['status']=='passed' and all(c['passed'] for c in verify['checks'])
    assert gate['verification_SHA256']==sha(OUT/'INDEPENDENT_VERIFICATION.json')
    assert gate['prepared_manifest_SHA256']==sha(OUT/'OUTPUT_MANIFEST.tsv')
    states={}
    for job in ['prepare','verify','geometry']:
        code=int((WORK/(job+'_exit_code.txt')).read_text().strip())
        pid=int((WORK/(job+'_pid.txt')).read_text().strip())
        exists=(Path('/proc')/str(pid)).exists()
        assert code==0 and not exists,(job,code,pid,exists)
        states[job]={'exit_code':code,'pid':pid,'process_exists':exists}
    baseline=json.loads((WORK/'SOURCE_BASELINE.json').read_text())
    records=[]
    for row in baseline:
        path=WORK/row['file']
        actual=sha(path)
        assert actual==row['SHA256'],path.name
        records.append({'file':path.name,'local_SHA256':row['SHA256'],'remote_SHA256':actual,'matches':True})
    for item in pd.read_csv(OUT/'VERIFIED_OUTPUT_MANIFEST.tsv',sep='\t').itertuples(index=False):
        path=Path(item.server_path)
        assert path.stat().st_size==item.bytes and sha(path,'md5')==item.MD5 and sha(path)==item.SHA256,path.name
    geometry_gate=json.loads((OUT/'GEOMETRY_COMPLETION_GATE.json').read_text())
    geometry=json.loads((OUT/'GEOMETRY_VIEWS_VERIFICATION.json').read_text())
    assert geometry_gate['status']==geometry['status']=='passed' and all(c['passed'] for c in geometry['checks'])
    assert geometry_gate['verification_SHA256']==sha(OUT/'GEOMETRY_VIEWS_VERIFICATION.json')
    assert geometry_gate['manifest_SHA256']==sha(OUT/'GEOMETRY_VIEWS_MANIFEST.tsv')
    assert geometry_gate['original_strict_gate_SHA256']==sha(OUT/'FINAL_COMPLETION_GATE.json')
    assert geometry['strict_verification_SHA256']==sha(OUT/'INDEPENDENT_VERIFICATION.json')
    for item in pd.read_csv(OUT/'GEOMETRY_VIEWS_MANIFEST.tsv',sep='\t').itertuples(index=False):
        path=Path(item.server_path)
        assert path.stat().st_size==item.bytes and sha(path,'md5')==item.MD5 and sha(path)==item.SHA256,path.name
    evidence={'status':'passed','checked_utc':datetime.now(timezone.utc).isoformat(),'jobs':states,
        'source_files_matched':len(records),'independent_checks':verify['checks_passed'],
        'verified_output_manifest_SHA256':sha(OUT/'VERIFIED_OUTPUT_MANIFEST.tsv'),'all_verified_outputs_current':True,
        'geometry_view_checks':geometry['checks_passed'],'voxel_grid_usable_pairs':geometry['voxel_grid_usable_pairs'],
        'strict_geometry_usable_pairs':geometry['strict_geometry_usable_pairs'],'geometry_conflict_pairs':geometry['source_geometry_conflict_pairs'],
        'geometry_view_manifest_SHA256':sha(OUT/'GEOMETRY_VIEWS_MANIFEST.tsv')}
    execution=OUT/'execution'
    execution.mkdir(exist_ok=True)
    (execution/'SOURCE_INTEGRITY.json').write_text(json.dumps(records,indent=2)+'\n')
    (execution/'FINAL_EXECUTION_EVIDENCE.json').write_text(json.dumps(evidence,indent=2)+'\n')
    print(json.dumps(evidence,indent=2))


if __name__=='__main__':
    main()

