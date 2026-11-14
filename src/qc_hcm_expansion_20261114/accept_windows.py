"""Extract to staging, verify every hash, load the three complete datasets, then promote."""
from pathlib import Path
from datetime import datetime,timezone
import hashlib,json,tarfile,shutil
import numpy as np,pandas as pd,anndata as ad
from datasets import load_from_disk
ROOT=Path('D:/hcm-geneformer-data')
INCOMING=ROOT/'_incoming/hcm_expansion_20261114'
PACKAGE=json.loads((INCOMING/'PACKAGE_STATUS.json').read_text())
ARCHIVE=INCOMING/Path(PACKAGE['archive']).name
def digest(p,algorithm='sha256'):
    h=hashlib.new(algorithm)
    with p.open('rb') as stream:
        for part in iter(lambda:stream.read(4*1024*1024),b''):h.update(part)
    return h.hexdigest()
assert ARCHIVE.stat().st_size==PACKAGE['archive_bytes'] and digest(ARCHIVE)==PACKAGE['archive_sha256']
stage=INCOMING/'unpacked';assert not stage.exists(),stage
stage.mkdir()
with tarfile.open(ARCHIVE) as tar:
    for m in tar.getmembers():
        assert not m.issym() and not m.islnk()
        resolved=(stage/m.name).resolve()
        assert resolved.is_relative_to(stage.resolve()),m.name
    tar.extractall(stage,filter='data')
payload=stage/'hcm_expansion_20261114'
manifest=pd.read_csv(payload/'DELIVERY_MANIFEST.tsv',sep='\t')
assert set(manifest.relative_path)=={p.relative_to(payload).as_posix() for p in payload.rglob('*') if p.is_file() and p.name!='DELIVERY_MANIFEST.tsv'}
for r in manifest.itertuples():
    p=payload/r.relative_path;assert p.stat().st_size==r.bytes and digest(p)==r.sha256 and digest(p,'md5')==r.md5,(r.relative_path,'hash mismatch')
token=payload/'tokenized/geneformer_v1_20261114'
prepared=json.loads((token/'PREPARATION_VERIFICATION.json').read_text())
accepted=[]
for r in prepared['datasets']:
    ds=r['dataset'];p=token/'datasets'/(ds+'.dataset')
    assert (p/'state.json').exists() and (p/'dataset_info.json').exists() and list(p.glob('*.arrow'))
    hf=load_from_disk(str(p));assert len(hf)==r['cells']
    source=payload/'processed'/Path(r['source']).parent.name/'cm_qc_full.h5ad'
    a=ad.read_h5ad(source)
    assert a.n_obs==len(hf) and a.n_vars==r['genes']
    assert a.obs.donor_id.nunique()==r['donors'] and a.obs_names.is_unique
    assert (a.X!=a.layers['counts']).nnz==0
    assert np.isfinite(a.X.data).all() and (a.X.data>=0).all() and np.equal(a.X.data,np.round(a.X.data)).all()
    offset=0
    for batch in hf.iter(batch_size=256):
        end=offset+len(batch['input_ids']);obs=a.obs.iloc[offset:end]
        assert batch['cell_id']==[ds+'::'+str(v) for v in obs.index]
        assert batch['source_cell_id']==list(obs.index.astype(str))
        for c in ['donor_id','sample_id','label']:
            assert list(map(str,batch[c]))==list(obs[c].astype(str)),(ds,c,offset)
        assert all(v=='unassigned' for v in batch['split'])
        assert np.array_equal(batch['n_counts'],np.asarray(a.X[offset:end].sum(1)).ravel())
        for seq,length in zip(batch['input_ids'],batch['length']):
            assert 0<length==len(seq)<=2048 and len(set(seq))==len(seq)
        offset=end
    assert offset==len(hf)
    accepted.append({'dataset':ds,'cells':len(hf),'genes':a.n_vars,'donors':r['donors'],'label_counts':a.obs.label.value_counts().to_dict(),'hf_files':[q.name for q in sorted(p.iterdir()) if q.is_file()],'full_metadata_and_lengths_checked':True})
assert sum(r['cells'] for r in accepted)==prepared['total_cells']
# Only verified new dated targets are exposed in the official tree.
dest_token=ROOT/'tokenized/geneformer_v1_20261114';assert not dest_token.exists()
for p in (payload/'processed').iterdir():assert not (ROOT/'processed'/p.name).exists()
(ROOT/'processed').mkdir(exist_ok=True);(ROOT/'tokenized').mkdir(exist_ok=True)
for p in (payload/'processed').iterdir():shutil.move(str(p),str(ROOT/'processed'/p.name))
shutil.move(str(token),str(dest_token))
audit=ROOT/'deliveries/hcm_expansion_20261114';assert not audit.exists();audit.mkdir(parents=True)
for name in ['DELIVERY_README.md','DELIVERY_MANIFEST.tsv']:shutil.copy2(payload/name,audit/name)
result={'status':'passed','host':'DESKTOP-74N1LUH','utc':datetime.now(timezone.utc).isoformat(),'root':str(ROOT),'files_verified':len(manifest),'all_size_and_sha256_match':True,'all_md5_match':True,'archive_sha256':PACKAGE['archive_sha256'],'datasets':accepted,'total_cells':sum(r['cells'] for r in accepted),'official_token_directory':str(dest_token),'training_split':'unassigned'}
(audit/'WINDOWS_ACCEPTANCE.json').write_text(json.dumps(result,indent=2)+'\n')
print(json.dumps(result),flush=True)
