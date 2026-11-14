"""Prepare separate V1 copies with original counts, feature masks and metadata intact."""
from pathlib import Path
from datetime import datetime,timezone
import hashlib,json,re
import anndata as ad
import numpy as np
ad.settings.allow_write_nullable_strings=True

ROOT=Path('/home/zhome/hcm-geneformer-data')
WORK=Path('/home/zhome/hcm-geneformer-analysis/hcm_expansion_20261114/nuclear_qc')
OUT=ROOT/'tokenized/geneformer_v1_20261114'
PREP=OUT/'prepared_inputs'
PREP.mkdir(parents=True,exist_ok=True)
CONFIG={ds:ROOT/f'processed/{ds}_qc_20261114/cm_qc_full.h5ad' for ds in ['GSE255296','GSE292754']}
CONFIG['GSE138262_whole_cell']=ROOT/'processed/GSE138262_whole_cell_qc_20261114/cm_qc_full.h5ad'

def digest(path):
    with Path(path).open('rb') as f:return hashlib.file_digest(f,'sha256').hexdigest()

def exact(a,b):
    return a.shape==b.shape and (a!=b).nnz==0

results=[];frozen={}
for ds,src in CONFIG.items():
    if ds.startswith('GSE138262'):
        gate=json.loads((src.parent/'FINAL_COMPLETION_GATE.json').read_text())
        assert gate['complete'] is True and gate['checks_passed']==71,gate
    else:
        assert json.loads((src.parent/'INDEPENDENT_VERIFICATION.json').read_text())['status']=='passed'
    record={'path':str(src),'bytes':src.stat().st_size,'sha256':digest(src)}
    frozen[ds]=record
    a=ad.read_h5ad(src)
    assert a.n_obs>0 and a.obs_names.is_unique and a.var_names.is_unique
    assert np.isfinite(a.X.data).all() and (a.X.data>=0).all() and np.equal(a.X.data,np.round(a.X.data)).all()
    assert exact(a.X,a.layers['counts'])
    a.X.sum_duplicates();a.X.eliminate_zeros();a.X.sort_indices()
    if 'ensembl_id' in a.var:ids=a.var.ensembl_id.astype(str);basis='source_var_ensembl_id'
    elif 'gene_id' in a.var:ids=a.var.gene_id.astype(str);basis='source_var_gene_id'
    else:ids=a.var_names.astype(str);basis='source_var_index'
    # Source IDs take priority over symbols and composite index labels.
    assert all(re.fullmatch(r'ENSG\d+(?:\.\d+)?',g) for g in ids)
    a.var['ensembl_id']=[g.split('.')[0] for g in ids]
    a.var['source_feature_id']=a.var_names.astype(str)
    original_obs_index_name=a.obs.index.name or ''
    original_var_index_name=a.var.index.name or ''
    a.obs.index.name='source_obs_index'
    a.var.index.name='source_var_index'
    a.obs['cell_id']=[ds+'::'+str(g) for g in a.obs_names]
    if 'source_cell_id' in a.obs:a.obs['original_source_cell_id']=a.obs.source_cell_id.astype(str)
    a.obs['source_cell_id']=a.obs_names.astype(str)
    if 'barcode_original' not in a.obs:
        a.obs['barcode_original']=a.obs['author_cell_id'].astype(str) if 'author_cell_id' in a.obs else a.obs_names.astype(str)
    for col in ['donor_id','sample_id']:
        assert col in a.obs and a.obs[col].notna().all()
        a.obs[col]=a.obs[col].astype(str)
    a.obs['source_disease']=a.obs.disease.astype(str)
    a.obs['label']=a.obs.label.astype(np.int64)
    assert all((int(l)==1 and d=='HCM') or (int(l)==0 and d in ['NF','Control']) for l,d in zip(a.obs.label,a.obs.source_disease))
    a.obs['disease']=np.where(a.obs.label.eq(0),'Control','HCM')
    a.obs['dataset']=ds;a.obs['study']=ds
    for col in ['split','split_status']:
        if col in a.obs:a.obs['source_'+col]=a.obs[col].astype(str)
    a.obs['split']='unassigned'
    a.obs['split_status']='pending_donor_overlap_and_model_history_review'
    a.obs['cell_type']='Cardiomyocyte'
    a.obs['n_counts']=np.asarray(a.X.sum(1)).ravel()
    a.obs['filter_pass']=np.ones(a.n_obs,dtype=np.int8)
    if 'age_stratum' not in a.obs:a.obs['age_stratum']='adult_whole_cell_resource'
    if 'training_use' not in a.obs:a.obs['training_use']='supplementary_pediatric_only' if ds=='GSE292754' else 'unassigned_pending_cohort_review'
    for col in ['age_stratum','age_years','sex','technology','tissue','clinical_stage','training_use','reference_domain','donor_independence','count_processing_history']:
        if col in a.obs:a.obs[col]=a.obs[col].astype(str)
    a.uns['tokenizer_preparation']={'source_path':str(src),'source_sha256':record['sha256'],'model_version':'V1','gene_mapping_basis':basis,'all_source_genes_preserved':True,'matrix_semantics':'unchanged_unscaled_source_counts','training_split_finalized':False,'source_obs_index_name':original_obs_index_name,'source_var_index_name':original_var_index_name}
    target=PREP/(ds+'.h5ad');assert not target.exists(),target
    partial=target.with_suffix('.partial.h5ad');assert not partial.exists(),partial
    a.write_h5ad(partial,compression='gzip')
    b=ad.read_h5ad(partial);original=ad.read_h5ad(src)
    assert exact(b.X,original.X) and exact(b.layers['counts'],original.layers['counts'])
    assert b.obs_names.equals(original.obs_names) and b.var_names.equals(original.var_names)
    for key in original.varm:assert np.array_equal(b.varm[key],original.varm[key])
    for col in original.obs:
        if col in ['source_cell_id','disease','split','split_status']:continue
        assert np.array_equal(b.obs[col].astype(str),original.obs[col].astype(str)),(ds,col)
    assert digest(src)==record['sha256']
    partial.rename(target)
    results.append({'dataset':ds,'source':str(src),'source_bytes':record['bytes'],'source_sha256':record['sha256'],'prepared_path':str(target),'prepared_bytes':target.stat().st_size,'prepared_sha256':digest(target),'cells':b.n_obs,'genes':b.n_vars,'donors':b.obs.donor_id.nunique(),'gene_mapping_basis':basis,'source_counts_and_feature_masks_preserved':True,'disease_counts':b.obs.disease.value_counts().to_dict(),'split':'unassigned'})
    print(json.dumps(results[-1]),flush=True)
(WORK/'INPUT_CONFIG.json').write_text(json.dumps(frozen,indent=2)+'\n')
report={'status':'passed','utc':datetime.now(timezone.utc).isoformat(),'datasets':results,'total_cells':sum(r['cells'] for r in results)}
(PREP/'PREPARATION_VERIFICATION.json').write_text(json.dumps(report,indent=2)+'\n')
