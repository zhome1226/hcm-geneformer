"""Audit whether weak cross-lineage marker detection reflects robust mixed programs."""
from pathlib import Path
import json
import anndata as ad
import numpy as np
import pandas as pd
from scipy.io import mmread
from threadpoolctl import threadpool_limits

ROOT=Path.home()/'xiaoe-internship-data'
WORK=ROOT/'_processing/all_data_cleaning_20261111/liu2023'
REF=ROOT/'processed/SCP1303_review_20261111/identity'
PANELS={'Cardiomyocyte':['TNNT2','TTN','TNNI3','ACTC1','MYH7','MYBPC3','MYL2','TNNC1','RYR2','PLN','ACTN2','MYH6'],'Fibroblast':['COL1A1','COL1A2','DCN','LUM','PDGFRA'],'Endothelial':['PECAM1','VWF','KDR','EMCN','CLDN5'],'Immune':['PTPRC','LST1','TYROBP','CD3D','CD79A'],'Mural':['RGS5','PDGFRB','CSPG4','MYH11','ACTA2','TAGLN']}
PAIRS=[('Cardiomyocyte','Fibroblast'),('Cardiomyocyte','Endothelial'),('Cardiomyocyte','Immune'),('Cardiomyocyte','Mural'),('Fibroblast','Endothelial'),('Fibroblast','Immune'),('Endothelial','Immune'),('Immune','Mural')]

def evidence(x,symbols):
 total=np.asarray(x.sum(axis=1)).ravel(); d={}
 for label,panel in PANELS.items():
  values=x[:,symbols.isin(panel).to_numpy()].toarray()
  d[label]={'n':(values>0).sum(axis=1),'n2':(values>=2).sum(axis=1),'umi':values.sum(axis=1),'program':np.log1p(values*(10000/np.maximum(total,1))[:,None]).mean(axis=1)}
 return d

genes=pd.read_csv(REF/'reference_genes.tsv',sep='\t')
meta=pd.read_csv(REF/'calibration_metadata.tsv',sep='\t')
with threadpool_limits(limits=2):
 x=mmread(REF/'calibration_counts.mtx.gz').T.tocsr()
ref=evidence(x,genes.gene_symbol)
truth=meta.broad_type.replace({'Pericyte':'Mural','Smooth_muscle':'Mural'})
bounds={}
for label,d in ref.items():
 positive=d['program'][truth.eq(label).to_numpy()];negative=d['program'][truth.ne(label).to_numpy()]
 bounds[label]={'positive_q10':float(np.quantile(positive,.1)),'negative_q99':float(np.quantile(negative,.99))}
 bounds[label]['threshold']=max(bounds[label].values())

def flags(d):
 weak=np.zeros(len(d['Cardiomyocyte']['n']),dtype=bool);strong=weak.copy()
 for left,right in PAIRS:
  l,r=d[left],d[right]
  ratio=np.minimum(l['program'],r['program'])/np.maximum(np.maximum(l['program'],r['program']),1e-12)
  w=(l['n']>=3)&(r['n']>=3)&(ratio>=.25)
  weak|=w
  strong|=w&(l['n2']>=3)&(r['n2']>=3)&(l['program']>=bounds[left]['threshold'])&(r['program']>=bounds[right]['threshold'])
 return weak,strong

rw,rs=flags(ref)
reference={'nuclei':len(meta),'weak_flags':int(rw.sum()),'strong_flags':int(rs.sum()),'by_type':pd.DataFrame({'truth':truth,'weak':rw,'strong':rs}).groupby('truth').agg(nuclei=('weak','size'),weak=('weak','sum'),strong=('strong','sum')).reset_index().to_dict('records')}
rows=[]
for row in pd.read_csv(WORK/'samples.tsv',sep='\t').itertuples():
 a=ad.read_h5ad(WORK/row.sample_id/'author_counts.h5ad')
 a=a[a.obs.pass_strict_numeric].copy()
 ev=evidence(a.X,a.var.gene_symbol);weak,strong=flags(ev)
 for author_type,idx in a.obs.groupby('author_cellType',observed=True).indices.items():
  rows.append({'sample_id':row.sample_id,'author_type':str(author_type),'nuclei':len(idx),'weak_mixed':int(weak[idx].sum()),'strong_mixed':int(strong[idx].sum())})
 pd.DataFrame({'barcode':a.obs_names,'weak_mixed':weak,'strong_mixed':strong}).to_csv(WORK/row.sample_id/'mixed_program_review.tsv',sep='\t',index=False)
pd.DataFrame(rows).to_csv(WORK/'MIXED_PROGRAM_REVIEW.tsv',sep='\t',index=False)
(WORK/'MIXED_PROGRAM_REVIEW.json').write_text(json.dumps({'bounds':bounds,'reference':reference,'threshold_source':'Fixed held-out adult LV nuclear reference; thresholds max(positive q10,negative q99), plus >=3 panel genes each with >=2 UMIs and original ratio>=.25. Exploratory project diagnostic, not experimentally calibrated doublet probability. Same holdout used to set/report thresholds.','samples_by_author_type':rows},indent=2)+'\n')
print(json.dumps({'reference':reference,'total_weak':sum(r['weak_mixed'] for r in rows),'total_strong':sum(r['strong_mixed'] for r in rows),'CM_weak':sum(r['weak_mixed'] for r in rows if r['author_type']=='CM'),'CM_strong':sum(r['strong_mixed'] for r in rows if r['author_type']=='CM')},indent=2))
