"""Freeze strong mixed-program thresholds using the existing adult LV reference."""
from pathlib import Path
import json
import numpy as np
import pandas as pd
from scipy.io import mmread
from threadpoolctl import threadpool_limits

ROOT=Path.home()/'xiaoe-internship-data'
WORK=Path('/home/zhome/hcm-geneformer-analysis/hcm_expansion_20261114/nuclear_qc')
REF=ROOT/'processed/SCP1303_review_20261111/identity'
panels={'Cardiomyocyte':['TNNT2','TTN','TNNI3','ACTC1','MYH7','MYBPC3','MYL2','TNNC1','RYR2','PLN','ACTN2','MYH6'],'Fibroblast':['COL1A1','COL1A2','DCN','LUM','PDGFRA'],'Endothelial':['PECAM1','VWF','KDR','EMCN','CLDN5'],'Immune':['PTPRC','LST1','TYROBP','CD3D','CD79A'],'Mural':['RGS5','PDGFRB','CSPG4','MYH11','ACTA2','TAGLN']}
pairs=[('Cardiomyocyte','Fibroblast'),('Cardiomyocyte','Endothelial'),('Cardiomyocyte','Immune'),('Cardiomyocyte','Mural'),('Fibroblast','Endothelial'),('Fibroblast','Immune'),('Endothelial','Immune'),('Immune','Mural')]
genes=pd.read_csv(REF/'reference_genes.tsv',sep='\t')
meta=pd.read_csv(REF/'calibration_metadata.tsv',sep='\t')
with threadpool_limits(limits=2):x=mmread(REF/'calibration_counts.mtx.gz').T.tocsr()
total=np.asarray(x.sum(axis=1)).ravel()
truth=meta.broad_type.replace({'Pericyte':'Mural','Smooth_muscle':'Mural'})
bounds={};ev={}
for label,panel in panels.items():
 values=x[:,genes.gene_symbol.isin(panel).to_numpy()].toarray()
 score=np.log1p(values*(10000/np.maximum(total,1))[:,None]).mean(axis=1)
 positive=score[truth.eq(label)];negative=score[truth.ne(label)]
 bounds[label]={'positive_q10':float(np.quantile(positive,.1)),'negative_q99':float(np.quantile(negative,.99))}
 bounds[label]['threshold']=max(bounds[label].values())
 ev[label]={'n':(values>0).sum(axis=1),'n2':(values>=2).sum(axis=1),'program':score}
weak=np.zeros(len(meta),dtype=bool);strong=weak.copy()
for left,right in pairs:
 l,r=ev[left],ev[right]
 ratio=np.minimum(l['program'],r['program'])/np.maximum(np.maximum(l['program'],r['program']),1e-12)
 valid=(l['n']>=3)&(r['n']>=3)&(ratio>=.25)
 weak|=valid
 strong|=valid&(l['n2']>=3)&(r['n2']>=3)&(l['program']>=bounds[left]['threshold'])&(r['program']>=bounds[right]['threshold'])
result={'bounds':bounds,'reference':{'nuclei':len(meta),'weak_flags':int(weak.sum()),'strong_flags':int(strong.sum())},
 'threshold_source':'Fixed adult LV reference holdout: max(positive q10,negative q99) panel score; >=3 genes each >=2 UMIs per lineage and smaller/larger program ratio >=.25. Same holdout defines and evaluates thresholds; not independent doublet calibration. Disease labels are never used.'}
(WORK/'MIXED_PROGRAM_REVIEW.json').write_text(json.dumps(result,indent=2)+'\n')
print(json.dumps(result,indent=2))
