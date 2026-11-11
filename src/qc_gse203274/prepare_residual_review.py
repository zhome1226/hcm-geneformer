"""Find coherent residual mixtures after author doublet removal; freeze batch triggers."""
from pathlib import Path
import json
import anndata as ad
import numpy as np
import pandas as pd
from scipy.io import mmwrite
from threadpoolctl import threadpool_limits
from qc_common import PANELS,PAIRS

ROOT=Path.home()/'xiaoe-internship-data'
WORK=ROOT/'_processing/all_data_cleaning_20261111/GSE203274'
raw=ad.read_h5ad(WORK/'author_full_counts.h5ad')
strict=raw.obs.pass_strict_numeric.to_numpy()
obs=raw.obs.loc[strict,['sample_id','donor_id','disease','author_MainCellType']].copy()
x=raw.X[strict];total=np.asarray(x.sum(axis=1)).ravel()
bounds=json.loads((WORK/'MIXED_PROGRAM_REVIEW.json').read_text())['bounds']
for label,panel in PANELS.items():
    loc=np.flatnonzero(raw.var.gene_symbol.isin(panel).to_numpy())
    assert len(loc)>=3,(label,len(loc))
    values=x[:,loc].toarray()
    obs[label+'_markers']=(values>0).sum(axis=1)
    obs[label+'_markers_ge2_UMI']=(values>=2).sum(axis=1)
    obs[label+'_program']=np.log1p(values*(10000/np.maximum(total,1))[:,None]).mean(axis=1)
mixed=np.zeros(len(obs),dtype=bool);weak=mixed.copy();reason=np.full(len(obs),'',dtype=object)
for left,right in PAIRS:
    l=obs[left+'_program'].to_numpy();r=obs[right+'_program'].to_numpy()
    ratio=np.minimum(l,r)/np.maximum(np.maximum(l,r),1e-12)
    weak|=(obs[left+'_markers'].to_numpy()>=3)&(obs[right+'_markers'].to_numpy()>=3)&(ratio>=.25)
    flag=(obs[left+'_markers_ge2_UMI'].to_numpy()>=3)&(obs[right+'_markers_ge2_UMI'].to_numpy()>=3)&(ratio>=.25)
    flag&=(l>=bounds[left]['threshold'])&(r>=bounds[right]['threshold'])
    reason[flag&~mixed]=left+'+'+right;mixed|=flag
obs['coherent_mixed_panel_flag']=mixed;obs['weak_coherent_mixed_panel_flag']=weak;obs['mixed_panel_pair']=reason
obs['cm_program_dominates']=obs.Cardiomyocyte_program.ge(obs[[p+'_program' for p in PANELS if p!='Cardiomyocyte']].max(axis=1))
obs.to_csv(WORK/'strict_marker_evidence.tsv.gz',sep='\t',index_label='cell_id')
summaries=[]
for sample,group in obs.groupby('sample_id',observed=True):
    count=int(group.coherent_mixed_panel_flag.sum());fraction=count/len(group)
    trigger=count>=10 and fraction>=.01
    all_sample=raw.obs.sample_id.eq(sample).to_numpy()
    model=all_sample&raw.obs.n_genes_by_counts.ge(500).to_numpy()
    directory=WORK/'captures'/sample;directory.mkdir(exist_ok=True,parents=True)
    if trigger:
        with threadpool_limits(limits=4):mmwrite(directory/'residual_model_counts.mtx',raw.X[model].T)
        pd.Series(raw.obs_names[model]).to_csv(directory/'residual_model_cells.txt',index=False,header=False)
    else:
        pd.DataFrame({'cell_id':group.index,'scDblFinder_score':np.nan,'scDblFinder_class':'not_rerun_author_mixed_cluster_QC'}).to_csv(directory/'doublets.tsv',sep='\t',index=False)
    summaries.append({'sample_id':sample,'strict_nuclei':len(group),'strong_mixed_nuclei':count,'strong_mixed_fraction':fraction,
                      'weak_mixed_nuclei':int(group.weak_coherent_mixed_panel_flag.sum()),'residual_model_triggered':trigger,
                      'model_nuclei':int(model.sum()),'model_min_genes':500,
                      'model_gene_count':raw.n_vars,'published_input_post_QC':True})
pd.DataFrame(summaries).to_csv(WORK/'RESIDUAL_REVIEW_SUMMARY.tsv',sep='\t',index=False)
summary={'strict_nuclei':len(obs),'strong_mixed_nuclei':int(mixed.sum()),'weak_mixed_nuclei':int(weak.sum()),
         'triggered_captures':[s['sample_id'] for s in summaries if s['residual_model_triggered']],
         'trigger_rule':'After inherited author CellBender and global/subtype mixed-cluster removal: supplemental scDblFinder only if >=10 AND >=1% strict nuclei have coherent strong mixed programs; model uses ALL published >=500-gene nuclei from that orig.ident, including >5000 and high-MT if present. Residual-screen heuristic, original loading and removal history unknown.',
         'exclusion_rule':'Only model doublet AND coherent strong distinct-lineage pair; mixtures without supplemental support remain pending, not proven doublets.'}
(WORK/'RESIDUAL_REVIEW_POLICY.json').write_text(json.dumps(summary,indent=2)+'\n')
print(json.dumps(summary,indent=2),flush=True)
