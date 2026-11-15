"""All-gene exploratory consensus; evidence tiers are distinct from V1 coverage."""
import hashlib
import json
import sys
from datetime import datetime, timezone
from pathlib import Path
import numpy as np
import pandas as pd

def sha(p):
    with Path(p).open('rb') as f:return hashlib.file_digest(f,'sha256').hexdigest()
def read(p):return pd.read_csv(p,sep='\t',keep_default_na=False,na_values=['NA'])
def save(df,p):df.to_csv(p,sep='\t',index=False,na_rep='NA',float_format='%.15g')
def boolean(s):
    m=s.map({True:True,False:False,'TRUE':True,'FALSE':False,'True':True,'False':False})
    assert m[s.notna()].notna().all()
    return m.fillna(False).astype(bool)
def direction(a,b):
    return np.select([~np.isfinite(a)|~np.isfinite(b),(a==0)|(b==0),a*b>0,a*b<0],
                     ['not_evaluable','zero','same','opposite'],default='not_evaluable')
def ranked(x):
    x=x.sort_values(['rank_score_worst_FDR','rank_score_min_abs_logFC','gene_id'],ascending=[True,False,True],kind='stable').copy()
    x.insert(0,'candidate_rank',np.arange(1,len(x)+1))
    return x

def main():
    w=Path(sys.argv[1]);out=w/'final'
    assert not out.exists()
    for tag in ['prepare','fit','coverage']:assert (w/f'logs/{tag}.exit').read_text().strip()=='0',tag
    prep=json.loads((w/'input/PREPARATION.json').read_text())
    assert prep['design_commit']=='aefa26c' and prep['design_sha256']==sha(w/'code/DESIGN.md')
    for x in prep['files']:assert sha(w/x['path'])==x['sha256'],x['path']
    assert sha(w/'input/GENE_UNIVERSE.tsv')==prep['gene_universe_sha256']
    previous=read(w/'input/previous/CANDIDATES.tsv')
    assert sha(w/'input/previous/CANDIDATES.tsv')=='d2fc63ed14d76c8d291afd09c067c0339584f38550369c8691248e0126bfa09d'
    d=read(w/'input/GENE_UNIVERSE.tsv');ids=d.gene_id.tolist()
    assert d.gene_id.is_unique
    for mod,path in [('scp','input/scp/discovery.tsv'),('liu','stats/liu_all.tsv.gz'),
                     ('chin','stats/chin_all.tsv.gz'),('whole','stats/whole_all.tsv.gz')]:
        s=read(w/path);assert s.gene_id.is_unique and set(s.gene_id)<=set(ids)
        d=d.merge(s.rename(columns={c:mod+'_'+c for c in s if c!='gene_id'}),on='gene_id',how='left',validate='1:1',sort=False)
        if mod=='scp':
            src=read(w/'input/pb/scp_raw_genes.tsv')
            d['scp_source_feature_present']=d.gene_id.isin(src.gene_id)
            d['scp_passes_expression_filter']=d.gene_id.isin(s.gene_id)
        else:
            for c in ['source_feature_present','passes_expression_filter']:d[mod+'_'+c]=boolean(d[mod+'_'+c])
        d[mod+'_availability']=np.select([~d[mod+'_source_feature_present'],~d[mod+'_passes_expression_filter']],
                                        ['absent_from_analysis_axis','below_expression_filter'],default='evaluable')
    assert d.gene_id.tolist()==ids
    cov=read(w/'coverage/ALL_COVERAGE.tsv.gz')
    assert not cov.duplicated(['dataset','disease','gene_id']).any()
    columns=['donors','cells','native_cells','expr_fraction','token_fraction','donors_native',
             'donors_token_ge10pct','median_donor_token_fraction']
    prefixes={'SCP1303':'scp','liu2023':'liu','chin_GSE174691':'chin_hcm','chin_GSE161921':'chin_control',
              'GSE138262_whole_cell':'whole','GSE255296':'gse255296','GSE292754':'gse292754'}
    for (ds,group),s in cov.groupby(['dataset','disease'],sort=True):
        assert set(s.gene_id)==set(ids)
        s=s.set_index('gene_id').loc[ids].reset_index()
        p='cov_'+prefixes[ds]+'_'+group+'_'
        d=d.merge(s[['gene_id']+columns].rename(columns={c:p+c for c in columns}),on='gene_id',validate='1:1',sort=False)
    source=cov[(cov.dataset=='SCP1303')&(cov.disease=='HCM')].set_index('gene_id').loc[ids]
    d['in_v1_vocab']=boolean(source.in_v1_vocab).to_numpy()
    d['token_id']=pd.array(source.token_id.to_numpy(),dtype='Int64')
    for ds,s in cov.groupby(['dataset','disease']):
        s=s.set_index('gene_id').loc[ids]
        assert np.array_equal(boolean(s.in_v1_vocab).to_numpy(),d.in_v1_vocab)
        assert np.allclose(s.token_id.to_numpy(dtype=float),d.token_id.to_numpy(dtype=float),equal_nan=True)
    d['scp_deg']=d.scp_FDR.lt(.05)&d.scp_logFC.abs().ge(.5)
    d['liu_deg']=d.liu_FDR.lt(.05)&d.liu_logFC.abs().ge(.5)
    d['comparison_status']=direction(d.scp_logFC,d.liu_logFC)
    d['same_direction']=d.comparison_status.eq('same')
    d['shared_deg']=d.scp_deg&d.liu_deg&d.same_direction
    d['scp_stable']=(d.scp_cb_logFC*d.scp_logFC).gt(0)&d.scp_loo_same_direction_fraction.eq(1)&d.scp_loo_min_abs_logFC.ge(.25)
    omitted=['liu_logFC_omit_Liu2023_CL0412','liu_logFC_omit_Liu2023_CL2']
    d['liu_stable']=(d[omitted].to_numpy()*d.liu_logFC.to_numpy()[:,None]>0).all(axis=1)
    a=d.shared_deg&d.scp_stable&d.liu_stable
    b=d.same_direction&(d.scp_deg^d.liu_deg)&d.scp_stable&d.liu_stable
    d['tier']=np.select([a,b],['A_shared_strong','B_single_study_directional'],default='none')
    d['main_support']=np.select([d.scp_deg&d.liu_deg,d.scp_deg,d.liu_deg],['both','SCP1303','Liu2023'],default='neither')
    d['scp_hcm_donors_token_ge10pct']=d.cov_scp_HCM_donors_token_ge10pct
    d['liu_hcm_donors_token_ge10pct']=d.cov_liu_HCM_donors_token_ge10pct
    d['token_eligible']=d.in_v1_vocab&d.token_id.notna()&d.scp_hcm_donors_token_ge10pct.ge(8)&d.liu_hcm_donors_token_ge10pct.ge(5)
    d['eligible_A']=a&d.token_eligible;d['eligible_B']=b&d.token_eligible
    d['evidence_category']=np.select([a,b,d.comparison_status.eq('not_evaluable'),d.comparison_status.eq('zero'),
        d.comparison_status.eq('opposite'),d.shared_deg,d.scp_deg^d.liu_deg],
        ['A_shared_strong','B_single_study_directional','not_evaluable','zero','opposite',
         'shared_direction_sensitive','single_study_direction_sensitive'],default='same_neither_DEG')
    d['rank_score_worst_FDR']=np.maximum(d.scp_FDR,d.liu_FDR)
    d['rank_score_min_abs_logFC']=np.minimum(d.scp_logFC.abs(),d.liu_logFC.abs())
    d['association_direction']=np.where(d.scp_logFC>0,'higher_in_HCM',np.where(d.scp_logFC<0,'lower_in_HCM','not_evaluable'))
    for mod in ['chin','whole']:
        d[mod+'_direction']=direction(d[mod+'_descriptive_log2_cpm_difference'],d.scp_logFC)
    d['descriptive_same_count']=d[['chin_direction','whole_direction']].eq('same').sum(axis=1)
    d['descriptive_opposite_count']=d[['chin_direction','whole_direction']].eq('opposite').sum(axis=1)
    d['descriptive_evaluable_count']=d[['chin_direction','whole_direction']].isin(['same','opposite','zero']).sum(axis=1)
    for ds,cut in [('gse255296',2),('gse292754',1)]:
        d[ds+'_half_donors_token_covered']=d['cov_'+ds+'_HCM_donors_token_ge10pct'].ge(cut)
    d['first_round_rank']=d.gene_id.map(dict(zip(previous.gene_id,previous['rank']))).astype('Int64')
    d['selection_notes']=[';'.join(note for condition,note in [
        (r.scp_availability!='evaluable','scp_'+r.scp_availability),(r.liu_availability!='evaluable','liu_'+r.liu_availability),
        (r.comparison_status=='opposite','opposite_direction'),(not r.scp_deg,'scp_DEG_threshold_not_met'),
        (not r.liu_deg,'liu_DEG_threshold_not_met'),(not r.scp_stable,'scp_stability_not_met'),
        (not r.liu_stable,'liu_stability_not_met'),(not r.token_eligible,'V1_coverage_not_met')] if condition) or 'all_A_requirements_met'
        for r in d.itertuples()]
    d['interpretation']='exploratory_HCM_association;beneficial_perturbation_undetermined'
    tables={'A_BIOLOGICAL':ranked(d[a]),'B_BIOLOGICAL':ranked(d[b]),
            'A_SIMULATABLE':ranked(d[d.eligible_A]),'B_SIMULATABLE':ranked(d[d.eligible_B])}
    tables['PRIORITY_TOP200']=tables['A_SIMULATABLE'].head(200).copy()
    d['A_simulatable_rank']=d.gene_id.map(dict(zip(tables['A_SIMULATABLE'].gene_id,tables['A_SIMULATABLE'].candidate_rank))).astype('Int64')
    d['B_simulatable_rank']=d.gene_id.map(dict(zip(tables['B_SIMULATABLE'].gene_id,tables['B_SIMULATABLE'].candidate_rank))).astype('Int64')
    d['priority_top200']=d.gene_id.isin(tables['PRIORITY_TOP200'].gene_id)
    old=previous[['rank','gene_id']].rename(columns={'rank':'old_rank'}).merge(d,on='gene_id',validate='1:1',sort=False)
    old['new_status']=np.select([old.eligible_A,old.tier.eq('A_shared_strong'),old.eligible_B,old.tier.eq('B_single_study_directional')],
                               ['A_simulatable','A_low_coverage','B_simulatable','B_low_coverage'],default='not_selected')
    assert old.gene_id.tolist()==previous.gene_id.tolist()
    out.mkdir()
    save(d,out/'ALL_GENE_EVIDENCE.tsv.gz')
    for name,t in tables.items():save(t,out/(name+'.tsv'))
    save(old,out/'PREVIOUS_200_COMPARISON.tsv')
    summary={'generated_utc':datetime.now(timezone.utc).isoformat(),'design_commit':'aefa26c','design_sha256':prep['design_sha256'],
        'gene_universe':len(d),'studies':{},'both_evaluable':int((d.scp_passes_expression_filter&d.liu_passes_expression_filter).sum()),
        'comparison_counts':d.comparison_status.value_counts().to_dict(),'shared_DEG':int(d.shared_deg.sum()),
        'shared_stability_failed':int((d.shared_deg&~a).sum()),'evidence_category_counts':d.evidence_category.value_counts().to_dict(),
        'table_rows':{k:len(v) for k,v in tables.items()},'prior_200_status':old.new_status.value_counts().to_dict(),
        'prior_200_sha256':sha(w/'input/previous/CANDIDATES.tsv'),'candidate_groups':{},
        'held_out_validation':False,'selection_scope':'full_gene_union_not_previous_200','joint_FDR_claim':False}
    for mod in ['scp','liu','chin','whole']:
        summary['studies'][mod]={'source_genes':int(d[mod+'_source_feature_present'].sum()),'expression_filtered':int(d[mod+'_passes_expression_filter'].sum())}
        if mod in ['scp','liu']:summary['studies'][mod]['DEG']=int(d[mod+'_deg'].sum())
    for name,t in tables.items():
        summary['candidate_groups'][name]={'higher_HCM':int(t.scp_logFC.gt(0).sum()),'lower_HCM':int(t.scp_logFC.lt(0).sum()),
            'in_previous_200':int(t.first_round_rank.notna().sum()),'new_vs_previous_200':int(t.first_round_rank.isna().sum()),
            'main_support':t.main_support.value_counts().to_dict(),
            'chin_directions':t.chin_direction.value_counts().to_dict(),'whole_directions':t.whole_direction.value_counts().to_dict(),
            'gse255296_half_donors_covered':int(t.gse255296_half_donors_token_covered.sum()),
            'gse292754_half_donors_covered':int(t.gse292754_half_donors_token_covered.sum())}
    (out/'SUMMARY.json').write_text(json.dumps(summary,indent=2,allow_nan=False)+'\n')
    (out/'OUTPUT_HASHES.json').write_text(json.dumps({p.name:sha(p) for p in sorted(out.iterdir()) if p.is_file()},indent=2)+'\n')
    print(json.dumps({k:summary[k] for k in ['gene_universe','studies','both_evaluable','shared_DEG','shared_stability_failed','table_rows','prior_200_status']},allow_nan=False))

if __name__=='__main__':main()
