"""Readable views preserve ranked candidates and separate disease evidence from coverage."""
import json
import sys
from pathlib import Path
import numpy as np
import pandas as pd

w=Path(sys.argv[1]);out=w/'review';out.mkdir(exist_ok=True)
def read(path):return pd.read_csv(path,sep='\t',keep_default_na=False,na_values=['NA'])
labels={'same':'同向','opposite':'反向','not_evaluable':'无法比较','zero':'零变化',
        'higher_in_HCM':'HCM较高','lower_in_HCM':'HCM较低','both':'SCP和Liu都达标',
        'SCP1303':'仅SCP达到差异门槛','Liu2023':'仅Liu达到差异门槛'}
columns={'candidate_rank':'本层顺序','gene_id':'基因编号','gene_symbol':'基因名称',
         'association_direction':'HCM表达方向','main_support':'达到差异门槛的研究',
         'scp_logFC':'SCP变化幅度_log2','scp_FDR':'SCP统计值_多重校正',
         'liu_logFC':'Liu变化幅度_log2','liu_FDR':'Liu统计值_多重校正',
         'scp_loo_min_abs_logFC':'SCP逐人剔除后最小绝对变化',
         'liu_logFC_omit_Liu2023_CL0412':'Liu仅保留CL2对照时变化',
         'liu_logFC_omit_Liu2023_CL2':'Liu仅保留CL0412对照时变化',
         'scp_hcm_donors_token_ge10pct':'SCP模型可见供体数_共15人',
         'liu_hcm_donors_token_ge10pct':'Liu模型可见供体数_共10人',
         'chin_direction':'Chin补充方向','whole_direction':'完整细胞补充方向',
         'cov_gse255296_HCM_donors_token_ge10pct':'GSE255296模型可见供体数_共4人',
         'cov_gse292754_HCM_donors_token_ge10pct':'儿童GSE292754模型可见供体数_共2人',
         'first_round_rank':'旧200名单中的序号_空白为新增'}
tables={}
for name,source in [('COMMON_CANDIDATES_REVIEW','A_SIMULATABLE'),('EXPLORATORY_CANDIDATES_REVIEW','B_SIMULATABLE')]:
    d=read(w/f'final/{source}.tsv');tables[source]=d
    simple=d[list(columns)].copy()
    for c in ['association_direction','main_support','chin_direction','whole_direction']:simple[c]=simple[c].map(labels)
    for c in ['candidate_rank','first_round_rank','scp_hcm_donors_token_ge10pct','liu_hcm_donors_token_ge10pct',
              'cov_gse255296_HCM_donors_token_ge10pct','cov_gse292754_HCM_donors_token_ge10pct']:simple[c]=simple[c].astype('Int64')
    simple.rename(columns=columns).to_csv(out/(name+'.csv'),index=False,encoding='utf-8-sig',na_rep='')
    assert simple.gene_id.tolist()==d.gene_id.tolist()
old=read(w/'final/PREVIOUS_200_COMPARISON.tsv')
status={'A_simulatable':'进入本轮共同候选','A_low_coverage':'共同证据通过但模型覆盖不足',
        'B_simulatable':'进入本轮探索候选','B_low_coverage':'探索证据通过但模型覆盖不足','not_selected':'未通过本轮分层要求'}
v=old[['old_rank','gene_id','gene_symbol','new_status','selection_notes']].copy()
v.new_status=v.new_status.map(status)
reason={'scp_DEG_threshold_not_met':'SCP未达到差异门槛','liu_DEG_threshold_not_met':'Liu未达到差异门槛',
        'opposite_direction':'两研究方向相反','scp_stability_not_met':'SCP稳定性未通过','liu_stability_not_met':'Liu对照敏感',
        'V1_coverage_not_met':'模型覆盖不足','all_A_requirements_met':'共同候选要求均通过',
        'scp_absent_from_analysis_axis':'SCP没有可比较特征','liu_absent_from_analysis_axis':'Liu没有可比较特征',
        'scp_below_expression_filter':'SCP信号不足','liu_below_expression_filter':'Liu信号不足'}
v.selection_notes=v.selection_notes.map(lambda s:'；'.join(reason.get(x,x) for x in s.split(';')))
v.rename(columns={'old_rank':'旧序号','gene_id':'基因编号','gene_symbol':'基因名称','new_status':'本轮去向','selection_notes':'门槛说明'}).to_csv(out/'PREVIOUS_200_REVIEW.csv',index=False,encoding='utf-8-sig')
patterns=[]
for name,t in tables.items():
    p=t.groupby(['chin_direction','whole_direction','gse255296_half_donors_token_covered','gse292754_half_donors_token_covered'],dropna=False).size().rename('genes').reset_index()
    p.insert(0,'candidate_table',name);patterns.append(p)
pd.concat(patterns).to_csv(out/'SUPPORT_PATTERNS.tsv',sep='\t',index=False)
full=pd.concat(tables.values(),ignore_index=True)
assert full.gene_id.is_unique
records=[]
study_roles=[('SCP1303','scp','joint_discovery'),('liu2023','liu','joint_discovery'),
             ('chin_GSE174691','chin_hcm','paired_description'),('chin_GSE161921','chin_control','paired_description'),
             ('GSE138262_whole_cell','whole','confounded_description'),
             ('GSE255296','gse255296','HCM_only_coverage'),('GSE292754','gse292754','pediatric_HCM_only_coverage')]
for r in full.itertuples():
    for dataset,prefix,role in study_roles:
        group='Control' if dataset=='chin_GSE161921' else 'HCM'
        cov='cov_'+prefix+'_'+group+'_'
        direction='not_applicable_no_control' if 'coverage' in role else 'same' if role=='joint_discovery' else r.chin_direction if role=='paired_description' else r.whole_direction
        records.append({'tier':r.tier,'candidate_rank':r.candidate_rank,'gene_id':r.gene_id,'gene_symbol':r.gene_symbol,
                        'dataset':dataset,'role':role,'direction_relative_to_SCP':direction,
                        'comparison_unit':'Chin_pair' if role=='paired_description' else dataset,
                        'coverage_group':group,'donors':getattr(r,cov+'donors'),'native_cells':getattr(r,cov+'native_cells'),
                        'expression_fraction_native':getattr(r,cov+'expr_fraction'),'token_fraction':getattr(r,cov+'token_fraction'),
                        'donors_token_ge10pct':getattr(r,cov+'donors_token_ge10pct')})
pd.DataFrame(records).to_csv(out/'SEVEN_DATASET_EVIDENCE.tsv',sep='\t',index=False,na_rep='NA',float_format='%.15g')
summary={}
for name,t in tables.items():
    summary[name]={'rows':len(t),'both_descriptive_same':int((t.chin_direction.eq('same')&t.whole_direction.eq('same')).sum()),
                  'at_least_one_descriptive_opposite':int((t.chin_direction.eq('opposite')|t.whole_direction.eq('opposite')).sum()),
                  'both_HCM_only_half_donors_covered':int((t.gse255296_half_donors_token_covered&t.gse292754_half_donors_token_covered).sum())}
(out/'REVIEW_SUMMARY.json').write_text(json.dumps(summary,indent=2)+'\n')
print(json.dumps(summary))
