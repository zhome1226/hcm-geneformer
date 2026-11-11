"""Report attrition against the published author cell types, before and after QC."""
from pathlib import Path
import hashlib
import json
import pandas as pd

root=Path.home()/'xiaoe-internship-data'
w=root/'_processing/all_data_cleaning_20261111/liu2023'
out=root/'processed/liu2023_qc_20261111'
verification=json.loads((out/'INDEPENDENT_VERIFICATION.json').read_text())
assert verification['status']=='passed'
meta=pd.read_csv(w/'author_metadata.tsv.gz',sep='\t',index_col='barcode',keep_default_na=False)
audit=pd.concat([pd.read_csv(out/f'{r.sample_id}_released_nuclei_audit.tsv.gz',sep='\t',index_col='barcode') for r in pd.read_csv(w/'samples.tsv',sep='\t').itertuples()])
assert len(audit)==55122 and audit.index.is_unique and audit.barcode_original.is_unique
original=meta.loc[audit.barcode_original]
assert set(audit.barcode_original)==set(meta.index)
assert (original.cellType.to_numpy()==audit.author_cellType.to_numpy()).all()
audit['author_type']=original.cellType.to_numpy()
audit['author_subcluster']=original.subcluster.to_numpy()
audit['core']=audit.final_action.eq('retain_CM')
audit['pending']=audit.final_action.eq('pending_CM')
audit['doublet']=audit.final_action.eq('exclude_suspected_doublet')
assert audit.core.sum()==verification['core_CM'] and audit.pending.sum()==verification['pending_CM']
table=audit.groupby(['disease','author_type','author_subcluster'],observed=True).agg(
    input_nuclei=('core','size'),strict_nuclei=('pass_strict_numeric','sum'),
    core_CM=('core','sum'),pending_CM=('pending','sum'),suspected_doublets=('doublet','sum')).reset_index()
table.to_csv(out/'author_QC_retention_by_subcluster.tsv',sep='\t',index=False)
cm=audit[audit.author_type.eq('CM')]
summary=cm.groupby('disease',observed=True).agg(author_CM_input=('core','size'),
    author_CM_strict=('pass_strict_numeric','sum'),author_CM_core=('core','sum'),
    author_CM_pending=('pending','sum'),author_CM_suspected_doublet=('doublet','sum')).reset_index()
summary['numeric_excluded']=summary.author_CM_input-summary.author_CM_strict
summary['core_fraction_of_author_input']=summary.author_CM_core/summary.author_CM_input
summary['core_fraction_of_author_strict']=summary.author_CM_core/summary.author_CM_strict
summary.to_csv(out/'author_CM_attrition_summary.tsv',sep='\t',index=False)
assert summary.author_CM_input.sum()==10221
assert summary.author_CM_strict.sum()==8927
assert summary.author_CM_core.sum()==8547
assert summary.author_CM_pending.sum()==380
assert summary.author_CM_suspected_doublet.sum()==0
result={'status':'passed_attrition_accounting','author_CM_input':int(len(cm)),
    'author_CM_strict':int(cm.pass_strict_numeric.sum()),'author_CM_core':int(cm.core.sum()),
    'author_CM_pending':int(cm.pending.sum()),'author_CM_suspected_doublet':int(cm.doublet.sum()),
    'pending_total':verification['pending_CM'],'pending_author_non_CM':int(audit.loc[audit.author_type.ne('CM'),'pending'].sum()),
    'disease_author_CM':summary.to_dict('records'),
    'script_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest()}
(out/'ATTRITION_SUMMARY.json').write_text(json.dumps(result,indent=2)+'\n')
lines=['\n## 作者心肌标签的质控损失\n',
       '作者原心肌10221核，经项目数值筛选剩8927核，联合身份核查后核心8547核、待定380核；原作者心肌没有被本轮强双细胞规则排除。待定总池7918核中，其余7538核的作者注释为非心肌，仅因参考/标记冲突进入广义待定；它们不能当确证心肌、也不能当7918个原心肌损失。\n',
       '| 来源组 | 作者原CM | 严格候选CM | 核心CM | 待定CM | 数值筛除CM |\n| --- | ---: | ---: | ---: | ---: | ---: |\n']
for r in summary.itertuples():
    lines.append(f'| {r.disease} | {r.author_CM_input} | {r.author_CM_strict} | {r.author_CM_core} | {r.author_CM_pending} | {r.numeric_excluded} |\n')
lines.append('\n逐作者亚群的输入→数值QC→最终归属见author_QC_retention_by_subcluster.tsv；该版本仍有数值筛选造成的选择偏差，后续不能把严格候选损失解释成真实疾病差异。\n')
report=out/'REPORT.md';text=report.read_text();marker='\n## 作者心肌标签的质控损失\n'
if marker in text:text=text.split(marker)[0]
report.write_text(text+''.join(lines))
print(json.dumps(result,indent=2))
