"""Publish verified SCP1303 QC results and explicit selection-loss accounting."""
from pathlib import Path
import json
import pandas as pd

root = Path.home() / 'xiaoe-internship-data'
out = root / 'processed/SCP1303_qc_reviewed_final_20261111'
verification = json.loads((out/'INDEPENDENT_VERIFICATION.json').read_text())
summary = json.loads((out/'EXPORT_SUCCESS.json').read_text())
assert verification['status']=='passed' and all(verification['checks'].values())
d = pd.read_csv(out/'review_decisions.tsv.gz',sep='\t',index_col='barcode',low_memory=False)
rows = []
for (disease, donor, subtype), part in d.groupby(['disease','donor_id','cell_type_original'],observed=True):
    baseline=part['eligible_strict_5000']
    rows.append({'disease':disease,'donor_id':donor,'subtype':subtype,
        'original_strict':int(baseline.sum()),'reviewed_core':int(part['kept_reviewed_core'].sum()),
        'pending':int((baseline & part['action'].eq('pending')).sum()),
        'suspected_doublet':int((baseline & part['action'].eq('exclude_suspected_doublet')).sum()),
        'noncm':int((baseline & part['action'].eq('exclude_noncm')).sum()),
        'above_5000':int((~baseline).sum())})
loss=pd.DataFrame(rows)
loss.to_csv(out/'SELECTION_LOSS.tsv',sep='\t',index=False)
priority=pd.read_csv(root/'processed/SCP1303_review_20261111/priority_26.tsv',sep='\t',index_col='barcode')
d.loc[priority.index].to_csv(out/'priority_26_decisions.tsv',sep='\t',index_label='barcode')
files='\n'.join(f"| {name} | {entry['bytes']} | {entry['md5']} |" for name,entry in verification['hashes'].items())
hcm_iii=loss.loc[loss['disease'].eq('HCM') & loss['subtype'].eq('Cardiomyocyte_III')]
report=f'''# SCP1303 清洗复核完成

验收时间：{verification['verified_utc']}。全部{len(verification['checks'])}项独立检查通过；没有训练、tokenization或疾病差异分析。

## 最终使用路径

`{out}/`

| 对象 | 核数 | 说明 |
| --- | ---: | --- |
| cm_qc_full.h5ad | {verification['core_nuclei']} | HCM 54344，NF 60367；31供体，36601基因 |
| cm_qc_light.h5ad | {verification['light_nuclei']} | HCM/NF各10000，31供体；44个原轻量核按供体/亚型补抽替换 |
| cm_review_pending.h5ad | {verification['pending_nuclei']} | 证据不足或冲突者单独保留，不进入核心版 |

| 文件 | 字节数 | 实测MD5 |
| --- | ---: | --- |
{files}

原严格基线114927核中，保留114711核；67核符合多证据疑似双细胞规则，5核符合非心肌规则，144核待定。疑似双细胞不是实验确证；供体疾病标签也不是逐核病理金标准。

## 规则及证据

- 主规则保持500≤CellBender检测基因≤5000，CellRanger计数重算MT<5%；NF=0、HCM=1，DCM不进入心肌二分类核心。
- 独立成人左心室单核参考使用SingleR；4个参考供体留出。单核级原始CM预测precision约0.899，因此不单独按参考标签删核，也不把参考分数称为概率。
- 为集中于CM_III的稀疏核补充整体表达程序：只用参考训练供体提取最多100个CM特异基因，并检查Leiden群体及至少3位供体的群体一致性。参考留出30个类型×供体聚合谱全部预测正确，包含4个CM谱；这只是参考适配检查，不代表HCM分类性能。
- 补充scDblFinder涵盖43个作者保留的完整10x文库、多种细胞类型，输入为CellRanger计数。残余率假设dbr=.01、sd=.015。只有工具调用、作者样本内高分和混合谱系等共同支持时排除；不声称恢复最初液滴检测。
- 其余无新信号核复用作者QC，逐核状态明确，不声称全部接受过独立参考确认。所有987个完整警示、185个轻量警示及5043个高基因数候选均有记录。
- 高于5000者保持项目上限排除，即使其心肌身份可信，也不记作双细胞。5000/6000覆盖比较见upper_limit_sensitivity.tsv；未悄悄放宽主规则。

## 群体与供体覆盖

HCM CM_III原严格基线{int(hcm_iii['original_strict'].sum())}核，复核核心{int(hcm_iii['reviewed_core'].sum())}核，待定{int(hcm_iii['pending'].sum())}核，疑似双细胞{int(hcm_iii['suspected_doublet'].sum())}核，非心肌{int(hcm_iii['noncm'].sum())}核。没有通过删除整个疾病状态群体来取得完成结果。

逐疾病、供体、作者亚群的损失见SELECTION_LOSS.tsv；26个重点核见priority_26_decisions.tsv。31位供体全部保留；现有供体split仍是暂定，不代表最终微调方案。跨研究供体独立性仍未证实。

## 验收与保留

重新打开三个导出对象，对全部保留核/基因的三套计数与原严格对象逐值比较；轻量版对完整对象也逐值比较。重算QC、核对条码/诊断/供体划分/完整基因空间，并独立核对输出大小与摘要。两个原严格H5AD的MD5保持不变。

完整未归一化计数及原CellRanger层保留，没有插补、标准化、整合或HVG裁切交付矩阵。规则、逐核记录、校验和、补抽记录均在本目录。原严格基线保留；`SCP1303_qc_reviewed_20261111/`为早期复核草稿，最终使用本目录。
'''
(out/'REPORT.md').write_text(report)
draft=root/'processed/SCP1303_qc_reviewed_20261111'
(draft/'DRAFT_NOT_FOR_TRAINING.md').write_text('早期复核草稿，尚未纳入群体稀疏性复核。最终使用 '+str(out)+'/；本目录不是最终交付。\n')
print(json.dumps({'status':'verified_complete','priority_26':d.loc[priority.index,'action'].value_counts().to_dict(),
    'HCM_CM_III':hcm_iii.sum(numeric_only=True).to_dict()},indent=2))
