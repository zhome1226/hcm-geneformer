"""Publish scope-adjusted text only after both final audits pass."""
import csv
import io
import json
import shutil
from datetime import datetime, timezone
from pathlib import Path

HERE=Path(__file__).resolve().parent
LOCAL=HERE.parent
REPO=Path(r'D:\VirtualCell\xiaoe-internship')
DEST=REPO/'results/qc/all_data_cleaning_20261111'
integrity=json.loads((HERE/'FINAL_INTEGRITY_AUDIT.json').read_text())
completion=json.loads((HERE/'COMPLETION_AUDIT.json').read_text())
assert integrity['status']=='passed' and integrity['fresh_checks_failed']==0
assert completion['status']=='passed' and completion['failed']==0
assert integrity['datasets_in_scope']==13
now=datetime.now(timezone.utc).isoformat()

f=LOCAL/'CLEANING_STATUS.tsv'
rows=list(csv.DictReader(f.open(encoding='utf-8-sig'),delimiter='\t'))
for row in rows:
    if row['dataset']=='cox7b_GSE227734':
        row['state']='excluded_by_user_missing_gene_annotation'
assert sum(r['state']=='complete_verified' for r in rows)==13
assert sum(r['state'].startswith('excluded_by_user') for r in rows)==1
with f.open('w',encoding='utf-8',newline='') as out:
    w=csv.DictWriter(out,fieldnames=list(rows[0]),delimiter='\t',quoting=csv.QUOTE_ALL,lineterminator='\n')
    w.writeheader();w.writerows(rows)

f=LOCAL/'PROTOCOL.md'
t=f.read_text(encoding='utf-8-sig')
t=t.replace('目标：完成下载清单中的14套生物/影像数据清洗，',
'''目标（按用户最新决定调整）：下载清单共14套，GSE227734因缺少权威基因映射由用户排除、不再清洗；完成其余13套生物/影像数据清洗，''')
t=t.replace('Liu2023、GSE227734、GSE138262','Liu2023、GSE138262')
t=t.replace('全部完成后逐项核验14套状态，再判定总目标完成。',
'''最终逐项核验13套在范围内数据，并单独记录1套用户排除；排除项未清洗、不可训练，不将它计为完成。2026-11-12最终复核见FINAL_REVIEW_20261112.md。''')
f.write_text(t,encoding='utf-8')

f=LOCAL/'GSE227734_REPORT.md'
t=f.read_text(encoding='utf-8-sig')
title,_,rest=t.partition('\n')
f.write_text('# GSE227734：用户排除，不再清洗\n\n'
             '2026-11-12 UTC：用户决定将此数据集排除。本套未完成分子质控、不可用于训练；原件及历史准备/来源审计保留。作者邮件未发送，发送授权已撤回，不继续联系作者或追查映射。以下为排除前的历史审计。\n\n'+rest.lstrip(),encoding='utf-8')
f=LOCAL/'gse227734/AUTHOR_CONTACT_STATUS.json'
j=json.loads(f.read_text(encoding='utf-8-sig'))
j.update({'updated_at':now,'timezone':'UTC','status':'not_sent_dataset_excluded_by_user',
          'current_authorization_to_send':False,'message_sent':False,
          'dataset_scope':'excluded_by_user; do not continue cleaning or author contact'})
f.write_text(json.dumps(j,ensure_ascii=False,indent=2)+'\n',encoding='utf-8')

report=f'''# 最终清洗检查：13套完成，1套用户排除

检查时间：{integrity['checked_utc']}。GSE227734按用户决定不再清洗；其余13套在约定模态与规则范围内完成清洗/整理及最终验收。邮件未发送，未联系作者。

服务器：`zhome@10.158.203.133`，数据根目录`/home/zhome/xiaoe-internship-data/`，最终产物在`processed/`下；逐套路径见[CLEANING_STATUS.tsv](CLEANING_STATUS.tsv)。

本轮重新读取{integrity['unique_files_rehashed']}个真实文件，共{integrity['total_bytes_rehashed']:,}字节（{integrity['total_bytes_rehashed']/1024**3:.3f} GiB）；与既有字节数及MD5/SHA256基线全部一致（基线仅有SHA256的文件按SHA256核验）。核验了HDF5结构和索引目标。主审计{integrity['fresh_checks_passed']}项通过、0项失败；终态退出/完成记录与ACDC几何补充审计{completion['passed']}项通过、0项失败。13份既有独立验收记录逐字节一致，确认{integrity['historical_independent_checks_confirmed']:,}项检查全部通过；另有ACDC体素复核2,100项通过。

本轮只作完整性及交付复查，没有重跑成功清洗、改变计数、训练模型或生成占位文件。历史失败日志和旧范围报告保留；当前范围与结论以本报告和最新状态表为准。

| 数据集 | 最终状态 | 使用限制 |
| --- | --- | --- |
| SCP1303 | 完成，reviewed_final版本 | 核心114711核、轻量20000；跨研究供体重叠及最终训练划分待定 |
| Chin GSE174691 | 完成 | HCM核心10150；报告待定核和逐供体损失 |
| Chin GSE161921 | 完成 | Control核心3534；与HCM来源分别保存 |
| Liu2023 | 完成 | 核心8547；CL0412A/B绑定同一供体 |
| GSE181764 | 完成 | 对照核心仅5核，不能视为充分二分类验证 |
| GSE138262 | 完成严格版 | 仅6个HCM核心、0对照，本版本不批准训练 |
| GSE203274 | 完成 | 儿童队列仅1位HCM供体，不证明跨供体泛化 |
| bulk GSE249925 | 完成发布矩阵整理/审计 | raw/TMM来源冲突保留，不认证未归一化计数，不用于单细胞训练 |
| VCC2025 | 完成发布数据验证和真实索引 | 扰动比赛，非HCM队列，保留官方划分 |
| NeurIPS2023 | 完成发布数据验证和真实索引 | 多组学/DE分模态，不赋HCM标签 |
| NeurIPS2021 | 完成发布数据验证和配对索引 | 骨髓数据，非HCM队列，ATAC派生异常坐标另记 |
| CausalBench2023 | 完成发布计数验证和扰动索引 | K562/RPE1，非HCM队列，未创建训练划分 |
| ACDC2017 | 完成影像审计和可用索引 | 300对可做体素分割；72对空间头冲突需另审 |
| GSE227734 | 用户排除；未完成清洗 | 缺少9760列权威基因映射；全部不可训练；原件保留 |

500≤检测基因≤5000、MT<5%是约定的HCM单细胞/单核严格基线，不是行业统一标准，也未套用到bulk或各类比赛。心肌身份和疑似双细胞按各研究已验收的多证据规则处理。清洗完成不表示13套都适合HCM微调；训练候选、供体重叠及训练/验证划分需下一阶段独立决定。

证据：[完整性审计](final_audit_20261112/FINAL_INTEGRITY_AUDIT.json)、[完成记录复核](final_audit_20261112/COMPLETION_AUDIT.json)、[逐文件复核记录](final_audit_20261112/PROGRESS.jsonl)。服务器审计执行目录：`_processing/final_cleaning_audit_20261112/`。GitHub仅发布文本和脚本。
'''
(LOCAL/'FINAL_REVIEW_20261112.md').write_text(report,encoding='utf-8')

f=DEST/'README.md'
t=f.read_text(encoding='utf-8-sig')
old='状态日期：2026-11-11 UTC。14套中13套完成各自模态的整理/质控验收，**GSE227734仍缺权威基因映射，全部数据清洗尚未完成**。逐套状态以[CLEANING_STATUS.tsv](CLEANING_STATUS.tsv)为准；来源准备审计通过不等于完成分子质控。'
new='状态日期：2026-11-12 UTC。**按用户最新决定，13套在范围内数据已完成各自模态的清洗/整理及最终验收；GSE227734排除，不再清洗。** 它仍缺权威基因映射、不可训练，不计为完成。最终复核见[FINAL_REVIEW_20261112.md](FINAL_REVIEW_20261112.md)，逐套路径见[CLEANING_STATUS.tsv](CLEANING_STATUS.tsv)。'
assert old in t
t=t.replace(old,new)
t=t.replace('**未完成：9760列匿名，无基因映射**','**用户排除：未完成清洗，不可训练**')
t=t.replace('权威顺序映射或可逐值核对的修正版仍是下一步必要输入。','用户已决定排除此套，不再追查或清洗。作者邮件未发送，授权已撤回。')
t+='\n2026-11-12最终审计重新核对真实产物、原始矩阵与索引目标；详见最终报告。`GIT_DELIVERY_AUDIT.json`及`GITHUB_DELIVERY_VERIFICATION.json`保留2026-11-11旧范围的历史交付证据，不代表当前排除后的范围。\n'
f.write_text(t,encoding='utf-8')
for name in ('PROTOCOL.md','CLEANING_STATUS.tsv','GSE227734_REPORT.md','FINAL_REVIEW_20261112.md','gse227734/AUTHOR_CONTACT_STATUS.json'):
    shutil.copyfile(LOCAL/name,DEST/name)
shutil.copytree(HERE,DEST/HERE.name,dirs_exist_ok=True)
print('Final report, scope and text delivery updated after passing both audits.')
