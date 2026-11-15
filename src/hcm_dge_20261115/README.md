# 2026-11-15 HCM 供体级 DGE

[设计](../../results/target_discovery/hcm_dge_20261115/DESIGN.md)和[结果报告](../../results/target_discovery/hcm_dge_20261115/REPORT.md)是本批次入口。脚本只写新的分析工作区，不重写输入 h5ad 或 dataset。Python 输入路径、样本身份/别名和冻结 SHA 在 pseudobulk.py、coverage.py 中显式限定；不能未经核对套用到别的对象。

## 已运行的阶段

|脚本|职责|
|---|---|
|pseudobulk.py|roster 冻结研究内供体；discovery 求 SCP raw、CellBender、HHATL 汇总；support 在锁定后生成 Liu、Chin、完整细胞汇总|
|coverage.py|直接读取原生计数与实际 V1 序列；核对细胞/供体一致，区分未测量与零表达|
|fit.R|edgeR TMM + robust QL；SCP 年龄/性别校正、31 次删一与两项敏感性；Liu 疾病主效应和两次单对照；两个混杂队列只给描述性值|
|lock_candidates.py|按预先设计排序，最多锁定 200；生成 SHA、锁定时间、漏斗，并把冻结名单设为只读|
|integrate.py|联接锁定后的证据、不改序号；分开保留同向、对照敏感、反向和无法评估|
|audit.py|从计数/覆盖明细独立重算过滤、SCP BH、删一指标、排名及表间对应；不重新拟合模型|
|export_evidence.py|验收通过后导出 80 项分析证据，排除原始细胞及全量供体计数矩阵|

环境：Linux R 4.5.3、edgeR 4.8.2；详见结果目录 session_discovery.txt、session_validation.txt。Python分别使用既有 QC 与 Geneformer V1 虚拟环境；不是本次重装。V1 2048-token 序列沿用已验收对象，不重新 tokenize。

## 重现顺序

下面展示本批次的逻辑命令。实际通过 nohup 式独立进程执行，各自有 logs/*.log 和 *.exit；SCP 31 次删一最多四个进程，各进程 BLAS 限一线程。完整日志在 Linux 和本机 evidence/ 中。

在**新的空工作目录**运行；不要覆盖当前冻结结果或复用已有锁文件。code/DESIGN.md 必须是提交 5019f70 的预先方案。

```bash
QC_PY=/home/zhome/xiaoe-internship-data/_processing/SCP1303_20261111/venv/bin/python
V1_PY=/home/zhome/xiaoe-internship-data/_processing/geneformer_v1_tokenization_20261112/venv/bin/python
R_RUN=/home/zhome/xiaoe-internship-data/_processing/all_data_cleaning_20261111/r_env/bin/Rscript
WORK=/path/to/new/empty/analysis_workspace
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 HCM_LODO_WORKERS=4

"$QC_PY" code/pseudobulk.py "$WORK" roster
"$QC_PY" code/pseudobulk.py "$WORK" discovery
"$V1_PY" code/coverage.py --work "$WORK" --phase discovery
"$R_RUN" code/fit.R "$WORK" discovery
# 确认前述成功，并按实际退出状态保存 logs/fit_d.exit；不可伪造成功标记。
"$QC_PY" code/lock_candidates.py "$WORK"

# 只有锁定成功且候选非空，才运行复核：
"$QC_PY" code/pseudobulk.py "$WORK" support
"$V1_PY" code/coverage.py --work "$WORK" --phase support
"$R_RUN" code/fit.R "$WORK" validation
# 为各阶段记录真实的 logs/{fit_d,lock,pb_s,coverage_s,fit_s}.exit。
"$QC_PY" code/integrate.py "$WORK"
"$QC_PY" code/audit.py "$WORK" --final final/CANDIDATE_EVIDENCE.tsv
# integrate、audit 的真实退出状态也须记录，方可导出。
"$QC_PY" code/export_evidence.py "$WORK"
```

audit.py 中 DESIGN_SHA 固定原始方案（LF规范化后），同时检查原始锁定输入的逐文件 SHA。复核不会基于结果回头修改候选；任何规则改动须新建分析版本，不称为本批次重现。

## 最终表怎么读

- CANDIDATES.tsv：完整冻结发现名单，SHA 与后续复核一致。
- CANDIDATE_EVIDENCE.tsv：200 行详细证据；所有 scp_* 列来自冻结发现结果，liu_* 为锁定后结果，chin_* / whole_* 为描述性对照。
- CANDIDATES_REVIEW.csv：UTF-8 BOM 中文简表，适合 Excel；没有改候选顺序。
- scp_logFC / liu_logFC：正值 HCM 较高，负值 HCM 较低。不是干预获益方向。
- scp_loo_min_abs_logFC：31 次删一效应绝对值最小值，不是置信区间下界。
- liu_logFC_omit_Liu2023_CL0412：剔除 CL0412，**只剩 CL2 对照**；另一列相反。
- liu_evidence_category：same_both_single_controls、same_control_sensitive、opposite、not_evaluable 分别表示“主分析和两次单对照都同向”“主同向但对照敏感”“主反向”“不可评估”。
- *_source_feature_present 与 *_passes_expression_filter：分别标记原生/共同可测特征和模型表达过滤；NA 不等于零。
- *_expr_fraction：计数检测率，分母为原生测量细胞；*_token_fraction：实际序列出现率，分母为该组全部核/细胞。二者不能混用。
- *_donors_token_ge10pct：其核/细胞中该 token 出现率至少 10% 的供体数；分母见 *_donors。不是研究间重复显著性。
- interpretation：统一保留“关联候选，干预获益方向未定”。

此阶段没有模型训练、模拟删除、疗效预测或 ToDesk 传输。
