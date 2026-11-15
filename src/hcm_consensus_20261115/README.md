# 第二版：七套数据的跨研究共同候选

本批是扩展探索：第一版SCP全量结果与其200候选的Liu结果已经揭盲，Liu在本轮改为共同发现，不再作为留出验证。[执行前固定方案](../../results/target_discovery/hcm_consensus_20261115/DESIGN.md)，提交aefa26c。只读既有QC/V1数据，不训练或模拟，不操作ToDesk。

## 程序与顺序

1. `prepare.py WORK`：将上一轮通过验收的供体矩阵、全部SCP结果、31次删一、原200和元数据复制到新工作目录，逐文件核对源/副本/源前后摘要。建立37,944基因并集，不从旧200出发。
2. `Rscript fit_all.R WORK`：Liu全基因疾病模型与两次单对照敏感性，Chin/完整细胞全可测轴描述比较；分别保留缺测、未通过表达过滤和可评估结果。Liu全体受检基因BH校正后才选行，与原200的FC/FDR等逐列一致才导出。
3. `all_coverage.py --work WORK`：逐研究调用原已验收的`coverage_base.py`，核对七套来源h5ad、准备对象、实际.dataset SHA及逐细胞元数据；输出37,944基因×74供体的2,807,856行明细和10组汇总。SCP表达检测使用raw，token覆盖仍来自实际既有序列。
4. `integrate.py WORK`：从全基因联接结果分层，保留A/B生物候选，再独立筛选模型覆盖；原200只用于事后对比。
5. `review_exports.py WORK`：中文简表、新旧比较和七套证据明细；Chin两个数据源共享一个成对比较，不计作两项独立方向证据。
6. `audit.py WORK`：独立重算表达过滤、SCP/Liu全基因BH、31次删一、每供体覆盖与汇总、所有联接/候选门槛/排序及新旧名单。不重复拟合模型，不将敏感性范围称为置信区间。
7. `export_evidence.py WORK [successful_audit_tag]`：核对所有阶段真实退出状态及最终验收，导出证据；原始细胞和全量供体计数矩阵留在Linux。

`coverage_base.py`从上一轮Linux `code/coverage.py`逐字节复制并以SHA固定，本仓源码在`../hcm_dge_20261115/coverage.py`，不改原算法。复用的SCP分析源码为`../hcm_dge_20261115/fit.R`。全部实际输入、版本和执行脚本随证据包保存。

## 运行环境

- 新工作目录：`/home/zhome/hcm-geneformer-analysis/hcm_consensus_20261115/`。
- QC Python：`/home/zhome/xiaoe-internship-data/_processing/SCP1303_20261111/venv/bin/python`；用于prepare/integrate/review/audit/export。
- V1 Python：`/home/zhome/xiaoe-internship-data/_processing/geneformer_v1_tokenization_20261112/venv/bin/python`；用于coverage。
- Rscript：`/home/zhome/xiaoe-internship-data/_processing/all_data_cleaning_20261111/r_env/bin/Rscript`，R4.5.3/edgeR4.8.2。
- 通过独立后台进程运行，BLAS单线程；日志和真实退出码在`logs/`。重新运行须新建目录，不覆盖成功输出；失败记录保留。

合成数据测试`test_fit_all.R`核对新R脚本与旧算法的一致性、两次单对照、NA及篡改旧结果时拒绝导出。真实数据的`PREVIOUS_LIU_RECONCILIATION.tsv`记录原200逐列对照；全部源码和结果凭据见本批结果目录及证据包。

## 字段及文件

- `ALL_GENE_EVIDENCE.tsv.gz`：全部37,944条；`*_source_feature_present`和`*_passes_expression_filter`分开记录来源可测与表达过滤；数值NA不填零。
- `A_BIOLOGICAL.tsv` / `B_BIOLOGICAL.tsv`：达到生物证据和稳定性规则的全部候选，尚未要求模型覆盖。
- `A_SIMULATABLE.tsv`：两研究均达差异门槛、同向、稳定，且两项研究的模型覆盖均达标。
- `B_SIMULATABLE.tsv`：只有一个研究达到差异门槛，另一研究同向；稳定性和模型覆盖达标。只作探索，不混作两研究共同差异。
- `PRIORITY_TOP200.tsv`：A可模拟层按预定顺序取最多200，实际少于200时保留真实数量。文件名表示上限。
- `rank_score_worst_FDR`：两研究校正值的较大者，仅排序，不是联合FDR；`rank_score_min_abs_logFC`是两研究绝对变化幅度的较小者。
- `cov_*_token_fraction`：实际模型序列出现率；`expr_fraction`是原生计数检出率，缺测时NA；不得混作同一指标。
- `PREVIOUS_200_COMPARISON.tsv`：旧顺序保留，列出本轮A、B、覆盖不足或不入选原因。
- `SEVEN_DATASET_EVIDENCE.tsv`：每候选七行；Chin两行的`comparison_unit=Chin_pair`，HCM-only方向标记`not_applicable_no_control`。
- `SUPPORT_PATTERNS.tsv`：描述研究方向和两HCM-only覆盖组合，不把不可比较基因算反向，也不把覆盖当成疾病证据。

最终验收只确认规则、计算与来源对应；没有证明跨研究供体完全独立、未测协变量已排除或敲除有益。
