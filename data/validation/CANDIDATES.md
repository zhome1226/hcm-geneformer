# HCM 补充数据检索与取舍

检索及下载决策日期：2026-11-10（Asia/Shanghai）。目标来自本仓库 [#1](https://github.com/zhome1226/xiaoe-internship/issues/1)、[#6](https://github.com/zhome1226/xiaoe-internship/issues/6)：用人心肌细胞的原始计数训练 Healthy/HCM Geneformer 分类器，随后验证虚拟敲除候选靶点。今晚仅检索、下载、完整性校验及元数据整理，不执行清洗、细胞筛选、差异表达或模型评估。

优先条件：人体心肌、真实 HCM、可追溯供体/样本、公开可下载的计数矩阵。分类验证优先具有对照的数据；疾病组单独数据及 bulk 数据只承担明确的补充用途。研究数、受试者数、组织取样数、技术重复数、细胞数分别记录。

## 决定下载

| 数据集及来源 | 已核实构成与文件 | 用途与限制 |
| --- | --- | --- |
| [GSE181764](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE181764)，Codden/Chin | 1 例梗阻性 HCM、6 例非梗阻性 HCM、2 例供体对照；9 份成人室间隔标本；原 tar 301056000 字节，27 个 MTX/TSV gzip 成员 | 扩展 HCM 亚型验证。属于现有 Chin 数据的同实验室扩展，不能把它和 GSE174691/GSE161921 算作三个完全独立研究。样本编号为 14–22；保留 GEO 元数据供后续按受试者核查。 |
| [GSE227734](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE227734)，Cox7B 研究 | 5 名成人 HCM 患者 P1–P5，室间隔；作者报告 338 个心肌细胞；`GSE227734_hcm_raw.h5ad.gz` 2398031 字节 | 独立疾病组心肌细胞特征及候选靶点复核。GEO 的单细胞提交未提供配套健康单细胞对照，不能单独计算 Healthy/HCM 分类 AUC；不能把细胞数当患者数。未载入 h5ad 检查表达矩阵。 |
| [GSE249925](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE249925)，Garmany et al. 2024 | 官方总体设计为 97 例 HCM、23 例对照；心肌组织 bulk RNA-seq；原始计数 TSV gzip 6968101 字节 | 后续在组织层面复核候选基因方向、疾病及性别相关性。不是单细胞/单核数据，不直接加入 Geneformer 单细胞分类训练或声称细胞特异性。 |

前两套组织/年龄信息与 bulk 的组别均来自官方 GSM 元数据，已原样保存在服务器；[SUPPLEMENTAL_SAMPLES.tsv](SUPPLEMENTAL_SAMPLES.tsv) 为 134 条 GEO 样本记录的索引，不是细胞注释或分析结果。GSE181764 与已下载的 GSE174691/GSE161921 的 GEO 样本编号及公开标本名称无重复，记录见服务器 `discovery20261110/chin_sample_identity_check.json`。正式划分数据前，仍须与 SCP1303 的供体标识及论文样本表交叉核对，才能确认受试者层面的训练/测试不重叠；本次没有声称已完成该核对。

## 暂不下载

| 候选及来源 | 原因 |
| --- | --- |
| [SCP1303 / Chaffin et al. 2022](https://www.nature.com/articles/s41586-022-04817-8) | 15 HCM、11 DCM、16 非衰竭心脏，主训练来源（#1）；用户后续已授权下载，范围与清单见 [training/README.md](../training/README.md)，不计为独立验证集。 |
| [GSE203274](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE203274) | 儿童先心病/心肌病混合研究。GEO 原始单核提交中 P86 为 11 岁 HCM 患者，GSM6165807/5808 为 rep1/rep2；供体亦含技术重复及左右心室取样。“2 HCM 样本”不能直接解释为两名独立患者。处理矩阵是否整合额外外部样本尚未核实，本次不作为成人独立验证集下载；未下载 scaled、PBMC、空间或其他模态。 |
| [GSE138262](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE138262) | 官方设计为 5 份 HCM 室间隔切除标本，没有同批健康单细胞对照；本轮已有体积更小的 GSE227734 承担疾病组补充用途，暂不继续扩充同类数据。 |
| [EGA EGAD50000002046](https://ega-archive.org/datasets/EGAD50000002046) | 较新的 HCM 单核资源，84 个样本来自 55 人；页面为受控访问 BAM，需要数据访问申请。不是今晚可直接下载的计数矩阵，不发起申请。 |
| [GSE145154](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE145154) | 不能用作明确的 HCM 队列；仓库 #1 已指出其无 HCM，旧 README 的推荐已过时。 |
| [GSE130652](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE130652) | 小鼠精子发生的 H3K27ac ChIP-seq，与人体 HCM 转录组任务无关。 |

## 完整性与复现

- 三个新增数据文件的字节数对照 NCBI HTTP `Content-Length`；GSE181764 tar 和 27 个解包成员另逐一对照官方 `filelist.txt`。两个 gzip 文件仅测试 gzip CRC，没有将表达矩阵载入分析软件。
- 每个新增服务器目录均有 `MD5SUMS`。本次未获得这三个数据文件的发布方参考 MD5，MD5 是本地实测并复核的完整性记录，不能称为“与官方参考 MD5 一致”。此前四个数据文件已有参考 MD5 匹配的结论保持不变。
- 原始元数据快照包括 GEO series 和 sample 文本。元数据字节数和 MD5 如实登记；没有独立发布方 MD5 或固定长度时，不声称完成了独立校验。
- 服务器路径、逐文件大小、哈希和用途见 [MANIFEST.tsv](MANIFEST.tsv)。下载脚本见 [download_supplemental.py](download_supplemental.py)，需传入已保存的、未经修改的 GEO 元数据缓存；本次在服务器运行 `python3 -u ~/xiaoe-internship-data/validation/_download/download_supplemental.py --metadata-cache ~/xiaoe-internship-data/validation/_download/discovery20261110`，使用 nohup 并记录 PID。
- 本次检索的官方页面文本、HTTP 响应头、样本索引保存在服务器 `validation/_download/discovery20261110/`，本地检索记录位于 `D:\VirtualCell\sources\`。详细下载日志为服务器 `validation/_download/supplemental_download.log`。
