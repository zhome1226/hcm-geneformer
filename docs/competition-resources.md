# 可复用的单细胞比赛资源

核查日期：2026-11-11（Asia/Shanghai）。围绕本项目的 HCM 分类、虚拟敲除与靶点复核，优先收集已结束比赛的数据入口、获奖代码及评估实现。本轮不执行第三方代码、不安装其环境、不清洗或分析表达矩阵。

这些比赛多数不以 HCM 为研究对象：VCC2025 为胚胎干细胞，NeurIPS2023 为外周血细胞，CausalBench 为细胞系，NeurIPS2021 为骨髓。ACDC 确实含 HCM，但属于 MRI。它们用于学习方法与评估，HCM 转录组训练仍采用 [SCP1303](../data/training/README.md) 等临床研究数据。

## 推荐顺序

| 优先级 | 比赛与状态 | 数据、代码现状 | 对我们最有用的部分 |
| --- | --- | --- | --- |
| 1 | [Arc Virtual Cell Challenge 2025](https://arcinstitute.org/news/virtual-cell-challenge-2025-wrap-up)，已结束 | H1 人胚胎干细胞，约 30 万单细胞、300 个 CRISPRi 靶基因；[官方 Atlas](https://github.com/ArcInstitute/arc-virtual-cell-atlas/tree/main/virtual-cell-challenge) 已公开 train/validation/test；[Cell-Eval](https://github.com/ArcInstitute/cell-eval) 为 MIT | 学习如何用真实扰动实验检验预测，并组合 DES、PDS、MAE 等指标；先看评估代码，再决定需要哪份表达矩阵。 |
| 2 | [NeurIPS 2023 / Open Problems – Single-Cell Perturbations](https://openproblems.bio/events/2023-08_neurips/)，已结束 | 人外周血细胞的小分子扰动；[冠军代码](https://github.com/Jean-KOUAGOU/1st-place-solution-single-cell-pbs) MIT；另有[亚军](https://github.com/Eliorkalfon/single_cell_pb)、[季军](https://github.com/okon2000/single_cell_perturbations) 公开仓库，但本轮只存档并核实冠军许可 | 学习训练集统计特征、交叉验证、简单模型与集成基线。比赛的 DE 预测表不是单细胞原始 UMI 计数，不能直接当 Geneformer Tokenizer 输入。 |
| 3 | [GSK.ai CausalBench 2023](https://gsk.ai/causalbench-challenge/)，已结束 | K562、RPE1 遗传扰动数据；[官方框架](https://github.com/causalbench/causalbench) Apache-2.0；[GuanLab 冠军](https://github.com/GuanLab/CausalBench2023_Challenge) MIT；数据下载地址见官方 `causalscbench/data_access/datasets/download_weissmann.py` | 学习干预数据下的基因调控网络推断与评估，可辅助审查靶点的机制依据。这里特指 GSK 的生物学 CausalBench，不能与同名 LLM/通用因果平台混淆。 |
| 4 | [NeurIPS 2021 Multimodal Single-Cell Integration](https://openproblems.bio/events/2021-09_neurips/)，已结束 | 人骨髓多组学；官方活动页仍链接 [GSE194122](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE194122)；[官方优胜方法合集](https://github.com/openproblems-bio/neurips2021_multimodal_topmethods) MIT | 学习跨供体、跨批次表示学习与评估，也可参考 RNA/ATAC/蛋白整合流程；目前单 RNA 的 HCM 任务不优先下载全量多组学。 |
| 备选 | [MICCAI 2017 ACDC](https://www.creatis.insa-lyon.fr/Challenge/acdc/databases.html)，已结束 | 官方说明比赛结束后数据和标注仍公开；100 人训练集含 20 HCM，50 人测试集含 10 HCM；属于心脏 cine MRI，访问按原平台规则 | 确实含 HCM，适合将来做影像分割/疾病分类学习；本项目没有配对的影像与转录组，因此不下载、不并入 Geneformer。 |

截至本次核查，主办方新闻给出的 [VCC 2026](https://arcinstitute.org/news/virtual-cell-challenge-2026) 最终提交截止为 11 月 5 日，获奖公布在 11 月中下旬。本轮未核实到完整获奖代码发布，故不把它列作已可复现的优胜方案。当前比赛域名已是 2026 版，2025 数据需使用上表归档入口；比赛结束也不意味着所有参赛者必须公开代码或权重。

## 已实际保存的资源

服务器目录：`/home/zhome/xiaoe-internship-data/benchmarks/competitions/`，与 HCM 临床验证目录分开。

- 5 份固定 commit 的源码 tar.gz：Cell-Eval、NeurIPS 2023 冠军、CausalBench 官方框架及冠军、NeurIPS 2021 官方优胜方法合集；保留原许可证，没有执行或修改第三方代码。
- VCC 2025 的 4 份 CSV：`gene_names.csv` 及 train/validation/test 的扰动细胞数清单。按官方 GCS 对象大小与 MD5 校验。
- 逐文件来源、commit、代码许可、字节数和本地 MD5 见 [competition_MANIFEST.tsv](../data/benchmarks/competition_MANIFEST.tsv)；服务器另有 `MD5SUMS`。
- VCC 三份矩阵及 CSV 的完整对象索引见 [competition_DATA_CATALOG.tsv](../data/benchmarks/competition_DATA_CATALOG.tsv)。本轮**未下载任何比赛表达矩阵或模型权重**：train 15482497461、validation 6928967541、test 11950739168 字节，合计约 34.36 GB；不下载约 5.28 TiB 的原始 FASTQ。
- 服务器不能连接 GitHub codeload，本轮代码先在 PC 下载，再通过 SCP 传入；官方 GCS CSV 可由服务器匿名下载。OpenProblems S3 列表探测在服务器超时、PC 返回 403，不能据此推断对象一定不可读取；本轮未取得 Kaggle 数据，也未登录或接受竞赛条款。
- 代码与数据的许可分别判断。上表的 MIT/Apache 是核查过的代码许可，不自动授权第三方数据或权重的任何用途。正式复制代码时保留相应 LICENSE、版权及必要 NOTICE；数据依原发布方条款与引用要求使用。

存档脚本为 [download_competition_resources.py](../data/benchmarks/download_competition_resources.py)，使用固定 commit 清单和 GCS 对象索引；校验报告、版本清单和日志保存在服务器 `_catalog/`。代码 tar.gz 的 MD5 为本地实测并复核，未声称具有发布方参考 MD5；4 个 CSV 则匹配发布方 GCS MD5。

## 如何用于当前项目

1. 保持现有 HCM 分类数据线。上述干细胞、血细胞和细胞系数据用于方法学习与实验扰动基准，不能替代 HCM 心肌独立验证。
2. 先阅读 NeurIPS 2023 冠军与 CausalBench 冠军，提炼简单基线、划分方式和误差评估。模型训练前再审查其特征构造是否严格仅使用训练折信息，不能照搬一个比赛的交叉验证设置到跨供体临床验证。
3. 虚拟敲除阶段再选择 VCC 2025 的适当数据分区作真实扰动基准。当前 Geneformer 分类器的 ΔP/PRS 不是预测表达矩阵，不能直接交给 Cell-Eval 的 DES/PDS/MAE；需要先定义适配的可检验输出，或另设表达预测基线。CRISPRi 是转录抑制，也不等同于完全基因敲除。

以上是资源筛选与后续建议，尚未完成任何优胜方案复现或证明 HCM 靶点有效。
