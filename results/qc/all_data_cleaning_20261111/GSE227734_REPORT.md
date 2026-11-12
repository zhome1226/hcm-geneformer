# GSE227734：用户排除，不再清洗

2026-11-12 UTC：用户决定将此数据集排除。本套未完成分子质控、不可用于训练；原件及历史准备/来源审计保留。作者邮件未发送，发送授权已撤回，不继续联系作者或追查映射。以下为排除前的历史审计。

2026-11-11T12:32:20.595350+00:00；60 项准备审计检查通过。服务器：`/home/zhome/xiaoe-internship-data/processed/cox7b_GSE227734_qc_20261111_pending_annotations/`。

338 个已发表成人原代心肌细胞，来自 5 位 HCM 切肌患者；所有 9760 列原始计数与作者元数据保留。逐值重读原始 gzip 内的 HDF5，计数、细胞顺序、匿名列顺序及作者计数/检出数全部一致，原件摘要未变化。这里的通过仅指准备及来源审计，不能解释为已完成分子质控。

## 必须解决的缺失

GEO 官方 HDF5 只有 X、obs、var；var 无字段，列名仅 0–9759，没有 gene symbol、Ensembl ID、raw 或其他基因映射。不能按字母、参考基因组顺序或其他队列的基因顺序猜配。因而线粒体百分比为 NaN（未测得），不是 0%；心肌标记、参考身份及合适的混合谱系复核也无法执行。240 个细胞符合匿名检出列数 500–5000，98 个高于 5000；这不等于 240 个细胞通过 MT<5% 和身份质控。全部338个细胞 training eligibility=false。

技术来自 STRT-seq 分离成人原代心肌细胞，不能套用 10x 液滴 EmptyDrops 或按装载量估计双细胞。8 碱基条码反复使用，338 个细胞只有 10 种条码；以原始患者前缀细胞 ID 保持唯一，不按条码序列去重。作者 Large/Small 是细胞形态分类，保留原值及面积、长宽、年龄和性别；不转成 Control/HCM 或扰动标签。可访问的元数据未确认形态数值单位，未臆定单位。

GEO 系列设计与 PMID37979444 摘要说明这338个原代心肌来自HCM患者，五个GSM与作者样本别名按原记录绑定；本发布对象没有健康对照。论文包含小鼠AAV9实验，不把人细胞标为COX7B干预。作者心肌标签仅作为来源注释，尚未独立验证；跨研究供体重叠和训练划分未确认。

产物 author_counts_annotation_pending.h5ad 是真实计数与元数据的待注释整理版，不是可用于基因模型的清洗对象；未生成虚假的核心或轻量对象。cell_preparation_audit.tsv.gz、feature_annotation_audit.tsv.gz、sample_preparation_summary.tsv、qc_policy.json 和 PREPARATION_VERIFICATION.json 保存证据。

下一步必须取得与原9760列顺序对应的权威基因表或计数可对应的修正版矩阵，核对来源、顺序与计数后补做MT及身份复核。若只能从原始测序重建，将属于另一个需明确文库/细胞映射、参考版本和计算流程的步骤，不能把新计数冒充原作者矩阵。

## 官方补充材料来源复核（2026-11-11）

新增取得[官方补充Word](https://ars.els-cdn.com/content/image/1-s2.0-S0022282823001876-mmc1.docx)，3,853,376字节，SHA256 `6b6ef2ca9c0737b0f300fa87244bf9ba9acd9bcd048839e8951b046bc2c5de54`。独立解析6张表及附件结构并与原HDF5核对，9项来源核查通过：作者过滤前390个、过滤后338个，各供体37/75/70/66/90及10种8碱基条码均一致。该发布输入已由作者过滤；不能恢复为最初390个细胞。

补充表没有9760列的顺序基因映射，也没有嵌入的数据文件。Table S5只有20个面积相关基因，不能对应匿名列。P5年龄原矩阵/GEO为63岁、补充Table S3为64岁，冲突未解决；保留两来源，不自动改写原元数据。原始文库到作者细胞的完整映射、UMI布局和具体注释版本仍不足以支持无猜测的FASTQ重建。

详见[gse227734/SOURCE_LOOKUP_SUMMARY.md](gse227734/SOURCE_LOOKUP_SUMMARY.md)与该目录新增`source_lookup_followup_20261111/SUPPLEMENT_SOURCE_AUDIT.json`、`SUPPLEMENT_DONOR_CROSSCHECK.tsv`。服务器对应证据在待注释目录的`source_evidence/`；原Word在`references/GSE227734_publisher_20261111/`。计数未改，全部338个细胞仍不可训练，清洗状态保持未完成。
