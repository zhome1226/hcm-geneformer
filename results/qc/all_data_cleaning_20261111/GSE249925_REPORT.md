# GSE249925：bulk 表达矩阵清洗与逐值验收

2026-11-11T14:28:33.349000+00:00；430 项检查通过。服务器产物：`/home/zhome/xiaoe-internship-data/processed/bulk_GSE249925_qc_20261111/`。

## 完成内容

原发布 64253 行特征×120 个样本，97 HCM、23 Control，全数通过样本对应与非零信息检查。可用矩阵保留 48629 行，移除全队列恒为零的 15624 行；结构无效特征 0 行。原件和完整发布矩阵全部保留，所有数值与原始 gzip 内表格逐值一致，未归一化、插补、重建或汇总任何表达值。清洗仅表示该发布矩阵、源标签与 QC 产物已验收；不是原始计数语义认证或训练许可。

两个真实 AnnData：all_published_expression.h5ad 保存所有源特征与样本；bulk_qc_expression.h5ad 保存通过特征 QC 的矩阵；bulk_qc_matrix.tsv.gz 保存相同可用矩阵及原始 Gene_Index/Length。X 与 published_values 层完全相同；没有伪称 raw counts 的 counts 层。sample_metadata_QC.tsv、feature_audit.tsv.gz 和排除索引可追溯每个样本及原特征行。sample_exclusions.tsv 的空结果表示真实没有被排除样本，不是占位数据。

## 来源与必须保留的语义限制

GEO GSE249925（https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE249925）；Garmany et al., Journal of Molecular and Cellular Cardiology 2024，PMID39255898，DOI10.1016/j.yjmcc.2024.09.005。原研究报告梗阻性 HCM 切肌病例和供体心脏对照。GEO 的 hg38、Ensembl v78、STAR/HT-Seq/edgeR 处理说明已存档。

GEO 文件名及补充文件描述称其为原始基因计数，但同一处理字段又明确提到所提供读数经过 TMM 标准化。全部数值虽为非负整数，这一点不足以排除已处理或取整来源。已检查 GEO、论文元数据与可访问的出版商摘要/方法片段，尚未取得能解决“该具体文件是否经过 TMM”的权威说明；出版商全文返回403、预印本全文在本次查询中不可访问。EuropePMC在服务器超时，本地成功取得期刊论文元数据；该期刊记录没有可直接读取的PMC全文。保留这些证据与检索限制，**certified_unnormalized_counts=False**。文件可用于受此限制的后续表达验证设计，不能冒充已确认的未归一化 HT-Seq 原计数直接投入计数模型。

## 基因与样本 QC

原始 chromosome@Ensembl@symbol@start@end、原 Length、顺序和所有行均保留；Length 是源注释，不按起止跨度重新推定。重复 Ensembl 行 0，重复 symbol 行 7064，只标记、保留原始特征，不按符号或 Ensembl 擅自求和。恒零行从可用矩阵移出，完整矩阵及 feature_exclusions.tsv.gz 仍含这些原始零行。未按低表达或疾病差异筛基因。

所有标题与 GSM/BioSample/SRA 精确对应，标签直接取 GEO group：HCM=1、Control=0。数字字符串和 R 等前后缀原样保留，不据此猜性别、供体或技术重复。每个 GSM 的源关系均核对；逐样本 donor、age、sex、batch 未提供，donor_id 留空并显式标为未知，不将120个GSM自动认证为120位可独立划分患者，也不将论文性别总数猜配到样本。同一病例跨研究或多标本重叠在训练前仍需审查。

发布表达总量 42135002–109606788；检测特征数 25390–32339。log总量/检出特征的5MAD标记分别 2/0 个，仅用于复核、不作删除规则。bulk MT仅按来源染色体计算描述性比例，未套用MT<5%。未执行液滴识别、双细胞检测、心肌提取、单细胞500–5000阈值、差异表达或模型训练。该矩阵是组织混合表达，不能把每个样本当作一枚心肌细胞。

## 验收与使用边界

三个原件大小和此前记录的本地 MD5 均一致；旧下载清单没有 SHA256，本轮新增测得值写入验收记录，不能说与不存在的旧 SHA256 比对通过。该 MD5 是下载基线复查，不宣称获得出版方独立 MD5。首轮验收误把空SHA256当作旧校验值导致退出1，失败日志、源码和原清单已归档；修正后只重跑独立验收，不重跑已经成功的清洗。独立脚本从 gzip 原表及 SOFT 重新解析，核对所有表达值、基因轴、样本轴、120份GEO分组和注册关系、全部元数据/QC选择、文件摘要及执行源码。FINAL_COMPLETION_GATE.json 只在这些检查全部通过后生成；还须核实 verify_exit_code.txt=0 才能记完成。

bulk清洗完成不消除 raw/TMM 语义疑点，也不证明测序质量、性别/批次无混杂或跨队列供体独立；FASTQ/比对级QC本次未执行。训练划分未分配。所有矩阵、逐特征审计保留服务器，Git仅收代码和文本汇总。
