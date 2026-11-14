# HCM数据资源盘点与模型来源核对

日期：2026-11-14（北京时间）。本轮完成已有QC对象盘点、官方模型配置核对及新增候选检索。没有重新清洗、拼接表达矩阵、划分训练集、运行DGE或训练。

## 现在有哪些可以用的资源

此前的“13套”是7套心肌单细胞/单核资源、1套心肌组织bulk和5套比赛资源，不是13个独立HCM病例对照队列。GSE227734按用户决定排除，不发作者邮件。

本轮连接原Linux数据服务器`ad2012`，以只读方式核对7个`cm_qc_full.h5ad`的供体、样本、疾病标签、基因轴及分散抽取的计数值。文件大小和修改时间未变。检查通过表示元数据及所抽计数一致，不等于重新全量验收或证明数据能支持疾病结论。

| 已有资源 | 现有心肌核心 | 核心中的研究内供体 | 当前用途与限制 |
|---|---:|---|---|
| SCP1303 | 114,711核 | 15 HCM + 16 NF | 可作发现/参考；使用官方心肌病微调权重时，不能自动当全新外部验证。NF为非衰竭供体。 |
| GSE174691 | 10,150核 | 9 HCM | 成人HCM支持；与GSE161921配对时，疾病与accession关联。 |
| GSE161921 | 3,534核 | 4对照 | 成人对照支持；批次校正不能自动解决疾病/accession混杂。 |
| Liu2023 | 8,547核 | 10 HCM + 2对照 | 建议优先保留作后续验证候选；CL0412A/B同供体，年龄/性别缺失，两个对照的效能有限。独立性尚未批准。 |
| GSE181764 | 8,227核 | 6 HCM + 2对照 | 8,222 HCM核、仅5对照核，适合作有限HCM补充；含终末期病例，保留缺测基因掩码。 |
| GSE138262 | 6个完整细胞 | 2 HCM + 0对照 | 旧严格版不足以进行疾病比较；优先复核完整细胞专用QC。 |
| GSE203274 | 34,532核 | 1儿童HCM + 4对照 | 儿童资源单独保留；P86的多个样本不能当多位患者。 |

GSE181764的实际保留供体是**6 HCM + 2对照**，缺失的是HCM供体Chin_2930。两个对照Chin_2900和Chin_3045分别仅1核和4核。任何此前“7 HCM + 1对照”的推算应以本次实际对象和[供体表](../results/target_discovery/hcm_catalog_20261114/DONOR_CATALOG.tsv)为准。

不含SCP1303的现有成人单核集合为Liu2023、GSE174691、GSE161921、GSE181764，共**30,458核：23,090 HCM、7,368对照**。合计33个研究内供体条目（25 HCM、8对照），不等于已证实33位跨研究独立的人，也不等于模型从未见过这些数据。

服务器已建立按研究指向原有QC对象的索引目录：

```text
/home/zhome/hcm-geneformer-analysis/hcm_cohort_catalog_20261114_v2/
  adult_non_scp_existing_qc/
    chin_GSE161921.h5ad -> 原有QC对象
    chin_GSE174691.h5ad -> 原有QC对象
    liu2023.h5ad        -> 原有QC对象
    chin_GSE181764.h5ad -> 原有QC对象
  QC_RESOURCE_CATALOG.tsv
  DONOR_CATALOG.tsv
  SAMPLE_CATALOG.tsv
  CATALOG_AUDIT.json
```

这里的`v2`只是本轮目录报告修订号，集合中的模型编码目标仍为Geneformer V1。目录不包含新拼出的“统一训练矩阵”。原数据真实路径仍在`/home/zhome/xiaoe-internship-data/processed/`，本轮没有更名或移动。ToDesk的`D:\hcm-geneformer-data`未操作。

SCP1303、Liu及Chin两套已经分别完成V1序列化；GSE181764、GSE138262、GSE203274没有因此自动产生V1 `.dataset`。本轮没有重启传输或检查用户已自行重传的Arrow。

## 你给的模型链接意味着什么

已核对官方仓库revision `1f7fbae4e469a5f4f1af8c111a529cfe1b3829f5`及其小型配置，没有下载权重。一个仓库包含不同模型，必须指定子目录：

| 路径 | 实际配置 | 对本项目的含义 |
|---|---|---|
| `Geneformer-V1-10M` | `BertForMaskedLM`，2048长度，25,426词表 | V1预训练基础模型；可在我们冻结供体划分后自行微调。 |
| `fine_tuned_models/Geneformer-V1-10M_CellClassifier_cardiomyopathies_220224` | `BertForSequenceClassification`，3个输出标签 | 官方心肌病微调分类器；原论文心肌病任务用过Chaffin来源材料。SCP1303须按已接触风险处理，具体训练/留出供体仍未核验。 |
| 仓库根目录当前默认权重 | V2-316M，4096长度，20,275词表 | 与我们现有V1编码不匹配，不能仅指定仓库根目录就默认得到V1。 |

来源：[官方模型卡](https://huggingface.co/ctheodoris/Geneformer/blob/1f7fbae4e469a5f4f1af8c111a529cfe1b3829f5/README.md)、[V1基础配置](https://huggingface.co/ctheodoris/Geneformer/blob/1f7fbae4e469a5f4f1af8c111a529cfe1b3829f5/Geneformer-V1-10M/config.json)、[心肌病分类器配置](https://huggingface.co/ctheodoris/Geneformer/blob/1f7fbae4e469a5f4f1af8c111a529cfe1b3829f5/fine_tuned_models/Geneformer-V1-10M_CellClassifier_cardiomyopathies_220224/config.json)、[原论文](https://pmc.ncbi.nlm.nih.gov/articles/PMC10949956/)。

分类器配置仅记录`LABEL_0/1/2`，不能把它们擅自对应成我们二分类的`0=Control、1=HCM`。使用该分类器前需核验原标签映射与训练记录。

Genecorpus-30M仓库中确实有心肌病**示例微调数据**。示例文件位于同一仓库，不证明其属于预训练语料。须分别记录“预训练是否接触”和“有监督微调是否接触”；不能仅凭论文年份、GEO公开日期或不同accession证明未见过。用户提供的是总仓地址，尚未确认实际加载的本地checkpoint。[官方语料仓库](https://huggingface.co/datasets/ctheodoris/Genecorpus-30M/tree/main)

模型接触史影响模型评估的独立性，但不直接否定SCP1303供体级DGE的统计用途。可以用它发现候选，不能再用同一材料宣称对模型和候选做了完全独立验证。最终需把“DGE发现”“模型微调/选参”“候选验证”分别记录，保留供体级留出。

## 为什么GSE138262要单独复核QC

完整来源对象仍有7,537个细胞：5位HCM患者的6,114细胞，以及另一研究/技术来源的7位健康供体的1,423伴随细胞；并未丢失。

既有QC指标显示：对照MT比例最低11.26%、中位59.26%；HCM中位42.52%。旧规则`500≤genes≤5000、MT<5%`会排除全部对照，最终仅余6个HCM核心。完整细胞与单核的线粒体信号不同；这说明统一严格基线不适用，不能反推所有被删细胞都质量差。

下一步读取作者原始纳入标记、方法和逐板/供体分布，联合复杂度、测序深度、细胞身份与混合信号制定完整细胞专用版本，保留旧版和全部来源对象，量化不同疾病的筛除偏差。不能为了增加数量临时套一个宽松MT阈值；即使重新QC成功，HCM与健康来源的技术/研究混杂仍然存在。依据：[既有逐值验收报告](../results/qc/all_data_cleaning_20261111/GSE138262_REPORT.md)、[作者更正仓库](https://github.com/vanrooij-lab/scRNAseq-HCM-human)。

## 新查到的HCM数据及重复风险

以下只下载了公开元数据，**没有下载表达矩阵，也没有清洗**。逐样本来源及下载URL见[新增候选表](../results/target_discovery/hcm_catalog_20261114/NEW_HCM_CANDIDATE_SAMPLES.tsv)。

- **GSE255296**：官方元数据共29样本，其中人4、猫10、鼠15。人样本为GSM8068510/5246、GSM8068511/5042、GSM8068512/2883、GSM8423481/2795。论文的人体比较为5 HCM + 4对照，其中4对照和1 HCM复用GSE161921/GSE174691；新增四个人样本可作HCM扩展候选，但不能把整篇研究算作全新病例对照队列，更不能纳入猫鼠数据。[GEO](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE255296)、[Ali et al., JAHA 2025](https://pmc.ncbi.nlm.nih.gov/articles/PMC12054467/)
- **GSE292754**：公开新增HCM2/P186与HCM3/P194两个儿童样本（GSM8864696、GSM8864697）。论文的第三例HCM1=P86此前已在GSE203275公开，并已包含在我们的GSE203274资源中；不能重复记一位新患者。儿童扩展单独保存，复用公共对照须去重。[GEO](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE292754)、[Tadros et al., 2025](https://pmc.ncbi.nlm.nih.gov/articles/PMC13067291/)

两套可扩展不同年龄层，均尚未获得“与模型训练完全独立”的批准；每个新增样本还需核对临床、供体别名及计数来源。

## 其余资源如何处理

GSE249925是97 HCM + 23对照的组织bulk；现有raw/TMM来源语义仍待解决，之后可用于组织层面方向复核，不能拼进心肌单细胞矩阵。VCC2025、NeurIPS2023、NeurIPS2021、CausalBench2023用于学习扰动/多模态/因果推断方法，不能赋成HCM病例对照。ACDC2017虽包含HCM MRI，但没有基因表达矩阵，不能作为Geneformer输入。本轮未重做这6套资源的下载或文件完整性验收。

## 执行顺序

1. 先复核GSE138262完整细胞QC，另存版本；GSE181764保留为受限HCM补充，儿童保持单独分层。
2. 固定具体V1 checkpoint、标签映射和供体级发现/验证角色；Liu优先保留验证，不把全体研究用于筛选后再称独立验证。若使用官方心肌病微调权重，SCP按已接触风险处理。
3. 完成上述边界后，按研究/供体建立可比较的pseudobulk，再做DGE及跨研究复核；保持原始或来源校正计数语义，不用整合表达、token序列或缺测补零代替计数。

完整清单、来源摘要和审计记录：[结果目录](../results/target_discovery/hcm_catalog_20261114/)。
