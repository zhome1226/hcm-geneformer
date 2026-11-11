# GSE138262：修正版计数清洗与独立验收

2026-11-11T13:19:30.643913+00:00；735项检查通过。输出：`/home/zhome/xiaoe-internship-data/processed/wehrens_GSE138262_qc_20261111/`。

已清洗验收，不等于适合训练。严格基线仅留下8个数值候选、6个心肌核心、2个身份待定；核心均为HCM，只有2位供体，无健康对照。本套不批准进入HCM/Control微调或独立疾病验证。没有为凑数量放宽规则、复制细胞或混入参考图谱。

## 输入更正与计数来源

作者仓库 https://github.com/vanrooij-lab/scRNAseq-HCM-human 明确更正：GEO旧TSV不是Wehrens et al., Cell Reports 2022实际使用的计数表。按作者链接取得概念DOI 10.5281/zenodo.18206553指向的版本DOI 10.5281/zenodo.18206554。5个文件全部符合Zenodo大小和MD5；32个完整修正UMI计数表保留。原下载清单的53个文件大小/MD5再次一致，旧GEO数据留作来源存档，没有覆盖。

16个CEL-seq2分选完整心肌细胞板来自5位HCM切肌患者，GSM及患者映射核对GEO和作者文件名；另外16个iCell8板/部分来自7位Wang健康供体，是作者按LV细胞群选取后重处理的对照伴随数据，保留其独立研究/技术来源。N2的多板、各板part1/2均绑定原供体；不把板当独立患者。Seurat ALL中另外的Litvinukova细胞单独列入来源审计，不重复加入当前队列，避免与现有参考图谱重用。

每个修正计数表均独立解析并与native、并集、最终各对象逐值核对；48组共享板/作者Seurat RNA计数全部一致。RNA assay已经删去线粒体基因；作者保存的线粒体计数元数据也与完整修正计数重算一致。MT百分比来自完整表，绝非从裁掉MT后的RNA层算出0%。原始命名、Ensembl ID及空gene symbol全部保留；并集28663个基因，native轴与逐样本存在掩码保留，不以共有基因替代完整计数。

## 固定严格规则与真实损失

500≤检出基因≤5000、MT<5%是既定项目基线，不是所有HCM完整细胞的统一质量定义。逐疾病统计如下：

```text
         source_cells  genes_500to5000  MT_lt5pct  strict_cells  CM_core  pending_identity
disease                                                                                   
Control          1423             1410          0             0        0                 0
HCM              6114             2635         53             8        6                 2
```

完整心肌细胞的线粒体RNA占比可远高于单核数据；本套强烈筛除，不能据此称其余细胞均为坏细胞。完整7537个来源细胞及所有计数保存在all_source_counts.h5ad；不满足严格数值规则者单独保存outside_strict_numeric.h5ad，逐细胞及32个样本/12位供体损失有表可查。未来若制定完整细胞专用阈值，应另发版本并与本严格版比较，本轮未改阈值。

## 身份、混合信号与适用边界

SingleR使用27664个共同Ensembl基因和固定供体分隔的成人左心室单核参考；参考域与本完整细胞并不完全一致。联合多基因程序、程序优势及参考分数尾部判断心肌身份；评分不是后验概率。只有8个严格候选，不使用小样本邻域去挽救身份。2个身份未通过心肌核心门槛的细胞保留待定，不强制冒充心肌。

本数据为板式分选CEL-seq2和已选iCell8细胞，不运行EmptyDrops，不套用10x装载量双细胞率；高基因数不等于双细胞。强跨谱系程序按固定参考阈值、多基因/多UMI支持复核，结果仅为可疑混合信号，不能排除同型双细胞或证明单细胞纯度。未经校准的单一工具标签不用于强删。诊断log1p(CP10k)副本有明确标识，计数对象X与counts逐值不变，未插补、整合、按HVG裁切或开展差异表达/训练。

此套HCM与健康伴随来源具有研究、技术、组织取材及选择规则差异；跨研究供体重叠仍未解决，任何后续划分须按供体绑定。本严格版缺乏可用对照，6个核心不足以代表原5位HCM供体，更不支持疾病建模结论。

## 产物与验收

all_source_counts.h5ad：全部来源计数及审计；all_cells_qc.h5ad：8个数值候选；cm_qc_full.h5ad与cm_qc_light.h5ad：6个真实心肌核心；cm_review_pending.h5ad：2个待定。mixed_lineage_review.h5ad为真实规则过滤结果，本轮0个，不是伪造候选。outside_strict_numeric.h5ad保留7529个数值规则外的细胞。所有计数对象保留完整并集基因轴与原值，诊断对象单独标识。

INDEPENDENT_VERIFICATION.json记录逐值验收范围、原件与源码摘要、供体分隔、数值/身份/混合规则和小样本限制；FINAL_COMPLETION_GATE.json只在全部检查通过后写入。参考校准表、sample_qc_summary.tsv、donor_qc_summary.tsv、sample_metadata.tsv、all_cell_audit.tsv.gz与cell_identity_audit.tsv.gz可追溯。只将代码与文本汇总放Git，矩阵及逐细胞表留服务器。
