# 虚拟细胞：Geneformer HCM 靶点挖掘

本仓库 [hcm-geneformer](https://github.com/zhome1226/hcm-geneformer) 用于 HCM 心肌细胞建模与 Geneformer 候选靶点研究。

**第一次了解项目，请先读[《项目工作总览》](PROJECT_OVERVIEW.md)**：用于说明数据为什么这样选、如何清洗、转化了多少，以及 200 个候选基因从哪里来。

项目计划用 [Geneformer](https://github.com/ctheodoris/Geneformer)（Theodoris et al., *Nature* 2023）模拟改变基因，研究肥厚型心肌病（HCM）心肌细胞的状态是否可能向对照靠近。具体模型和评价方案仍需确定，不能把表达差异直接当作敲除有效的证据。

截至 2026-11-15，本工作线已完成七套心肌核心数据的清洗和 V1 格式转化，共 151,359 条细胞或细胞核记录；按人比较患者与对照后，固定了 200 个候选基因。Liu2023 中 158 个总体同向，其中 141 个在两次单对照检查中仍同向；42 个反向，均保留原名单和顺序。见[中文候选表](results/target_discovery/hcm_dge_20261115/CANDIDATES_REVIEW.csv)和[结果报告](results/target_discovery/hcm_dge_20261115/REPORT.md)。

本流程尚未训练或运行模拟敲除。下一步按[候选与模拟策略](docs/hcm-target-strategy.md)核对具体模型、已接触的数据、比较对象和评价办法，再制定逐基因模拟方案。ToDesk 的新增三套数据及本次结果包仍待用户手动传输后确认，数据制作完成不等于目标端收全。

硬件：RTX 4090。主要发现来源为 [SCP1303 / Chaffin et al., Nature 2022](https://singlecell.broadinstitute.org/single_cell/study/SCP1303)。各阶段的详细报告及历史状态可从总览末尾的链接查阅。


时间节点见 Milestones：11.13 数据与流程跑通 → 11.20 Tokenization / 分类器 / 靶点池 → 11.24 模型交付 → 12.05 虚拟敲除与验证。

