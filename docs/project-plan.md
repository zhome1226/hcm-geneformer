# 项目计划

核心目标：针对 HCM（及文档中的扩张型心肌病背景），微调 Geneformer 疾病分类器，再对 DEG 得到的约 100 个候选基因做 In Silico Knockout，用分类概率逆转度（ΔP）和向量流形漂移度（PRS）排序逆转靶点。

## 时间线

| 截止 | 任务 | 负责人 |
| --- | --- | --- |
| 11.13（纪要：11.14 前完成下载与格式对齐） | 公开数据质控，产出带标签 AnnData；对齐官方微调数据格式 | 张硕（格式对齐协同驰南） |
| 11.13 | 官方示例数据演示并跑通微调，拿到正向结果 | 驰南演示，宇栋跑通，张硕一起跟流程 |
| 11.20 | 公开数据 Rank-value Tokenization，产出 Hugging Face `.dataset` | 宇栋、张硕 |
| 11.20 | 新数据疾病分类器微调，验证集 AUC-ROC ≥ 0.90 | 宇栋（张硕参与） |
| 11.20 | DEG（p ≤ 0.01 且 log2FC > 1）筛前 100 个转录因子、激酶、病理分子 | 张硕 |
| 11.24 | 对比微调前 baseline 与微调后模型，交付权重、ROC-AUC、混淆矩阵 | 宇栋 |
| 12.05 | 对 100 个候选基因虚拟敲除，计算 ΔP、PRS，出 Top 5 与 UMAP | 张硕 |
| 12.05 | GO/KEGG 富集，与 HDAC、CALM 等已知靶点对齐，写生物学验证 | 张硕 |

## 关键参数

- 微调参考：`batch_size` 8 或 16，Gradient Accumulation = 4，`learning_rate` 5e-5，Cross-Entropy。文档划分为 70/15/15，纪要为 8:2，待南珏确认。
- 敲除：`InSilicoPerturber`，`perturb_type="delete"`，单基因敲除，`model_type="CellClassifier"`。
- ΔP = P_Disease(原始病理细胞) − P_Disease(敲除后细胞)，越接近 1 越好。
- PRS 为敲除后相对健康质心的距离缩短率，> 0 表示向健康流形靠拢。
