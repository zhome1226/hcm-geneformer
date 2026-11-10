# 代码

建议按阶段分目录，例如：

- `src/qc/` 质控与心肌亚群筛选
- `src/tokenize/` Geneformer Rank-value Tokenization
- `src/finetune/` 疾病分类器微调与评估
- `src/deg/` 差异表达与候选靶点筛选
- `src/knockout/` InSilicoPerturber、ΔP、PRS
- `src/enrichment/` GO/KEGG 与已知靶点对齐

预训练权重：`ctheodoris/Geneformer`（优先 `Geneformer-6L-30M`）。
