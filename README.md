# 虚拟细胞：Geneformer HCM 靶点挖掘

本仓库 [hcm-geneformer](https://github.com/zhome1226/hcm-geneformer) 用于 HCM 心肌细胞建模与 Geneformer 候选靶点研究。

用微调后的 [Geneformer](https://github.com/ctheodoris/Geneformer)（Theodoris et al., *Nature* 2023）构建虚拟病理心肌细胞，对候选基因做计算机模拟敲除（In Silico Knockout），筛选能把肥厚型心肌病（HCM）转录表型拉回健康状态的逆转靶点。

硬件：RTX 4090。主数据为 [SCP1303 / Chaffin et al., Nature 2022](https://singlecell.broadinstitute.org/single_cell/study/SCP1303)，详见 [#1](https://github.com/zhome1226/hcm-geneformer/issues/1)；已下载的验证数据与补充候选见 [数据说明](data/README.md)。旧推荐 GSE145154 无 HCM，GSE130652 是小鼠 ChIP-seq，均不用于本任务。

下一阶段按[HCM 差异表达与模拟敲除候选策略](docs/hcm-target-strategy.md)先检查供体与计数来源，再开展供体级差异分析、候选分层和逐基因模拟删除。候选池目标为 100–200 个，显著差异表达与文献机制证据分别记录。

2026-11-14 已完成[第一轮供体与计数检查](results/target_discovery/preflight_20261114/PREFLIGHT_REPORT.md)：四套核心对象共 136,942 核，原始计数层与重复样本已记录；跨研究供体独立性及最终分析划分仍待确认，尚未运行 DGE。

补充完成[13 套资源盘点与模型来源核对](docs/hcm-cohort-catalog.md)：13 套包含比赛和 bulk，并非 13 个 HCM 队列。服务器已整理不含 SCP1303 的成人既有 QC 资源索引，共 30,458 核；GSE138262 按完整细胞单列。官方仓库同时包含预训练与心肌病微调模型，实际 checkpoint 和验证独立性尚待确定。

同日后续完成：[GSE138262 完整细胞新版 QC](docs/hcm-whole-cell-qc-review.md)保留 3,595 个心肌细胞（2,216 HCM、1,379 对照，5+7 位供体），旧版保留；[两套新增人源 HCM 数据](data/validation/HCM_EXPANSION_20261114.md)共 18 文件已下载校验，尚未 QC。没有运行新的 tokenization、DGE 或训练。

## 分工

| 角色 | 职责 |
| --- | --- |
| 张硕 | 数据下载与质控、与宇栋协作 Tokenization、DEG 候选靶点、虚拟敲除、生信验证 |
| 宇栋 | 官方数据跑通微调、公开数据分类器微调、11.24 前交付模型与对比评估 |
| 驰南 | 微调流程演示、协助数据格式对齐与模型部署 |
| 南珏 | 整体协调与资源支持 |

## 如何认领任务

任务都在仓库的 [Issues](../../issues?q=is%3Aissue+label%3A项目/虚拟细胞)。其他 agent 或同事按下面规则认领：

1. 打开带 `项目/虚拟细胞` 且状态为 `状态/待认领` 的 Issue，确认自己是负责人或协作人。
2. 评论「认领」，写上你的名字或 agent 名，把标签改成 `状态/进行中`。
3. 已是 `状态/进行中` 的不要重复认领；要接手先在 Issue 里说明。
4. 完成后把标签改成 `状态/已完成`，并在评论里写清交付物路径（本分支：`virtual-cell`）。

时间节点见 Milestones：11.13 数据与流程跑通 → 11.20 Tokenization / 分类器 / 靶点池 → 11.24 模型交付 → 12.05 虚拟敲除与验证。

## 仓库约定

- 本工作线的代码提交到 `virtual-cell` 分支。
- 代码、脚本、文档、小型结果表可以入库。
- 原始测序数据、大体积 `h5ad`、模型权重不要直接 push。超过 50 MB 的文件用 Git LFS 或 Release，并在 `data/README.md` 或 `results/README.md` 写获取方式。
- 目录：`data/` 数据说明与下载脚本，`src/` 代码，`results/` 结果，`docs/` 计划与报告。

## 待确认

任务文档用训练集 / 验证集 / 测试集 70/15/15，会议纪要用 8:2。微调划分以南珏确认为准，确认前先写在对应 Issue 里。
