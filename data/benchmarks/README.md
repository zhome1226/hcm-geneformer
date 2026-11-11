# 比赛资源与方法基准

这些资源用于学习单细胞扰动预测、调控网络推断和多组学评估，与 `data/validation/` 的 HCM 队列分别管理。

核查日期：2026-11-11（Asia/Shanghai）。详细取舍、官方数据及获奖代码入口见 [比赛资源报告](../../docs/competition-resources.md)。

用户随后已授权补齐完整比赛数据；新批次范围、许可、实时状态与已验证文件清单见 [expansion/README.md](../expansion/README.md)。下列 9 文件清单记录首轮代码/元数据存档。

- 服务器：`zhome@10.158.203.133:/home/zhome/xiaoe-internship-data/benchmarks/competitions/`。
- 已下载资源：[competition_MANIFEST.tsv](competition_MANIFEST.tsv)，包括 5 份固定 commit 的代码存档和 VCC 2025 的 4 份数据索引 CSV。
- 已核实可用的 VCC 对象：[competition_DATA_CATALOG.tsv](competition_DATA_CATALOG.tsv)，明确区分已下载元数据和未下载表达矩阵。
- 下载与校验脚本：[download_competition_resources.py](download_competition_resources.py)，输入为 [固定版本清单](competition_REPOSITORIES.json) 和 [GCS 对象快照](competition_VCC_OBJECTS.json)。原始记录亦保存在服务器 `_catalog/`；代码由 PC 下载并 SCP 后，可使用 `--verify-existing-code` 复核已有存档。
- 数据、源码存档、权重不直接提交 Git；本目录只存说明、脚本及小型清单。本轮没有清洗、分析、模型训练或运行第三方代码。
