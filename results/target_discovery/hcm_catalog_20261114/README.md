# HCM资源目录审计

报告：[HCM数据资源盘点与模型来源核对](../../../docs/hcm-cohort-catalog.md)。日期：2026-11-14（北京时间）。

- `ALL_13_RESOURCE_ROLES.tsv`：此前13套资源的模态、现状和用途，不等于13个HCM队列。
- `QC_RESOURCE_CATALOG.tsv`：服务器7个既有QC核心对象的路径、大小、数量、计数语义及限制。
- `DONOR_CATALOG.tsv`、`SAMPLE_CATALOG.tsv`：聚合供体/样本表，不含逐细胞表达或token。
- `CATALOG_AUDIT.json`：只读检查范围及GSE138262完整来源QC指标。计数值为抽样核查，没有重新做全矩阵验收。
- `COLLECTION_RECEIPT.json`：服务器正式目录及4个真实软链接的存在、目标路径与大小回执。
- `CATALOG_INPUT_CONFIG.json`：脚本配置；SHA256沿用既有验收，本轮没有重算全文件哈希。
- `MODEL_PROVENANCE.json`和三个`*CONFIG.json`：官方模型配置及本地checkpoint仍待确认的边界。
- `NEW_HCM_CANDIDATE_SAMPLES.tsv`：6个新增人HCM候选样本；仅元数据已取回，表达矩阵未下载。
- `SOURCES.tsv`：本机小型来源证据及哈希；论文全文和压缩元数据不提交Git。

生成代码：`src/target_discovery/catalog_hcm_resources.py`，与同目录`preflight_inputs.py`共同运行。

```bash
python catalog_hcm_resources.py --config CATALOG_INPUT_CONFIG.json --output /path/to/new_catalog
```

输出目录必须不存在。成功后创建指向4套既有非SCP成人QC对象的软链接，不拷贝、改写或合并矩阵。默认不分配训练/验证角色。

本轮正式服务器目录：`/home/zhome/hcm-geneformer-analysis/hcm_cohort_catalog_20261114_v2/`。较早无后缀目录是初次盘点，正式目录修正了GSE181764对照供体数的说明；两个目录均未修改源对象。
