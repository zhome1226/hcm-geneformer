# 新增 HCM 单核 QC、V1 编码与交付

本批仅处理 GSE255296 四个人源样本、GSE292754 的 P186/P194；将通过复核的核心与此前验收的 GSE138262 完整心肌细胞版分别序列化。复用已有环境与已审计的单核流程，保留全基因计数、原生基因轴及缺测掩码。疾病标签来自供体诊断。

运行目录：`/home/zhome/hcm-geneformer-analysis/hcm_expansion_20261114/nuclear_qc/`。正式 Linux 输出根：`/home/zhome/hcm-geneformer-data/`；ToDesk 根：`D:/hcm-geneformer-data/`。

顺序：`run_qc.sh` → 核实 `qc_exit_code.txt=0` 与两套 `INDEPENDENT_VERIFICATION.json` → `run_v1.sh` → `v1_exit_code.txt=0` 与 smoke/full 验收 → `package_delivery.py` → 传输归档与 PACKAGE_STATUS → ToDesk 执行 `accept_windows.py`。

- `stage.py`：按样本读原矩阵并计算 QC；`SAMPLES.json` 和选定 GEO SOFT 快照固定来源。两位儿童的年龄和性别直接核对 GEO。
- `identity.R`：固定成人 LV 单核参考、供体分离校准、SingleR、逐捕获 scDblFinder；没有低 UMI 背景，未运行 EmptyDrops。
- `calibrate_mixed.py`、`finalize.py`：固定混合谱系证据，心肌程序与参考/邻域支持联合筛选，待复核单列。单核数值基线 500≤基因≤5000、MT<5%；不将算法单独调用当作实证双细胞。
- `verify.py`：重新读原 MTX，核对全部输出计数、元数据、过滤分区、参考阈值、混合程序和邻域规则；保留全部来源条码损失审计。
- `prepare_v1.py`：在副本添加 Ensembl 和编码元数据；GSE138262 使用真正的 `var.gene_id`，不把复合索引当基因符号。保留源计数、layers、varm 和原有注释。
- `tokenize_v1.py`：原版官方 TranscriptomeTokenizer、V1、2048、无 CLS/EOS；先 smoke 再全量。来源 revision `1f7fbae4e469a5f4f1af8c111a529cfe1b3829f5`，运行时核对源码摘要。
- `verify_v1.py`：不用官方排名函数，独立计算词表映射、median-scaled 排名及截断边界，逐细胞核对所有 token 与元数据。
- `package_delivery.py`、`accept_windows.py`：交付三个核心 `.h5ad`、完整 `.dataset`（Arrow 和两个 JSON）及报告；ToDesk 全文件大小/SHA256、实际加载、全部行元数据和长度验收后放入新日期目录。

GSE255296 年龄未逐人确认、跨研究独立性未知；两套新增子集均只有 HCM，不能补造对照。GSE292754 儿童和 GSE138262 完整细胞分别保存；后者疾病与研究/技术混杂。参考尾部和混合程序阈值是项目启发式，未验证疾病域的敏感性/特异性。计数守恒及编码验收通过不代表同意训练或确证独立验证；本批无 DGE、模型训练和模拟敲除。
