# VCC2025 发布数据检查与索引

对象是 Arc Virtual Cell Challenge 2025 的 H1 人胚胎干细胞 CRISPRi 扰动数据。它不是 HCM 患者队列；non-targeting 是实验阴性对照，不能标成健康心肌细胞。

- 原件：`/home/zhome/xiaoe-internship-data/benchmarks/competitions/vcc2025/`，7 文件，包括 train/validation/test 的 H5AD、3 张 perturbation counts 表和 gene_names.csv。
- 执行目录：`/home/zhome/xiaoe-internship-data/_processing/all_data_cleaning_20261111/VCC2025/`。
- 交付目录：`/home/zhome/xiaoe-internship-data/processed/VCC2025_qc_20261111/`。
- 使用已安装的 `_processing/SCP1303_20261111/venv/bin/python`，无需模型权重。

先核对作业与退出码，避免重复运行。`inspect_vcc2025.py` 仅做元数据检查。`run_prepare.sh` 在后台运行 `prepare_vcc2025.py`，复查 7 个原件的大小、MD5、可用的历史 SHA256，再逐块检查所有计数和 CSR 结构，生成逐细胞、逐基因、扰动及官方划分索引。原始 H5AD 保持只读，不制作重复大矩阵，不修改作者标准化结果，也不增加 HCM 标签。保留发布方已有 QC；不套用心肌细胞 500–5000/MT<5% 规则。

三个官方文件重复包含 38176 个 non-targeting 细胞；检查它们的标识、元数据和规范化稀疏计数行 SHA256，唯一对照索引指向 train 副本。扰动与官方划分保持不变，不能把重复阴性对照当成独立测试结果。

准备退出 0 后，`run_verify.sh` 重新计算原件摘要、扫描每个计数、用不同分块和行列汇总方法核对所有 QC 数值，并检查细胞/基因轴、原始元数据、逐行索引、扰动集合和对照重复。全部检查通过才写 `FINAL_COMPLETION_GATE.json`；还需确认 `verify_exit_code.txt=0`、进程已退出、源码 SHA256 与本地一致，才能更新总状态。数据验收不证明扰动效力、单细胞身份或 HCM 适用性。

首轮准备已完成 train 扫描，随后在 `valid_meta &= ...` 处因 Pandas 返回只读数组而退出 1，未导出完成产物。旧日志、源码及 PID/退出记录保存在 `attempt1_readonly_mask/`。修复仅为布尔掩码显式 `copy=True`，`check_recovery.py` 验证可写性及缺失值处理。原件未改动；恢复需重新完成准备与独立验收，不能把首轮失败当作完成。

Git 仅收源码、规则、汇总和摘要；逐细胞/逐基因表、完整检查记录及原始矩阵保留服务器。最终结果见本地 `deliverables/all-data-cleaning-20261111/VCC2025_REPORT.md` 与 `VCC2025_OUTPUT_MANIFEST.tsv`。
