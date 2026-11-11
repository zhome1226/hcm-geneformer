# VCC2025：发布计数检查与官方划分索引

2026-11-11 15:08 UTC 独立验收完成，1141 项检查全部通过，准备与验收均退出 0，相关进程已结束。服务器交付目录：`/home/zhome/xiaoe-internship-data/processed/VCC2025_qc_20261111/`。

## 数据与结果

来源为 [Arc Virtual Cell Challenge 官方说明](https://github.com/ArcInstitute/arc-virtual-cell-atlas/blob/main/virtual-cell-challenge/README.md) 和 [实验设计说明](https://arcinstitute.org/news/behind-the-data-virtual-cell-challenge)。原件位于服务器 `benchmarks/competitions/vcc2025/`，7 文件：三份 H5AD、三份扰动细胞数 CSV 和 gene_names.csv。H1 人胚胎干细胞 CRISPRi 扰动数据用于虚拟细胞方法学习与基准测试，不能作为 HCM 患者/健康心肌验证队列。

| 官方划分 | 发布细胞记录 | 通过技术检查 | 扰动靶基因 | 重复包含的阴性对照 | 非对照细胞 |
|---|---:|---:|---:|---:|---:|
| train | 221273 | 221273 | 150 | 38176 | 183097 |
| validation | 98927 | 98927 | 50 | 38176 | 60751 |
| test | 170846 | 170846 | 100 | 38176 | 132670 |

三个文件共 491046 条记录，38176 个 non-targeting 细胞在每个划分中重复发布。按标识去除这些重复后，共 414694 个唯一细胞，其中 376518 个是扰动细胞；阴性对照唯一索引指向 train 副本。扰动靶基因及非对照细胞标识跨划分不重复。不能把重复对照计为三套独立观测，也不能将 non-targeting 标成健康患者。

全部 18080 个发布基因保留，顺序与 gene_names.csv 逐列相同；训练文件原生 gene_id 仅在完全相同的基因轴上提供可追溯映射，不猜补基因。train/validation/test 各有 4/48/12 个全零特征，保留并标记，避免按测试表达选择特征。零信息细胞和缺少关键元数据的细胞均为 0，因此没有为了“清洗”制造删除。

## 如何使用产物

`DATASET_INDEX.tsv` 指向经过完整验证的原始 H5AD 及各划分 usable_index；`*_cell_QC.tsv.gz` 保存原始行位置、细胞标识、target_gene、guide_id、batch、UMI、检出基因、面板内 MT 描述值、技术检查与用途标记。`*_held_out_perturbation_index.tsv.gz` 精确索引 validation/test 的扰动细胞，排除共享阴性对照。`shared_controls_unique_index.tsv.gz` 提供唯一参考对照；`feature_index_QC.tsv.gz` 保存完整特征轴及各划分统计。

这是“保留并验证发布方已清洗矩阵 + 建立可用索引”的交付，不复制约 34 GB 的矩阵，不另造 H5AD 占位。读取原始矩阵时按索引的 source_row 提取，并保留官方划分。逐细胞、逐基因表与完整检查记录仅在服务器保存，本地仅留汇总、规则、摘要与执行记录。

## 检查与适用边界

准备和独立验收分别扫描每个存储值，检查 CSR 结构、索引范围、有限性、非负与整数性，复算每个细胞/基因的统计；原文件 X 虽以 float32 存储，逐值检查的读数为非负整数，未凭 dtype 直接判作标准化表达。三套最大存储读数为 2365/2083/2406，存储零值及非规范稀疏块均为 0；无 raw 或额外 layers。全部 300 个扰动靶点的细胞数、UMI 中位数与官方 CSV 对照。共享对照需逐标识、原元数据与排序后的稀疏计数行 SHA256 一致，而不只是名字相同。

保留发布方已有 QC，不重建 FASTQ/液滴识别，不额外运行心肌身份判断或双细胞过滤；不使用 HCM 500–5000 基因、MT<5% 规则。面板内 MT 仅描述。没有根据敲低程度、疾病信号或测试标签调阈值，未做归一化、插补、批次整合、DE、tokenization 或模型训练。结构和元数据检查通过不等于独立认证每个细胞的单细胞身份或生物学敲低效力。

7 个原件大小和 MD5 均与下载基线一致；3 个 H5AD 已有历史 SHA256，可再次比对；另 4 个元数据文件本轮补记 SHA256，不能声称与缺失的历史值比对。失败首轮及恢复记录保留在执行目录 `attempt1_readonly_mask/`：首轮 train 扫描完成后，Pandas 返回只读布尔数组使 `valid_meta &= ...` 报错退出 1。修复为显式可写副本；回归检查通过后重新准备。原件未修改，不复用失败作业充当完成记录。

执行源码在 `D:\VirtualCell\sources\qc_vcc2025\`，7 文件的本地/远端 SHA256 已核对相同。准备退出 0、1141 项独立检查通过、verify 退出 0、进程结束和 `FINAL_COMPLETION_GATE.json` 均已核实；验收记录和产物清单下载后再次核对 SHA256。总表可将 VCC2025 记为 complete_verified；全部数据仍有其他待完成项。本轮未提交 Git，按交接约定在全部数据逐套验收后统一核查文本/源码交付。
