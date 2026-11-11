# NeurIPS2023：分模态检查与可用索引

数据为健康人 PBMC 小分子扰动及同一批 3 位供体的 RNA/ATAC 基线。它不是 HCM 或心肌细胞队列。范围为下载清单中全部 8 个原件，原始文件保持只读；不下载训练权重，不训练模型。

- 原件根目录：`/home/zhome/xiaoe-internship-data/benchmarks/competitions/neurips2023/`。
- 执行目录：`/home/zhome/xiaoe-internship-data/_processing/all_data_cleaning_20261111/NeurIPS2023/`。
- 产物目录：`/home/zhome/xiaoe-internship-data/processed/NeurIPS2023_qc_20261111/`。
- Python：`/home/zhome/xiaoe-internship-data/_processing/SCP1303_20261111/venv/bin/python`。
- 最终状态以独立验收、`FINAL_COMPLETION_GATE.json`、退出码和进程实况共同确认，不能以目录存在或作业启动代替。

| 数据单元 | 原生主矩阵 | 含义 |
|---|---|---|
| sc_original | sc_counts.h5ad 的 layers/counts | 301785 个发布细胞的未归一化 RNA 计数 |
| sc_reannotated | sc_counts_reannotated_with_counts.h5ad 的 raw/X | 同批细胞中 298087 个的作者重注释版本 |
| baseline_rna | multiome_counts.h5mu 的 mod/rna/raw/X | 26335 个基线细胞的 RNA 计数 |
| baseline_atac | multiome_counts.h5mu 的 mod/atac/layers/binarized | 同 26335 个细胞的二值峰可及性，不是原始片段计数 |
| de_train / de_test | 两份发布 DE H5AD 的 X | 614 / 255 个细胞类型与药物条件的有符号差异表达统计，不是细胞计数 |

RNA 和 ATAC 的 X 是作者处理后的数值，不能直接叫 raw counts。DE 的负数有效；X 在正 P 值位置逐值符合 `sign(logFC) * -log10(P.Value)`。原五个 DE 层全部保留。`obsm/protein_counts` 的 12 列实际是 CMO_301–CMO_312 细胞标签计数，HTO_clr 是处理后的标签信号，不标成生物蛋白 ADT。

`inspect_neurips2023.py` 先检查来源结构。`run_prepare.sh` 经 nohup 运行 `prepare_neurips2023.py`：核对 8 原件大小、下载基线 MD5/SHA256，完整扫描 37 个原生 assay/obsm 矩阵，检查 CSR、数值域、标识和元数据，建立真实 observation/feature/usable/quarantine 索引。只添加技术统计与用途，不重新归一化、整合、做 DE、判断心肌身份、重做双细胞或套 HCM 500–5000/MT<5% 阈值。全部基因/峰保留，零 DE 响应不剔除。

版本对应表检查共享细胞的每个 RNA 计数行及 CMO 标签计数相同，保留两版标签与作者移出的 3698 个细胞，不把两版相加当新细胞。独立细胞索引以原版六类细胞注释对应官方 DE 条件；重注释版的 T 类合并等变化仍单独保留。官方 control 包含 Dabrafenib 1µM、Dimethyl Sulfoxide 14.1µM 和 Belinostat 0.1µM，药物名、剂量和原 split 全部原样保留。

多组学 RNA/ATAC 的 obsmap、varmap、每个细胞与峰/基因顺序和配对均检查。官方 train/control/public_test/private_test 不改；该划分不能描述为供体独立测试，因为基线和扰动都来自同 3 位供体。

准备退出 0 后才运行 `run_verify.sh`。独立脚本用 AnnData 的 backed sparse reader、不同分块和直接行归约重新读取每个矩阵，复核所有计数/统计、原元数据、逐行/逐基因索引、发布版本对应、多组学映射、官方测试 ID 与药物注释。全部通过才写完成门禁；还须核对 verify 退出 0、进程结束、源码和本地 SHA256 一致。验收不重新认证作者生物学注释、扰动效力或单细胞身份。

首轮独立验收在重注释版 `dispersions_norm` 元数据字段处退出：源 float32 写出为可还原到 float32 的十进制串，旧验收错误地要求它等于扩展后的 float64 精确值。`diagnose_metadata.py` 核实全部字段按源 dtype 还原完全相同；修复后的比较仍逐值精确，不放宽容差。原始矩阵及准备产物未改，仅重跑独立验收，旧源码/日志/失败记录保存在 `verification_attempt1_dtype/`。最终另交 `SOURCE_METADATA_DTYPES.json` 说明源字段类型。

第二轮已通过两版单细胞及 RNA 基线，在 ATAC X 的存储顺序标记比较处退出。原始 CSR 索引未排序但有效；AnnData/SciPy 转换会自动排序，使旧验收读取后的 canonical 标志丢失原始存储状态。`diagnose_matrix.py` 证实所有数值统计一致，修复后直接检查磁盘索引及行边界，再与排序后的稀疏数值扫描交叉检查。不把未排序当坏细胞，不修改原文件。旧记录保存在 `verification_attempt2_matrix_summary/`。

Git 仅收源码和文本汇总。原件、逐细胞/逐特征表和完整检查记录留服务器；本地交付为 `deliverables/all-data-cleaning-20261111/NEURIPS2023_REPORT.md`、`NEURIPS2023_OUTPUT_MANIFEST.tsv` 及 `neurips2023/` 汇总/执行记录。
