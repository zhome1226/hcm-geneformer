# NeurIPS2023：PBMC 扰动、基线多组学和差异表达结果

2026-11-11 15:52 UTC 独立验收完成，383 项检查全部通过；准备和验收均退出 0，相关进程已结束。服务器产物目录为 `/home/zhome/xiaoe-internship-data/processed/NeurIPS2023_qc_20261111/`。本套可记为 complete_verified；全部数据清洗仍有其他未完成项。

## 数据及交付

来源：[OpenProblems 官方任务](https://github.com/openproblems-bio/task_perturbation_prediction)、[NeurIPS2023 Kaggle 比赛](https://www.kaggle.com/competitions/open-problems-single-cell-perturbations)。8 个下载原件共约 22.49 GB，全部大小、MD5、SHA256 与下载基线一致。研究在 3 位健康供体的 PBMC 中测量小分子处理 24 小时后的反应，并提供同供体 RNA/ATAC 基线；不作为 HCM 或心肌细胞二分类数据。

| 数据单元 | 发布行数 | 特征数 | 通过技术检查 | 主矩阵 |
|---|---:|---:|---:|---|
| 原版单细胞 RNA | 301785 | 21265 | 301785 | layers/counts |
| 作者重注释版 RNA | 298087 | 21265 | 298087 | raw/X |
| 配对基线 RNA | 26335 | 17438 | 26335 | mod/rna/raw/X |
| 配对基线 ATAC | 26335 | 135418 | 26335 | mod/atac/layers/binarized |
| DE train | 614 个条件 | 18211 | 614 | X：有符号 DE 分数 |
| DE test | 255 个条件 | 18211 | 255 | X：有符号 DE 分数 |

重注释版是原版的子集，不能算新增 298087 个细胞；RNA/ATAC 基线的 26335 个细胞是配对测量，也不能相加成 52670 个细胞。DE 行是细胞类型、药物及剂量条件，不能当单个细胞。原版 301785 个扰动实验细胞和 26335 个基线细胞使用分别标明实验来源的 UID；最终训练划分仍需按任务另行设计。

交付保留所有经过验证的原始 H5AD/H5MU，`DATASET_INDEX.tsv` 指向真实原件、主矩阵位置及可用行/特征索引。`*_observation_QC.tsv.gz` 保留全部作者元数据及原始行位置、新增技术统计/用途；`*_usable_index.tsv.gz`、`*_quarantine_index.tsv` 给出可用与隔离记录。原件不复制、不覆盖，不制作占位 H5AD。

## 检查与整理规则

完整检查 37 个原生矩阵，覆盖 X、RNA 计数、ATAC 二值层、五个 DE 层、CMO 标签信号、作者评分和降维坐标。分别检查结构、有限性、负数/小数/二值/P 值范围与原始轴，重算每个主矩阵的逐行和逐特征统计。保留全部发布基因与峰；单细胞两版各有 6 个全零 RNA 特征，其余单元无全零特征。所有单元零信息行与新增技术隔离均为 0；这表示沿用已通过检查的发布方清洗结果，不证明每个细胞已获独立生物学认证。

RNA 计数从指定 counts/raw 槽读取；小数表达 X 继续按作者处理后表达保存。ATAC X 为处理后可及性，binarized 是二值峰，不冒充未提供的原始片段计数。`protein_counts` 实为 CMO_301–CMO_312 的 12 个细胞标签，HTO_clr 可有负数，不写成生物学蛋白面板或疾病标签。

DE 的 15825359 个 X 值分别与原 P.Value、logFC 对照，所有 P>0，无须猜补下限；公式为 `sign(logFC)*-log10(P.Value)`。train/test 最大绝对公式误差约 2.84e-14/1.42e-14。原 is_de / is_de_adj 与各自 P<0.05 的逐值差异均为 0；只是核对原发布结果，未重新计算 DE。保留官方 train/control/public_test/private_test，869 个 DE 条件全部可匹配原版六类细胞的条件记录，官方测试 255 行 ID 映射保持对应。

## 注释、版本和配对

共享的 298087 个细胞要求 RNA 计数行 SHA256 以及全部 12 个 CMO 计数逐值相同；原始供体、药物、剂量、时点、板/孔和 split 均相同。作者移出的 3698 个细胞仍保留在原版及对应表中，不将其自动断言为坏细胞。两版注释的变化单独保存，例如 38865 个原 NK 标签在作者重注释版成为 T cells；这属于发布版本差异，本轮不以某版标签作重新训练的判决。原版六类标签用于对应官方 DE，重注释版四类标签不强行覆盖它。

`multiome_paired_index.tsv.gz` 逐行记录基线 RNA/ATAC 配对，原 obsmap/varmap 以及全部峰/基因位置均已核对一致。作者原 QC 指标继续保留；新增指标清楚标为按发布基因/峰重算。26167 个基线细胞的原 total_counts、n_genes 与当前发布基因轴重算不同，而原 n_counts 与本轮全部重算总量一致；这些差异单独记录，不覆盖来源统计或据此删除细胞。ATAC 的原 n_features_per_cell 与二值层重算逐行一致。

control 是发布方的实验划分：Dabrafenib 1µM（26346 细胞）、Dimethyl Sulfoxide 14.1µM（28462）和 Belinostat 0.1µM（24059），合计 78867。不会统一改成 DMSO 或健康心肌细胞。基线与扰动共享 Donor 1–3，官方划分不是独立供体测试。MOA 表全部 147 行保留，缺失机制/引文及重复标识只标记，不猜补化学或生物信息。

## 验收与边界

独立验收重新读取 8 原件、扫描全部矩阵数值、复核原元数据、每个细胞/条件/特征统计、版本关系、多组学配对、DE 语义、药物注释和导出摘要。383 项检查通过，prepare/verify 退出 0、进程终止、8 个源码文件本地/服务器 SHA256 一致及完成门禁均已核实；验收记录与产物清单下载后再次核对 SHA256。重注释版 cell_type_orig 与原版六类注释逐行一致，两版本标签变化的 185505 个共享细胞完整可追溯。

首轮验收在作者 float32 的 dispersions_norm 字段比较处退出。CSV 十进制串可以逐值还原到原 float32，旧脚本却要求它精确等于扩展后的 float64 数值；诊断确认全部作者字段按源类型恢复相同。修复为按原 dtype 精确比较，未放宽数值容差、未修改数据；只重跑独立验收，首轮源码/失败日志归档于 `verification_attempt1_dtype/`，源字段类型另记 `SOURCE_METADATA_DTYPES.json`。

第二轮在 ATAC X 的存储顺序标记处退出；所有数值统计一致。原始 CSR 列索引未排序，而 AnnData/SciPy 读取转换会排序，旧脚本错误地把转换后的状态与磁盘状态直接比较。修复为另行逐块检查原生索引的行内顺序，仍完整复算数值；未排序 CSR 是有效可读格式，不因此删除细胞或改写矩阵。第二轮源码/失败记录保存在 `verification_attempt2_matrix_summary/`。

不套 HCM 500–5000/MT<5% 或心肌身份标准，不重做 EmptyDrops、双细胞筛除、归一化、插补、整合、HVG 裁切、DE、tokenization 或模型训练。数据验收不重新认证作者细胞身份、扰动效力、供体生物学独立性或 HCM 适用性。原件、逐细胞/逐特征结果留服务器，Git 只交代码和文本汇总；按交接在全部数据逐套验收后统一核查 Git 交付。
