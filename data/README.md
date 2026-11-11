# 数据

- 主训练来源：[SCP1303 / Chaffin et al., Nature 2022](https://singlecell.broadinstitute.org/single_cell/study/SCP1303)，15 HCM、11 DCM、16 非衰竭心脏，下载范围、原始计数层和完整性清单见 [training/README.md](training/README.md)，任务背景见 [#1](https://github.com/zhome1226/xiaoe-internship/issues/1)。旧推荐 GSE145154 无 HCM，GSE130652 是小鼠 ChIP-seq，均不用于本任务。
- 已下载的验证数据见 [validation/README.md](validation/README.md)，新增候选的适用范围和排除原因见 [validation/CANDIDATES.md](validation/CANDIDATES.md)。
- 可学习的比赛数据、获奖代码及评估工具见 [benchmarks/README.md](benchmarks/README.md)；与 HCM 临床验证数据分开管理。
- 后续公开数据补充批次见 [expansion/README.md](expansion/README.md)：约 83.3 GB 的完整比赛数据、ACDC、GSE203274/GSE138262；是否完成以其中 STATUS.json 为准，计划清单不等于已下载。
- SCP1303 严格基线清洗已完成：[结果、服务器路径和独立验收](../results/qc/SCP1303_20261111/README.md)。完整 114,927 核，轻量 20,000 核，全部保留 36,601 基因。
- 本轮规则为 500≤检测基因数≤5000、原计数线粒体百分数<5%；这是本项目基线。复用作者背景/双细胞QC和心肌类型，以多标记核验；标签 `NF=0, HCM=1` 来自供体诊断。划分暂定，跨研究供体重叠未核实。
- 实际原件位于服务器 `/home/zhome/xiaoe-internship-data/`，质控对象位于其中 `processed/SCP1303_qc_20261111/`。矩阵与逐核审计表不入库，下载/清洗脚本、元数据和小型报告可以提交。
