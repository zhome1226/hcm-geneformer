# Chin 两套队列清洗

范围：GSE174691 的9位HCM供体和GSE161921的4位正常供体，室间隔snRNA-seq。固定数据根为服务器 `/home/zhome/xiaoe-internship-data/`；原始文件保持不变。

正式入口为 `run_chin.sh`，依次完成原始计数检查与暂存、EmptyDrops、各捕获文库scDblFinder、参考注释、联合心肌身份判断、导出及独立验收。长任务使用nohup与PID文件；恢复前先核实具体进程是否仍运行，不能只看旧退出码。

依赖：既有SCP1303 Python环境、all_data_cleaning_20261111/r_env、GEO官方样本与系列元数据、Litvinukova2020成人左心室单核参考的固定训练/供体留出矩阵。本队列按共同Ensembl基因重新训练SingleR，不把疾病标签用于身份模型。版本与参考适配结果记录在输出中。

## 步骤与文件

- `stage_chin.py`：逐原始文件SHA256、完整液滴QC；10x H5暂存保持原始计数。`prefilter.h5ad`仅作内存友好暂存，≥200基因不等于正式通过QC。
- `run_chin.R`：EmptyDrops lower=100、FDR≤.001，必要时增加Monte Carlo次数；双细胞模型使用已识别且≥500基因的完整多类型捕获文库，保留高MT和>5000基因输入；SingleR注释严格QC候选。
- `finalize_chin.py`：严格500–5000基因、MT<5%；多证据疑似双细胞排除，心肌核心与待定分开；每研究独立保存全基因、未归一化X/counts及≤20k唯一轻量核。
- `label_diagnostic_metadata.py`：仅修正归一化诊断副本的元数据说明，并核对矩阵摘要不变；不修改原始计数对象。
- `verify_chin.py`：重新读取13个原始GEO矩阵，逐值比较全部交付计数，核对元数据、条码、完整捕获模型输入、QC、归一化副本的声明和值、轻量/完整关系及摘要。

输出：`processed/chin_GSE161921_qc_20261111/`、`processed/chin_GSE174691_qc_20261111/`。`diagnostic_normalized.h5ad`是明确标注的log1p(CP10k)分析副本；正式计数对象是`all_cells_qc.h5ad`、`cm_qc_full.h5ad`、`cm_qc_light.h5ad`、`cm_review_pending.h5ad`及`suspected_doublets.h5ad`。

本轮旧原型曾用≥200基因进入双细胞模型，在低复杂度样本中抬高了按捕获量推算的预期率。该评分仅保留在`doublets.tsv`中作历史记录，不用于任何最终过滤。`revise_doublets.R`和`run_chin_v2.sh`用于本轮从旧原型恢复；最终仅使用`doublets_v2.tsv`，严格筛选阈值未改变。

大队列诊断邻域明确使用sklearn精确近邻，以避开当前UMAP/Numba对近似后端非连续数组的兼容问题。已完成的较小对照队列保留其有效邻域结果；病例/对照未合并表达矩阵。近邻与标记证据相关，不作为独立软件投票。

以`INDEPENDENT_VERIFICATION.json`的当前passed状态及真实文件为验收依据；启动进程、导出成功标记或历史通过记录都不能代替最终验收。待定身份不声称被实验确证。跨研究供体重叠、最终训练划分及模型性能不在本次清洗中解决。
