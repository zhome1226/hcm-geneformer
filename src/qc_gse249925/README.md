# GSE249925 bulk 表达清洗

服务器根 `/home/zhome/xiaoe-internship-data/`；原件 `validation/bulk_GSE249925/`；执行目录 `_processing/all_data_cleaning_20261111/GSE249925/`；交付 `processed/bulk_GSE249925_qc_20261111/`。

1. `inspect_gse249925.py` 复查三原件摘要、完整矩阵数值/轴、SOFT样本映射和来源处理语义，保存 INPUT_INSPECTION.json 及来源检索记录。
2. `prepare_gse249925.py` 精确绑定120个样本的GEO标题、GSM、BioSample、SRA与疾病标签。保存所有原发布表达值和原生复合基因轴；可用矩阵仅移除恒零或结构无效特征，不聚合重复Ensembl/symbol，不按病例组筛基因。缺失donor/age/sex/batch明确留空或标未知。
3. `verify_gse249925.py` 用CSV流与独立SOFT块解析重读原件，逐值比较完整/可用H5AD与TSV，核对每个GEO分组和注册关系、全部QC规则、文件摘要及源码SHA256。验收生成最终报告/门禁；只有退出码0且门禁与实际文件一致后才将总表记完成。

`run_prepare.sh`、`run_verify.sh` 使用文件锁、PID、日志与退出码，通过nohup启动。恢复前先核实PID对应命令、退出码和产物，不重复运行成功的准备阶段。

旧下载基线只有大小和MD5，SHA256字段为空；首轮验收将空值当校验值而失败，源码、清单和日志留在 `verification_attempt1_missing_SHA/`。修正后只重跑验收，原清洗产物未改；新测SHA256与“无旧SHA基线”的标记一并进入最终记录。

GEO原始计数标签与TMM说明冲突尚未解决：整数值不能证明是未归一化计数。输出用 `published_values` 层而非 `counts`，`certified_unnormalized_counts=False`。来源证据在交付的 source_lookup/；这项来源限制不能因数值验收成功而省略。

本套是组织bulk表达，不是单细胞心肌；不套用单细胞基因数/MT阈值、EmptyDrops、双细胞检测或心肌提取。不运行标准化、差异表达、tokenization或训练。训练前仍须解决计数来源与供体/临床对应问题。

原件、大矩阵及逐特征表仅保留服务器；Git仅保存源码与文本汇总。Python复用既有venv，无需新增模型或工具。
