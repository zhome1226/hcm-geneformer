# GSE227734 准备与阻塞审计

本入口不宣称完成清洗。官方发布矩阵有338个细胞、9760个匿名列，缺少真实基因标识。所有列名仅为0–9759，原HDF5无其他映射。未取得可靠映射前，MT比例、标记/参考身份及基因模型输入无法建立。

`prepare_gse227734.py` 保留全部原计数与作者元数据，按GEO把5个样本别名绑定供体与样本；补充数值审计并明确所有细胞不可用于训练。当前AnnData API的None键是X的别名，不能误认为源文件有额外计数或基因注释；直接读取原始HDF5的X/obs/var，只额外存储准确counts层。

`verify_gse227734.py` 重新解压原gzip到内存，逐值比较真实计数、细胞/列顺序与元数据，重算数值审计并检查未虚构MT、基因名、健康组或完成状态。准备审计通过不等于数据清洗完成。

服务器：`/home/zhome/xiaoe-internship-data/`。原件 `validation/cox7b_GSE227734/`；执行目录 `_processing/all_data_cleaning_20261111/GSE227734/`；待注释产物 `processed/cox7b_GSE227734_qc_20261111_pending_annotations/`。长任务使用nohup运行 `run_preparation.sh`，保存PID、日志和退出码；恢复前核实进程。

新增 `audit_publisher_supplement.py`：只读取官方补充Word及原HDF5，核对条码、供体细胞数和临床信息；9项来源核查通过。补充Table S4的390→338及逐供体数量与原矩阵一致；P5年龄GEO/矩阵63与补充表64的冲突保留在`source_evidence/`。补充6张小表没有9760列顺序映射，清洗状态仍为缺基因注释；审计不改变计数或训练资格。补充原件在`references/GSE227734_publisher_20261111/`，详细查找见交付`gse227734/SOURCE_LOOKUP_SUMMARY.md`。

STRT-seq成人分离心肌细胞不套10x液滴流程。重复8碱基条码不代表重复细胞，Large/Small形态不代表Control/HCM，原论文小鼠Cox7b干预不代表这些人细胞受干预。计数/逐细胞表仅在服务器保存，Git仅保留代码与文本汇总。
