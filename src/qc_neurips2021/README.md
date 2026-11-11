# NeurIPS2021：多模态整理、可用索引与独立验收

来源为 [GEO GSE194122](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE194122) 发布的两个经过作者质控的健康骨髓 BMMC 对象。它们不是 HCM 队列；保留作者已有双细胞去除，不套用 HCM 的 500–5000/MT<5% 规则。

- 原件：`/home/zhome/xiaoe-internship-data/benchmarks/competitions/neurips2021/`。
- 执行与展开件：`/home/zhome/xiaoe-internship-data/_processing/all_data_cleaning_20261111/NeurIPS2021/`。
- 可用索引及全部审计：`/home/zhome/xiaoe-internship-data/processed/NeurIPS2021_qc_20261111/`。
- Python：`/home/zhome/xiaoe-internship-data/_processing/SCP1303_20261111/venv/bin/python`。

`stage_inspect.py` 已完成2原件大小、MD5、SHA256与下载基线核对，读到gzip结尾确认CRC后解包并保留压缩原件。stage成功退出0后不重新解包。

`run_prepare.sh` 在nohup下执行准备：复查压缩原件和展开件摘要，完整扫描15个主矩阵及附加矩阵，以原件/展开件和真实行列索引交付，避免重复保存数GB矩阵。`ASSAY_INDEX.tsv` 给出4种测量的文件、计数位置、行索引与特征选择；逐细胞及逐特征统计均在服务器保留。

| 发布对象 | 配对细胞 | 测量与特征数 | 计数位置 |
|---|---:|---|---|
| CITE-seq | 90261 | RNA 13953；ADT 134 | `layers/counts` 的各自模态列 |
| Multiome | 69249 | RNA 13431；ATAC 116490 | `layers/counts` 的各自模态列 |

全部159510个发布细胞通过计数可用性及身份字段结构检查，两对象共10个不同DonorID。GEO总体设计描述12供体，不等于本次这两个发布对象包含12个供体。X、ATAC gene activity、isotype controls及嵌入是作者处理后的表示，不叫原始计数。重算MT指标仅基于已发布的基因轴，不据此声称恢复最初完整转录组的MT百分数。

源细胞、全部特征、作者标签和原始缺失值保留。CITE里36个名称同时用于RNA和蛋白，涉及72列；以模态加名称及原列号分别标识，不合并。98个ADT和所有ATAC峰的gene_id缺失不猜填；GEX的Ensembl ID完整且唯一。伪时间缺失不插补、不作为删除细胞的理由。

Multiome的`obsm/ATAC_lsi_red`有44170个非有限值，影响36147行，只位于原第37、38列（从0开始计数）。`DERIVED_COORDINATE_QC.tsv.gz`排除包含非有限值的这两列，余下0–36列可按索引取用；`ATAC_lsi_full`的50列均有限。原件原坐标不改，不因这个派生表示问题删除计数细胞。全部异常行在`DERIVED_NONFINITE_ROWS.tsv.gz`可追溯。

配对仅指一个对象同一原始行的RNA/ADT或RNA/ATAC。两个对象来自不同细胞，不能按相同条码文本建立跨对象细胞配对；同DonorID跨站点/模态绑定为一个供体。CITE的原`is_train`含train、test、iid_holdout，逐行保留；本次Multiome对象没有该列，不虚构官方划分，不把它描述成供体独立测试。

准备成功后运行`run_verify.sh`。`verify_neurips2021.py`不导入准备代码，使用AnnData磁盘稀疏读取器、不同块大小和原生CSC列方向重新扫描15矩阵，逐值复算每细胞/每特征计数、配对、身份字段、划分、特征碰撞和全部派生坐标隔离，并复查全部摘要。通过后写`INDEPENDENT_VERIFICATION.json`、`FINAL_COMPLETION_GATE.json`及`VERIFIED_OUTPUT_MANIFEST.tsv`；最终还核查退出码、进程结束与源码一致。

未做新的归一化、插补、整合、差异表达、HCM疾病标注或模型训练。Git只保存源码和文本汇总；逐细胞/逐特征表、矩阵和完整日志留服务器。
