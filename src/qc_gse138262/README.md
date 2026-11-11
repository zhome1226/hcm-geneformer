# GSE138262：作者修正版完整心肌细胞计数清洗

作者：Wehrens et al., Cell Reports 2022，DOI 10.1016/j.celrep.2022.110809。作者仓库明确纠正GEO旧TSV并非论文计数表；本流程使用作者指向的Zenodo概念DOI 10.5281/zenodo.18206553对应版本10.5281/zenodo.18206554，保留GEO旧下载原件。

服务器输入：`~/xiaoe-internship-data/validation/hcm_GSE138262/author_corrected_zenodo18206554/`；执行：`_processing/all_data_cleaning_20261111/gse138262/`；输出：`processed/wehrens_GSE138262_qc_20261111/`。

## 顺序

1. `download_corrected.py`通过curl续传5个作者文件，核对Zenodo字节数/MD5、ZIP CRC与解包路径。`run_download.sh`使用nohup、flock、PID和退出码；MacOS资源叉不是计数表。已通过校验的文件复用。
2. `inspect_author.R`使用R原生读取两个Seurat RDS，保留RNA counts、全元数据和作者降维坐标；可用的额外counts assay也导出。发布对象实际只保留RNA assay，线粒体计数只能从完整修正表或作者已保存元数据检查，不能对裁掉MT后的RNA层重算MT%。不以data/scale.data替代counts。
3. `stage_gse138262.py`解析作者使用空格分隔但扩展名为TSV的32个计数表，核实作者filename/GSM/plate/patient/prefix关系。JE03/JE04对应的文件次序沿用作者实际ID定义，不能按编号猜配。作者未提供symbol的Ensembl行保留；原生轴及并集逐样本存在掩码保留。每表与两个Seurat对象的共享RNA计数逐值核对。导出完整未归一化计数、逐细胞数值审计和参考适配输入。
4. `run_identity.R`使用现有成人左心室单核参考的供体分隔训练/校准集，对严格候选进行SingleR注释。分选CEL-seq2和预选iCell8完整细胞不套用10x装载量双细胞模型，不运行EmptyDrops。
5. `calibrate_mixed.py`复算既有固定参考的强混合程序阈值；`finalize_gse138262.py`保留500≤检出基因≤5000、MT<5%的项目严格版，联合参考尾部和多标记定义核心。来源研究分别建立诊断图，小于50个严格候选不以邻域挽救身份；混合信号留待定，不冒充实验确证双细胞。所有数值规则外细胞另存，计数原值不改。
6. `verify_gse138262.py`用独立逐行解析器重读全部32个原表，对native、并集、所有最终计数对象逐值核验；复查53个旧原件、5个修正下载、48组作者RNA比较、MT元数据、疾病/供体关系、参考供体隔离、标记/阈值/身份规则、诊断变换与本地/远端源码SHA256。全部通过才写FINAL_COMPLETION_GATE.json。

阶段驱动为`run_preflight.sh`、`run_stage.sh`、`run_identity.sh`、`run_finalize.sh`、`run_verify.sh`。恢复前核实PID和实际命令；不能因SSH观察超时重启。预查中发现的元数据编码、R脚本解析及SciPy索引兼容问题已修正，历史失败日志保留；完成的R导出和下载不要重复运行。每次执行验收前更新`CODE_EXPECTED_SHA256SUMS.txt`并同步实际源码。

## 适用限制

HCM：5供体、16板，CEL-seq2分选完整细胞。健康伴随：Wang 7供体、16个iCell8板/part，作者预先选择LV群并重新处理；绑定同供体多板和part，保留其研究/技术来源。Seurat ALL额外含Litvinukova参考细胞，仅单独审计，不作为新增独立健康队列并入。原始7537细胞、28663基因完整保留。

严格版仅8个数值候选，6个心肌核心（2位HCM供体）、2个待定，没有健康对照。不得因清洗验收完成就投入HCM/Control微调或疾病验证。统一MT<5%规则对完整心肌细胞的筛除极重，不等于所有排除者都是坏细胞；本轮不自行放宽。完整计数可供后续另定完整细胞方案时比较。诊断log1p(CP10k)对象与原始计数对象分开；未插补、整合、裁切交付基因、差异表达、分配训练测试或进行模型训练。矩阵和逐细胞表仅保留服务器，Git仅存源码及文本汇总。
