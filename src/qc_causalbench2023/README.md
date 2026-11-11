# CausalBench2023：K562 / RPE1 原始计数整理

仅整理下载清单中的两份Replogle Perturb-seq原件，不进行网络推断或训练。官方仓库固定提交为`1a2143cffdc85f835b41ce8d52034be1bf903e71`，其下载代码引用的Figshare文件ID为35773219、35775606，与本地下载基线一致。原官方`create_dataset.py`、`preprocessing.py`、`splitting.py`和`download_weissmann.py`从已下载代码归档读取并保留来源摘要。

- 原件：`/home/zhome/xiaoe-internship-data/benchmarks/competitions/causalbench2023/data/`。
- 执行目录：`/home/zhome/xiaoe-internship-data/_processing/all_data_cleaning_20261111/CausalBench2023/`。
- 产物：`/home/zhome/xiaoe-internship-data/processed/CausalBench2023_qc_20261111/`。
- Python：`/home/zhome/xiaoe-internship-data/_processing/SCP1303_20261111/venv/bin/python`。

官方流程先读X再归一化，本次X完整扫描要求有限、非负、整数值。原件只读，通过原文件加真实可用行列索引交付，全部已发布基因保留，不按基因symbol合并Ensembl列。完整细胞/基因统计、扰动和捕获关系留服务器，Git只存代码和文本汇总。

K562包含310385细胞×8563基因，RPE1包含247914细胞×8749基因。两者为不同细胞系，不伪造供体身份或HCM标签；同名条码及gem_group按细胞系命名。non-targeting必须同时按gene和gene_id的精确值识别；NANOG是实际扰动靶基因，不因包含“nan”字样变成对照或缺失。sgID_AB的两个sgRNA保留，不当作双细胞证据。目标未出现在发布表达基因轴的细胞保留并明确标记，不能把目标硬配到其他基因列。

不套HCM基因数/MT阈值，不重跑已有实验质控，不归一化、插补或做新的扰动强度选择。原mitopercent字段单位未在当前来源中独立证实，保留原值；重算MT仅基于发布基因轴并注明范围。

两原件均没有固定行级训练/测试划分。已归档官方参数化配方：在任务预处理后按interventions分层，test_size=0.2、random_state=0。官方还可选择strong-perturbation过滤、目标表达分位筛选及>100细胞的网络基因选择；这些属于后续基准任务，不套到通用计数保存。本次不虚构固定官方split，不产生新的模型训练划分。

`run_prepare.sh`经nohup运行`prepare.py`，核对原件大小/MD5/SHA256，完整读数值、逐细胞/逐基因统计并输出可用/隔离行索引。`run_verify.sh`只在准备退出0后运行`verify.py`，用HDF5 read_direct、不同块大小独立重读两矩阵；逐值复算所有统计和元数据，复查扰动/捕获映射、官方来源及产物摘要。最后由`collect_completion.py`要求prepare/verify退出0、对应PID已消失、本地/远端源码SHA256一致及最终产物摘要一致，才形成最终执行证据。
