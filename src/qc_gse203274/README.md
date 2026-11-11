# GSE203274 儿童心脏单核清洗

来源：Hill et al., *Nature* 2022，*Integrated multi-omic characterization of congenital heart disease*，DOI `10.1038/s41586-022-04989-3`，GEO GSE203274。

服务器数据根 `/home/zhome/xiaoe-internship-data/`；原件 `validation/pediatric_GSE203274/`；执行目录 `_processing/all_data_cleaning_20261111/GSE203274/`；产物 `processed/pediatric_GSE203274_qc_20261111/`。

## 处理与恢复

1. `stage_gse203274.py` 从两个已下载原件核对大小/MD5，暂存全部 157293 核、29266 个源基因及作者元数据。原件是 CellBender 与作者 QC 后发布的未归一化计数，不是完整原始液滴。按 labID 绑定 13 位供体、30 个捕获；保留 LV/RV、年龄冲突及候选 GSM 集合。其他心脏病 label=-1，不转换成健康对照。
2. `repair_stage_metadata.py` 修复暂存布尔元数据类型，逐个稀疏计数组件 SHA256 证明计数未变；原暂存摘要留在 `stage_attempt1_metadata_dtype/`。后续步骤使用修复后的暂存对象，不重跑已成功的暂存任务。
3. `run_identity.R` 在源与成人 LV 参考均唯一、完全同名的 18393 个符号上运行 SingleR；参考训练/留出供体分隔。派生参考 ID 与 hg19 原生基因身份区分，所有源基因保留；成人参考用于儿童存在域差异。
4. `calibrate_mixed.py`、`prepare_residual_review.py` 和 `run_residual.R` 复核多谱系信号。作者已有多轮混合簇移除，只有每捕获强混合≥10 且≥1%才触发补测 scDblFinder；未触发者明确记为未重跑，不虚构分数或零双细胞率。
5. `finalize_gse203274.py` 保持项目严格基线 500≤检测基因≤5000、MT<5%，结合作者 CM、SingleR、多标记与邻域保存核心、待定及原计数。二分类核心只含 HCM/Control，其他疾病心肌另存。诊断标准化副本与计数对象分开。
6. `verify_gse203274.py` 独立重读原 tar 的三个成员及原始矩阵，逐值核对全部计数与九个输出，重算 QC、疾病/供体/年龄、身份与残余模型决策；检查供体/临床元数据和完整基因轴，核对源码 SHA256。只有验收成功、`verify_exit_code.txt=0` 且最终门禁存在才能记为完成。

各入口 `run_*.sh` 使用文件锁、PID、日志及退出码；长作业通过 nohup 运行。重启前同时查看 PID 对应命令与退出记录，不能仅凭旧日志、PID 文件或观测超时重跑。

## 使用边界

只有一位 HCM 供体 P86，两次捕获不能当作两个患者。儿童/成人、部位、年龄、跨研究重复及供体划分须在训练前复核。训练与 tokenization 尚未进行；清洗验收不等于批准直接训练或证实疾病泛化能力。

所有矩阵和逐核表仅保留服务器；本地及 Git 保存代码、来源证据和文本汇总。环境复用现有 Python venv、R 环境与心脏参考。
