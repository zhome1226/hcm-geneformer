# SCP1303 定向复核代码

本目录保留严格基线之后的身份、混合谱系/残余双细胞复核、导出及独立验收脚本。最终报告与规则在`results/qc/all_data_cleaning_20261111/scp1303_reviewed/`。最终服务器对象是`/home/zhome/xiaoe-internship-data/processed/SCP1303_qc_reviewed_final_20261111/`，114711核心核、20000轻量核、144待定核；82项独立验收通过。

原计数与所有基因保留；500–5000及MT<5%的项目基线不变。`prepare_identity_reference.py`及SingleR脚本使用成人心脏参考并保留训练/留出供体；`prepare_secondary_doublets.py`与R脚本在完整作者捕获构成中补测；`export_reviewed_scp1303.py`、`verify_reviewed_scp1303.py`负责最终导出与逐项验收。参考分数不作为概率，单个参考标签/单基因阴性不直接删核。

`prepare_external_inputs.py`、`inspect_export_liu.R`、`inspect_all_data_inputs.py`等是后续多数据集流程的输入暂存/来源检查依赖。实际数据集清洗与验收入口见相邻`qc_*/README.md`，不将输入暂存当成清洗完成。

服务器既有Python虚拟环境为`_processing/SCP1303_20261111/venv/`，R环境为`_processing/all_data_cleaning_20261111/r_env/`，均相对数据根`/home/zhome/xiaoe-internship-data/`。脚本中的固定路径对应本轮实际运行位置。长任务使用已有nohup启动器、PID及退出码记录，运行前核查活进程；已成功阶段不重跑。Geneformer权重、模型训练和tokenization不是清洗前置步骤。
