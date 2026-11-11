# ACDC2017：心脏MRI与分割标注检查

范围是下载清单全部904个原件：150位患者的750个NIfTI、150份Info.cfg及4份许可/引用文件。训练患者001–100、测试患者101–150，原分组DCM/HCM/MINF/NOR/RV和官方划分保留，不转换成单细胞疾病标签。

- 原件：`/home/zhome/xiaoe-internship-data/benchmarks/competitions/acdc2017/`。
- 执行目录：`/home/zhome/xiaoe-internship-data/_processing/all_data_cleaning_20261111/ACDC2017/`。
- 交付目录：`/home/zhome/xiaoe-internship-data/processed/ACDC2017_qc_20261111/`。
- Python：`/home/zhome/xiaoe-internship-data/_processing/SCP1303_20261111/venv/bin/python`；本次增加NiBabel 5.3.2，安装日志保留在执行目录。

`run_prepare.sh`经nohup运行`prepare_acdc.py`：904原件分别对照下载大小/MD5/SHA256；每个gzip完整解压验证CRC后通过NiBabel从内存字节读取，检查全部像素、尺寸、仿射、正体素间距、数据类型和分割标签。分割0=背景、1=右心室腔、2=心肌、3=左心室腔。ED/ES和文件帧编号从1开始，4D数组索引减1；每个相位图像与原4D帧逐值比较，图像/标注/序列必须几何一致。

输出真实病例、全文件、300个相位图像-标注配对、可用与隔离索引；每个像素数组记录SHA256，原文件保持只读。异常对象不填占位值，不重采样或改动原强度及标注。不同诊断组仍按原始组保留，影像不冒充RNA计数。

准备成功后`run_verify.sh`执行独立`verify_acdc.py`，采用从原始路径打开gzip文件的另一条NiBabel读取路径，重新读取所有750个完整体积并比较像素SHA256、几何、每类标注体素数、原相位、300个配对与全部可用/隔离规则；904个文件摘要再次核对。通过后再以`collect_completion.py`核查prepare/verify退出0、PID结束、源码本地/远端SHA256及产物摘要一致。

交付为已验证的原影像/标注和可直接读取的索引，不复制大体积或训练模型。原CC BY-NC-SA 4.0许可与强制引用一并保留。Git仅存代码和文本汇总，不提交NIfTI影像或分割文件。

实际检查发现72对（36位患者）的原3D图像/标注彼此对齐、像素精确等于原4D对应帧，但3D头与4D头仿射不同。原严格空间几何索引228对和72条隔离记录保留；`run_geometry.sh`在严格独立验收后运行`geometry_views.py`，再从原路径读取全部600个3D文件，增加300对体素网格分割可用索引、2100项检查及逐对几何冲突证据。两种用途分开：体素网格可用不等于物理空间信息已修复。原头信息不改，额外几何门禁、摘要与退出码也纳入最终收集器。
