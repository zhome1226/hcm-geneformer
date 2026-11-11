# SCP1303 严格基线清洗

实际结果和参数见 [结果报告](../../results/qc/SCP1303_20261111/README.md)。本流程只筛选/抽样作者已处理的计数，不做训练、DEG、整合或 tokenization。

## 文件

- `prepare_input.py`：对照下载清单核验源 gzip 的大小和 SHA-256，原子写入独立解压副本。
- `bootstrap_scp1303.sh`：建立独立环境并检查输入；首次使用要求工作目录已有本目录脚本及 Scanpy skill 的 `skill/inspect_data.py`、`skill/_common.py`。
- `inspect_metadata.py`：形状、元数据与计数层概览。
- `qc_core.py` / `test_qc_core.py`：阈值、供体分层划分、无重复抽样及核心测试。
- `run_scp1303_qc.py`：成对计数层质控、作者注释复用、多标记诊断、完整及轻量 H5AD、逐核审计和哈希。
- `verify_scp1303_qc.py`：独立重建源筛选集合，重新打开输出核对计数、标签、供体划分和哈希；由验收表重绘供体保留率图。
- `launch_scp1303_{qc,verification}.sh`：本次固定路径的 nohup 启动器，保存日志和 PID。

## 本次环境与复现

工作目录：`/home/zhome/xiaoe-internship-data/_processing/SCP1303_20261111/`。Python 3.12 / Scanpy 1.12.5，精确依赖锁定见结果目录 `reports/environment.freeze.txt`；已有环境在 `venv/`。首次 bootstrap 的宽版本约束仅用于安装，复现应使用实际 freeze 版本。

本次将脚本同步到工作目录后运行 `bash launch_scp1303_qc.sh`，导出完成后运行 `bash launch_scp1303_verification.sh`。两次验收均通过，最后一次增加了供体保留率表的独立绘图；H5AD 未改写。

使用已验证环境在一个**新的输出目录**复现：

```bash
work=/home/zhome/xiaoe-internship-data/_processing/SCP1303_20261111
export MPLBACKEND=Agg OPENBLAS_NUM_THREADS=4 OMP_NUM_THREADS=4 NUMBA_NUM_THREADS=4
"$work/venv/bin/python" "$work/run_scp1303_qc.py" \
  --input "$work/author_input.h5ad" \
  --output /home/zhome/xiaoe-internship-data/processed/SCP1303_qc_reproduction \
  --seed 20261111 --light-size 20000
"$work/venv/bin/python" "$work/verify_scp1303_qc.py" \
  --input "$work/author_input.h5ad" \
  --output /home/zhome/xiaoe-internship-data/processed/SCP1303_qc_reproduction
```

长任务应使用 nohup 或 tmux；清洗脚本拒绝覆盖已有成功对象。输出的文件级哈希受 H5AD 中生成时间影响，复现比较应以成员、数值、元数据和参数为准。

核心测试命令（在本目录执行）：

```bash
python -m unittest -v test_qc_core.py
```

仅代码及文本汇总入库；源/输出 H5AD、逐核审计 gzip 和图留在服务器，图亦有本地交付副本。暂定 split 和混合标记 flag 的限制必须随对象一起传递。
