# 执行与恢复记录

日期：2026-11-14。错误被保留，未以占位文件补齐数据；原计数和已完成的旧四套交付未修改。

- QC 首次 R 启动因复用脚本时参考目录日期错误中止，修正为真实既有 `SCP1303_review_20261111`；复用已完成 staging，最终 `qc_exit_code.txt=0`。GSE255296 336 项、GSE292754 228 项均通过。
- 编码准备先遇到旧版 AnnData 默认不允许 nullable string 写入，随后遇到原完整细胞索引名 `cell_id` 与新增编码列重名。明确启用兼容写入，并只在编码副本改索引名称、保留原名称记录；失败产物独立隔离，来源原件摘要不变。
- 96 细胞 smoke 和全部官方 V1 编码成功。第一次 full 验收用 CSR 内部 data/index 数组顺序比较，因准备副本规范化排序而中止。三套全矩阵重新逐坐标比较，计数和 counts layer 均完全相同；校验改为逐坐标比较，再执行所有 full 排名/metadata 检查，最终全部通过，`v1_exit_code.txt=0`。没有重新编码以掩盖差异。
- 首次上传在 `reput` 发现远端归档尚不存在时退出，未传输数据。传输程序改为首次 `put`、存在时 `reput` 续传；正式目录只在目标端全部验收后放置。
- 随后直连 Tailscale 经 sfo 中继传至 4,403,200 字节后因慢速停止；经 Linux ProxyJump 的续传亦超时，未达到目标端验收条件。GitHub 私有草稿 Release 作为传输备选，正在上传归档；已询问用户确认 ToDesk 开机联网。此阶段不能称 ToDesk 交付完成。
- 后续 GitHub 上传通过本机既有 `127.0.0.1:7890` 进程代理成功，耗时约413秒，GitHub 服务端 SHA256 与源归档一致。此前本机直连慢速、Linux 连接 uploads.github.com 超时均未成功，未把失败说成上传完成。
- ToDesk 的 HTTPS 续传在三个官方 CDN 地址均发生读取/连接超时；最后保留9,810,793字节，低于197,790,119预期。小范围读取测试也超时。所有有界尝试已退出；没有运行目标端接受脚本或生成虚假的 WINDOWS_ACCEPTANCE.json，正式目录未放置。已请用户提供远端可用代理。

长日志和失败产物索引保存在本机 `D:/VirtualCell/deliverables/hcm-qc-token-transfer-20261114/` 与 Linux `/home/zhome/hcm-geneformer-analysis/hcm_expansion_20261114/nuclear_qc/`；Git 仅保存小型脚本、报告和验收摘要。
