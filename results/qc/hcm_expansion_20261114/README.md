# HCM扩展下载及完整细胞QC结果

日期：2026-11-14（北京时间）。

- [下载报告](../../../data/validation/HCM_EXPANSION_20261114.md)：GSE255296四个人源HCM样本、GSE292754两位儿童HCM样本，18文件共438,817,642字节；下载校验完成，尚未QC。
- `DOWNLOAD_MANIFEST.tsv`：逐文件服务器路径、字节数、本地MD5/SHA256、来源URL及校验语义。
- `DOWNLOAD_STATUS.json`、`SAMPLE_MATRIX_CHECKS.tsv`：来源大小、gzip CRC、6个矩阵全量数值/维度检查及状态。
- `DOWNLOAD_RECHECK.json`：18文件的本地MD5SUMS再校验；来源未提供MD5参考值，不声称上游MD5已比对。
- [完整细胞QC报告](../../../docs/hcm-whole-cell-qc-review.md)：GSE138262新版保留3,595心肌细胞，97候选留待复核，71项计数/元数据验收通过。
- `whole_cell/`：QC政策、供体/样本筛除表、非互斥复核原因、输出大小/哈希及完成凭据。

逐细胞审计、表达矩阵和参考模型保留服务器，未提交Git。新数据没有做训练、DGE、扰动或序列化。
