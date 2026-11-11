# SCP1303 人左心室单核转录组

来源为 [Chaffin et al., Nature 2022](https://www.nature.com/articles/s41586-022-04817-8)，研究标题 *Single-nucleus profiling of human dilated and hypertrophic cardiomyopathy*；[Broad SCP1303](https://singlecell.broadinstitute.org/single_cell/study/SCP1303/) 为原发布入口。[Theodoris et al., Nature 2023](https://www.nature.com/articles/s41586-023-06139-9) 为本项目采用的 Geneformer 方法论文，数据来源与方法引用分别记录。

## 队列与用途

- 42 位供体：15 HCM、11 DCM、16 非衰竭对照（NF）；作者报告 592,689 个细胞核、36,601 个基因，包含多种细胞类型，并非全为心肌细胞。
- 组织为人左心室；作者 README 说明疾病样本来自接受心脏移植的终末期心衰患者离体心脏，NF 来自心脏移植供体，不能当作早期 HCM 或室间隔切除队列。
- 计划作为 HCM/NF 心肌分类的主训练来源；DCM 保留但不预先合并为 HCM。供体划分、细胞筛选与现有验证集的供体重叠检查留到后续明确授权的处理阶段。

## 实际下载范围

服务器目录：`/home/zhome/xiaoe-internship-data/training/SCP1303/`。文件来源为官方 [Human Cell Atlas 同研究项目](https://explore.data.humancellatlas.org/projects/85c0d6fa-f117-4d76-b01a-5d5e8f5f9188)，该页面明确链接 SCP1303，列出相同研究、42 位供体及左心室组织。下载日期：2026-11-11（Asia/Shanghai）。

保存作者提供的 `human_dcm_hcm_scportal_03.17.2022.h5ad.gz` 及 6 份配套文件：作者 README、细胞注释、基因索引、barcode 索引、预先计算的 UMAP、HCA 样本元数据 XLSX。这里的 UMAP 是发布方已有文件，未运行新的降维。

作者 README 描述完整 H5AD：

- `.X`：CellBender 调整后的计数；
- `.layers['cellranger_raw']`：Cell Ranger 原始计数；
- `.obs`：`donor_id`、`biosample_id`、`disease`（DCM/HCM/NF）、细胞类型、年龄、性别及作者已有 QC 字段；
- `.var`：基因符号、Ensembl `gene_ids` 等。

仅保存压缩原件，不解压或载入表达矩阵。没有保留重复的 11.74 GB `raw_counts_V1.mtx`、18.19 GB 表达矩阵或规范化索引；下载范围足以保留两种计数和供体/疾病注释。原始 FASTQ 属于 dbGaP 受控数据，不在本次范围内。

## 完整性与复现

逐文件字节数、实测 MD5、SHA-256 与来源 URL 见 [MANIFEST.tsv](MANIFEST.tsv)。服务器保存 `MD5SUMS`、`SHA256SUMS` 和 `_catalog/download_summary.json`。SHA-256 与大小对照 HCA 官方对象索引；MD5 为本地实测校验值，发布方目录未提供参考 MD5。压缩 H5AD 另做 gzip CRC 检查。

2026-11-11 09:41（Asia/Shanghai）全部完成：7 个文件合计 **3,779,317,638 字节**，大小与官方 SHA-256 全部匹配，H5AD gzip CRC 检查通过。完整对象大小为 3,533,964,484 字节，SHA-256 为 `088da9ed54b3d45e039a00b7f3db58e9a032496fb1e45b1ac32019c2280c9cb1`；其余逐文件结果及实测 MD5 见清单与 [MD5SUMS](MD5SUMS)、[SHA256SUMS](SHA256SUMS)。仅下载、校验与整理；没有额外清洗、表达分析或训练。

官方原始索引快照为 [SCP1303_HCA_FILES.json](SCP1303_HCA_FILES.json)，包含全部 11 个发布对象；被索引不等于已下载，实际下载以 MANIFEST 为准。其 API 为 `https://service.azul.data.humancellatlas.org/index/files`，查询过滤项目 ID `85c0d6fa-f117-4d76-b01a-5d5e8f5f9188`。文件使用索引中的 `azul_mirror_uri` 所指官方公共 S3 对象。

下载脚本 [download_scp1303.py](download_scp1303.py) 仅依赖 Python 标准库和 curl，使用 `curl -C -` 续传、最多两个并发传输；断线时重新启动 curl，以现有文件长度恢复，避免内部重试回退到旧偏移量。完成后校验大小、SHA-256、MD5 和 gzip CRC；不进行生物学处理。

```bash
cd /home/zhome/xiaoe-internship-data/training/SCP1303
nohup python3 -u _catalog/download_scp1303.py \
  --catalog _catalog/hca_files.json --base . \
  > _catalog/download.log 2>&1 < /dev/null &
echo $! > _catalog/download.pid
```

本次服务器直连 Broad/HCA S3 超时，实际通过仅监听服务器回环地址的 SSH 转发使用 PC 已配置的网络出口，下载进程仍在服务器 `nohup` 下运行。按需可传 `--proxy http://127.0.0.1:PORT`，转发停止后续传需重建出口或恢复服务器直连。

HCA 项目页面注明数据按 [CC BY 4.0](https://creativecommons.org/licenses/by/4.0/) 发布，并适用其 [Data Use Agreement](https://data.humancellatlas.org/about/consent-and-ethics/data-use-agreement)。使用时引用原研究及 HCA 项目。Git 只保存本目录的文本说明、索引、清单和下载脚本，数据原件留在服务器。
