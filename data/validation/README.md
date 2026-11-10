# HCM 独立验证数据

下载日期：2026-11-10（Asia/Shanghai，UTC+8）。服务器：`zhome@10.158.203.133`。

数据根目录：`/home/zhome/xiaoe-internship-data/validation/`。本目录提交说明、清单和下载脚本；数据保存在服务器，按下述用途分别用于独立验证、亚型补充或靶点复核。

第一批验证数据：

| 数据集 | 样本构成 | 组织部位 | 服务器子目录 |
| --- | --- | --- | --- |
| Chin / GSE174691 | 9 例梗阻性 HCM，9 个样本 | 室间隔，心肌切除标本 | `chin_GSE174691/` |
| Chin / GSE161921 | 4 例未用于移植的供体对照，4 个样本 | 室间隔 | `chin_GSE161921/` |
| Liu et al., Cell Discovery 2023 | 10 例 HCM（10 个样本）+ 2 例健康供体（3 个样本），合计 13 个样本 | 两组均为室间隔 | `liu2023/` |

Liu 数据中的 HEALTHY1A 和 HEALTHY1B 来自同一供体，因此 12 名受试者对应 13 个样本。此信息来自论文正文及图 1，未对下载的数据进行分析。

来源：

- [GSE174691](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE174691)：[RAW.tar](https://ftp.ncbi.nlm.nih.gov/geo/series/GSE174nnn/GSE174691/suppl/GSE174691_RAW.tar)、[filelist.txt](https://ftp.ncbi.nlm.nih.gov/geo/series/GSE174nnn/GSE174691/suppl/filelist.txt)。
- [GSE161921](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE161921)：[RAW.tar](https://ftp.ncbi.nlm.nih.gov/geo/series/GSE161nnn/GSE161921/suppl/GSE161921_RAW.tar)、[filelist.txt](https://ftp.ncbi.nlm.nih.gov/geo/series/GSE161nnn/GSE161921/suppl/filelist.txt)。
- Liu et al.：[论文](https://www.nature.com/articles/s41421-022-00490-3)、[开放全文](https://pmc.ncbi.nlm.nih.gov/articles/PMC9842679/)、[Figshare 合集](https://doi.org/10.6084/m9.figshare.c.5777948.v2)。[Seurat 对象页面](https://figshare.com/articles/dataset/Seurat_object_for_the_HCM_snRNAseq_dataset/17869103)、[下载](https://ndownloader.figshare.com/files/32582447)；[Cell Ranger 输出页面](https://figshare.com/articles/dataset/single-nucleus_RNA-seq_data_of_cardiac_tissues_from_HCM_patients_and_healthy_donors/19960811)、[下载](https://ndownloader.figshare.com/files/35521496)。

下载与校验：

- 开始前 `df -h ~` 显示约 349 GB 可用，满足至少 10 GB 的要求。下载由服务器 `nohup` 后台任务运行，使用 `wget -c` / `curl -C -`。Figshare 初始服务器直链返回 403，最终由 PC 获取其公开下载重定向，服务器从该有效存储链接下载。
- 第一批四个数据文件的大小和参考 MD5 全部匹配，详情见 [MANIFEST.tsv](MANIFEST.tsv)。GEO tar 及其 39 个解包文件的字节数逐一对照官方 `filelist.txt`；两个 `filelist.txt` 的大小对照来源 HTTP Content-Length。Liu 两个文件的字节数对照 Figshare 元数据。
- 三个服务器子目录各有 `MD5SUMS`，含文件字节数注释及实测 MD5，可在对应目录运行 `md5sum -c MD5SUMS`。四个下载数据文件具有独立参考 MD5；`filelist.txt` 和 GEO 解包文件记录本地 MD5，解包文件的来源完整性由已匹配参考 MD5 的原 tar 保证，不声称存在额外的发布方逐文件 MD5。
- 只将两套 GEO 的 `RAW.tar` 解包到各自的 `raw/`，保留原 tar，内部 `.gz` 不解压。Liu 的 `.Rdata` 和 `.tar.gz` 原样保留，未载入 Seurat、未清洗、未分析；未下载空间转录组数据文件。
- 第一批贡献 45 条清单记录：4 个下载数据文件、2 个 GEO filelist、39 个 GEO 解包文件。相对路径均以服务器数据根目录为基准；解包文件的来源 URL 指向包含它的原 tar。操作脚本、运行日志、PID 和校验报告保存在服务器 `_download/`，不属于数据清单。

四个数据文件：

| 相对路径 | 字节数 | MD5 |
| --- | ---: | --- |
| `chin_GSE174691/GSE174691_RAW.tar` | 649707520 | `8b3f5a48add3dd0a1d4d1feb769b716a` |
| `chin_GSE161921/GSE161921_RAW.tar` | 225945600 | `5daa30b0d89827bc1618fb47fc7c878a` |
| `liu2023/HCM_snRNAseq_seuratObj.Rdata` | 2376598237 | `22153f776c3e667d87feead247272799` |
| `liu2023/snRNA-seq.tar.gz` | 450043610 | `3d850cd3017a9c8623febdc4269e71fa` |

## 同日补充下载

| 数据集 | 构成 | 服务器子目录 | 用途 |
| --- | --- | --- | --- |
| GSE181764 | 1 梗阻性 HCM + 6 非梗阻性 HCM + 2 供体对照；成人室间隔 | `chin_GSE181764/` | Chin 同实验室扩展，补充 HCM 亚型验证 |
| GSE227734 | 5 名成人 HCM 患者，作者报告 338 个心肌细胞；室间隔 | `cox7b_GSE227734/` | 疾病组心肌细胞特征复核；无配套健康单细胞对照，不单独评估二分类 |
| GSE249925 | 97 例 HCM + 23 例对照；心肌组织 bulk RNA-seq | `bulk_GSE249925/` | 靶点的组织层面复核，不直接用于单细胞分类器 |

来源、取舍和适用限制见 [CANDIDATES.md](CANDIDATES.md)。三套数据的 GEO 原始样本文本已保存，[SUPPLEMENTAL_SAMPLES.tsv](SUPPLEMENTAL_SAMPLES.tsv) 登记 134 条样本记录供后续核对，不作为额外患者或细胞计数。

新增数据文件的来源字节数分别为 301056000、2398031、6968101。均核对 NCBI `Content-Length`；GSE181764 仅解包到 `raw/`，保留原 tar，成员大小核对官方 filelist。另两个 gzip 仅校验压缩完整性，原样保留，未载入矩阵、未清洗、未分析。各新增目录均有实测 `MD5SUMS`；未获得发布方参考 MD5，故不声称新增文件与官方 MD5 匹配。

总清单现为 82 条：第一批 45 条 + 补充批次 37 条（3 个数据文件、1 个 GEO filelist、27 个解包文件、6 个原始 GEO 元数据文本）。下载及校验脚本为 [download_supplemental.py](download_supplemental.py)，运行日志和报告位于服务器 `_download/`。
