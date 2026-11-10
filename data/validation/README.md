# HCM 独立验证数据

下载日期：2026-11-10（Asia/Shanghai，UTC+8）。服务器：`zhome@10.158.203.133`。

数据根目录：`/home/zhome/xiaoe-internship-data/validation/`。本目录仅提交来源说明和清单；数据保存在服务器，用途均为独立验证。

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
- 四个数据文件的大小和参考 MD5 全部匹配，详情见 [MANIFEST.tsv](MANIFEST.tsv)。GEO tar 及其 39 个解包文件的字节数逐一对照官方 `filelist.txt`；两个 `filelist.txt` 的大小对照来源 HTTP Content-Length。Liu 两个文件的字节数对照 Figshare 元数据。
- 三个服务器子目录各有 `MD5SUMS`，含文件字节数注释及实测 MD5，可在对应目录运行 `md5sum -c MD5SUMS`。四个下载数据文件具有独立参考 MD5；`filelist.txt` 和 GEO 解包文件记录本地 MD5，解包文件的来源完整性由已匹配参考 MD5 的原 tar 保证，不声称存在额外的发布方逐文件 MD5。
- 只将两套 GEO 的 `RAW.tar` 解包到各自的 `raw/`，保留原 tar，内部 `.gz` 不解压。Liu 的 `.Rdata` 和 `.tar.gz` 原样保留，未载入 Seurat、未清洗、未分析；未下载空间转录组数据文件。
- `MANIFEST.tsv` 含 45 条记录：4 个下载数据文件、2 个 GEO filelist、39 个 GEO 解包文件。相对路径均以服务器数据根目录为基准；解包文件的来源 URL 指向包含它的原 tar。操作脚本、运行日志、PID 和校验报告保存在服务器 `_download/`，不属于数据清单。

四个数据文件：

| 相对路径 | 字节数 | MD5 |
| --- | ---: | --- |
| `chin_GSE174691/GSE174691_RAW.tar` | 649707520 | `8b3f5a48add3dd0a1d4d1feb769b716a` |
| `chin_GSE161921/GSE161921_RAW.tar` | 225945600 | `5daa30b0d89827bc1618fb47fc7c878a` |
| `liu2023/HCM_snRNAseq_seuratObj.Rdata` | 2376598237 | `22153f776c3e667d87feead247272799` |
| `liu2023/snRNA-seq.tar.gz` | 450043610 | `3d850cd3017a9c8623febdc4269e71fa` |
