# GSE227734 基因注释查找与恢复条件

本轮确认：服务器上的官方 gzip 解压后，真实 HDF5 仅含 X、obs、var；var/_index 是 0–9759，var 没有其他列。没有 raw、uns 或 layers 可供恢复基因名。直接读取 gzip 内的数据也得到同样结果，排除了本轮解压或转换丢失基因名的可能。

已访问的来源与记录：

- [GEO 官方补充目录](https://ftp.ncbi.nlm.nih.gov/geo/series/GSE227nnn/GSE227734/suppl/)：HTTP200，只有原来的 GSE227734_hcm_raw.h5ad.gz；当前目录没有独立 features 文件。目录快照在 source_lookup/geo_index.txt。
- [GEO 系列](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE227734) 与5个GSM：给出STRT-seq、患者别名和原始数据关系，未提供原9760列的映射。
- [原论文 PubMed](https://pubmed.ncbi.nlm.nih.gov/37979444/)：明确338个原代心肌细胞来自HCM切肌患者。DOI 10.1016/j.yjmcc.2023.11.005。出版商页面和公开内容API本轮均返回403，未取得完整补充材料，不能声称已排除其中存在映射的可能。
- 网络搜索找到使用本数据的第三方 [AdiVM/SLMC_single-cell](https://github.com/AdiVM/SLMC_single-cell)；递归文件清单没有专门的9760列基因映射文件。本轮只把它作为查找线索，未将其当作原作者数据，也未执行其代码；清单快照在 source_lookup/github_repository_tree.txt。
- [ENA 原始读段记录](https://www.ebi.ac.uk/ena/portal/api/filereport?accession=PRJNA946727&result=read_run&fields=run_accession,sample_alias,fastq_bytes,fastq_ftp&format=tsv)：39个双端运行，压缩FASTQ合计65,959,532,565字节（61.43GiB）。完整索引在 source_lookup/ena_reads.txt；本轮未下载FASTQ。它不是现成的基因映射。

各请求状态和URL记录在 source_lookup/LOOKUP_RESULTS.json。GEO目录的filelist.txt返回404，不等于下载失败；该系列目录确实没有该附属文件。

## 2026-11-11 17:15 UTC：官方补充材料新增证据

已取得[出版商官方补充Word文件](https://ars.els-cdn.com/content/image/1-s2.0-S0022282823001876-mmc1.docx)：3,853,376字节，SHA256 `6b6ef2ca9c0737b0f300fa87244bf9ba9acd9bcd048839e8951b046bc2c5de54`。原件另存服务器 `references/GSE227734_publisher_20261111/mmc1.docx`，不放入Git。正文/内容API的403不能用于推断补充材料也不可获取。

独立解析全部6张表、297段文本和DOCX附件/外链：没有嵌入文件或外部数据链接；最大表21行，没有原9760列的顺序映射。Table S5仅20个与面积相关的基因及相关系数，不能对应匿名列，也不能充作完整基因表。

- Table S1给出10个8碱基细胞条码，集合与原矩阵完全一致；并未给出原始文库到338个作者细胞的完整映射、完整UMI/读段布局或确切基因注释版本，不能按39个run的顺序推断。
- Table S4：过滤前P1–P5为40/80/90/90/90，共390；过滤后37/75/70/66/90，共338，与原矩阵逐供体一致。该发布输入已经过作者过滤。
- Table S3：五名供体均男性；P1–P4年龄与GEO一致。P5原矩阵/GEO为63，补充表为64，保持来源冲突未解决，不自动覆盖任一值。
- 新增远端审计 `processed/cox7b_GSE227734_qc_20261111_pending_annotations/source_evidence/`，9项来源核查通过；`SUPPLEMENT_SOURCE_AUDIT.json`及`SUPPLEMENT_DONOR_CROSSCHECK.tsv`记录证据，计数未改，清洗仍未完成。
- 第三方SLMC仓库实际读取本机`GSE227734_hcm_raw.h5ad`，记录同样338×9760，仅数值诊断；未发现基因映射或修正版入口。Europe PMC和OpenAlex未提供开放正文镜像；GitHub公开仓库名检索GSE227734/COX7B及作者邮箱前缀无结果。检索不构成不存在其他来源的证明。

本轮检索原始结果、官方来源快照和Word解析保存在本地 `sources/research_gse227734_followup_20261111/`，以及 `sources/research_gse227734_*20261111*.json`；远端历史16:55 UTC超时记录保持归档，不覆盖。未联系作者，未下载61.43GiB FASTQ。

清洗仍未完成。恢复条件是有权威的、按原9760列顺序排列的基因标识，或计数和细胞可核对的修正版矩阵。不能猜用全基因组顺序、其他队列的轴或算法推测标签。若从原始FASTQ重建计数，还需先查清39个文库/运行与338个作者细胞的对应、条码/UMI结构及参考注释版本，并独立保存新计数及差异；不能假装它们就是原作者的计数。
