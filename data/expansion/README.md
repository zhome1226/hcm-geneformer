# 公开数据补充下载批次

用户于 2026-11-11 明确授权将此前列出的公开比赛和补充队列先下载保存。本批启动于当日，**是否完成以 [STATUS.json](STATUS.json) 为准**；`downloading` 表示仍在传输，不能视为文件齐全。没有进行清洗、细胞筛选、表达分析、模型训练或执行比赛代码。

共计划 **922 个文件、83,328,342,290 字节**。其中 ACDC 占 904 个文件（150 位患者的影像、标注、Info.cfg 及公共许可/引用文件），文件数不能解释为患者数。

| 资源 | 选择范围 | 服务器相对目录（根为 `~/xiaoe-internship-data/`） |
| --- | --- | --- |
| VCC2025 | train / validation / test 三份完整 H5AD，共 34,362,204,170 字节 | `benchmarks/competitions/vcc2025/` |
| NeurIPS2023 | 完整 sc_counts、作者重注释计数版、基线 multiome、公开 DE train/test、药物注释、id_map、dataset_info；8 文件 | `benchmarks/competitions/neurips2023/` |
| CausalBench2023 | 官方框架指定的 K562 essential / RPE1 原始单细胞 H5AD；2 文件 | `benchmarks/competitions/causalbench2023/data/` |
| NeurIPS2021 / GSE194122 | 官方 GEO 的 CITE / Multiome 两份处理后 H5AD gzip，包含 RNA/蛋白或 RNA/ATAC | `benchmarks/competitions/neurips2021/` |
| ACDC2017 | 官方平台 training 100 人 / testing 50 人；150 人中含 30 HCM；MRI、标注和说明 | `benchmarks/competitions/acdc2017/` |
| GSE203274 | AllNuclei snRNA 原始计数 tar.gz 及 metadata gzip；2 文件 | `validation/pediatric_GSE203274/` |
| GSE138262 | 原 RAW.tar，25,886,720 字节；已核对 52 个成员与官方 filelist，解包到 raw/ 并通过内部 gzip CRC | `validation/hcm_GSE138262/` |

## 入口、许可与适用范围

- VCC：原 [官方 Atlas](https://github.com/ArcInstitute/arc-virtual-cell-atlas/tree/main/virtual-cell-challenge) 的 GCS 公共对象；逐对象 generation、大小和 CRC32C 已存档。三份矩阵并非 HCM 队列。
- NeurIPS2023：使用 [Open Problems 官方基准仓库](https://github.com/openproblems-bio/task_perturbation_prediction) 所属项目的公开 S3 `openproblems-bio/public/neurips-2023-competition/`。直接获取主办方公开对象，无需登录 Kaggle；未接受 Kaggle 条款。未下载重复历史训练矩阵、预测/评分结果或 LINCS 扩展。
- CausalBench：[官方固定版本下载实现](https://github.com/causalbench/causalbench/blob/1a2143cffdc85f835b41ce8d52034be1bf903e71/causalscbench/data_access/datasets/download_weissmann.py) 指定 Figshare 35773219 / 35775606；大小和 MD5 来自 [原研究 Figshare 元数据](https://api.figshare.com/v2/articles/20029387)。RPE1 使用其 `supplied_md5`，该对象的 `computed_md5` 为空。没有下载 65.8 GB 的 GWPS 全基因组对象，选定范围与 CausalBench 官方框架一致。
- NeurIPS2021：[活动页](https://openproblems.bio/events/2021-09_neurips/) 链接 [GSE194122](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE194122)。保存公开完整处理对象，未另外下载 ATAC fragments、FASTQ 或每个赛题的重复切分副本。
- ACDC：[官方入口](https://www.creatis.insa-lyon.fr/Challenge/acdc/databases.html) 指向公开 Girder 平台，匿名列出和下载文件。保留 [许可条款](ACDC_LICENSE_TERMS.md)：CC BY-NC-SA 4.0，限定非商业科研用途；使用时引用 Bernard et al., IEEE TMI 2018，DOI `10.1109/TMI.2018.2837502`。MRI 与当前 Geneformer 转录组分开保存，不当作同模态验证集。
- [GSE203274](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE203274)：儿童混合疾病队列，HCM 技术重复及来源关系仍待核实；保存数据不意味着已认定其为独立成人 HCM 队列。不下载 PBMC、scaled 矩阵或各亚群的重复副本。
- [GSE138262](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE138262)：HCM 单组室间隔数据，缺同批健康对照，作补充资料保存。
- [EGAD50000002046](https://ega-archive.org/datasets/EGAD50000002046) 仍未下载：官方页面明确受控访问，需数据访问委员会批准和 DAA；本批不申请或绕过权限。

代码许可不自动授权数据用途；各数据按原发布方许可、访问及引用要求使用。

## 状态与完整性

- [ASSETS.json](ASSETS.json)：全部计划对象，包含预期字节数、公开来源 URL 和可取得的发布方校验信息；被列入计划不代表已下载完成。
- [STATUS.tsv](STATUS.tsv)：逐文件 `not_started` / `downloading` / `verified` 状态及已写入字节数。
- [MANIFEST.tsv](MANIFEST.tsv)：**只登记已完成大小与完整性检查的文件**；[MD5SUMS](MD5SUMS) 为对应实测 MD5。
- VCC 对照官方 CRC32C；CausalBench 对照 Figshare MD5；NeurIPS2023 对照 S3 ETag，multipart 对象核对 8 MiB 分片复合 MD5。未取得发布方哈希的 GEO/ACDC 文件记录实测 MD5/SHA-256、核对官方大小，gzip 另校验 CRC，不声称有官方参考 MD5。
- 服务器日志、官方对象快照、每文件 verified.json、各批 summary 和实时清单均在 `~/xiaoe-internship-data/_download_expansion/`。下载源故障会记录错误，不生成占位文件。

## 后台运行与复现

下载器 [download_expansion.py](../download_expansion.py) 使用 `curl -C -`；每次重连重新读取实际长度，保持至少 10 GiB 空间余量。使用 nohup，在服务器运行；[discover_acdc.py](../discover_acdc.py) 仅列举官方公开目录，[monitor_expansion.py](../monitor_expansion.py) 持续写状态，所有批次退出后结束。

部分来源需经 PC 网络出口，临时 SSH 转发仅监听服务器 `127.0.0.1:39843`，对应 PC 已配置的 `127.0.0.1:7890`。[keep_download_tunnel.ps1](../keep_download_tunnel.ps1) 在 Windows 隐藏运行并自动重建转发；需要 PC 保持开机、联网和 VPN 可用。所有批次结束后，它收取最终报告到本地 `D:\VirtualCell\sources\expansion\final\`，停止转发并退出。

[finalize_download_batch.ps1](../finalize_download_batch.ps1) 随后只发布四份文本报告：MANIFEST、MD5SUMS、STATUS.tsv 和 STATUS.json。它检查工作区干净、分支仍为 virtual-cell、HEAD 与启动记录一致且 origin 未变；满足时提交并正常 push，核对远端 SHA。仓库有后续改动、网络失败或 push 冲突时保留报告并在 `final/SYNC_STATUS.json` 记录需人工处理，不覆盖其他工作或强推。下载失败时发布状态也如实为 partial_failure。

初始提交中的状态是带时间戳的快照；实时状态在服务器 `_download_expansion/STATUS.json` 和 PC `D:\VirtualCell\sources\expansion\LIVE_STATUS.json`。如需停止 PC 辅助进程，可创建本地 `sources/expansion/STOP_DOWNLOAD_MONITOR`；这会中止依赖该出口的传输，服务器已有文件可续传。
