# 新增HCM人源样本下载

后续进展（2026-11-14）：本报告记录的两套发布矩阵已完成 QC，分别得到 7,516 和 3,306 个心肌核核心，并已分别完成 Geneformer V1 编码。ToDesk 交付状态与规则见[后续报告](../../docs/hcm-expansion-qc-v1-delivery.md)。以下下载阶段的“未QC”描述是当时的来源状态。

下载日期：2026-11-14（北京时间）。在原Linux数据服务器`ad2012`直接以`nohup`运行`wget -c`，开始前剩余183 GB。正式目录：

```text
/home/zhome/hcm-geneformer-data/validation/hcm_expansion_20261114/
  GSE255296/  # 四个人源HCM样本，各自保留10x矩阵、features、barcodes
  GSE292754/  # 两位儿童HCM样本，各自保留上述三个文件
```

| 数据 | 此次下载样本 | 表达文件 | 字节数 | 发布条码数（未QC） |
|---|---|---:|---:|---:|
| GSE255296 | GSM8068510/5246、GSM8068511/5042、GSM8068512/2883、GSM8423481/2795 | 12 | 297,299,364 | 48,873 |
| GSE292754 | GSM8864696/P186、GSM8864697/P194 | 6 | 141,518,278 | 30,187 |
| 合计 | 6个人源HCM样本 | 18 | 438,817,642 | 79,060 |

这里的条码数是发布矩阵的列数，不能直接称为通过QC的细胞数或心肌核数。全部保留来源计数，未进行过滤、归一化、心肌筛选、tokenization、DGE或模型训练。

## 样本边界

- GSE255296是跨物种研究；本次仅下载4个人源HCM样本，猫鼠数据不纳入。论文整体人体比较使用了已有GSE161921四对照和GSE174691一HCM，不把复用材料算新患者。四个新增记录的年龄及跨研究供体独立性仍需核实，暂不能只凭不同accession称全新外部验证。[官方GEO](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE255296)、[原论文](https://pmc.ncbi.nlm.nih.gov/articles/PMC12054467/)
- GSE292754新增P186、P194，两者为儿童HCM，单独存放。论文HCM1=P86此前已在GSE203275发布，包含于我们GSE203274资源，此次不重复下载或新增计数；公共对照同样须去重。[官方GEO](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE292754)、[原论文](https://pmc.ncbi.nlm.nih.gov/articles/PMC13067291/)

## 校验结果

18/18文件大小与NCBI响应头`Content-Length`一致；18/18完整解压并通过gzip CRC。六个样本均为36,601个features；矩阵维度与features/barcodes完全匹配，两轴ID唯一，全部存储计数为有限、非负整数。

每个文件已计算本地MD5与SHA256；两个数据集目录各有`MD5SUMS`，另一次`md5sum -c`复核12/12和6/6通过。检查到的来源响应头没有提供MD5参考值，因此不能声称“与作者MD5全部一致”。此处明确区分来源大小/CRC/结构检查与本地哈希追踪。

完整[文件清单](../../results/qc/hcm_expansion_20261114/DOWNLOAD_MANIFEST.tsv)、[矩阵检查](../../results/qc/hcm_expansion_20261114/SAMPLE_MATRIX_CHECKS.tsv)、[完成状态](../../results/qc/hcm_expansion_20261114/DOWNLOAD_STATUS.json)、[本地MD5复核](../../results/qc/hcm_expansion_20261114/DOWNLOAD_RECHECK.json)。

## 复现下载

脚本：`src/data_expansion/download_hcm_geo.py`。样本配置沿用`results/target_discovery/hcm_catalog_20261114/NEW_HCM_CANDIDATE_SAMPLES.tsv`，该文件是下载前元数据快照；当前下载状态以本报告和上述完成状态为准。

```bash
nohup python download_hcm_geo.py \
  --samples NEW_HCM_CANDIDATE_SAMPLES.tsv \
  --output /home/zhome/hcm-geneformer-data/validation/hcm_expansion_20261114 \
  > download.log 2>&1 &
```

下一步如开展这两套QC，需按样本核对临床、去重关系、原始发布矩阵是否已进行细胞调用；儿童与成人保持分层，不能直接把所有条码加入现有训练集。
