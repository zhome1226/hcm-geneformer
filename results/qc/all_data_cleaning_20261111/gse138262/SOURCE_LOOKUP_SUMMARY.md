# GSE138262 输入来源更正

作者仓库`vanrooij-lab/scRNAseq-HCM-human`在根readme和`Projects/howtousedata_HCMSCS/data_notes/readme.md`明确说明：GEO旧TSV并非Wehrens et al., Cell Reports 2022论文使用的计数表。本地证据在`source_lookup/`，仓库树SHA为`621c874b6e1f98f6eee83084c6424a52e24d4f97`。

作者链接https://doi.org/10.5281/zenodo.18206553；公开API解析为版本记录18206554，DOI https://doi.org/10.5281/zenodo.18206554 ，CC-BY-4.0。下载5个文件共327015800字节，全部符合该版本记录的大小和MD5：两份计数表ZIP、两个Seurat对象、FACS_INDEX_data.zip。详见上级GSE138262_CORRECTED_DOWNLOAD_MANIFEST.tsv与本目录DOWNLOAD_VERIFICATION.json。

32个修正计数表实际用空格分隔，不能机械按扩展名指定tab。ZIP带macOS资源叉，须识别而非作为计数文件读取。部分基因只有Ensembl ID、symbol为空；保留原行身份，不猜配或删行。作者脚本JE03/JE04的实际文件映射必须与GSM/plate表联合核对，不能按编号自动对应。

作者Seurat对象的RNA assay已移除MT，两个对象实际只提供RNA assay；保存的nCount_rnaMitoCounts/percent.mt元数据不代表MT仍在RNA矩阵中。独立从完整修正计数表重算MT，与作者MT计数元数据逐细胞一致。48组板/作者对象共享RNA计数逐值一致。

HCM为5供体16个分选CEL-seq2完整心肌细胞板；健康伴随为Wang7供体16个iCell8板/part，已按原研究LV群预选后由Wehrens重处理。作者ALL对象额外引用Litvinukova2020图谱，不将其重复作为新独立对照。原53个GEO下载文件大小/MD5再次核对一致，继续保留用于历史溯源；修正文件不覆盖旧文件。

最后清洗结果及适用限制见REPORT.md：735项验收通过，但严格MT<5%等规则只保留8数值候选、6心肌核心、2待定，无健康对照，不批准进入HCM/Control训练。7537个来源细胞及28663个并集基因全部保留。
