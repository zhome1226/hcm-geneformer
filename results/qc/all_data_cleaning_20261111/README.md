# 全部已下载数据：清洗与验收交付

状态日期：2026-11-11 UTC。14套中13套完成各自模态的整理/质控验收，**GSE227734仍缺权威基因映射，全部数据清洗尚未完成**。逐套状态以[CLEANING_STATUS.tsv](CLEANING_STATUS.tsv)为准；来源准备审计通过不等于完成分子质控。

服务器：`zhome@10.158.203.133`；数据根目录：`/home/zhome/xiaoe-internship-data/`。原件、表达矩阵、逐细胞/逐特征记录和影像留在服务器；本目录只交付代码入口、规则、摘要、文件清单和验收证据。文件路径及摘要见各套MANIFEST。原始计数保留；不做模型训练、tokenization、差异表达或靶点排序。

| 数据 | 模态 | 本轮状态 | 报告 |
| --- | --- | --- | --- |
| SCP1303 | 成人心室snRNA | 完成定向复核，使用reviewed_final版本 | [最终报告](scp1303_reviewed/REPORT.md) |
| Chin GSE174691、GSE161921 | HCM/供体对照snRNA，两套独立保存 | 两套完成 | [Chin报告](CHIN_REPORT.md) |
| Liu2023 | HCM/健康snRNA | 完成 | [报告](LIU_REPORT.md) |
| GSE181764 | 不同临床阶段snRNA | 完成；对照核心仅5核 | [报告](GSE181764_REPORT.md) |
| GSE138262 | 不同平台完整心肌细胞 | 完成严格版；仅6个HCM核心，无对照，不批准本版本训练 | [报告](GSE138262_REPORT.md) |
| GSE203274 | 儿童多疾病snRNA | 完成；HCM仅1供体，其他疾病单独标识 | [报告](GSE203274_REPORT.md) |
| GSE227734 | 成人原代心肌scRNA | **未完成：9760列匿名，无基因映射** | [阻塞报告](GSE227734_REPORT.md)、[来源追查](gse227734/SOURCE_LOOKUP_SUMMARY.md) |
| GSE249925 | bulk RNA | 完成整理/审计；raw/TMM来源冲突未解决 | [报告](GSE249925_REPORT.md) |
| VCC2025 | 扰动scRNA比赛 | 完成发布输入验证与真实索引 | [报告](VCC2025_REPORT.md) |
| NeurIPS2023 | 扰动/多组学/DE比赛 | 完成，保留官方划分与数据类型 | [报告](NEURIPS2023_REPORT.md) |
| NeurIPS2021 | CITE-seq/Multiome比赛 | 完成，保持配对和官方划分 | [报告](NEURIPS2021_REPORT.md) |
| CausalBench2023 | K562/RPE1扰动scRNA比赛 | 完成 | [报告](CAUSALBENCH2023_REPORT.md) |
| ACDC2017 | 心脏MRI比赛 | 完成审计；72对空间头冲突另审 | [报告](ACDC2017_REPORT.md) |

## 规则与适用范围

[PROTOCOL.md](PROTOCOL.md)定义各模态流程。HCM单细胞/单核的500≤检测基因≤5000、MT<5%是用户确认的项目严格基线，不是行业统一硬标准。心肌身份结合参考、多基因程序和群体证据；疑似双细胞需多证据支持。作者已有去背景/过滤处理如实保留，不对预过滤对象虚称重新完成最初液滴调用。

完整细胞、单核、bulk、多组学、扰动和MRI按各自技术处理；比赛数据不赋HCM疾病标签。训练资格另行评估：完成清洗不能证明足以进行HCM二分类或微调。跨研究供体重叠、技术重复和最终训练/验证划分尚未定稿。

SCP1303当前核心114711核、轻量20000核、待定144核，路径为`processed/SCP1303_qc_reviewed_final_20261111/`。[旧严格基线](../SCP1303_20261111/README.md)保留追溯，不作为当前最终版。Chin/Liu等最终目录见各报告，历史失败尝试有记录且未冒充成功。

## 代码及证据

代码在仓库`src/qc_chin/`、`src/qc_liu/`、`src/qc_gse*/`、`src/qc_scp1303_review/`和各比赛`src/qc_*/`入口README。现有服务器Python/R环境已通过导入核查，见[环境实测](ENVIRONMENT_RECHECK_20261111.txt)。已完成作业不应因旧日志或观察超时而重启；先确认真实进程及退出记录。

各数据子目录保存独立验收、实际规则、样本/供体损失汇总及终态退出证据。验收项数量是覆盖范围记录，不能当作生物学准确率。`TEXT_DELIVERY_COPY_AUDIT.json`将归档内容逐文件映射至本地源文件与SHA256；本次Git交付门禁另见`GIT_DELIVERY_AUDIT.json`。

本批代码和源报告逐字节归档；`GIT_DIFF_CHECK.txt`如实保留旧代码/报告末尾空行、Markdown空格及TSV末尾空字段的格式诊断。这些不是质控失败，也未为消除格式提示而改动已验收代码或裁掉TSV字段。

GSE227734新增官方补充材料检查：390→338、五供体数量及10种细胞条码与发布矩阵一致；6张补充小表没有9760列映射。P5年龄GEO/原矩阵63与补充表64的冲突保留。9项来源检查通过不改变其未完成状态，也未生成占位清洗对象；权威顺序映射或可逐值核对的修正版仍是下一步必要输入。
