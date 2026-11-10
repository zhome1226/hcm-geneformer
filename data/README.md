# 数据

- 主训练来源：[SCP1303 / Chaffin et al., Nature 2022](https://singlecell.broadinstitute.org/single_cell/study/SCP1303)，15 HCM、11 DCM、16 非衰竭心脏，详见 [#1](https://github.com/zhome1226/xiaoe-internship/issues/1)。旧推荐 GSE145154 无 HCM，GSE130652 是小鼠 ChIP-seq，均不用于本任务。
- 已下载的验证数据见 [validation/README.md](validation/README.md)，新增候选的适用范围和排除原因见 [validation/CANDIDATES.md](validation/CANDIDATES.md)。
- 质控：线粒体基因比例 < 5%，检测到的基因数 500–5000。
- 亚群：提取心肌细胞（标记基因 TNNT2、MYH7、TTN），细胞病理标签为 `0: Control`、`1: HCM`。规模约 1–2 万细胞即可。
- 存放：原始数据放 `data/raw/`（已 gitignore），质控后 AnnData 放 `data/processed/`（默认不入库）。下载脚本和元数据可以提交。
