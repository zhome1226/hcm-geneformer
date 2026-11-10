# 数据

- 来源：GEO [GSE145154](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE145154) 或 [GSE130652](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE130652)，含健康对照（Donor Control / NF）与病理组（HCM / DCM）。
- 质控：线粒体基因比例 < 5%，检测到的基因数 500–5000。
- 亚群：提取心肌细胞（标记基因 TNNT2、MYH7、TTN），细胞病理标签为 `0: Control`、`1: HCM`。规模约 1–2 万细胞即可。
- 存放：原始数据放 `data/raw/`（已 gitignore），质控后 AnnData 放 `data/processed/`（默认不入库）。下载脚本和元数据可以提交。
