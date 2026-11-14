# HCM Geneformer V1 delivery

Packaged UTC: 2026-11-14T12:51:59.504814+00:00

Destination root: D:/hcm-geneformer-data; each dataset is separate; no training/validation roles assigned.
QC cores are unchanged source-count AnnData files, with all native genes and any native missingness masks.
Geneformer V1: 2048-token maximum, no CLS/EOS, official median-scaled ranking. All tokens and metadata independently verified before packaging.
New nuclear cohorts use 500<=detected genes<=5000 and MT<5%; whole-cell GSE138262 uses author-specific non-MT UMI>1000 and identity consensus.
GSE255296 HCM-only (age unresolved); GSE292754 HCM-only pediatric females 14/16; whole-cell disease/technology/study confounded.
No ambient-based independent droplet calls, no verified ambient correction history for new matrices, no independent validation/training approval.

| Dataset | CM cells/nuclei | Donor identifiers | Disease counts |
|---|---:|---:|---|
| GSE255296 | 7516 | 4 | {'HCM': 7516} |
| GSE292754 | 3306 | 2 | {'HCM': 3306} |
| GSE138262_whole_cell | 3595 | 12 | {'HCM': 2216, 'Control': 1379} |
