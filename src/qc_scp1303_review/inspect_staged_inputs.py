"""Inspect staged cohort metadata and count identities before cleaning."""
from pathlib import Path
import gzip
import json
import anndata as ad
import h5py
import pandas as pd

root=Path.home()/'xiaoe-internship-data'
work=root/'_processing/all_data_cleaning_20261111'
out=root/'processed/all_data_cleaning_20261111'
report={}
path=work/'GSE227734/GSE227734_hcm_raw.h5ad'
with h5py.File(path) as f:
    obs=ad.io.read_elem(f['obs']); var=ad.io.read_elem(f['var'])
    report['GSE227734']={'shape':[len(obs),len(var)],'obs_columns':list(obs.columns),
        'var_columns':list(var.columns),'var_index':list(var.index[:5]),
        'obs_examples':obs.head(3).astype(str).to_dict('records'),
        'categories':{c:obs[c].astype(str).value_counts().to_dict() for c in obs if obs[c].nunique()<30},
        'layers':list(f['layers']) if 'layers' in f else [],'has_raw':'raw' in f}
meta_path=root/'validation/pediatric_GSE203274/GSE203274_AllNuclei_snRNA_metadata.csv.gz'
meta=pd.read_csv(meta_path,low_memory=False)
report['GSE203274']={'shape':list(meta.shape),'columns':list(meta.columns),
    'examples':meta.head(3).fillna('').astype(str).to_dict('records'),
    'categories':{c:meta[c].astype(str).value_counts().to_dict() for c in meta if meta[c].nunique()<80}}
bulk=root/'validation/bulk_GSE249925/GSE249925_20231208-For-GEO-Raw-Counts.tsv.gz'
head=pd.read_csv(bulk,sep='\t',nrows=3)
report['GSE249925']={'columns':list(head.columns),'examples':head.astype(str).to_dict('records')}
for directory in ['NeurIPS2021_CITE','NeurIPS2021_multiome']:
    path=next((work/directory).glob('*.h5ad'))
    with h5py.File(path) as f:
        obs=ad.io.read_elem(f['obs']); var=ad.io.read_elem(f['var'])
        report[directory]={'path':str(path),'shape':[len(obs),len(var)],'obs_columns':list(obs.columns),
            'var_columns':list(var.columns),'layers':list(f['layers']) if 'layers' in f else [],'has_raw':'raw' in f,
            'categories':{c:var[c].astype(str).value_counts().to_dict() for c in ['feature_types','feature_type'] if c in var},
            'obs_categories':{c:obs[c].astype(str).value_counts().to_dict() for c in ['batch','donor','site','cell_type','celltype'] if c in obs}}
(out/'STAGED_SCHEMAS.json').write_text(json.dumps(report,indent=2,ensure_ascii=False)+'\n')
print('STAGED_METADATA_INSPECTED',flush=True)
