"""Repair a stage metadata dtype conversion; prove all sparse count components unchanged."""
from datetime import datetime, timezone
from pathlib import Path
import hashlib
import json
import shutil

import anndata as ad
import h5py
import numpy as np
import pandas as pd

ROOT=Path.home()/'xiaoe-internship-data'
WORK=ROOT/'_processing/all_data_cleaning_20261111/GSE203274'


def age(text):
    return str(text).replace('/','_')


def component_hashes(handle):
    result={}
    for name in ['data','indices','indptr']:
        dataset=handle['X'][name];hasher=hashlib.sha256()
        for start in range(0,len(dataset),10_000_000):hasher.update(dataset[start:start+10_000_000].tobytes())
        result[name]=hasher.hexdigest()
    return result


path=WORK/'author_full_counts.h5ad'
assert (WORK/'stage_exit_code.txt').read_text().strip()=='0'
stage=json.loads((WORK/'STAGE_SUCCESS.json').read_text())
archive=WORK/'stage_attempt1_metadata_dtype'
archive.mkdir(exist_ok=True)
if not (archive/'STAGE_SUCCESS.json').exists():shutil.copy2(WORK/'STAGE_SUCCESS.json',archive/'STAGE_SUCCESS.json')
with h5py.File(path,'r+') as f:
    before=component_hashes(f)
    obs=ad.io.read_elem(f['obs'])
    author=pd.read_csv(ROOT/'validation/pediatric_GSE203274/GSE203274_AllNuclei_snRNA_metadata.csv.gz',index_col=0,keep_default_na=False)
    assert np.array_equal(obs.age_original.to_numpy(),author.loc[obs.index,'age'].astype(str).to_numpy())
    # Compare full structured dates; text separators alone differ between source formats.
    conflict=obs.age_original.astype(str).map(age).ne(obs.geo_age_original.astype(str).map(age)).to_numpy(dtype=bool)
    obs['age_source_conflict']=conflict
    ad.io.write_elem(f,'obs',obs)
    uns=ad.io.read_elem(f['uns'])
    uns['source']='GSE203274 author released unnormalized integer RNA counts after CellBender v0.1.0 and author QC; rawCount is the source filename label, not uncorrected CellRanger counts or full droplets'
    ad.io.write_elem(f,'uns',uns)
    after=component_hashes(f)
    assert before==after,'Count matrix changed during metadata-only repair'
obs.to_csv(WORK/'cell_numeric_audit.tsv.gz',sep='\t')
stage['age_conflicting_donors']=sorted(obs.loc[conflict,'donor_original'].unique())
stage['metadata_repair']='Restored bool age_source_conflict from source text dates; sparse matrix components SHA256 unchanged'
(WORK/'STAGE_SUCCESS.json').write_text(json.dumps(stage,indent=2)+'\n')
result={'repaired_utc':datetime.now(timezone.utc).isoformat(),'before_X_SHA256':before,'after_X_SHA256':after,
        'age_conflicting_donors':stage['age_conflicting_donors'],'counts_unchanged':True}
(WORK/'METADATA_REPAIR_VERIFICATION.json').write_text(json.dumps(result,indent=2)+'\n')
print(json.dumps(result,indent=2),flush=True)
