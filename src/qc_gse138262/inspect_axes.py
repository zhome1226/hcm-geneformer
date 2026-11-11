from pathlib import Path
import json
import pandas as pd
from stage_gse138262 import sample_definitions

W = Path.home()/'xiaoe-internship-data/_processing/all_data_cleaning_20261111/gse138262'
result = []
for sample in sample_definitions():
    table = pd.read_csv(sample['source_path'],sep=r'\s+',index_col=0,keep_default_na=False)
    axis = table.index.astype(str)
    invalid = axis[~axis.str.match(r'^ENSG\d+_.+$')]
    record = {'sample_id':sample['sample_id'],'shape':list(table.shape),
              'unexpected_gene_count':len(invalid),'unexpected_gene_examples':invalid[:20].tolist(),
              'first_genes':axis[:3].tolist(),'first_columns':table.columns[:3].tolist()}
    result.append(record)
    print(json.dumps(record),flush=True)
(W/'NATIVE_AXIS_INSPECTION.json').write_text(json.dumps(result,indent=2)+'\n')
