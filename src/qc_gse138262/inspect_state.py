from pathlib import Path
import json
import re
import pandas as pd

W = Path.home() / 'xiaoe-internship-data/_processing/all_data_cleaning_20261111/gse138262'
S = Path.home() / 'xiaoe-internship-data/validation/hcm_GSE138262/author_corrected_zenodo18206554'
main = (W/'AUTHOR_MAIN.R').read_text()
result = {}
for name, end in [('HCM_SCS_ids','dataset_list_paths_HCM'),('WANG_ids','dataset_list_paths_WANG')]:
    block = main.split(name + ' =',1)[1].split(end,1)[0]
    result[name] = re.findall(r'([\w.]+)\s*=\s*\x27([^\x27]+)\x27',block)
for kind in ['ROOIJonly','ALL']:
    result[kind + '_fields'] = pd.read_csv(W/f'author_{kind}/METADATA_FIELDS.tsv',sep='\t').to_dict('records')
result['objects'] = pd.read_csv(W/'AUTHOR_OBJECT_INSPECTION.tsv',sep='\t').to_dict('records')
result['count_files'] = []
for path in S.rglob('*counts.table.tsv'):
    if path.name.startswith('._') or '__MACOSX' in path.parts:
        continue
    with path.open() as handle:
        result['count_files'].append({'path':str(path),'header':next(handle)[:160],'first':next(handle)[:160]})
(W/'SOURCE_INSPECTION.json').write_text(json.dumps(result,indent=2)+'\n')
print(json.dumps({k:v for k,v in result.items() if k not in ['count_files','ROOIJonly_fields','ALL_fields']},indent=2))
print(pd.read_csv(W/'author_ALL/METADATA_FIELDS.tsv',sep='\t').to_string(index=False))
print(result['count_files'][0])
