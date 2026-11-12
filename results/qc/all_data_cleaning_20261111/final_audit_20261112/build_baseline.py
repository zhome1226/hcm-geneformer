import csv
import hashlib
import json
from pathlib import Path

LOCAL = Path(__file__).resolve().parent.parent
OUT = Path(__file__).resolve().parent
REMOTE = '/home/zhome/xiaoe-internship-data/'
configs = [
    ('SCP1303', 'SCP1303_qc_reviewed_final_20261111', 'SCP1303_20261111'),
    ('chin_GSE174691', 'chin_GSE174691_qc_20261111', 'all_data_cleaning_20261111/chin_pair'),
    ('chin_GSE161921', 'chin_GSE161921_qc_20261111', 'all_data_cleaning_20261111/chin_pair'),
    ('liu2023', 'liu2023_qc_20261111', 'all_data_cleaning_20261111/liu2023'),
    ('chin_GSE181764', 'chin_GSE181764_qc_20261111', 'all_data_cleaning_20261111/gse181764'),
    ('GSE138262', 'wehrens_GSE138262_qc_20261111', 'all_data_cleaning_20261111/gse138262'),
    ('GSE203274', 'pediatric_GSE203274_qc_20261111', 'all_data_cleaning_20261111/GSE203274'),
    ('bulk_GSE249925', 'bulk_GSE249925_qc_20261111', 'all_data_cleaning_20261111/GSE249925'),
    ('VCC2025', 'VCC2025_qc_20261111', 'all_data_cleaning_20261111/VCC2025'),
    ('NeurIPS2023', 'NeurIPS2023_qc_20261111', 'all_data_cleaning_20261111/NeurIPS2023'),
    ('NeurIPS2021', 'NeurIPS2021_qc_20261111', 'all_data_cleaning_20261111/NeurIPS2021'),
    ('CausalBench2023', 'CausalBench2023_qc_20261111', 'all_data_cleaning_20261111/CausalBench2023'),
    ('ACDC2017', 'ACDC2017_qc_20261111', 'all_data_cleaning_20261111/ACDC2017'),
]
datasets = {k: {'output': REMOTE+'processed/'+p, 'execution': REMOTE+'_processing/'+e} for k,p,e in configs}
artifacts = []
for f in sorted(LOCAL.glob('*OUTPUT_MANIFEST.tsv')):
    for row in csv.DictReader(f.open(encoding='utf-8-sig'), delimiter='\t'):
        row = {k.lower(): v for k,v in row.items()}
        path = row['server_path']
        ds = next((k for k,d in datasets.items() if path.startswith(d['output']+'/')), None)
        ds = ds or row.get('dataset')
        assert ds in datasets, (f,path,ds)
        artifacts.append({'dataset': ds, 'path': path, 'bytes': int(row['bytes']),
                          'md5': row.get('md5'), 'sha256': row.get('sha256'), 'baseline': f.name})

scp = LOCAL.parent/'hcm-qc-scp1303-20261111/review/final/INDEPENDENT_VERIFICATION.json'
j=json.loads(scp.read_text())
for name,h in j['hashes'].items():
    h={k.lower():v for k,v in h.items()}
    artifacts.append({'dataset':'SCP1303','path':datasets['SCP1303']['output']+'/'+name,
                      'bytes':h['bytes'],'md5':h.get('md5'),'sha256':h.get('sha256'),
                      'baseline':str(scp)})

local_dirs = {'SCP1303':scp.parent, 'chin_GSE174691':LOCAL/'chin_pair/chin_GSE174691',
              'chin_GSE161921':LOCAL/'chin_pair/chin_GSE161921', 'liu2023':LOCAL/'liu2023',
              'chin_GSE181764':LOCAL/'gse181764', 'GSE138262':LOCAL/'gse138262',
              'GSE203274':LOCAL/'gse203274','bulk_GSE249925':LOCAL/'gse249925',
              'VCC2025':LOCAL/'vcc2025','NeurIPS2023':LOCAL/'neurips2023',
              'NeurIPS2021':LOCAL/'neurips2021','CausalBench2023':LOCAL/'causalbench2023',
              'ACDC2017':LOCAL/'acdc2017'}
for ds,folder in local_dirs.items():
    f=folder/'INDEPENDENT_VERIFICATION.json'
    v=json.loads(f.read_text(encoding='utf-8-sig'))
    checks=v.get('checks',v.get('all_checks'))
    flags=list(checks.values()) if isinstance(checks,dict) else [x['passed'] for x in checks]
    assert flags and all(x is True for x in flags),(ds,'failed_checks')
    assert v['status'] in ['passed','complete_verified'],ds
    assert v.get('checks_passed',len(flags))==len(flags),ds
    datasets[ds].update({'verification_sha256':hashlib.sha256(f.read_bytes()).hexdigest(),
                         'independent_checks':len(flags)})

payload={'datasets':datasets,'artifacts':artifacts,
         'excluded':{'dataset':'cox7b_GSE227734','reason':'User excluded from cleaning; missing authoritative gene mapping; no author contact authorized.'}}
(OUT/'BASELINE.json').write_text(json.dumps(payload,indent=2)+'\n',encoding='utf-8')
print(json.dumps({'datasets':len(datasets),'artifacts':len(artifacts),'bytes':sum(a['bytes'] for a in artifacts),
                  'independent_checks':sum(d['independent_checks'] for d in datasets.values())}))
