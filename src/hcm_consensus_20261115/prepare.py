"""Snapshot verified donor matrices and full discovery evidence into a new workspace."""
import hashlib
import json
import shutil
import sys
from datetime import datetime, timezone
from pathlib import Path
import pandas as pd

w=Path(sys.argv[1]).resolve()
old=Path('/home/zhome/hcm-geneformer-analysis/hcm_dge_20261115')
assert w.name=='hcm_consensus_20261115' and w.parent==old.parent and not (w/'input').exists()
def sha(p):
    with p.open('rb') as f:return hashlib.file_digest(f,'sha256').hexdigest()
audit=json.loads((old/'audit/FINAL_AUDIT.json').read_text())
assert audit['status']=='passed'
assert sha(old/'lock/CANDIDATES.tsv')=='d2fc63ed14d76c8d291afd09c067c0339584f38550369c8691248e0126bfa09d'
assert (w/'code/DESIGN.md').is_file()
records=[]
def cp(source,dest):
    before=sha(source); dest.parent.mkdir(exist_ok=True,parents=True)
    shutil.copyfile(source,dest)
    assert sha(dest)==before==sha(source)
    records.append(dict(source=str(source),path=str(dest.relative_to(w)),bytes=dest.stat().st_size,sha256=before))
for source in sorted((old/'pb').iterdir()):
    if source.is_file():cp(source,w/'input/pb'/source.name)
for name in ['discovery.tsv','scp_loo_logFC.tsv.gz','scp_design.tsv','scp_raw_libraries.tsv',
             'scp_cb_libraries.tsv','scp_hhatl_libraries.tsv','scp_hhatl_status.tsv','session_discovery.txt']:
    cp(old/'stats'/name,w/'input/scp'/name)
for folder,name in [('stats','liu.tsv'),('lock','CANDIDATES.tsv'),('lock','LOCK.json'),('audit','FINAL_AUDIT.json')]:
    cp(old/folder/name,w/'input/previous'/name)
cp(old/'coverage/discovery.tsv',w/'input/scp/coverage.tsv')
cp(old/'code/coverage.py',w/'code/coverage_base.py')
annotations=[]
for priority,mod in enumerate(['scp_raw','liu','chin','whole']):
    x=pd.read_csv(w/f'input/pb/{mod}_genes.tsv',sep='\t',keep_default_na=False)
    assert x.gene_id.is_unique and x.gene_id.str.fullmatch(r'ENSG[0-9]+').all()
    x['source']=mod;x['priority']=priority;annotations.append(x)
full=pd.concat(annotations,ignore_index=True)
rows=[]
for gene,x in full.groupby('gene_id',sort=True):
    x=x.sort_values('priority')
    symbols=list(dict.fromkeys(s for value in x.gene_symbol for s in value.split(';') if s))
    rows.append(dict(gene_id=gene,gene_symbol=symbols[0] if symbols else gene,
                     gene_symbol_aliases=';'.join(symbols),source_axes=';'.join(x.source)))
genes=pd.DataFrame(rows)
genes.to_csv(w/'input/GENE_UNIVERSE.tsv',sep='\t',index=False)
record={'created_utc':datetime.now(timezone.utc).isoformat(),'design_commit':'aefa26c',
        'design_sha256':sha(w/'code/DESIGN.md'),'previous_workspace':str(old),
        'prior_outcomes_seen':True,'Liu_role':'joint_exploratory_discovery_not_validation',
        'gene_universe':len(genes),'gene_universe_sha256':sha(w/'input/GENE_UNIVERSE.tsv'),'files':records}
(w/'input/PREPARATION.json').write_text(json.dumps(record,indent=2)+'\n')
print(json.dumps({'files':len(records),'gene_universe':len(genes),'design_commit':'aefa26c'}))
