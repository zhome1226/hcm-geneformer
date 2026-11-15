"""Seal only the predeclared SCP results, before any held-out contrasts are read."""
import hashlib,json,sys
from datetime import datetime,timezone
from pathlib import Path
import numpy as np
import pandas as pd
w=Path(sys.argv[1]);out=w/'lock'
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def tab(p):return pd.read_csv(p,sep='\t',keep_default_na=False,na_values=['NA'])
def save(df,p):df.to_csv(p,sep='\t',index=False,na_rep='NA',float_format='%.15g')
assert (w/'logs/fit_d.exit').read_text().strip()=='0'
assert not (w/'pb/SUPPORT_PSEUDOBULK.json').exists() and not (w/'stats/liu.tsv').exists()
assert not (w/'coverage/support.tsv').exists() and not out.exists()
files=[w/p for p in ['stats/discovery.tsv','pb/scp_raw_genes.tsv','coverage/discovery.tsv','pb/DONOR_ROSTER.tsv','pb/ROSTER_GATE.json','code/DESIGN.md']]
main,annotations,cover=map(tab,files[:3])
assert main.gene_id.is_unique and cover.gene_id.is_unique and annotations.gene_id.is_unique
d=main.merge(annotations,on='gene_id',how='left',validate='1:1').merge(cover,on='gene_id',how='left',validate='1:1')
assert d.gene_symbol.notna().all() and d.hcm_donors.eq(15).all() and d.hcm_cells.eq(54344).all()
assert d.in_v1_vocab.isin([False,True]).all() and d.FDR.between(0,1).all()
d['primary_DEG']=d.FDR.lt(.05) & d.logFC.abs().ge(.5)
d['cellbender_same_direction']=np.sign(d.cb_logFC)==np.sign(d.logFC)
d['loo_stable']=d.loo_same_direction_fraction.eq(1)&d.loo_min_abs_logFC.ge(.25)
d['robust_DEG']=d.primary_DEG & d.cellbender_same_direction & d.loo_stable
d['token_eligible']=d.in_v1_vocab & d.token_id.notna() & d.hcm_donors_token_ge10pct.ge(8)
d['eligible']=d.robust_DEG & d.token_eligible
d['association_direction']=np.where(d.logFC>0,'higher_in_HCM',np.where(d.logFC<0,'lower_in_HCM','zero'))
ordered=d.loc[d.eligible].assign(absolute_effect=lambda x:x.logFC.abs()).sort_values(['FDR','absolute_effect','gene_id'],ascending=[True,False,True],kind='mergesort').drop(columns='absolute_effect')
cand=ordered.head(200).copy();cand.insert(0,'rank',np.arange(1,len(cand)+1))
first=['gene_id','gene_symbol','association_direction','logFC','FDR','cb_logFC','cb_FDR','loo_min_abs_logFC','loo_same_direction_fraction','hcm_donors_token_ge10pct','hcm_token_fraction','hcm_expr_fraction','token_id','hhatl_logFC','hhatl_FDR']
cand=cand[['rank']+first+[x for x in cand if x not in ['rank',*first]]]
d['locked_rank']=d.gene_id.map(dict(zip(cand.gene_id,cand['rank']))).astype('Int64')
out.mkdir();save(cand,out/'CANDIDATES.tsv');save(d,out/'DISCOVERY_ALL_TESTED.tsv.gz')
digest=sha(out/'CANDIDATES.tsv');(out/'CANDIDATES.tsv.sha256').write_text(digest+'  CANDIDATES.tsv\n',encoding='ascii')
gates=[('tested',np.ones(len(d),dtype=bool)),('primary_DEG',d.primary_DEG),('same_direction_CB_among_DEG',d.primary_DEG&d.cellbender_same_direction),('robust_after_LODO',d.robust_DEG),('robust_and_V1_eligible',d.eligible),('locked_top200',d.locked_rank.notna())]
funnel=pd.DataFrame([dict(stage=k,genes=int(m.sum()),higher_HCM=int((m & d.logFC.gt(0)).sum()),lower_HCM=int((m & d.logFC.lt(0)).sum())) for k,m in gates]);save(funnel,out/'FUNNEL.tsv')
receipt=dict(locked_utc=datetime.now(timezone.utc).isoformat(),candidate_sha256=digest,selected=len(cand),genes_tested=len(d),cutoff_commit='5019f70',inputs={str(x.relative_to(w)):sha(x) for x in files},no_validation_outcomes_at_lock=True)
(out/'LOCK.json').write_text(json.dumps(receipt,indent=2)+'\n')
for n in ['CANDIDATES.tsv','CANDIDATES.tsv.sha256','LOCK.json']: (out/n).chmod(0o444)
print(funnel.to_string(index=False),flush=True)
