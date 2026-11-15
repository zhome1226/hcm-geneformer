"""Independent post-lock reconciliation from donor matrices, LODO folds and raw token tallies.

Usage: python audit.py WORK [--final RELATIVE_TSV] [--out RELATIVE_JSON]
Never call before the downstream jobs have finished. No model refits and no changes to inputs.
"""
import argparse
import hashlib
import json
import sys
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import pandas as pd


# From the pre-outcome panel/inventory, keyed by study (n, donors, HCM donors, role).
PANEL = {
    'SCP1303': (114711, 31, 15, 'primary_discovery'),
    'liu2023': (8547, 12, 10, 'held_out_directional_replication'),
    'chin_GSE174691': (10150, 9, 9, 'accession_confounded_support'),
    'chin_GSE161921': (3534, 4, 0, 'accession_confounded_support'),
    'GSE138262_whole_cell': (3595, 12, 5, 'study_technology_confounded_support'),
    'GSE255296': (7516, 4, 4, 'HCM_only_case_coverage'),
    'GSE292754': (3306, 2, 2, 'pediatric_HCM_only_case_coverage'),
}
DESIGN_SHA = '72fca3947e2ad32ce3a246d79b48ac09ace31247237f1a807d57f4517d4b5b16'
SOURCE_FILES = {
    'stats/discovery.tsv', 'pb/scp_raw_genes.tsv', 'coverage/discovery.tsv',
    'pb/DONOR_ROSTER.tsv', 'pb/ROSTER_GATE.json', 'code/DESIGN.md',
}
PK = ['dataset', 'disease', 'gene_id']


def require(test, msg):
    if not bool(test):
        raise ValueError(msg)


def same(actual, expected, msg):
    require(np.array_equal(np.asarray(actual), np.asarray(expected)), msg)


def near(actual, expected, msg, rtol=3e-9, atol=1e-12):
    a, b = map(lambda x: np.asarray(x, dtype=float), (actual, expected))
    require(a.shape == b.shape and np.allclose(a, b, rtol=rtol, atol=atol, equal_nan=True), msg)


def opts(s):
    # Never use astype(bool) on strings ("FALSE" would become True).
    b = s.map({True: True, False: False, 'TRUE': True, 'FALSE': False, 'True': True, 'False': False})
    require(b.notna().all(), 'Unrecognized boolean in ' + str(s.name))
    return b.to_numpy(dtype=bool)


def sha(path):
    with path.open('rb') as fp:
        return hashlib.file_digest(fp, 'sha256').hexdigest()


class Audit:
    def __init__(self, w):
        self.w = w
        self.done = []
        self.phase = ''
        self.record = {}

    def tab(self, rel):
        return pd.read_csv(self.w/rel, sep='\t', keep_default_na=False, na_values=['NA',''])

    def json(self, rel):
        return json.loads((self.w/rel).read_text(encoding='utf8'))

    def run(self, name, fn, *args):
        self.phase = name
        result = fn(*args)
        self.done.append(name)
        return result

    def seal(self):
        lock = self.tab('lock/CANDIDATES.tsv')
        self.ids = list(lock.gene_id)
        self.lock = lock.set_index('gene_id', verify_integrity=True)
        digest = sha(self.w/'lock/CANDIDATES.tsv')
        for f in ('lock/CANDIDATES.tsv.sha256', 'coverage/support.lock.sha256'):
            same((self.w/f).read_text().split(), [digest,'CANDIDATES.tsv'], f)
        receipt = self.json('lock/LOCK.json')
        v = self.tab('stats/validation_design.tsv')
        require(len(v)==1 and v.lock_sha256[0]==digest and v.locked_genes[0]==len(lock), 'Validation seal')
        require(receipt['candidate_sha256']==digest and receipt['selected']==len(lock)
                and receipt['no_validation_outcomes_at_lock'] is True and receipt['cutoff_commit']=='5019f70',
                'Discovery receipt')
        require(datetime.fromisoformat(receipt['locked_utc']).utcoffset().total_seconds()==0, 'Lock timestamp')
        require(set(receipt['inputs'])==SOURCE_FILES and
                all(sha(self.w/f)==h for f,h in receipt['inputs'].items()), 'Lock provenance hashes')
        require(hashlib.sha256((self.w/'code/DESIGN.md').read_bytes().replace(b'\r\n', b'\n')).hexdigest()==DESIGN_SHA,
                'Predeclared analysis rules changed')
        same(lock['rank'], np.arange(1, len(lock)+1), 'Candidate ranks')
        require(len(lock)<=200 and lock.gene_id.str.fullmatch(r'ENSG[0-9]+').all() and
                lock.token_id.notna().all() and lock.token_id.is_unique, 'Lock identity')
        self.record.update(lock_sha256=digest, locked=len(lock), locked_utc=receipt['locked_utc'])

    def roster(self):
        ro = self.tab('pb/DONOR_ROSTER.tsv')
        k = ['dataset','donor_id']
        require(not ro.duplicated(k).any() and len(ro)==74 and ro.nuclei.sum()==151359, 'Panel totals')
        for ds,(n,d,h,role) in PANEL.items():
            x=ro[ro.dataset.eq(ds)]
            require(len(x)==d and x.disease.eq('HCM').sum()==h and x.nuclei.sum()==n
                    and set(x.role)=={role} and set(x.disease)<={'HCM','Control'}, ds+' roster')
        same(ro.model_train_test_split.unique(),['unassigned'],'Model split not assigned')
        same(ro.cross_study_person_match.unique(),['unresolved'],'Person-level cross-study status')
        gate=self.json('pb/ROSTER_GATE.json')
        require(gate['roster_sha256']==sha(self.w/'pb/DONOR_ROSTER.tsv')
                and gate['donor_rows']==len(ro) and gate['cells']==int(ro.nuclei.sum()),'Roster receipt')
        s=self.tab('pb/SAMPLE_ROSTER.tsv')
        require(not s.duplicated(['dataset','sample_id']).any(), 'A capture occurs in multiple donor rows')
        r=ro.set_index(k).sort_index()
        a=s.groupby(k).agg(nuclei=('nuclei','sum'),n_samples=('sample_id','nunique'))
        require(a.index.equals(r.index), 'Sample donor manifest')
        for c in a: same(a[c],r[c],c+' sample-to-donor reconciliation')
        lab=s.merge(ro[k+['disease']],on=k,validate='m:1',suffixes=('','_donor'))
        same(lab.disease,lab.disease_donor,'Per-capture diagnosis')
        alias=s[(s.dataset=='liu2023') & s.sample_id.isin(['CL0412A','CL0412B'])]
        require(len(alias)==2 and alias.donor_id.nunique()==1 and alias.disease.eq('Control').all(),'Liu repeated control capture')
        selected=ro[ro.dataset.eq('SCP1303')].set_index('donor_id')
        require(selected.query("disease=='Control'").source_disease.eq('NF').all(), 'SCP controls should be non-failing')
        require(selected.query("disease=='HCM'").nuclei.sum()==54344, 'SCP HCM cells')
        design=self.tab('stats/scp_design.tsv').set_index('donor_id')
        require(set(design.index)==set(selected.index),'SCP discovery donor identities')
        m=selected.loc[design.index]
        same(list(design),['(Intercept)','age_centered','sexmale','diseaseHCM'],'Covariate design')
        age=m.age.astype(float)
        sex=m.sex.str.lower()
        require(set(sex)<= {'f','female','m','male'}, 'Unknown sex')
        near(design.values,np.column_stack([np.ones(len(m)),age-age.mean(),sex.isin(['male','m']),m.disease.eq('HCM')]),'Rebuild design')
        require(np.linalg.matrix_rank(design)==4,'Rank-deficient design')
        self.ro=ro
        self.record.update(donor_study_rows=len(ro),total_cells=int(ro.nuclei.sum()),datasets=len(PANEL))

    def counts(self):
        self.matrices={}
        self.filtered={}
        for mod in ('scp_raw','scp_cb','scp_hhatl','liu','chin','whole'):
            x=self.tab(f'pb/{mod}.tsv.gz').set_index('gene_id',verify_integrity=True)
            m=self.tab(f'pb/{mod}_donors.tsv').set_index('donor_id',verify_integrity=True)
            same(x.columns,m.index,mod+' donors')
            same(self.tab(f'pb/{mod}_genes.tsv').gene_id,x.index,mod+' feature axis')
            require(np.issubdtype(x.to_numpy().dtype,np.integer) and x.to_numpy().min()>=0,mod+' integral counts')
            lib=x.sum()
            same(lib,m.pseudobulk_library_size,mod+' sum of columns')
            if mod=='scp_hhatl':
                sub=self.tab('pb/SCP_SUBTYPE_RETENTION.tsv')
                eligible=sub[(sub.SubCluster=='CM-HHATL')&(sub.nuclei>=50)].set_index('donor_id')
                require(set(eligible.index)==set(m.index) and m.groupby('disease').size().min()>=5, 'Subtype donor support')
                same(eligible.loc[m.index,'nuclei'],m.nuclei,'Subtype eligible nuclei')
            if mod in ('whole','liu'):
                mapping=self.tab(f'pb/{mod}_feature_audit.tsv.gz')
                valid=mapping.gene_id.fillna('').str.fullmatch(r'ENSG[0-9]+').to_numpy() & opts(mapping.measured_in_all_included_samples)
                same(opts(mapping.included_in_model),valid,mod+' native and mapping mask')
                require(set(mapping.loc[valid,'gene_id'])==set(x.index),mod+' pseudobulk mapping')
            if mod=='chin':
                balance=self.tab('pb/chin_common_axis_count_balance.tsv').set_index('donor_id').reindex(lib.index)
                same(balance.after_intersection,lib,'Chin shared gene axis library total')
                same(balance.before_intersection-balance.after_intersection,balance.not_shared_count,'Chin removed counts')
            n=m.disease.value_counts().min()
            min_n=n if n<=10 else 10+(n-10)*.7
            # edgeR filterByExpr.default on unnormalized counts; 1e-14 is the package boundary tolerance.
            threshold=10/np.median(lib)*1e6
            keep=(x.div(lib,axis=1).mul(1e6).ge(threshold).sum(axis=1)>=min_n-1e-14)&(x.sum(axis=1)>=15-1e-14)
            self.filtered[mod]=list(x.index[keep])
            self.matrices[mod]=x
        self.record['source_genes']={k:len(v) for k,v in self.matrices.items()}

    def discovery(self):
        raw=self.tab('stats/discovery.tsv').set_index('gene_id',verify_integrity=True)
        same(raw.index,self.filtered['scp_raw'],'Recomputed discovery expression filter')
        p=raw.PValue.to_numpy()
        require(np.isfinite(p).all() and ((p>=0)&(p<=1)).all(),'P values')
        ix=np.argsort(p); bh=np.minimum.accumulate((p[ix]*len(p)/np.arange(1,len(p)+1))[::-1])[::-1]
        near(raw.FDR.iloc[ix],np.minimum(bh,1),'Independent genome-wide BH')
        folds=self.tab('stats/scp_loo_logFC.tsv.gz').set_index('gene_id',verify_integrity=True)
        same(folds.index,raw.index,'LODO gene axis')
        same(folds.columns,self.matrices['scp_raw'].columns,'31 single donor omissions')
        a=folds.to_numpy(); effect=raw.logFC.to_numpy(); delta=np.abs(a-effect[:,None])
        require(np.isfinite(a).all() and len(folds.columns)==31,'Finite LODO matrix')
        share=(np.sign(a)==np.sign(effect[:,None])).mean(axis=1)
        floor=np.abs(a).min(axis=1)
        for col,expected in [('loo_same_direction_fraction',share),('loo_min_abs_logFC',floor),
                             ('loo_min_logFC',a.min(axis=1)),('loo_max_logFC',a.max(axis=1)),
                             ('loo_max_abs_delta',delta.max(axis=1))]:
            near(raw[col],expected,col)
        # R may break a near tie differently after serializing; a claimed maximum must attain it.
        claimed=folds.columns.get_indexer(raw.loo_max_influence_omitted_donor)
        require((claimed>=0).all(),'Influential donor not in design')
        near(delta[np.arange(len(raw)),claimed],delta.max(axis=1),'Influential donor delta')
        cov=self.tab('coverage/discovery.tsv').set_index('gene_id',verify_integrity=True).loc[raw.index]
        detail=self.tab('coverage/discovery_donors.tsv.gz')
        self.cover(detail,self.tab('coverage/discovery.tsv').rename(columns=lambda s:s.removeprefix('hcm_')),True)
        eligible_tech=opts(cov.in_v1_vocab)&np.isfinite(cov.token_id)&(cov.hcm_donors_token_ge10pct>=8)
        deg=(raw.FDR.to_numpy()<.05)&(np.abs(effect)>=.5)
        cb=np.sign(raw.cb_logFC.to_numpy())==np.sign(effect)
        robust=deg&cb&(share==1)&(floor>=.25)
        eligible=np.flatnonzero(robust&eligible_tech)
        sorted_pool=sorted(eligible,key=lambda j:(raw.FDR.iloc[j],-abs(effect[j]),raw.index[j]))
        predicted=list(raw.index[sorted_pool[:200]])
        same(self.ids,predicted,'Independent candidate membership and order')
        near(self.lock.logFC,raw.loc[predicted,'logFC'],'Locked disease coefficient')
        same(self.lock.association_direction,np.where(self.lock.logFC>0,'higher_in_HCM','lower_in_HCM'),'Association labels')
        anno=self.tab('pb/scp_raw_genes.tsv').set_index('gene_id')
        same(self.lock.gene_symbol,anno.loc[predicted,'gene_symbol'],'Gene nomenclature')
        near(self.lock.token_id,cov.loc[predicted,'token_id'],'Frozen token IDs')
        funnel=self.tab('lock/FUNNEL.tsv').set_index('stage')
        tests={'tested':np.repeat(True,len(raw)), 'primary_DEG':deg, 'same_direction_CB_among_DEG':deg&cb,
               'robust_after_LODO':robust,'robust_and_V1_eligible':robust&eligible_tech.to_numpy(),
               'locked_top200':raw.index.isin(predicted)}
        same(list(funnel.index),list(tests),'Funnel order')
        for name,mask in tests.items():
            same(funnel.loc[name,['genes','higher_HCM','lower_HCM']], [np.sum(mask), np.sum(mask&(effect>0)), np.sum(mask&(effect<0))],name)
        self.record.update(genomewide_tests=len(raw),robust=int(robust.sum()),eligible=len(sorted_pool),LODO_folds=31)

    def cover(self, detail, supplied, discovery=False):
        key=PK+['donor_id']
        require(not detail.duplicated(key).any() and not supplied.duplicated(PK).any(),'Duplicate coverage keys')
        ro=self.ro
        if discovery:
            ro=ro[(ro.dataset=='SCP1303')&(ro.disease=='HCM')]
            loci=set(self.matrices['scp_raw'].index)
        else:
            loci=set(self.ids)
        # Each study-by-diagnosis contains every expected gene at every expected donor.
        donors=ro.set_index(['dataset','donor_id','disease']).nuclei.sort_index()
        g=detail.groupby(['dataset','donor_id','disease'])
        require(g.size().index.equals(donors.index) and g.size().eq(len(loci)).all() and set(detail.gene_id)==loci,
                'Full coverage grid')
        same(g.donor_cells.min(),donors,'Donor coverage denominator')
        same(g.donor_cells.max(),donors,'Donor coverage denominator')
        n=detail.donor_cells; t=detail.token_cells; h=detail.expressing_cells; nc=detail.native_cells
        require(((0<=t)&(t<=nc)&(0<=h)&(h<=nc)&(nc<=n)&(h<=detail.count_sum)).all(), 'Coverage arithmetic bounds')
        require((t[~detail.dataset.eq('SCP1303')]<=h[~detail.dataset.eq('SCP1303')]).all(),'Native non-SCP token without expression')
        near(detail.expr_fraction_native,h/nc.replace(0,np.nan),'Missing vs native per donor')
        near(detail.token_fraction_all_cells,t/n,'Per-donor token fraction')
        added=detail.assign(native=nc.gt(0),expressed=h.gt(0),positive=t.gt(0),pct5=20*t>=n,
                            pct10=10*t>=n,exact_token_fraction=t/n)
        gr=added.groupby(PK,sort=True)
        sums=gr.agg(donors=('donor_id','size'),cells=('donor_cells','sum'),native_cells=('native_cells','sum'),
                    donors_native=('native','sum'),expressing_cells=('expressing_cells','sum'),count_sum=('count_sum','sum'),
                    donors_expressing=('expressed','sum'),token_cells=('token_cells','sum'),
                    donors_token_positive=('positive','sum'),donors_token_ge5pct=('pct5','sum'),
                    donors_token_ge10pct=('pct10','sum'),median_donor_token_fraction=('exact_token_fraction','median'))
        sums['expr_fraction']=sums.expressing_cells/sums.native_cells.replace(0,np.nan)
        sums['token_fraction']=sums.token_cells/sums.cells
        reported=supplied.set_index(PK).sort_index()
        require(sums.index.equals(reported.index),'Coverage aggregate rows')
        for col in sums: near(reported[col],sums[col],'Coverage aggregate '+col)
        mods={'SCP1303':'scp_raw','liu2023':'liu','GSE138262_whole_cell':'whole',
              'chin_GSE161921':'chin','chin_GSE174691':'chin'}
        for ds,sub in detail.groupby('dataset'):
            if ds not in mods: continue
            x=self.matrices[mods[ds]]
            pivot=sub.pivot(index='gene_id',columns='donor_id',values='count_sum')
            axis=pivot.index.intersection(x.index)
            require(len(axis)>0,ds+' common expression axis')
            same(pivot.loc[axis].values,x.loc[axis,pivot.columns].values,ds+' coverage counts vs pseudobulk')
        if not discovery:
            same(supplied.token_id.values,self.lock.token_id.loc[supplied.gene_id].values,'Support token IDs')

    def final(self, path):
        f=self.tab(path).set_index('gene_id',verify_integrity=True)
        same(f.index,self.ids,'Final gene order')
        near(f['rank'],self.lock['rank'],'Final frozen ranks')
        for c in ('gene_symbol','association_direction','token_id'):
            same(f[c],self.lock[c],'Final '+c)
        for c in f:
            if c.startswith('scp_') and c!='scp_hhatl_direction':
                original=c.removeprefix('scp_')
                require(original in self.lock,'Unexpected discovery field '+c)
                a,b=f[c],self.lock[original]
                if pd.api.types.is_numeric_dtype(b): near(a,b,'Final '+c)
                else: same(a,b,'Final '+c)
        for mod in ('liu','chin','whole'):
            s=self.tab('stats/'+mod+'.tsv').set_index('gene_id').loc[self.ids]
            for c in s:
                a,b=f[mod+'_'+c],s[c]
                if pd.api.types.is_numeric_dtype(b): near(a,b,'Final '+mod+' '+c)
                else: same(a.fillna('NA'),b.fillna('NA'),'Final '+mod+' '+c)
        cov=self.tab('coverage/support.tsv')
        for ds,prefix in [('liu2023','liu'),('GSE255296','gse255296'),('GSE292754','gse292754')]:
            s=cov[(cov.dataset==ds)&(cov.disease=='HCM')].set_index('gene_id').loc[self.ids]
            for c in f:
                if c.startswith(prefix+'_HCM_') and c.removeprefix(prefix+'_HCM_') in s:
                    near(f[c],s[c.removeprefix(prefix+'_HCM_')],'Final coverage '+c)
        def direction(a,b):
            if not np.isfinite(a): return 'not_evaluable'
            if a==0: return 'zero'
            return 'same' if a*b>0 else 'opposite'
        expected=[]
        for _,r in f.iterrows():
            base=direction(r.liu_logFC,r.scp_logFC)
            one=direction(r.liu_logFC_omit_Liu2023_CL0412,r.scp_logFC)
            two=direction(r.liu_logFC_omit_Liu2023_CL2,r.scp_logFC)
            require((base,one,two)==(r.liu_direction,r.liu_direction_omit_Liu2023_CL0412,r.liu_direction_omit_Liu2023_CL2),'Liu direction labels')
            category=base
            if base=='same':
                category='same_both_single_controls' if one==two=='same' else 'same_sensitivity_not_evaluable' if 'not_evaluable' in (one,two) else 'same_control_sensitive'
            expected.append(category)
            require(direction(r.scp_hhatl_logFC,r.scp_logFC)==r.scp_hhatl_direction,'Subtype direction')
            for mod in ('chin','whole'):
                require(direction(r[mod+'_descriptive_log2_cpm_difference'],r.scp_logFC)==r[mod+'_descriptive_direction'],mod+' final direction')
            reason='feature_absent_or_unmapped' if not r.liu_source_feature_present else 'below_expression_filter' if not r.liu_passes_expression_filter else 'evaluated'
            require(reason==r.liu_not_evaluable_reason,'Missingness explanation')
        same(f.liu_evidence_category,expected,'Evidence categories')
        flag=f.liu_direction.eq('same')&f.liu_FDR.lt(.05)
        same(f.liu_FDR_lt05_same_direction,flag,'Liu significance with direction')
        s=self.json('final/SUMMARY.json')
        require(s['candidate_sha256']==self.record['lock_sha256'] and s['locked_candidates']==len(f),'Final summary seal')
        require(s['liu_evidence_category_counts']==pd.Series(expected).value_counts().to_dict(),'Final category totals')
        require(s['liu_direction_counts']==f.liu_direction.value_counts().to_dict(),'Final direction totals')
        require(s['higher_in_HCM']==int(f.scp_logFC.gt(0).sum()) and s['lower_in_HCM']==int(f.scp_logFC.lt(0).sum()),'Final disease direction totals')
        require(s['liu_FDR_lt05_same_direction']==int(flag.sum()) and s['liu_FDR_lt05_opposite']==int((f.liu_direction.eq('opposite')&f.liu_FDR.lt(.05)).sum()),'Final FDR summary')
        require(s['liu_evaluable']==int(f.liu_logFC.notna().sum()),'Final estimability summary')
        for ds,sm in s['coverage'].items():
            x=cov[(cov.dataset==ds)&(cov.disease=='HCM')]
            require(sm['HCM_donors']==int(x.donors.iloc[0]) and sm['HCM_cells']==int(x.cells.iloc[0]),ds+' final donor totals')
            for key,mask in [('genes_measured_any_cells',x.native_cells.gt(0)),('genes_token_in_any_cells',x.token_fraction.gt(0)),('genes_token_ge10pct_in_at_least_half_donors',2*x.donors_token_ge10pct>=x.donors)]:
                require(sm[key]==int(mask.sum()),ds+' coverage summary '+key)
            near(sm['median_gene_token_fraction'],x.token_fraction.median(),ds+' coverage median')
        self.record['final_evidence_sha256']=sha(self.w/path)
        self.record['liu_categories']=s['liu_evidence_category_counts']

    def support(self):
        d=self.tab('coverage/support_donors.tsv.gz')
        self.cover(d,self.tab('coverage/support.tsv'))
        validation=self.tab('stats/validation_design.tsv')
        require(validation.liu_tested_genomewide[0]==len(self.filtered['liu']),'Liu genomewide denominator')
        for mod in ('liu','chin','whole'):
            s=self.tab('stats/'+mod+'.tsv').set_index('gene_id',verify_integrity=True)
            same(s.index,self.ids,mod+' frozen scope and order')
            exists=s.index.isin(self.matrices[mod].index); keep=s.index.isin(self.filtered[mod])
            same(opts(s.source_feature_present),exists,mod+' feature available vs missing')
            same(opts(s.passes_expression_filter),keep,mod+' expression filter')
            if mod=='liu':
                require(len(s.filter(like='logFC_omit_').columns)==2,'Liu exactly two sensitivity fits')
                omit=s.filter(like='logFC_omit_')
                same(sorted(omit.columns),['logFC_omit_Liu2023_CL0412','logFC_omit_Liu2023_CL2'],'Liu omitted control identities')
                for c in ['logFC','PValue','FDR',*omit]:
                    same(np.isfinite(s[c]),keep,'Liu estimability '+c)
                require(s.loc[keep,'FDR'].between(0,1).all() and (s.loc[keep,'FDR']+1e-14 >= s.loc[keep,'PValue']).all(),
                        'Liu adjusted P values')
            else:
                require(not set(s)&{'FDR','PValue','logFC','pvalue','qvalue'},mod+' no disease-significance claim')
                same(np.isfinite(s.descriptive_log2_cpm_difference),keep,mod+' descriptive availability')
                near(s.descriptive_log2_cpm_difference,s.median_log2CPM_HCM-s.median_log2CPM_Control,mod+' descriptive contrast')
        self.record['support_donor_gene_rows']=len(d)


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('work',type=Path)
    parser.add_argument('--final',type=Path)
    parser.add_argument('--out',type=Path,default=Path('audit/FINAL_AUDIT.json'))
    args=parser.parse_args()
    w=args.work.resolve()
    for f in ('lock/CANDIDATES.tsv','stats/validation_design.tsv','coverage/support.lock.sha256'):
        if not (w/f).is_file(): sys.exit(f'Downstream output missing; audit has not started: {f}')
    a=Audit(w)
    dest=w/args.out
    dest.parent.mkdir(parents=True,exist_ok=True)
    try:
        for label,fn in [('verify frozen inputs, study design and validation seal',a.seal),
                         ('reconcile seven study roles, samples and donor covariates',a.roster),
                         ('count conservation, native masks and independent expression filtering',a.counts),
                         ('recalculate discovery coverage, BH, 31 omissions, thresholds and rankings',a.discovery),
                         ('reconcile held-out contrasts and seven-source actual sequence coverage',a.support)]:
            a.run(label,fn)
        if args.final:
            a.run('reconcile final evidence table',a.final,args.final)
        error=None
    except Exception as exc:
        error=f'{a.phase}: {type(exc).__name__}: {exc}'
    report={'status':'failed' if error else 'passed','mode':'post-lock independent reconciliation; no model refit',
            'completed_checks':a.done,'metrics':a.record,'failure':error,'utc':datetime.now(timezone.utc).isoformat()}
    tmp=dest.with_suffix('.partial')
    tmp.write_text(json.dumps(report,indent=2,allow_nan=False)+'\n',encoding='utf8');tmp.replace(dest)
    print(json.dumps(report,allow_nan=False))
    if error:sys.exit(1)


if __name__=='__main__':
    main()
