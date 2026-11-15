"""Donor pseudobulk on QC core counts. Roster/discovery never read held-out counts.

Usage: python pseudobulk.py WORK roster|discovery|support
Inputs are frozen in previously verified V1-preparation manifests. All data opens are read-only.
"""
import argparse
import hashlib
import json
import re
from pathlib import Path

import h5py
import numpy as np
import pandas as pd
import scipy.sparse as sp
from anndata.io import read_elem

HOME = Path('/home/zhome')
OLD = HOME/'xiaoe-internship-data/tokenized/geneformer_v1_20261112/prepared_inputs'
NEW = HOME/'hcm-geneformer-data/tokenized/geneformer_v1_20261114/prepared_inputs'
ROLES = dict(SCP1303='primary_discovery', liu2023='held_out_directional_replication',
             chin_GSE174691='accession_confounded_support', chin_GSE161921='accession_confounded_support',
             GSE138262_whole_cell='study_technology_confounded_support',
             GSE255296='HCM_only_case_coverage', GSE292754='pediatric_HCM_only_case_coverage')
N_CELLS = dict(SCP1303=114711, liu2023=8547, chin_GSE174691=10150, chin_GSE161921=3534,
               GSE138262_whole_cell=3595, GSE255296=7516, GSE292754=3306)
N_DONORS = dict(SCP1303=31, liu2023=12, chin_GSE174691=9, chin_GSE161921=4,
                GSE138262_whole_cell=12, GSE255296=4, GSE292754=2)
LAYERS = dict(scp_raw='layers/cellranger_raw', scp_cb='layers/counts', scp_hhatl='layers/cellranger_raw',
              liu='layers/counts', chin='layers/counts', whole='layers/counts')


def sha(path):
    with path.open('rb') as stream:
        return hashlib.file_digest(stream, 'sha256').hexdigest()


def txt(out, name, frame):
    frame.to_csv(out/name, sep='\t', index=False, na_rep='NA')


def save_json(out, name, obj):
    (out/name).write_text(json.dumps(obj, indent=2, ensure_ascii=False)+'\n')


def inputs():
    sources = {}
    for root in (OLD, NEW):
        manifest = json.loads((root/'PREPARATION_VERIFICATION.json').read_text())
        assert manifest['status'] == 'passed'
        for r in manifest['datasets']:
            ds = r['dataset']
            assert ds not in sources and r['cells'] == N_CELLS[ds]
            path = Path(r['source'])
            assert path.stat().st_size == r['source_bytes']
            sources[ds] = (path, r['source_sha256'])
    assert set(sources) == set(ROLES)
    return sources


def obs_frame(f, ds):
    obs = read_elem(f['obs'])
    assert obs.index.is_unique and len(obs) == N_CELLS[ds]
    for col in ('donor_id', 'sample_id', 'disease', 'label'):
        assert obs[col].notna().all()
    obs['donor_id'] = obs.donor_id.astype(str)
    obs['sample_id'] = obs.sample_id.astype(str)
    obs['source_disease'] = obs.disease.astype(str)
    obs['disease'] = obs.source_disease.replace({'NF':'Control'})
    assert obs.disease.isin(['Control','HCM']).all()
    assert np.array_equal(obs.label.astype(int), obs.disease.map({'Control':0,'HCM':1}))
    assert obs.groupby('donor_id').disease.nunique().max() == 1
    assert obs.groupby('sample_id').donor_id.nunique().max() == 1
    assert obs.donor_id.nunique() == N_DONORS[ds]
    if ds == 'liu2023':
        aliases = obs.loc[obs.sample_id.isin(['CL0412A','CL0412B'])]
        assert set(aliases.sample_id)=={'CL0412A','CL0412B'} and aliases.donor_id.nunique()==1
    return obs


def values(rows, cols, default='unknown'):
    for col in cols:
        if col in rows:
            val = sorted({str(v) for v in rows[col].dropna() if str(v).strip() not in ['', 'nan','NA','unknown']})
            if val:
                return ';'.join(val)
    return default


def donors(ds, obs):
    rec=[]
    for donor, rows in obs.groupby('donor_id', sort=True):
        disease=rows.disease.iloc[0]
        age=values(rows,['age','age_years'])
        sex=values(rows,['sex']).lower()
        assert ';' not in age and ';' not in sex
        rec.append(dict(dataset=ds, donor_id=donor, donor_original=values(rows,['donor_original'],''),
                        disease=disease, source_disease=values(rows,['source_disease']), nuclei=len(rows),
                        n_samples=rows.sample_id.nunique(), sample_ids=values(rows,['sample_id']),
                        age=age, sex=sex, tissue=values(rows,['tissue','tissue_site'],'left ventricle' if ds=='SCP1303' else 'unknown'),
                        technology=values(rows,['technology','assay'],'snRNA-seq' if ds=='SCP1303' else 'unknown'),
                        study_origin=values(rows,['study_origin'],ds), role=ROLES[ds],
                        source_clinical_stage=values(rows,['clinical_stage']),
                        age_stratum=values(rows,['age_stratum','author_age_group']),
                        cross_study_person_match='unresolved', model_train_test_split='unassigned'))
    return pd.DataFrame(rec)


def roster(sources, out):
    blocks=[]; samples=[]; subs=[]; manifests=[]
    for ds,(path,expected) in sources.items():
        with h5py.File(path, 'r') as f:
            obs=obs_frame(f,ds)  # metadata only; no layers or X access
            blocks.append(donors(ds,obs))
            sf=obs.groupby(['donor_id','sample_id','disease'], sort=True).size().rename('nuclei').reset_index()
            sf.insert(0,'dataset',ds); samples.append(sf)
            if ds=='SCP1303':
                st=obs.groupby(['donor_id','disease','SubCluster'],observed=True).size().rename('nuclei').reset_index()
                subs.append(st)
            manifests.append(dict(dataset=ds, path=str(path), bytes=path.stat().st_size, expected_sha256=expected,
                                  role=ROLES[ds], n_cells=len(obs), n_donors=obs.donor_id.nunique()))
    df=pd.concat(blocks,ignore_index=True)
    assert len(df)==sum(N_DONORS.values()) and df.nuclei.sum()==sum(N_CELLS.values())
    txt(out,'DONOR_ROSTER.tsv',df);txt(out,'SAMPLE_ROSTER.tsv',pd.concat(samples,ignore_index=True))
    txt(out,'SCP_SUBTYPE_RETENTION.tsv',pd.concat(subs,ignore_index=True))
    txt(out,'INPUT_SOURCES.tsv',pd.DataFrame(manifests))
    aliases={}
    for row in df.itertuples():
        for label, val in [('donor_id',row.donor_id),('donor_original',row.donor_original)]:
            for name in val.split(';'):
                if name:aliases.setdefault(name,[]).append((row.dataset,row.donor_id,label))
    conflicts=[dict(alias=k,matches=';'.join('|'.join(x) for x in v))
               for k,v in aliases.items() if len(set(x[0] for x in v))>1]
    txt(out,'CROSS_STUDY_ALIAS_FLAGS.tsv',pd.DataFrame(conflicts,columns=['alias','matches']))
    save_json(out,'ROSTER_GATE.json',{'role':'DGE_only_not_model_split','donor_rows':len(df),
                                      'cells':int(df.nuclei.sum()),'roster_sha256':sha(out/'DONOR_ROSTER.tsv'),
                                      'alias_scope':'exact available source IDs/aliases only; absence does not establish cross-study independence',
                                      'exact_cross_study_aliases_to_review':len(conflicts)})


def get_axis(f,ds,obs):
    var=read_elem(f['var'])
    assert var.index.is_unique
    if ds=='liu2023':
        lookup=pd.read_csv(OLD/'liu2023_GENE_MAPPING.tsv.gz',sep='\t',dtype=str,keep_default_na=False)
        assert list(lookup['source_feature'])==list(var.index.astype(str))
        with h5py.File(OLD/'liu2023.h5ad','r') as prepared:
            prepared_var=read_elem(prepared['var'])
            assert np.array_equal(lookup.ensembl_id, prepared_var.ensembl_id.astype(str))
            assert prepared_var.index.equals(var.index)
        ids=pd.Index(lookup.ensembl_id)
        status=lookup.mapping_status.to_numpy()
    else:
        ids=pd.Index(var['ensembl_id'] if 'ensembl_id' in var else var['gene_id'] if 'gene_id' in var else var.index).astype(str)
        status=np.repeat('native_id',len(ids))
    symbols=var.gene_symbol.astype(str).to_numpy()
    native=np.ones(len(ids),bool)
    if ds=='GSE138262_whole_cell':
        mask=np.asarray(read_elem(f['varm/native_feature_presence_by_sample']),dtype=bool)
        order=list(read_elem(f['uns/feature_presence_sample_order']))
        selected=sorted(set(obs.sample_id))
        assert set(selected)<=set(order) and mask.shape==(len(ids),len(order))
        native=mask[:,[order.index(s) for s in selected]].all(axis=1)
    ids=np.asarray([re.sub(r'\.\d+$','',x) if re.fullmatch(r'ENSG\d+(?:\.\d+)?',x) else '' for x in ids])
    return ids,symbols,native,status


def aggregate(group, obs, hhatl=False):
    assert group.attrs['encoding-type'] in ['csr_matrix', b'csr_matrix']
    n,g=map(int,group.attrs['shape'])
    assert n==len(obs)
    ids=pd.Index(sorted(set(obs.donor_id)));d=ids.get_indexer(obs.donor_id)
    assert (d>=0).all()
    use=obs.SubCluster.astype(str).eq('CM-HHATL').to_numpy() if hhatl else np.zeros(n,dtype=bool)
    ptr=np.asarray(group['indptr'],dtype=np.int64)
    assert len(ptr)==n+1 and ptr[0]==0 and ptr[-1]==len(group['data'])==len(group['indices']) and np.all(np.diff(ptr)>=0)
    sums=np.zeros((len(ids),g),dtype=np.int64);sub=np.zeros_like(sums) if hhatl else None
    totals=np.zeros(len(ids),dtype=np.int64);subtotal=np.zeros(len(ids),dtype=np.int64)
    for i in range(0,n,1024):
        j=min(n,i+1024);b,e=ptr[i],ptr[j]
        counts=group['data'][b:e]
        assert np.isfinite(counts).all() and (counts>=0).all() and np.equal(counts,np.floor(counts)).all()
        cols=group['indices'][b:e]
        assert (cols>=0).all() and (cols<g).all()
        x=sp.csr_matrix((counts.astype(np.int64),cols,ptr[i:j+1]-b),shape=(j-i,g))
        weights=sp.csr_matrix((np.ones(j-i,dtype=np.int64),(d[i:j],np.arange(j-i))),shape=(len(ids),j-i))
        sums+=(weights@x).toarray()
        # A second, independent row-total reduction checks every donor; avoids float32/lib-size rounding.
        rowtotal=np.asarray(x.sum(axis=1)).ravel()
        np.add.at(totals,d[i:j],rowtotal)
        if hhatl:
            keep=use[i:j]
            sub+=(weights[:,keep]@x[keep,:]).toarray()
            np.add.at(subtotal,d[i:j][keep],rowtotal[keep])
    assert np.array_equal(totals,sums.sum(axis=1))
    if hhatl: assert np.array_equal(subtotal,sub.sum(axis=1))
    return ids,sums,sub


def collapse(ids, symbols, valid, matrix):
    keep=(ids!='') & valid
    unq, inv=np.unique(ids[keep],return_inverse=True)
    agg=np.zeros((len(unq),len(matrix)),dtype=np.int64)
    np.add.at(agg,inv,matrix[:,keep].T)
    assert np.array_equal(agg.sum(axis=0),matrix[:,keep].sum(axis=1))
    src=pd.DataFrame({'gene_id':ids[keep],'gene_symbol':symbols[keep]})
    unique=src.drop_duplicates().groupby('gene_id',sort=True)['gene_symbol'].agg(';'.join)
    out=pd.DataFrame(agg,index=pd.Index(unq,name='gene_id'))
    return out,unique.reindex(unq),keep


def emit(out,name,counts,df,symbols):
    assert list(counts.columns)==list(df.donor_id) and counts.index.is_unique
    assert len(counts)>0 and set(df.disease)=={'Control','HCM'} and df.donor_id.is_unique
    df=df.copy()
    df['pseudobulk_library_size']=counts.sum(axis=0).to_numpy()
    assert df.pseudobulk_library_size.gt(0).all()
    counts.to_csv(out/(name+'.tsv.gz'),sep='\t',index=True,index_label='gene_id',compression={'method':'gzip','compresslevel':1,'mtime':0})
    txt(out,name+'_donors.tsv',df)
    txt(out,name+'_genes.tsv',pd.DataFrame({'gene_id':counts.index,'gene_symbol':symbols.reindex(counts.index).to_numpy()}))


def run(sources,out,discover):
    audits=[]
    groups={'SCP1303':['scp_raw','scp_cb']} if discover else {'liu2023':['liu'],'chin_GSE161921':['chin'],'chin_GSE174691':['chin'],'GSE138262_whole_cell':['whole']}
    chin=[]
    for ds,names in groups.items():
        path,expected=sources[ds]
        assert sha(path)==expected, f'{ds} source sha256'
        with h5py.File(path,'r') as f:
            obs=obs_frame(f,ds);df=donors(ds,obs)
            ids,symbols,native,status=get_axis(f,ds,obs)
            for name in names:
                donorids,joint,sub=aggregate(f[LAYERS[name]],obs,hhatl=name=='scp_raw')
                assert list(donorids)==list(df.donor_id)
                counts,gmap,keep=collapse(ids,symbols,native,joint)
                counts.columns=donorids
                excluded=joint[:,~keep].sum(axis=1)
                totals=joint.sum(axis=1)
                assert np.array_equal(counts.sum(axis=0).to_numpy()+excluded,totals)
                audits.append(dict(dataset=ds,modality=name,layer=LAYERS[name],source_sha256=expected,source_cells=len(obs),
                                   source_genes=len(ids),native_across_all_samples=int(native.sum()),unmapped_features=int((ids=='').sum()),
                                   included_features=int(keep.sum()),collapsed_ensembl_ids=len(counts),source_total=int(totals.sum()),
                                   excluded_count=int(excluded.sum()),included_count=int(counts.to_numpy().sum()),direct_donor_conservation=True))
                txt(out,ds+'_'+name+'_count_balance.tsv',pd.DataFrame({'donor_id':donorids,'source_total':totals,
                    'excluded_unmapped_or_unmeasured':excluded,'included_before_cross_accession_intersection':counts.sum(axis=0).to_numpy()}))
                if name=='liu' or name=='whole':
                    lookup=pd.DataFrame({'source_feature_id':list(read_elem(f['var']).index),'gene_symbol':symbols,'gene_id':ids,
                                         'measured_in_all_included_samples':native,'mapping_status':status,'included_in_model':keep})
                    lookup.to_csv(out/(name+'_feature_audit.tsv.gz'),sep='\t',index=False)
                if name=='chin':
                    chin.append((counts,df,gmap))
                else:
                    emit(out,name,counts,df,gmap)
                    if name=='scp_raw':
                        subframe=donors(ds,obs.loc[obs.SubCluster.astype(str).eq('CM-HHATL')])
                        retain=subframe.loc[subframe.nuclei.ge(50)].copy()
                        assert retain.groupby('disease').donor_id.size().min()>=5
                        sc,sg,_=collapse(ids,symbols,native,sub)
                        sc.columns=donorids
                        emit(out,'scp_hhatl',sc[retain.donor_id],retain,sg)
                        subset_totals=sub.sum(axis=1)
                        subset_excluded=sub[:,~keep].sum(axis=1)
                        assert np.array_equal(sc.sum(axis=0).to_numpy()+subset_excluded,subset_totals)
                        audits.append(dict(dataset=ds,modality='scp_hhatl',source_cells=int(retain.nuclei.sum()),included_donors=len(retain),
                                           layer=LAYERS[name],min_per_donor=50,source_sha256=expected,direct_donor_conservation=True))
                print('aggregated',name,ds,len(counts),len(donorids),flush=True)
    if chin:
        (a,da,ga),(b,db,gb)=chin
        both=a.index.intersection(b.index,sort=False)
        allcounts=pd.concat([a.loc[both],b.loc[both]],axis=1)
        md=pd.concat([da,db],ignore_index=True)
        assert allcounts.columns.is_unique and md.donor_id.is_unique
        gm=pd.concat([ga.loc[both],gb.loc[both]],axis=1)
        symbols=gm.iloc[:,0].where(gm.iloc[:,0].eq(gm.iloc[:,1]),gm.iloc[:,0]+';'+gm.iloc[:,1])
        emit(out,'chin',allcounts,md,symbols)
        audit=pd.concat([a.sum().rename('before_intersection'),b.sum().rename('before_intersection')]).to_frame().join(allcounts.sum().rename('after_intersection'))
        audit['not_shared_count']=audit.before_intersection-audit.after_intersection
        assert audit.not_shared_count.ge(0).all()
        txt(out,'chin_common_axis_count_balance.tsv',audit.rename_axis('donor_id').reset_index())
        audits.append(dict(dataset='Chin_pair',modality='chin_shared_axis',shared_genes=len(both),n_donors=len(md),removed_accession_specific_count=int(audit.not_shared_count.sum())))
    save_json(out,'DISCOVERY_PSEUDOBULK.json' if discover else 'SUPPORT_PSEUDOBULK.json',audits)


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('work',type=Path);p.add_argument('phase',choices=['roster','discovery','support'])
    args=p.parse_args()
    out=args.work/'pb';out.mkdir(parents=True,exist_ok=True)
    sources=inputs()
    if args.phase=='roster':roster(sources,out);return
    gate=json.loads((out/'ROSTER_GATE.json').read_text())
    assert gate['roster_sha256']==sha(out/'DONOR_ROSTER.tsv'), 'role roster changed after freeze'
    if args.phase=='support':
        lock=args.work/'lock/CANDIDATES.tsv'
        assert lock.is_file() and lock.stat().st_size>0,'freeze candidates before held-out outcomes'
    run(sources,out,args.phase=='discovery')

if __name__=='__main__':main()
