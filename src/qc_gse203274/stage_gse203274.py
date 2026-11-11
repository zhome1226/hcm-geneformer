"""Stage published pediatric cardiac nuclei with immutable counts and donor binding."""
from datetime import datetime, timezone
from pathlib import Path
import hashlib
import json
import re

import anndata as ad
import numpy as np
import pandas as pd
from scipy.io import mmread, mmwrite
from threadpoolctl import threadpool_limits

ROOT=Path.home()/'xiaoe-internship-data'
WORK=ROOT/'_processing/all_data_cleaning_20261111/GSE203274'
INPUT=WORK/'AllNuclei_snRNA_rawCount'
SOURCE=ROOT/'validation/pediatric_GSE203274'
REF=ROOT/'processed/SCP1303_review_20261111/identity'
GROUPS={'Donor':('Control',0),'HCM':('HCM',1),'DCM':('DCM',-1),
        'TOF':('TOF',-1),'Neo_HLHS':('HLHS',-1),'HF_HLHS':('HLHS',-1)}


def digest(path,algorithm='sha256'):
    with Path(path).open('rb') as handle:
        return hashlib.file_digest(handle,algorithm).hexdigest()


def age_parts(text):
    match=re.fullmatch(r'(\d+)y(?:[_/](\d+)m)?(?:[_/](\d+)d)?',str(text))
    if not match:return None
    return [int(p or 0) for p in match.groups()]


def main():
    assert not (WORK/'STAGE_SUCCESS.json').exists(),'Inspect existing successful stage before rerunning'
    baseline=json.loads((WORK/'ORIGINAL_DOWNLOAD_BASELINE.json').read_text(encoding='utf-8-sig'))
    assert len(baseline)==2
    for entry in baseline:
        path=Path(entry['server_path'])
        assert path.stat().st_size==int(entry['bytes']) and digest(path,'md5')==entry['md5'],path
    print('Original archive and metadata bytes/MD5 match baseline',flush=True)
    genes=pd.read_csv(INPUT/'genes.tsv.gz',sep='\t',header=None,keep_default_na=False)[0].astype(str)
    cells=pd.read_csv(INPUT/'barcodes.tsv.gz',sep='\t',header=None,keep_default_na=False)[0].astype(str)
    meta=pd.read_csv(SOURCE/'GSE203274_AllNuclei_snRNA_metadata.csv.gz',index_col=0,keep_default_na=False)
    assert cells.is_unique and meta.index.is_unique and set(cells)==set(meta.index)
    meta=meta.loc[cells].copy();meta.index.name='cell_id'
    print('Reading full published integer Matrix Market counts',flush=True)
    with threadpool_limits(limits=4):x=mmread(INPUT/'expression_matrix.mtx.gz').T.tocsr()
    x.sum_duplicates();x.eliminate_zeros();x.sort_indices()
    assert x.shape==(len(cells),len(genes))
    assert np.isfinite(x.data).all() and (x.data>=0).all() and np.equal(x.data,np.rint(x.data)).all()
    assert np.array_equal(meta['orig.ident'].astype(str),cells.str.rsplit('_',n=1).str[0])
    feature_ids=pd.Index(genes)
    if not feature_ids.is_unique:
        feature_ids=pd.Index([symbol+'__source_row'+str(i) for i,symbol in enumerate(genes)])
    var=pd.DataFrame({'gene_symbol':genes.to_numpy(),'source_row':np.arange(len(genes)),
                      'source_symbol_is_unique':~genes.duplicated(keep=False).to_numpy()},index=feature_ids)
    var.index.name='source_feature_id';var['mt']=var.gene_symbol.str.startswith('MT-')
    total=np.asarray(x.sum(axis=1)).ravel();detected=np.diff(x.indptr)
    mt_counts=np.asarray(x[:,var.mt.to_numpy()].sum(axis=1)).ravel()
    mt_percent=np.divide(mt_counts*100.,total,out=np.full(len(cells),np.nan),where=total>0)
    assert np.array_equal(total,pd.to_numeric(meta.nCount_RNA).to_numpy())
    assert np.array_equal(detected,pd.to_numeric(meta.nFeature_RNA).to_numpy())
    published_mt=pd.to_numeric(meta['percent.mt']).to_numpy()
    mt_units={'published_equals_percentage':bool(np.allclose(mt_percent,published_mt,rtol=0,atol=1e-8)),
              'published_equals_fraction':bool(np.allclose(mt_percent/100,published_mt,rtol=0,atol=1e-8)),
              'max_difference_percentage':float(np.max(np.abs(mt_percent-published_mt)))}
    assert mt_units['published_equals_percentage'] or mt_units['published_equals_fraction'],mt_units
    obs=meta.add_prefix('author_')
    obs['sample_id']=meta['orig.ident'].astype(str)
    obs['donor_original']=meta.labID.astype(str)
    # The author patientID field splits 13_198_LV/RV; labID and GEO bind one donor.
    obs['donor_id']='GSE203274_'+obs.donor_original
    obs['source_patientID']=meta.patientID.astype(str)
    obs['source_patientID_is_donor']=meta.patientID.astype(str).eq(meta.labID.astype(str))
    obs['clinical_group']=meta.ClinicalRank.astype(str)
    assert set(obs.clinical_group)==set(GROUPS)
    assert meta['Diagnosis'].astype(str).eq(meta.ClinicalRank.astype(str)).all()
    obs['disease']=obs.clinical_group.map({k:v[0] for k,v in GROUPS.items()})
    obs['label']=obs.clinical_group.map({k:v[1] for k,v in GROUPS.items()}).astype(np.int8)
    obs['binary_disease_eligible']=obs.label.ge(0)
    obs['tissue_region']=meta.region.astype(str)
    obs['technology']='10x v3 snRNA-seq, published preselected nuclei'
    obs['age_original']=meta.age.astype(str)
    ages=obs.age_original.map(age_parts)
    assert ages.notna().all()
    obs['age_year_component']=[p[0] for p in ages]
    obs['age_month_component']=[p[1] for p in ages]
    obs['age_day_component']=[p[2] for p in ages]
    obs['age_years_approx_from_author_text']=[p[0]+p[1]/12+p[2]/365.25 for p in ages]
    obs['total_counts']=total;obs['n_genes_by_counts']=detected
    obs['mt_counts']=mt_counts;obs['pct_counts_mt']=mt_percent
    obs['pass_strict_numeric']=(detected>=500)&(detected<=5000)&(mt_percent<5)
    obs['numeric_exclusion_reason']=np.select([total==0,detected<500,detected>5000,mt_percent>=5],
                                             ['zero_counts','genes_lt500','genes_gt5000','MT_ge5pct'],default='pass')
    obs['split_status']='unassigned_bind_donor_across_replicates_and_LV_RV'
    obs['age_source_conflict']=False
    geo=json.loads((WORK/'GEO_SAMPLES_GSE203274.json').read_text(encoding='utf-8-sig'))
    geo=[r for r in geo if r['source_name_ch1']==['Cardiac']]
    assert len(geo)==30
    per_donor={}
    for r in geo:
        name=r['title'][0].split(',')[0].strip()
        characteristics=dict(v.split(': ',1) for v in r['characteristics_ch1'])
        per_donor.setdefault(name,[]).append((r['accession'],characteristics))
    samples=[]
    for sample,indices in obs.groupby('sample_id',observed=True).groups.items():
        part=obs.loc[indices]
        for column in ['donor_original','clinical_group','disease','label','tissue_region','age_original']:
            assert part[column].nunique()==1,(sample,column)
        donor=part.donor_original.iloc[0]
        alias=(('WU198LV' if part.tissue_region.iloc[0]=='LV' else 'WU198RV') if donor=='13_198'
               else 'WU13235' if donor=='13_235' else donor)
        candidates=per_donor[alias]
        expected={'HCM':'HCM','DCM':'DCM','TOF':'TOF','Control':'Donor'}
        if part.disease.iloc[0] in expected:
            assert all(c[1]['diagnosis']==expected[part.disease.iloc[0]] for c in candidates)
        else:
            assert all('HLH' in c[1]['diagnosis'] for c in candidates)
        geo_ages=sorted({c[1]['age'] for c in candidates})
        assert len(geo_ages)==1
        geo_age=geo_ages[0];author_age=part.age_original.iloc[0]
        same_age=age_parts(geo_age)==age_parts(author_age)
        obs.loc[indices,'geo_age_original']=geo_age
        obs.loc[indices,'age_source_conflict']=not same_age
        obs.loc[indices,'geo_diagnosis_original']=' | '.join(sorted({c[1]['diagnosis'] for c in candidates}))
        obs.loc[indices,'geo_donor_alias']=alias
        # All candidate accessions remain explicit; don't guess 9k/10k to rep1/rep2.
        obs.loc[indices,'candidate_geo_accessions']=';'.join(c[0] for c in candidates)
        obs.loc[indices,'capture_GSM_mapping_status']='candidate_set_within_same_donor_and_region_not_unique'
        row={'sample_id':sample,'donor_id':part.donor_id.iloc[0],'donor_original':donor,
             'source_patientID':'|'.join(sorted(part.source_patientID.unique())),
             'clinical_group':part.clinical_group.iloc[0],'disease':part.disease.iloc[0],'label':int(part.label.iloc[0]),
             'tissue_region':part.tissue_region.iloc[0],'author_age':author_age,'geo_age':geo_age,'age_source_conflict':not same_age,
             'geo_donor_alias':alias,'candidate_geo_accessions':';'.join(c[0] for c in candidates),
             'capture_GSM_mapping_status':'candidate_set_within_same_donor_and_region_not_unique',
             'source_nuclei':len(part),'strict_nuclei':int(part.pass_strict_numeric.sum()),
             'author_CM':int(part.author_MainCellType.eq('CM').sum()),
             'author_CM_strict':int((part.author_MainCellType.eq('CM')&part.pass_strict_numeric).sum()),
             'author_median_MT_percent':float(part.pct_counts_mt.median())}
        samples.append(row)
    assert len(samples)==30 and obs.donor_id.nunique()==13
    assert obs.loc[obs.disease.eq('HCM'),'donor_id'].nunique()==1
    # Fill missing symbols only from the existing authoritative reference, never infer order.
    reference=pd.read_csv(REF/'reference_genes.tsv',sep='\t',keep_default_na=False)
    unique_ref=reference.loc[~reference.gene_symbol.duplicated(keep=False)&reference.gene_symbol.ne('')]
    unique_source=var.loc[var.source_symbol_is_unique]
    common=unique_ref.loc[unique_ref.gene_symbol.isin(unique_source.gene_symbol)].copy()
    assert len(common)>10000
    source_lookup=pd.Series(np.arange(len(var)),index=var.gene_symbol)
    positions=source_lookup.loc[common.gene_symbol].to_numpy(dtype=int)
    var['reference_gene_id']=''
    var.iloc[positions,var.columns.get_loc('reference_gene_id')]=common.gene_id.to_numpy()
    var['gene_id_mapping_status']=np.where(var.reference_gene_id.ne(''),'unique_exact_symbol_match_to_existing_GRCh38_reference','unmapped_or_ambiguous_symbol_preserved')
    for column in obs.select_dtypes(include=['object','string']).columns:
        obs[column]=obs[column].fillna('').astype(str)
    a=ad.AnnData(x,obs=obs,var=var)
    a.uns['source']='GSE203274 author released unnormalized integer RNA counts after CellBender v0.1.0 and author QC; rawCount is the source filename label, not uncorrected CellRanger counts or full droplets'
    a.uns['count_source_assembly']='GEO metadata reports hg19; original genes are symbols, not native Ensembl IDs'
    a.uns['reference_mapping_policy']='Unique exact symbol matches to existing adult LV reference used only for identity; all source genes retained, mapped GRCh38 IDs are derived annotations not claimed native hg19 IDs'
    a.uns['author_MT_unit_check_json']=json.dumps(mt_units)
    a.write_h5ad(WORK/'author_full_counts.h5ad',compression='gzip')
    obs.to_csv(WORK/'cell_numeric_audit.tsv.gz',sep='\t')
    pd.DataFrame(samples).to_csv(WORK/'samples.tsv',sep='\t',index=False)
    common.assign(source_column=positions).to_csv(WORK/'reference_gene_mapping.tsv',sep='\t',index=False)
    common.gene_id.to_csv(WORK/'reference_common_genes.txt',index=False,header=False)
    strict=obs.pass_strict_numeric.to_numpy()
    with threadpool_limits(limits=4):mmwrite(WORK/'strict_reference_query.mtx',x[strict][:,positions].T)
    pd.Series(obs.index[strict]).to_csv(WORK/'strict_query_cells.txt',index=False,header=False)
    fields=[{'field':c,'unique_values':int(meta[c].nunique()),'examples':[str(v) for v in meta[c].unique()[:10]]} for c in meta]
    result={'checked_utc':datetime.now(timezone.utc).isoformat(),'source_nuclei':a.n_obs,'genes':a.n_vars,'nnz':int(x.nnz),
            'samples':len(samples),'donors':13,'HCM_donors':1,'Control_donors':4,'strict_nuclei':int(strict.sum()),
            'common_reference_genes':len(common),'mt_genes':int(var.mt.sum()),'mt_unit_check':mt_units,
            'age_conflicting_donors':sorted(obs.loc[obs.age_source_conflict.astype(bool),'donor_original'].unique()),
            'source_hashes':{str(p):{'bytes':p.stat().st_size,'sha256':digest(p)} for p in INPUT.iterdir() if p.is_file()},
            'metadata_fields':fields}
    (WORK/'STAGE_SUCCESS.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps({k:v for k,v in result.items() if k not in ['source_hashes','metadata_fields']},indent=2),flush=True)


if __name__=='__main__':main()
