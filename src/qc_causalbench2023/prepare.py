"""Curated Perturb-seq releases: audit all raw counts without benchmark selection."""
from datetime import datetime,timezone
from pathlib import Path
import hashlib
import json
import shutil
import anndata as ad
import h5py
import numpy as np
import pandas as pd

ROOT=Path.home()/'xiaoe-internship-data'
WORK=ROOT/'_processing/all_data_cleaning_20261111/CausalBench2023'
OUT=ROOT/'processed/CausalBench2023_qc_20261111'


def digest(path,kind='sha256'):
    with path.open('rb') as f:
        return hashlib.file_digest(f,kind).hexdigest()


def dump(name,value):
    (OUT/name).write_text(json.dumps(value,indent=2,allow_nan=False)+'\n')


def main():
    OUT.mkdir(parents=True,exist_ok=True)
    assert not (OUT/'EXPORT_SUCCESS.json').exists(),'Successful preparation exists; do not repeat'
    baseline=json.loads((WORK/'ORIGINAL_DOWNLOAD_BASELINE.json').read_text())
    assert len(baseline)==2
    summary,originals,types,matrices=[],[],{},[]
    for item in baseline:
        path=Path(item['server_path'])
        cell_line='K562' if path.name.startswith('K562') else 'RPE1'
        md5,sha=digest(path,'md5'),digest(path)
        assert path.stat().st_size==int(item['bytes']) and md5==item['md5'] and sha==item['sha256']
        originals.append({'file':path.name,'path':str(path),'bytes':path.stat().st_size,'MD5':md5,'SHA256':sha,'source_url':item['source_url']})
        with h5py.File(path,'r') as f:
            obs,var=ad.io.read_elem(f['obs']),ad.io.read_elem(f['var'])
            x=f['X']
            assert isinstance(x,h5py.Dataset) and x.shape==(len(obs),len(var))
            assert obs.index.is_unique and var.index.is_unique
            assert (obs.index.astype(str).str.strip()!='').all() and var.index.astype(str).str.startswith('ENSG').all()
            assert set(f)=={'X','obs','var'}
            types[cell_line]={'obs':{k:str(v) for k,v in obs.dtypes.items()},'var':{k:str(v) for k,v in var.dtypes.items()}}
            n,p=x.shape
            sums,n_gene=np.zeros(n),np.zeros(n,dtype='int64')
            feature_sum,feature_n=np.zeros(p),np.zeros(p,dtype='int64')
            mt_mask=var.gene_name.astype(str).str.startswith('MT-').to_numpy()
            mt_sum=np.zeros(n)
            minimum,maximum=np.inf,-np.inf
            nonzero=0
            for start in range(0,n,768):
                stop=min(start+768,n)
                a=x[start:stop,:].astype(np.float64)
                assert np.isfinite(a).all() and (a>=0).all() and (a==np.rint(a)).all() and (a<np.iinfo(np.int64).max).all()
                sums[start:stop]=a.sum(axis=1)
                hit=a>0
                n_gene[start:stop]=hit.sum(axis=1)
                feature_sum+=a.sum(axis=0)
                feature_n+=hit.sum(axis=0)
                mt_sum[start:stop]=a[:,mt_mask].sum(axis=1)
                minimum,maximum=min(minimum,float(a.min())),max(maximum,float(a.max()))
                nonzero+=int(hit.sum())
            qc=obs.copy()
            qc.index.name='source_cell_id'
            qc['source_row']=np.arange(n)
            qc['source_file']=str(path)
            qc['counts_slot']='X'
            qc['cell_line']=cell_line
            qc['cell_uid']='CausalBench2023:'+cell_line+':'+obs.index.astype(str)
            qc['capture_uid']='CausalBench2023:'+cell_line+':gem:'+obs.gem_group.astype(str).to_numpy()
            qc['perturbation_gene_id']=obs.gene_id.astype(str).to_numpy()
            control=obs.gene.astype(str).eq('non-targeting').to_numpy()
            assert np.array_equal(control,obs.gene_id.astype(str).eq('non-targeting'))
            qc['is_non_targeting_control']=control
            qc['target_in_published_gene_axis']=obs.gene_id.astype(str).isin(var.index.astype(str)).to_numpy()
            qc['target_expression_missing_from_axis']=~control & ~qc.target_in_published_gene_axis
            qc['published_count_sum']=sums
            qc['published_detected_genes']=n_gene
            qc['published_mt_count_sum']=mt_sum
            qc['published_axis_mt_pct']=np.divide(100*mt_sum,sums,out=np.full(n,np.nan),where=sums>0) if mt_mask.any() else np.full(n,np.nan)
            qc['mt_percentage_scope']='published_feature_axis_only' if mt_mask.any() else 'no_mitochondrial_genes_in_published_axis'
            qc['author_mitopercent_unit']='unverified_preserved_without_conversion'
            required=['gene','gene_id','gene_transcript','sgID_AB','gem_group']
            valid=obs[required].notna().all(axis=1).to_numpy() & obs[required].astype(str).apply(lambda col:col.str.strip().ne('')).all(axis=1).to_numpy()
            qc['required_metadata_valid']=valid
            valid&=sums>0
            qc['usable_counts']=valid
            qc['quarantine_reason']=np.where(~qc.required_metadata_valid,'missing_identity_metadata',np.where(sums<=0,'zero_count_cell',''))
            qc['official_split']='not_embedded_in_raw_release_reference_code_defines_parameterized_split'
            qc['donor_identity']='not_provided_cell_line_data'
            qc['HCM_label_assigned']=False
            qc.to_csv(OUT/f'{cell_line}_observation_QC.tsv.gz',sep='\t')
            qc.loc[valid,['source_row','source_file','counts_slot','cell_uid','capture_uid','perturbation_gene_id','is_non_targeting_control','target_in_published_gene_axis','official_split']].to_csv(OUT/f'{cell_line}_usable_rows.tsv.gz',sep='\t')
            qc.loc[~valid].to_csv(OUT/f'{cell_line}_quarantine.tsv.gz',sep='\t')
            feature=var.copy()
            feature.index.name='source_gene_id'
            feature['source_column']=np.arange(p)
            feature['feature_uid']='CausalBench2023:'+cell_line+':'+var.index.astype(str)
            feature['duplicate_symbol']=var.gene_name.duplicated(keep=False).to_numpy()
            feature['published_count_sum']=feature_sum
            feature['cells_with_nonzero_counts']=feature_n
            feature['usable_feature']=feature_n>0
            feature.to_csv(OUT/f'{cell_line}_feature_QC.tsv.gz',sep='\t')
            conditions=qc.groupby(['gene','gene_id','gene_transcript','sgID_AB','is_non_targeting_control','target_in_published_gene_axis'],observed=True,dropna=False).agg(cells=('usable_counts','size'),usable_cells=('usable_counts','sum')).reset_index()
            conditions.to_csv(OUT/f'{cell_line}_PERTURBATION_SUMMARY.tsv',sep='\t',index=False)
            batches=qc.groupby(['gem_group','capture_uid'],observed=True).agg(cells=('usable_counts','size'),usable_cells=('usable_counts','sum'),control_cells=('is_non_targeting_control','sum')).reset_index()
            batches.to_csv(OUT/f'{cell_line}_CAPTURE_SUMMARY.tsv',sep='\t',index=False)
            matrix={'cell_line':cell_line,'rows':n,'features':p,'minimum':minimum,'maximum':maximum,'nonzero_values':nonzero,'nonfinite_values':0,'negative_values':0,'fractional_values':0,'zero_count_cells':int((sums==0).sum()),'all_zero_features':int((feature_n==0).sum()),'count_sum_total':float(sums.sum())}
            matrices.append(matrix)
            summary.append({'cell_line':cell_line,'source_cells':n,'source_genes':p,'usable_cells':int(valid.sum()),'quarantined_cells':int((~valid).sum()),'control_cells':int(control.sum()),'perturbation_labels_including_control':int(obs.gene_id.nunique()),'captures':int(obs.gem_group.nunique()),'target_missing_axis_cells':int(qc.target_expression_missing_from_axis.sum()),'duplicate_symbol_rows':int(feature.duplicate_symbol.sum()),'published_mt_genes':int(mt_mask.sum()),'source_file':str(path),'counts_slot':'X'})
            print('Full counts and metadata prepared:',cell_line,n,p,int(valid.sum()),flush=True)
    pd.DataFrame(summary).to_csv(OUT/'DATASET_QC_SUMMARY.tsv',sep='\t',index=False)
    pd.DataFrame(matrices).to_csv(OUT/'MATRIX_QC_SUMMARY.tsv',sep='\t',index=False)
    pd.DataFrame(originals).to_csv(OUT/'ORIGINAL_FILES.tsv',sep='\t',index=False)
    dump('SOURCE_METADATA_DTYPES.json',types)
    policy={'dataset':'CausalBench2023','source_data':'Replogle Perturb-seq; exact two Figshare files used by pinned CausalBench code',
        'counts_location':'X; full finite/nonnegative/integer scan plus official preprocessing code before normalization',
        'originals_read_only':True,'primary_delivery':'verified original H5AD plus exact usable row and feature indices',
        'count_quality':'nonzero cell counts and complete source perturbation/capture identity; no HCM numerical filters',
        'control_definition':'exact gene and gene_id == non-targeting; NANOG remains a targeted gene',
        'guides':'sgID_AB retained verbatim; two guides for same target are part of CRISPRi design, not doublet evidence',
        'targets_outside_gene_axis':'preserved as known interventions; explicit flag, not guessed into a gene column or deleted',
        'duplicate_gene_symbols':'retain distinct Ensembl IDs; no aggregation',
        'author_mitopercent':'unit not independently established from downloaded metadata; preserve field unchanged, no assumed fraction-to-percent conversion',
        'MT_computation':'only available published gene axis, no filtering',
        'split':'no fixed row split embedded in the two raw files; official reference code applies stratified test_size=0.2/random_state=0 after task-specific preprocessing; preserve recipe, create no new training split',
        'strong_perturbation_filter':'optional benchmark filtering needs separate summary XLSX and cell knockdown selection; not run during count cleaning',
        'feature_selection':'official graph-inference selection of targets with >100 cells not applied to general-purpose count archive',
        'HCM_rules_applied':False,'HCM_labels_assigned':False,'normalization':False,'imputation':False,'doublet_recall':False,'network_inference':False,'training':False}
    dump('qc_policy.json',policy)
    evidence={'official_code_commit':'1a2143cffdc85f835b41ce8d52034be1bf903e71','repository':'https://github.com/causalbench/causalbench',
        'download_file_IDs':[35773219,35775606],'official_split_parameters':{'test_size':0.2,'random_state':0,'stratify':'post_preprocessing_interventions'},
        'actual_fixed_row_split_in_downloaded_files':False,'new_split_generated':False,'source_code':[]}
    lookup=OUT/'source_lookup'
    lookup.mkdir(exist_ok=True)
    for source in sorted((WORK/'source_lookup').glob('*.py')):
        shutil.copyfile(source,lookup/source.name)
        evidence['source_code'].append({'file':source.name,'SHA256':digest(source)})
    dump('SOURCE_EVIDENCE.json',evidence)
    (OUT/'REPORT.md').write_text('# CausalBench2023：原始计数与扰动索引\n\n来源为Replogle Perturb-seq的K562/RPE1细胞系，官方CausalBench固定代码指向同两个Figshare文件；不是HCM队列。X在官方归一化步骤前读取，本次完整核查数值域并保留全部发布基因。\n\n'+'\n'.join(f"- {s['cell_line']}：{s['source_cells']}细胞、{s['source_genes']}基因；可用{s['usable_cells']}，隔离{s['quarantined_cells']}，non-targeting对照{s['control_cells']}，{s['captures']}捕获。" for s in summary)+'\n\n计数原件只读，以真实行列索引交付。扰动基因Ensembl ID、gene/transcript、sgID_AB和gem_group完整保留；双sgRNA不作双细胞，NANOG不作对照。Ensembl列唯一，重复symbol不合并；缺少目标基因表达列的细胞单独标记，但不伪造列或删除其有效表达。\n\n不套HCM阈值，不做归一化/插补/整合或网络推断；原mitopercent字段单位未独立核实而不擅自换算。当前对象无固定split，原官方参数化分层划分及可选strong-perturbation过滤配方存档，本次不新建训练划分或按扰动强度筛选。所有逐细胞/逐基因文件只留服务器。\n\n准备后须独立验收并核查退出码、进程和源码一致性。\n')
    records=[]
    for p in sorted(OUT.rglob('*')):
        if p.is_file():
            records.append({'file':str(p.relative_to(OUT)),'server_path':str(p),'bytes':p.stat().st_size,'MD5':digest(p,'md5'),'SHA256':digest(p)})
    pd.DataFrame(records).to_csv(OUT/'OUTPUT_MANIFEST.tsv',sep='\t',index=False)
    dump('EXPORT_SUCCESS.json',{'status':'prepared_not_independently_verified','completed_utc':datetime.now(timezone.utc).isoformat(),'units':summary})
    print('Preparation complete; independent verification required.',flush=True)


if __name__=='__main__':
    main()
