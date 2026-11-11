"""Independently validate the two original dense count matrices and all indices."""
from datetime import datetime,timezone
from pathlib import Path
import hashlib
import json
import tarfile
import traceback
import anndata as ad
import h5py
import numpy as np
import pandas as pd

ROOT=Path.home()/'xiaoe-internship-data'
WORK=ROOT/'_processing/all_data_cleaning_20261111/CausalBench2023'
OUT=ROOT/'processed/CausalBench2023_qc_20261111'
CHECKS=[]


def check(name,result):
    CHECKS.append({'name':name,'passed':bool(result)})
    if not result:
        raise AssertionError(name)


def digest(path,kind='sha256'):
    h=hashlib.new(kind)
    with path.open('rb') as f:
        while block:=f.read(16*1024*1024):
            h.update(block)
    return h.hexdigest()


def table(name,index=False):
    return pd.read_csv(OUT/name,sep='\t',index_col=0 if index else None,keep_default_na=False,float_precision='round_trip')


def source_metadata_equal(source,saved):
    if source.index.astype(str).tolist()!=saved.index.astype(str).tolist():
        return False
    for name in source:
        if name not in saved:
            return False
        a=source[name]
        if pd.api.types.is_numeric_dtype(a) and not pd.api.types.is_bool_dtype(a):
            b=pd.to_numeric(saved[name],errors='coerce').to_numpy().astype(a.to_numpy().dtype)
            if not np.array_equal(a.to_numpy(),b,equal_nan=True):
                return False
        else:
            expected=a.astype(object).where(a.notna(),'').astype(str)
            if not np.array_equal(expected,saved[name].astype(str)):
                return False
    return True


def same_summary_rows(actual,expected):
    return actual.columns.tolist()==expected.columns.tolist() and sorted(map(tuple,actual.astype(str).to_numpy()))==sorted(map(tuple,expected.astype(str).to_numpy()))


def main():
    assert not (OUT/'FINAL_COMPLETION_GATE.json').exists(),'Successful verification exists; do not rerun'
    baseline=json.loads((WORK/'ORIGINAL_DOWNLOAD_BASELINE.json').read_text())
    originals=table('ORIGINAL_FILES.tsv').set_index('path')
    units=table('DATASET_QC_SUMMARY.tsv').set_index('cell_line')
    matrices=table('MATRIX_QC_SUMMARY.tsv').set_index('cell_line')
    types=json.loads((OUT/'SOURCE_METADATA_DTYPES.json').read_text())
    check('exact_two_original_inputs',len(baseline)==len(originals)==len(units)==len(matrices)==2)
    all_uids=[]
    total_cells=0
    total_usable=0
    for item in baseline:
        path=Path(item['server_path'])
        line='K562' if path.name.startswith('K562') else 'RPE1'
        original=originals.loc[str(path)]
        check(line+'_source_bytes',path.stat().st_size==int(item['bytes'])==original.bytes)
        check(line+'_source_MD5',digest(path,'md5')==item['md5']==original.MD5)
        check(line+'_source_SHA256',digest(path)==item['sha256']==original.SHA256)
        with h5py.File(path,'r') as f:
            obs,var=ad.io.read_elem(f['obs']),ad.io.read_elem(f['var'])
            qc,feature=table(line+'_observation_QC.tsv.gz',True),table(line+'_feature_QC.tsv.gz',True)
            usable=table(line+'_usable_rows.tsv.gz',True)
            quarantine=table(line+'_quarantine.tsv.gz',True)
            check(line+'_native_only_raw_matrix_and_axes',set(f)=={'X','obs','var'} and isinstance(f['X'],h5py.Dataset))
            check(line+'_all_author_cell_metadata_exact',source_metadata_equal(obs,qc))
            check(line+'_all_author_gene_metadata_exact',source_metadata_equal(var,feature))
            check(line+'_source_metadata_dtypes',types[line]=={'obs':{k:str(v) for k,v in obs.dtypes.items()},'var':{k:str(v) for k,v in var.dtypes.items()}})
            n,p=f['X'].shape
            check(line+'_source_shape_and_unique_axes',(n,p)==(len(obs),len(var)) and obs.index.is_unique and var.index.is_unique)
            check(line+'_source_row_and_column_indices',np.array_equal(qc.source_row,np.arange(n)) and np.array_equal(feature.source_column,np.arange(p)))
            check(line+'_source_file_and_count_pointers',qc.source_file.eq(str(path)).all() and qc.counts_slot.eq('X').all())
            check(line+'_unique_cell_and_capture_namespaces',np.array_equal(qc.cell_uid,'CausalBench2023:'+line+':'+obs.index.astype(str)) and np.array_equal(qc.capture_uid,'CausalBench2023:'+line+':gem:'+obs.gem_group.astype(str)) and qc.cell_uid.is_unique)
            check(line+'_feature_Ensembl_axis_unchanged',var.index.astype(str).str.startswith('ENSG').all() and np.array_equal(feature.feature_uid,'CausalBench2023:'+line+':'+var.index.astype(str)))
            check(line+'_duplicate_symbols_retained',np.array_equal(feature.duplicate_symbol,var.gene_name.duplicated(keep=False)))
            mt=var.gene_name.astype(str).str.startswith('MT-').to_numpy()
            row_sum,row_n,row_mt=np.zeros(n),np.zeros(n,dtype='int64'),np.zeros(n)
            feature_sum,feature_n=np.zeros(p),np.zeros(p,dtype='int64')
            lo,hi=np.inf,-np.inf
            nonzero=0
            for start in range(0,n,511):
                stop=min(n,start+511)
                # HDF5 read_direct uses a separate I/O path and different block sizes.
                values=np.empty((stop-start,p),dtype=f['X'].dtype)
                f['X'].read_direct(values,np.s_[start:stop,:])
                check(line+'_finite_nonnegative_integers_rows_'+str(start),np.isfinite(values).all() and (values>=0).all() and np.equal(values,np.floor(values)).all())
                a=values.astype(np.float64)
                row_sum[start:stop]=np.add.reduce(a,axis=1)
                row_n[start:stop]=np.count_nonzero(a,axis=1)
                row_mt[start:stop]=np.add.reduce(a[:,mt],axis=1)
                feature_sum+=np.add.reduce(a,axis=0)
                feature_n+=np.count_nonzero(a,axis=0)
                nonzero+=int(np.count_nonzero(a))
                lo,hi=min(lo,float(a.min())),max(hi,float(a.max()))
            check(line+'_all_cell_count_sums',np.array_equal(qc.published_count_sum,row_sum))
            check(line+'_all_cell_detected_genes',np.array_equal(qc.published_detected_genes,row_n))
            check(line+'_all_cell_MT_count_sums',np.array_equal(qc.published_mt_count_sum,row_mt))
            expected_mt=np.divide(100*row_mt,row_sum,out=np.full(n,np.nan),where=row_sum>0) if mt.any() else np.full(n,np.nan)
            check(line+'_available_gene_axis_MT_percentage',np.array_equal(pd.to_numeric(qc.published_axis_mt_pct,errors='coerce'),expected_mt,equal_nan=True))
            check(line+'_author_MT_not_silently_rescaled',qc.author_mitopercent_unit.eq('unverified_preserved_without_conversion').all())
            check(line+'_all_gene_count_sums',np.array_equal(feature.published_count_sum,feature_sum))
            check(line+'_all_gene_detection_counts',np.array_equal(feature.cells_with_nonzero_counts,feature_n))
            check(line+'_usable_gene_rule',np.array_equal(feature.usable_feature,feature_n>0))
            control=obs.gene.astype(str).eq('non-targeting').to_numpy()
            check(line+'_control_label_consistency',np.array_equal(control,obs.gene_id.astype(str).eq('non-targeting')))
            check(line+'_exact_non_targeting_controls',np.array_equal(qc.is_non_targeting_control,control))
            check(line+'_NANOG_not_misclassified_as_control',not qc.loc[obs.gene.astype(str).eq('NANOG'),'is_non_targeting_control'].any())
            target=obs.gene_id.astype(str).isin(var.index.astype(str)).to_numpy()
            check(line+'_target_gene_IDs_preserved',np.array_equal(qc.perturbation_gene_id,obs.gene_id.astype(str)))
            check(line+'_target_axis_flags',np.array_equal(qc.target_in_published_gene_axis,target) and np.array_equal(qc.target_expression_missing_from_axis,~target&~control))
            required=['gene','gene_id','gene_transcript','sgID_AB','gem_group']
            meta=obs[required].notna().all(axis=1).to_numpy()&obs[required].astype(str).apply(lambda col:col.str.strip().ne('')).all(axis=1).to_numpy()
            valid=meta&(row_sum>0)
            check(line+'_required_metadata_rule',np.array_equal(qc.required_metadata_valid,meta))
            check(line+'_usable_cell_rule',np.array_equal(qc.usable_counts,valid))
            check(line+'_usable_rows_exact',usable.index.tolist()==obs.index[valid].tolist() and np.array_equal(usable.source_row,np.flatnonzero(valid)))
            check(line+'_quarantine_rows_exact',quarantine.index.tolist()==obs.index[~valid].tolist() and np.array_equal(quarantine.source_row,np.flatnonzero(~valid)))
            check(line+'_quarantine_reasons',np.array_equal(qc.quarantine_reason,np.where(~meta,'missing_identity_metadata',np.where(row_sum<=0,'zero_count_cell',''))))
            check(line+'_no_HCM_labels_or_invented_donors',not qc.HCM_label_assigned.any() and qc.donor_identity.eq('not_provided_cell_line_data').all() and qc.cell_line.eq(line).all())
            check(line+'_missing_official_split_not_invented',qc.official_split.eq('not_embedded_in_raw_release_reference_code_defines_parameterized_split').all())
            s=units.loc[line]
            expected={'source_cells':n,'source_genes':p,'usable_cells':int(valid.sum()),'quarantined_cells':int((~valid).sum()),'control_cells':int(control.sum()),'perturbation_labels_including_control':int(obs.gene_id.nunique()),'captures':int(obs.gem_group.nunique()),'target_missing_axis_cells':int((~control&~target).sum()),'duplicate_symbol_rows':int(var.gene_name.duplicated(keep=False).sum()),'published_mt_genes':int(mt.sum())}
            for key,val in expected.items():
                check(line+'_dataset_summary_'+key,s[key]==val)
            expected_matrix={'rows':n,'features':p,'minimum':lo,'maximum':hi,'nonzero_values':nonzero,'nonfinite_values':0,'negative_values':0,'fractional_values':0,'zero_count_cells':int((row_sum==0).sum()),'all_zero_features':int((feature_n==0).sum()),'count_sum_total':float(row_sum.sum())}
            for key,val in expected_matrix.items():
                check(line+'_matrix_summary_'+key,matrices.loc[line,key]==val)
            expected_cond=qc.groupby(['gene','gene_id','gene_transcript','sgID_AB','is_non_targeting_control','target_in_published_gene_axis'],observed=True,dropna=False).agg(cells=('usable_counts','size'),usable_cells=('usable_counts','sum')).reset_index()
            actual_cond=table(line+'_PERTURBATION_SUMMARY.tsv')
            check(line+'_every_perturbation_and_guide_pair_count',same_summary_rows(actual_cond,expected_cond))
            expected_cap=qc.groupby(['gem_group','capture_uid'],observed=True).agg(cells=('usable_counts','size'),usable_cells=('usable_counts','sum'),control_cells=('is_non_targeting_control','sum')).reset_index()
            actual_cap=table(line+'_CAPTURE_SUMMARY.tsv')
            check(line+'_every_capture_count',same_summary_rows(actual_cap,expected_cap))
            all_uids.extend(qc.cell_uid.tolist())
            total_cells+=n
            total_usable+=int(valid.sum())
            print('Full independent matrix and index checks passed:',line,n,p,flush=True)
    check('global_cell_UIDs_unique_across_cell_lines',len(all_uids)==len(set(all_uids)))
    evidence=json.loads((OUT/'SOURCE_EVIDENCE.json').read_text())
    check('official_source_and_split_recipe',evidence['official_code_commit']=='1a2143cffdc85f835b41ce8d52034be1bf903e71' and evidence['download_file_IDs']==[35773219,35775606] and evidence['official_split_parameters']=={'test_size':0.2,'random_state':0,'stratify':'post_preprocessing_interventions'} and not evidence['actual_fixed_row_split_in_downloaded_files'] and not evidence['new_split_generated'])
    archive=ROOT/'benchmarks/competitions/code/causalbench__causalbench/1a2143cffdc85f835b41ce8d52034be1bf903e71.tar.gz'
    with tarfile.open(archive,'r:gz') as tar:
        official_paths={'create_dataset.py':'causalscbench/data_access/create_dataset.py','preprocessing.py':'causalscbench/data_access/utils/preprocessing.py','splitting.py':'causalscbench/data_access/utils/splitting.py','download_weissmann.py':'causalscbench/data_access/datasets/download_weissmann.py'}
        for source in evidence['source_code']:
            members=[m for m in tar.getmembers() if m.name.endswith('/'+official_paths[source['file']]) and m.isfile()]
            check('one_pinned_official_source_'+source['file'],len(members)==1)
            expected_sha=hashlib.sha256(tar.extractfile(members[0]).read()).hexdigest()
            check('official_source_SHA_'+source['file'],digest(OUT/'source_lookup'/source['file'])==expected_sha==source['SHA256'])
    policy=json.loads((OUT/'qc_policy.json').read_text())
    check('no_unrequested_HCM_filter_or_model_operations',all(policy[k] is False for k in ['HCM_rules_applied','HCM_labels_assigned','normalization','imputation','doublet_recall','network_inference','training']) and policy['originals_read_only'])
    for row in table('OUTPUT_MANIFEST.tsv').itertuples(index=False):
        path=Path(row.server_path)
        check('output_size_'+row.file,path.stat().st_size==row.bytes)
        check('output_MD5_'+row.file,digest(path,'md5')==row.MD5)
        check('output_SHA256_'+row.file,digest(path)==row.SHA256)
    now=datetime.now(timezone.utc).isoformat()
    verification={'status':'passed','checked_utc':now,'checks_passed':len(CHECKS),'matrices_fully_recomputed':2,'source_cells':total_cells,'usable_cells':total_usable,'checks':CHECKS}
    (OUT/'INDEPENDENT_VERIFICATION.json').write_text(json.dumps(verification,indent=2)+'\n')
    gate={'status':'passed','checked_utc':now,'independent_checks':len(CHECKS),'verification_SHA256':digest(OUT/'INDEPENDENT_VERIFICATION.json'),'prepared_manifest_SHA256':digest(OUT/'OUTPUT_MANIFEST.tsv'),'source_immutable':True,'scope':'both downloaded K562/RPE1 raw count releases and complete usable/perturbation/capture indices','does_not_certify_HCM_applicability_or_fixed_training_split':True}
    (OUT/'FINAL_COMPLETION_GATE.json').write_text(json.dumps(gate,indent=2)+'\n')
    records=[]
    for path in sorted(OUT.rglob('*')):
        if path.is_file() and path.name!='VERIFIED_OUTPUT_MANIFEST.tsv':
            records.append({'file':str(path.relative_to(OUT)),'server_path':str(path),'bytes':path.stat().st_size,'MD5':digest(path,'md5'),'SHA256':digest(path)})
    pd.DataFrame(records).to_csv(OUT/'VERIFIED_OUTPUT_MANIFEST.tsv',sep='\t',index=False)
    print('All',len(CHECKS),'independent checks passed.',flush=True)


if __name__=='__main__':
    try:
        main()
    except Exception as error:
        OUT.mkdir(parents=True,exist_ok=True)
        (OUT/'VERIFICATION_FAILURE.json').write_text(json.dumps({'error':repr(error),'traceback':traceback.format_exc(),'checks':CHECKS},indent=2)+'\n')
        raise
