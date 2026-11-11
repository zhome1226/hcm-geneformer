"""Re-read every source matrix independently and validate all published indices."""
from datetime import datetime, timezone
from pathlib import Path
import hashlib
import json
import traceback
import anndata as ad
import h5py
import numpy as np
import pandas as pd
from scipy import sparse

ROOT=Path.home()/'xiaoe-internship-data'
WORK=ROOT/'_processing/all_data_cleaning_20261111/NeurIPS2021'
OUT=ROOT/'processed/NeurIPS2021_qc_20261111'
CHECKS=[]


def check(name,passed):
    CHECKS.append({'name':name,'passed':bool(passed)})
    if not passed:
        raise AssertionError(name)


def digest(path,algorithm='sha256'):
    hasher=hashlib.new(algorithm)
    with path.open('rb') as f:
        while block:=f.read(8*1024*1024):
            hasher.update(block)
    return hasher.hexdigest()


def table(name,index=False):
    return pd.read_csv(OUT/name,sep='\t',index_col=0 if index else None,keep_default_na=False,float_precision='round_trip')


def metadata_equal(native,saved):
    if native.index.astype(str).tolist()!=saved.index.astype(str).tolist():
        return False
    for name in native:
        if name not in saved:
            return False
        source=native[name]
        if pd.api.types.is_numeric_dtype(source) and not pd.api.types.is_bool_dtype(source):
            restored=pd.to_numeric(saved[name],errors='coerce').to_numpy().astype(source.to_numpy().dtype)
            if not np.array_equal(restored,source.to_numpy(),equal_nan=True):
                return False
        else:
            expected=source.astype(object).where(source.notna(),'').astype(str).to_numpy()
            if not np.array_equal(saved[name].astype(str).to_numpy(),expected):
                return False
    return True


def independent_scan(node,role,mods=None,mt=None):
    enc=node.attrs.get('encoding-type','array')
    is_sparse=enc in {'csr_matrix','csc_matrix'}
    if is_sparse:
        shape=tuple(node.attrs['shape'])
        axis=0 if enc=='csr_matrix' else 1
        pointer=node['indptr'][:]
        check(node.name+'_sparse_pointer',len(pointer)==shape[axis]+1 and pointer[0]==0 and pointer[-1]==len(node['data'])==len(node['indices']) and (np.diff(pointer)>=0).all())
        for left in range(0,len(node['indices']),4_000_000):
            indices=node['indices'][left:left+4_000_000]
            check(node.name+'_index_bounds_'+str(left),not len(indices) or (indices.min()>=0 and indices.max()<shape[1-axis]))
        reader=ad.io.sparse_dataset(node)
    elif enc=='dataframe':
        reader=ad.io.read_elem(node)
        shape=reader.shape
        axis=0
    else:
        reader=node
        shape=node.shape
        axis=0
    n,p=shape
    rs,rn,rb=np.zeros(n),np.zeros(n,dtype='int64'),np.zeros(n,dtype='int64')
    cs,cn,cb=np.zeros(p),np.zeros(p,dtype='int64'),np.zeros(p,dtype='int64')
    groups={m:[np.zeros(n),np.zeros(n,dtype='int64')] for m in sorted(set(mods))} if mods is not None else {}
    rmt=np.zeros(n)
    low,high=np.inf,-np.inf
    negative,fractional=0,0
    for start in range(0,shape[axis],1024 if axis==0 else 128):
        stop=min(shape[axis],start+(1024 if axis==0 else 128))
        if is_sparse:
            part=reader[start:stop,:] if axis==0 else reader[:,start:stop]
            part=part.astype(np.float64)
        else:
            values=reader.iloc[start:stop].to_numpy() if isinstance(reader,pd.DataFrame) else reader[start:stop,:]
            part=sparse.csr_matrix(values.astype(np.float64))
        part.sum_duplicates()
        part.sort_indices()
        data=part.data
        finite=np.isfinite(data)
        vals=data[finite]
        if len(vals):
            low,high=min(low,float(vals.min())),max(high,float(vals.max()))
        negative+=int((vals<0).sum())
        fractional+=int((vals!=np.floor(vals)).sum())
        if role=='counts':
            assert finite.all() and (vals>=0).all() and (vals==np.floor(vals)).all() and (vals<np.iinfo(np.int64).max).all()
        elif role=='published_processed_nonnegative':
            assert finite.all() and (vals>=0).all()
        elif role!='author_embedding':
            assert finite.all()
        bad=part.copy()
        bad.data=(~finite).astype('int64')
        bad.eliminate_zeros()
        part.data[~finite]=0
        part.eliminate_zeros()
        # Direct sparse reductions differ from the preparation reader and block sizes.
        total_r=np.asarray(part.sum(axis=1)).ravel()
        total_c=np.asarray(part.sum(axis=0)).ravel()
        hit_r,hit_c=part.getnnz(axis=1),part.getnnz(axis=0)
        bad_r,bad_c=np.asarray(bad.sum(axis=1)).ravel(),np.asarray(bad.sum(axis=0)).ravel()
        if axis==0:
            rs[start:stop],rn[start:stop],rb[start:stop]=total_r,hit_r,bad_r
            cs+=total_c
            cn+=hit_c
            cb+=bad_c
            for m,(tot,det) in groups.items():
                chosen=part[:,mods==m]
                tot[start:stop]=np.asarray(chosen.sum(axis=1)).ravel()
                det[start:stop]=chosen.getnnz(axis=1)
            if mt is not None:
                rmt[start:stop]=np.asarray(part[:,mt].sum(axis=1)).ravel()
        else:
            cs[start:stop],cn[start:stop],cb[start:stop]=total_c,hit_c,bad_c
            rs+=total_r
            rn+=hit_r
            rb+=bad_r
    summary={'rows':n,'features':p,'finite_nonzero_values':int(rn.sum()),'nonfinite_values':int(rb.sum()),
        'nonfinite_rows':int((rb>0).sum()),'nonfinite_columns':int((cb>0).sum()),'negative_values':negative,'fractional_values':fractional,
        'minimum_finite_stored_value':float(low) if np.isfinite(low) else 0.0,'maximum_finite_stored_value':float(high) if np.isfinite(high) else 0.0,
        'zero_information_finite_rows':int((rn==0).sum()),'all_zero_finite_columns':int((cn==0).sum()),'finite_row_sum_total':float(rs.sum())}
    return summary,{'row_sum':rs,'row_n':rn,'row_bad':rb,'col_sum':cs,'col_n':cn,'col_bad':cb,'groups':groups,'mt':rmt}


def main():
    assert not (OUT/'FINAL_COMPLETION_GATE.json').exists(),'Successful independent gate exists; do not overwrite'
    stage=json.loads((WORK/'STAGED_INPUTS.json').read_text())
    baseline=json.loads((WORK/'ORIGINAL_DOWNLOAD_BASELINE.json').read_text())
    original=table('ORIGINAL_FILES.tsv')
    check('exact_two_originals_and_two_expanded_files',len(stage)==len(baseline)==2 and len(original)==4)
    for item,base in zip(stage,baseline):
        check('original_baseline_identity_'+item['original_file'],item['original_path']==base['server_path'] and item['original_MD5']==base['md5'] and item['original_SHA256']==base['sha256'] and item['original_bytes']==int(base['bytes']))
        for kind in ['original','expanded']:
            path=Path(item[kind+'_path'])
            row=original.loc[original.path.eq(str(path))]
            check(kind+'_one_manifest_row_'+path.name,len(row)==1)
            check(kind+'_bytes_'+path.name,path.stat().st_size==item[kind+'_bytes']==int(row.iloc[0].bytes))
            check(kind+'_MD5_'+path.name,digest(path,'md5')==item[kind+'_MD5']==row.iloc[0].MD5)
            check(kind+'_SHA256_'+path.name,digest(path)==item[kind+'_SHA256']==row.iloc[0].SHA256)
    print('All compressed and expanded inputs independently rehashed.',flush=True)
    matrix_summary=table('MATRIX_QC_SUMMARY.tsv').set_index(['unit','slot'])
    axes=table('DERIVED_COORDINATE_QC.tsv.gz')
    bad_aux=table('DERIVED_NONFINITE_ROWS.tsv.gz')
    units=table('UNIT_QC_SUMMARY.tsv').set_index('unit')
    assays=table('ASSAY_INDEX.tsv')
    missing=table('SOURCE_METADATA_MISSING.tsv')
    types=json.loads((OUT/'SOURCE_METADATA_DTYPES.json').read_text())
    ids=json.loads((OUT/'IDENTIFIER_AUDIT.json').read_text())
    obs_all,seen_slots,all_donor=[],[],[]
    for unit,item in zip(['cite','multiome'],stage):
        path=Path(item['expanded_path'])
        with h5py.File(path,'r') as f:
            obs,var=ad.io.read_elem(f['obs']),ad.io.read_elem(f['var'])
            qc,feature=table(unit+'_observation_QC.tsv.gz',True),table(unit+'_feature_QC.tsv.gz',True)
            usable=table(unit+'_usable_pairs.tsv.gz',True)
            quarantine=table(unit+'_quarantine.tsv.gz',True)
            check(unit+'_all_author_observation_fields_exact',metadata_equal(obs,qc))
            check(unit+'_all_author_feature_fields_exact',metadata_equal(var,feature))
            check(unit+'_source_metadata_dtypes_recorded',types[unit]=={'obs':{k:str(v) for k,v in obs.dtypes.items()},'var':{k:str(v) for k,v in var.dtypes.items()}})
            check(unit+'_original_row_and_column_positions',np.array_equal(qc.source_row,np.arange(len(obs))) and np.array_equal(feature.source_column,np.arange(len(var))))
            check(unit+'_unique_source_cell_axis',obs.index.is_unique and (obs.index.astype(str).str.strip()!='').all())
            check(unit+'_all_source_feature_names_nonempty',(var.index.astype(str).str.strip()!='').all())
            check(unit+'_source_pointers',qc.source_file.eq(str(path)).all() and qc.counts_slot.eq('layers/counts').all())
            check(unit+'_exact_cell_UIDs',np.array_equal(qc.cell_uid,'NeurIPS2021:'+unit+':'+obs.index.astype(str)) and qc.cell_uid.is_unique)
            check(unit+'_exact_donor_UIDs',np.array_equal(qc.donor_uid,'NeurIPS2021:donor:'+obs.DonorID.astype(str)))
            check(unit+'_exact_feature_UIDs',feature.feature_uid.tolist()==[f'NeurIPS2021:{unit}:feature:{i}' for i in range(len(var))])
            mods=var.feature_types.astype(str).to_numpy()
            expected_mods={'GEX','ADT'} if unit=='cite' else {'GEX','ATAC'}
            check(unit+'_correct_modality_set',set(mods)==expected_mods)
            check(unit+'_unique_assay_specific_feature_keys',np.array_equal(feature.assay_name_key,mods+':'+var.index.astype(str)) and feature.assay_name_key.is_unique)
            check(unit+'_missing_gene_IDs_preserved',np.array_equal(feature.gene_id_missing_original,var.gene_id.isna()))
            check(unit+'_GEX_gene_IDs_complete_and_unique',var.loc[mods=='GEX','gene_id'].notna().all() and var.loc[mods=='GEX','gene_id'].is_unique)
            check(unit+'_duplicate_names_not_collapsed',np.array_equal(feature.duplicate_name_across_native_axis,var.index.duplicated(keep=False)) and ids[unit]['duplicate_name_rows']==int(var.index.duplicated(keep=False).sum()))
            for axis,source in [('obs',obs),('var',var)]:
                for name in source:
                    recorded=missing.loc[missing.unit.eq(unit)&missing.axis.eq(axis)&missing.field.eq(name)]
                    check(unit+'_'+axis+'_'+name+'_missingness',len(recorded)==1 and int(recorded.iloc[0]['missing'])==int(source[name].isna().sum()) and int(recorded.iloc[0].total)==len(source))
            expected_aux={'ADT_X_pca','ADT_X_umap','ADT_isotype_controls','GEX_X_pca','GEX_X_umap'} if unit=='cite' else {'ATAC_gene_activity','ATAC_lsi_full','ATAC_lsi_red','ATAC_umap','GEX_X_pca','GEX_X_umap'}
            check(unit+'_every_auxiliary_assay_covered',set(f['obsm'])==expected_aux)
            count_metrics=None
            for slot in ['X','layers/counts']+['obsm/'+key for key in sorted(expected_aux)]:
                role='counts' if slot=='layers/counts' else 'published_processed_nonnegative' if slot in {'X','obsm/ATAC_gene_activity'} else 'author_isotype_processed' if slot=='obsm/ADT_isotype_controls' else 'author_embedding'
                mt=(mods=='GEX')&var.index.astype(str).str.startswith('MT-')
                info,m=independent_scan(f[slot],role,mods if slot=='layers/counts' else None,mt if slot=='layers/counts' else None)
                recorded=matrix_summary.loc[(unit,slot)]
                check(unit+'_'+slot+'_role',recorded.role==role)
                for field,value in info.items():
                    check(unit+'_'+slot+'_'+field,np.isclose(float(recorded[field]),float(value),rtol=1e-12,atol=1e-8) if field=='finite_row_sum_total' and role!='counts' else recorded[field]==value)
                check(unit+'_'+slot+'_source_cell_axis',info['rows']==len(obs))
                if slot in {'X','layers/counts'}:
                    check(unit+'_'+slot+'_source_feature_axis',info['features']==len(var))
                if slot=='layers/counts':
                    count_metrics=m
                if slot.startswith('obsm/'):
                    ax=axes.loc[axes.unit.eq(unit)&axes.slot.eq(slot)].sort_values('source_column')
                    check(unit+'_'+slot+'_all_derived_columns_indexed',np.array_equal(ax.source_column,np.arange(info['features'])))
                    check(unit+'_'+slot+'_nonfinite_per_column',np.array_equal(ax.nonfinite_values,m['col_bad']))
                    check(unit+'_'+slot+'_nonzero_per_column',np.array_equal(ax.finite_nonzero_rows,m['col_n']))
                    check(unit+'_'+slot+'_usable_coordinate_rule',np.array_equal(ax.usable_derived_coordinate,(m['col_bad']==0)&(m['col_n']>0)))
                    expected_reasons=np.where(m['col_bad']>0,'contains_nonfinite',np.where(m['col_n']==0,'all_zero','finite_nonzero'))
                    check(unit+'_'+slot+'_coordinate_reasons',np.array_equal(ax.reason,expected_reasons))
                    bad=bad_aux.loc[bad_aux.unit.eq(unit)&bad_aux.slot.eq(slot)].sort_values('source_row')
                    rows=np.flatnonzero(m['row_bad'])
                    check(unit+'_'+slot+'_exact_nonfinite_rows',np.array_equal(bad.source_row,rows) and np.array_equal(bad.nonfinite_values,m['row_bad'][rows]) and bad.source_cell_id.tolist()==obs.index[rows].astype(str).tolist())
                    if slot=='obsm/ATAC_gene_activity':
                        expected_labels=list(map(str,ad.io.read_elem(f['uns/ATAC_gene_activity_var_names'])))
                    elif isinstance(f[slot],h5py.Group) and f[slot].attrs.get('encoding-type')=='dataframe':
                        original_table=ad.io.read_elem(f[slot])
                        check(unit+'_'+slot+'_table_row_alignment',original_table.index.equals(obs.index))
                        expected_labels=original_table.columns.astype(str).tolist()
                    else:
                        expected_labels=list(map(str,range(info['features'])))
                    check(unit+'_'+slot+'_exact_derived_feature_axis',ax.source_label.astype(str).tolist()==expected_labels)
                seen_slots.append((unit,slot))
                print('Independently verified full matrix:',unit,slot,flush=True)
            required=['batch','DonorID','Samplename','Site','DonorNumber','cell_type','Modality']
            valid=obs[required].notna().all(axis=1).to_numpy()&obs[required].astype(str).apply(lambda x:x.str.strip().ne('')).all(axis=1).to_numpy()
            check(unit+'_required_metadata_rule',np.array_equal(qc.metadata_valid,valid))
            for mod,(tot,det) in count_metrics['groups'].items():
                check(unit+'_'+mod+'_all_row_totals_exact',np.array_equal(qc[mod+'_published_count_sum'],tot))
                check(unit+'_'+mod+'_all_row_detected_counts_exact',np.array_equal(qc[mod+'_published_detected_features'],det))
                check(unit+'_'+mod+'_nonzero_flag',np.array_equal(qc[mod+'_nonzero_counts'],det>0))
                valid&=det>0
                idx=assays.loc[assays.unit.eq(unit)&assays.modality.eq(mod)]
                check(unit+'_'+mod+'_assay_loader_index',len(idx)==1 and idx.iloc[0].source_file==str(path) and idx.iloc[0].counts_slot=='layers/counts' and int(idx.iloc[0].rows)==len(obs) and int(idx.iloc[0].features)==int((mods==mod).sum()) and int(idx.iloc[0].usable_pairs)==int(valid.sum()))
            check(unit+'_all_feature_totals_exact',np.array_equal(feature.published_count_sum,count_metrics['col_sum']))
            check(unit+'_all_feature_detection_counts_exact',np.array_equal(feature.cells_with_nonzero_counts,count_metrics['col_n']))
            check(unit+'_usable_feature_rule',np.array_equal(feature.usable_feature,count_metrics['col_n']>0))
            check(unit+'_all_MT_sums_exact',np.array_equal(qc.GEX_published_mt_count_sum,count_metrics['mt']))
            total=count_metrics['groups']['GEX'][0]
            expected_mt=np.divide(100*count_metrics['mt'],total,out=np.full(len(obs),np.nan),where=total>0) if mt.any() else np.full(len(obs),np.nan)
            check(unit+'_published_axis_MT_percentage',np.array_equal(pd.to_numeric(qc.GEX_published_axis_mt_pct,errors='coerce'),expected_mt,equal_nan=True))
            check(unit+'_usable_pair_rule',np.array_equal(qc.usable_paired_counts,valid))
            check(unit+'_exact_usable_rows',usable.index.tolist()==obs.index[valid].tolist() and np.array_equal(usable.source_row,np.flatnonzero(valid)))
            check(unit+'_exact_quarantine_rows',quarantine.index.tolist()==obs.index[~valid].tolist() and np.array_equal(quarantine.source_row,np.flatnonzero(~valid)))
            reasons=np.where(~qc.metadata_valid,'missing_required_metadata',np.where(~valid,'zero_counts_in_one_or_more_modalities',''))
            check(unit+'_exact_quarantine_reasons',np.array_equal(qc.quarantine_reason,reasons))
            check(unit+'_paired_cells_have_both_modalities',qc.pairing_basis.eq('same_native_row_two_modalities').all() and usable.pairing_basis.eq('same_native_row_two_modalities').all())
            check(unit+'_no_HCM_labels_created',not qc.disease_label_assigned.any() and 'disease_label' not in qc)
            if unit=='cite':
                check(unit+'_official_split_preserved',np.array_equal(qc.official_split,obs.is_train.astype(str)) and qc.official_split_available.all() and set(qc.official_split)=={'train','test','iid_holdout'})
            else:
                check(unit+'_absent_official_split_not_invented','is_train' not in obs and not qc.official_split_available.any() and qc.official_split.eq('not_supplied_in_published_multiome_object').all())
            unitrow=units.loc[unit]
            check(unit+'_unit_summary',unitrow.source_rows==len(obs) and unitrow.source_features==len(var) and unitrow.usable_paired_cells==valid.sum() and unitrow.quarantined_cells==(~valid).sum() and unitrow.donors==obs.DonorID.nunique() and unitrow.batches==obs.batch.nunique() and unitrow.mt_genes_in_published_axis==mt.sum())
            obs_all.append(qc)
            donor=obs[['DonorID','DonorNumber','Site','batch','Samplename','DonorAge','DonorGender','VendorLot']].copy()
            donor.insert(0,'unit',unit)
            donor['official_split']=qc.official_split.to_numpy()
            donor['usable_paired_counts']=valid
            all_donor.append(donor)
    check('exact_all_fifteen_native_matrices',len(seen_slots)==15 and set(matrix_summary.index)==set(seen_slots) and matrix_summary.index.is_unique)
    check('exact_four_assay_indices',len(assays)==4 and not assays[['unit','modality']].duplicated().any())
    combined=pd.concat(obs_all)
    check('global_cell_UIDs_do_not_pair_separate_assays',combined.cell_uid.is_unique)
    donors=pd.concat(all_donor,ignore_index=True)
    expected=donors.groupby(['unit','DonorID','DonorNumber','Site','batch','Samplename','official_split'],observed=True,dropna=False).agg(cells=('usable_paired_counts','size'),usable_pairs=('usable_paired_counts','sum')).reset_index()
    recorded=table('DONOR_SITE_SPLIT_SUMMARY.tsv')
    check('all_donor_site_split_group_rows',len(recorded)==len(expected) and all(np.array_equal(recorded[c].astype(str),expected[c].astype(str)) for c in expected))
    expected=donors.groupby('DonorID',observed=True).agg(cells=('unit','size'),units=('unit','nunique'),sites=('Site','nunique'),age_values=('DonorAge','nunique'),gender_values=('DonorGender','nunique'),donor_number_values=('DonorNumber','nunique'),vendor_lot_values=('VendorLot','nunique')).reset_index()
    recorded=table('DONOR_AUDIT.tsv')
    check('cross_site_and_modality_donors_bound',len(recorded)==len(expected) and all(np.array_equal(recorded[c].astype(str),expected[c].astype(str)) for c in expected))
    check('source_donor_demographics_consistent',expected[['age_values','gender_values','donor_number_values','vendor_lot_values']].eq(1).all().all())
    policy=json.loads((OUT/'qc_policy.json').read_text())
    check('scope_preserves_counts_and_avoids_HCM_rules',all(policy[k] is False for k in ['HCM_rules_applied','HCM_labels_assigned','new_normalization','new_doublet_filter','imputation','model_training']) and policy['originals_read_only'])
    manifest=table('OUTPUT_MANIFEST.tsv')
    for row in manifest.itertuples(index=False):
        path=Path(row.server_path)
        check('output_bytes_'+row.file,path.stat().st_size==row.bytes)
        check('output_MD5_'+row.file,digest(path,'md5')==row.MD5)
        check('output_SHA256_'+row.file,digest(path)==row.SHA256)
    now=datetime.now(timezone.utc).isoformat()
    verification={'status':'passed','checked_utc':now,'checks_passed':len(CHECKS),'matrices_fully_recomputed':15,
        'compressed_and_expanded_files_verified':4,'source_cells':len(combined),'usable_paired_cells':int(combined.usable_paired_counts.sum()),'unique_donors':int(donors.DonorID.nunique()),'checks':CHECKS}
    (OUT/'INDEPENDENT_VERIFICATION.json').write_text(json.dumps(verification,indent=2)+'\n')
    gate={'status':'passed','checked_utc':now,'independent_checks':len(CHECKS),'matrices':15,'verification_SHA256':digest(OUT/'INDEPENDENT_VERIFICATION.json'),
        'prepared_manifest_SHA256':digest(OUT/'OUTPUT_MANIFEST.tsv'),'source_immutable':True,'scope':'all two downloaded GSE194122 release objects; original counts plus usable modality/pairing indices',
        'not_a_certificate_of_HCM_applicability_or_donor_independent_training_split':True}
    (OUT/'FINAL_COMPLETION_GATE.json').write_text(json.dumps(gate,indent=2)+'\n')
    (OUT/'VERIFICATION_REPORT.md').write_text(f"# NeurIPS2021 独立验收\n\n{now}；{len(CHECKS)} 项检查全部通过。2原件及2展开件摘要与下载/解包基线一致；15个完整矩阵重新扫描，逐细胞、逐特征及配对/划分索引复算通过。共 {len(combined)} 个细胞、{int(donors.DonorID.nunique())} 位供体、{int(combined.usable_paired_counts.sum())} 个可用双模态细胞。\n\n缺失的ATAC降维坐标按维度隔离，不填值；RNA与蛋白同名特征不合并。所有计数/作者元数据仍可溯源，未赋HCM疾病标签或生成训练划分。最终交付还应核对脚本退出码、进程结束及本地/远端源码一致。\n")
    records=[]
    for p in sorted(OUT.iterdir()):
        if p.is_file() and p.name!='VERIFIED_OUTPUT_MANIFEST.tsv':
            records.append({'file':p.name,'server_path':str(p),'bytes':p.stat().st_size,'MD5':digest(p,'md5'),'SHA256':digest(p)})
    pd.DataFrame(records).to_csv(OUT/'VERIFIED_OUTPUT_MANIFEST.tsv',sep='\t',index=False)
    print('All',len(CHECKS),'independent checks passed.',flush=True)


if __name__=='__main__':
    try:
        main()
    except Exception as error:
        OUT.mkdir(parents=True,exist_ok=True)
        (OUT/'VERIFICATION_FAILURE.json').write_text(json.dumps({'error':repr(error),'traceback':traceback.format_exc(),'checks':CHECKS},indent=2)+'\n')
        raise
