"""Audit native multimodal counts and publish reversible usable row/column indices."""
from datetime import datetime, timezone
from pathlib import Path
import hashlib
import json
import shutil
import anndata as ad
import h5py
import numpy as np
import pandas as pd
from scipy import sparse

ROOT = Path.home() / 'xiaoe-internship-data'
WORK = ROOT / '_processing/all_data_cleaning_20261111/NeurIPS2021'
OUT = ROOT / 'processed/NeurIPS2021_qc_20261111'
UNITS = ['cite', 'multiome']


def digest(path, algorithm='sha256'):
    with path.open('rb') as handle:
        return hashlib.file_digest(handle, algorithm).hexdigest()


def save_json(path, value):
    path.write_text(json.dumps(value, indent=2, allow_nan=False) + '\n')


def native_blocks(node):
    """Iterate in native sparse orientation; never repeatedly slice CSC by rows."""
    if isinstance(node, h5py.Group) and node.attrs.get('encoding-type') in {'csr_matrix', 'csc_matrix'}:
        is_csr = node.attrs['encoding-type'] == 'csr_matrix'
        shape = tuple(node.attrs['shape'])
        pointer = node['indptr'][:]
        n_axis, other = shape if is_csr else shape[::-1]
        assert len(pointer) == n_axis + 1 and pointer[0] == 0
        assert pointer[-1] == len(node['data']) == len(node['indices']) and (np.diff(pointer) >= 0).all()
        step = 1536 if is_csr else 192
        for start in range(0, n_axis, step):
            stop = min(n_axis, start + step)
            low, high = int(pointer[start]), int(pointer[stop])
            indexes = node['indices'][low:high]
            assert not len(indexes) or (indexes.min() >= 0 and indexes.max() < other)
            constructor = sparse.csr_matrix if is_csr else sparse.csc_matrix
            block_shape = (stop-start, shape[1]) if is_csr else (shape[0], stop-start)
            matrix = constructor((node['data'][low:high].astype(np.float64), indexes,
                pointer[start:stop+1]-low), shape=block_shape)
            # Canonicalize only this in-memory copy; disk storage remains read-only.
            matrix.sum_duplicates()
            matrix.sort_indices()
            yield start, matrix, is_csr
    else:
        table = ad.io.read_elem(node) if isinstance(node, h5py.Group) else node
        for start in range(0, table.shape[0], 1536):
            values = table.iloc[start:start+1536].to_numpy() if isinstance(table, pd.DataFrame) else table[start:start+1536]
            yield start, sparse.csr_matrix(values.astype(np.float64)), True


def scan(node, role, modalities=None, mt=None):
    if isinstance(node, h5py.Group) and node.attrs.get('encoding-type') in {'csr_matrix','csc_matrix'}:
        shape = tuple(node.attrs['shape'])
    elif isinstance(node, h5py.Group):
        shape = ad.io.read_elem(node).shape
    else:
        shape = node.shape
    n, p = shape
    row_sum, row_n, row_bad = np.zeros(n), np.zeros(n, dtype='int64'), np.zeros(n,dtype='int64')
    col_sum, col_n, col_bad = np.zeros(p), np.zeros(p,dtype='int64'), np.zeros(p,dtype='int64')
    by_modality = {m: [np.zeros(n),np.zeros(n,dtype='int64')] for m in np.unique(modalities)} if modalities is not None else {}
    row_mt = np.zeros(n)
    lo, hi, neg, frac, zeros = np.inf, -np.inf, 0, 0, 0
    for start, matrix, csr in native_blocks(node):
        data = matrix.data
        finite = np.isfinite(data)
        local_bad = matrix.copy()
        local_bad.data = (~finite).astype('int64')
        local_bad.eliminate_zeros()
        clean = matrix.copy()
        clean.data[~finite] = 0
        clean.eliminate_zeros()
        vals = data[finite]
        if len(vals):
            lo,hi=min(lo,float(vals.min())),max(hi,float(vals.max()))
        neg += int((vals<0).sum())
        frac += int((vals!=np.rint(vals)).sum())
        zeros += int((vals==0).sum())
        if role == 'counts':
            assert finite.all() and (vals >= 0).all() and (vals==np.rint(vals)).all()
            assert (vals < np.iinfo(np.int64).max).all()
        elif role == 'published_processed_nonnegative':
            assert finite.all() and (vals >= 0).all()
        elif role != 'author_embedding':
            assert finite.all(), node.name
        sums_r, sums_c = np.asarray(clean.sum(axis=1)).ravel(), np.asarray(clean.sum(axis=0)).ravel()
        n_r, n_c = clean.getnnz(axis=1), clean.getnnz(axis=0)
        b_r, b_c = np.asarray(local_bad.sum(axis=1)).ravel(), np.asarray(local_bad.sum(axis=0)).ravel()
        if csr:
            stop=start+matrix.shape[0]
            row_sum[start:stop],row_n[start:stop],row_bad[start:stop]=sums_r,n_r,b_r
            col_sum+=sums_c
            col_n+=n_c
            col_bad+=b_c
            for mod, (total, detected) in by_modality.items():
                chosen=clean[:,modalities==mod]
                total[start:stop]=np.asarray(chosen.sum(axis=1)).ravel()
                detected[start:stop]=chosen.getnnz(axis=1)
            if mt is not None:
                row_mt[start:stop]=np.asarray(clean[:,mt].sum(axis=1)).ravel()
        else:
            stop=start+matrix.shape[1]
            col_sum[start:stop],col_n[start:stop],col_bad[start:stop]=sums_c,n_c,b_c
            row_sum+=sums_r
            row_n+=n_r
            row_bad+=b_r
    info={'rows':n,'features':p,'role':role,'finite_nonzero_values':int(row_n.sum()),
        'nonfinite_values':int(row_bad.sum()),'nonfinite_rows':int((row_bad>0).sum()),
        'nonfinite_columns':int((col_bad>0).sum()),'negative_values':neg,'fractional_values':frac,
        'minimum_finite_stored_value':lo if np.isfinite(lo) else 0.0,
        'maximum_finite_stored_value':hi if np.isfinite(hi) else 0.0,
        'zero_information_finite_rows':int((row_n==0).sum()),'all_zero_finite_columns':int((col_n==0).sum()),
        'finite_row_sum_total':float(row_sum.sum())}
    return info, {'row_sum':row_sum,'row_n':row_n,'row_bad':row_bad,'col_sum':col_sum,'col_n':col_n,
        'col_bad':col_bad,'modality':by_modality,'mt':row_mt}


def required_valid(obs):
    names=['batch','DonorID','Samplename','Site','DonorNumber','cell_type','Modality']
    return obs[names].notna().all(axis=1).to_numpy() & obs[names].astype(str).apply(lambda x:x.str.strip().ne('')).all(axis=1).to_numpy()


def main():
    OUT.mkdir(parents=True,exist_ok=True)
    assert not (OUT/'EXPORT_SUCCESS.json').exists(), 'Successful preparation exists; verify before rerun'
    staged=json.loads((WORK/'STAGED_INPUTS.json').read_text())
    original_baseline=json.loads((WORK/'ORIGINAL_DOWNLOAD_BASELINE.json').read_text())
    assert len(staged)==len(original_baseline)==2
    original_records=[]
    for stage,baseline in zip(staged,original_baseline):
        assert stage['original_path']==baseline['server_path']
        for prefix,path in [('original',Path(stage['original_path'])),('expanded',Path(stage['expanded_path']))]:
            md5,sha=digest(path,'md5'),digest(path)
            assert path.stat().st_size==stage[prefix+'_bytes'] and md5==stage[prefix+'_MD5'] and sha==stage[prefix+'_SHA256']
            if prefix=='original':
                assert md5==baseline['md5'] and sha==baseline['sha256'] and path.stat().st_size==int(baseline['bytes'])
            original_records.append({'kind':prefix,'path':str(path),'bytes':path.stat().st_size,'MD5':md5,'SHA256':sha,'source_url':stage['source_url']})
    pd.DataFrame(original_records).to_csv(OUT/'ORIGINAL_FILES.tsv',sep='\t',index=False)
    matrices,axes,bad_aux,units,assays,missings,donor_tables=[],[],[],[],[],[],[]
    metadata_types,collisions={},{}
    for unit,stage in zip(UNITS,staged):
        path=Path(stage['expanded_path'])
        with h5py.File(path,'r') as f:
            obs,var=ad.io.read_elem(f['obs']),ad.io.read_elem(f['var'])
            assert obs.index.is_unique and (obs.index.astype(str).str.strip()!='').all()
            assert (var.index.astype(str).str.strip()!='').all()
            mods=var.feature_types.astype(str).to_numpy()
            assert set(mods)==({'GEX','ADT'} if unit=='cite' else {'GEX','ATAC'})
            assert var.loc[mods=='GEX','gene_id'].notna().all() and var.loc[mods=='GEX','gene_id'].is_unique
            metadata_types[unit]={'obs':{k:str(v) for k,v in obs.dtypes.items()},'var':{k:str(v) for k,v in var.dtypes.items()}}
            for axis,table in [('obs',obs),('var',var)]:
                for column in table:
                    missings.append({'unit':unit,'axis':axis,'field':column,'missing':int(table[column].isna().sum()),'total':len(table)})
            slots=['X','layers/counts']+['obsm/'+key for key in f['obsm']]
            count_metrics=None
            mt=(mods=='GEX') & var.index.astype(str).str.startswith('MT-')
            for slot in slots:
                role='counts' if slot=='layers/counts' else 'published_processed_nonnegative' if slot in {'X','obsm/ATAC_gene_activity'} else 'author_isotype_processed' if slot=='obsm/ADT_isotype_controls' else 'author_embedding'
                node=f[slot]
                if isinstance(node,h5py.Group) and node.attrs.get('encoding-type')=='dataframe':
                    table=ad.io.read_elem(node)
                    assert table.index.equals(obs.index) and table.columns.is_unique
                    labels=table.columns.astype(str).tolist()
                elif slot=='obsm/ATAC_gene_activity':
                    labels=list(map(str,ad.io.read_elem(f['uns/ATAC_gene_activity_var_names'])))
                    assert len(set(labels))==len(labels)
                else:
                    labels=None
                info,metrics=scan(node,role,mods if slot=='layers/counts' else None,mt if slot=='layers/counts' else None)
                assert info['rows']==len(obs)
                if slot in {'X','layers/counts'}:
                    assert info['features']==len(var)
                matrices.append({'unit':unit,'slot':slot,**info})
                if slot=='layers/counts':
                    count_metrics=metrics
                if slot.startswith('obsm/'):
                    if labels is not None:
                        assert len(labels)==info['features']
                    for col in range(info['features']):
                        axes.append({'unit':unit,'slot':slot,'source_column':col,'source_label':labels[col] if labels else str(col),
                            'nonfinite_values':int(metrics['col_bad'][col]),'finite_nonzero_rows':int(metrics['col_n'][col]),
                            'usable_derived_coordinate':bool(metrics['col_bad'][col]==0 and metrics['col_n'][col]>0),
                            'reason':'contains_nonfinite' if metrics['col_bad'][col] else 'all_zero' if not metrics['col_n'][col] else 'finite_nonzero'})
                    for row in np.flatnonzero(metrics['row_bad']):
                        bad_aux.append({'unit':unit,'slot':slot,'source_row':int(row),'source_cell_id':str(obs.index[row]),'nonfinite_values':int(metrics['row_bad'][row])})
                print('Fully scanned',unit,slot,info['rows'],info['features'],'nonfinite',info['nonfinite_values'],flush=True)
            qc=obs.copy()
            qc.index.name='source_cell_id'
            qc['source_row']=np.arange(len(obs))
            qc['source_file']=str(path)
            qc['counts_slot']='layers/counts'
            qc['cell_uid']='NeurIPS2021:'+unit+':'+obs.index.astype(str)
            qc['donor_uid']='NeurIPS2021:donor:'+obs.DonorID.astype(str).to_numpy()
            qc['official_split']=obs.is_train.astype(str).to_numpy() if 'is_train' in obs else 'not_supplied_in_published_multiome_object'
            qc['official_split_available']='is_train' in obs
            qc['metadata_valid']=required_valid(obs)
            qc['pairing_basis']='same_native_row_two_modalities'
            qc['disease_label_assigned']=False
            valid=qc.metadata_valid.to_numpy(copy=True)
            for mod,(total,detected) in count_metrics['modality'].items():
                qc[mod+'_published_count_sum']=total
                qc[mod+'_published_detected_features']=detected
                qc[mod+'_nonzero_counts']=detected>0
                valid &= detected>0
            qc['GEX_published_mt_count_sum']=count_metrics['mt']
            gex_total=count_metrics['modality']['GEX'][0]
            qc['GEX_published_axis_mt_pct']=np.divide(100*count_metrics['mt'],gex_total,out=np.full(len(obs),np.nan),where=(gex_total>0)) if mt.any() else np.full(len(obs),np.nan)
            qc['mt_percentage_scope']='published_feature_axis_only' if mt.any() else 'no_mitochondrial_genes_in_published_feature_axis'
            qc['usable_paired_counts']=valid
            qc['quarantine_reason']=np.where(~qc.metadata_valid,'missing_required_metadata',np.where(~valid,'zero_counts_in_one_or_more_modalities',''))
            qc.to_csv(OUT/f'{unit}_observation_QC.tsv.gz',sep='\t')
            columns=['source_row','source_file','counts_slot','cell_uid','donor_uid','official_split','pairing_basis']
            qc.loc[valid,columns].to_csv(OUT/f'{unit}_usable_pairs.tsv.gz',sep='\t')
            qc.loc[~valid].to_csv(OUT/f'{unit}_quarantine.tsv.gz',sep='\t')
            feature=var.copy()
            feature.index.name='source_feature_name'
            feature['source_column']=np.arange(len(var))
            feature['feature_uid']=[f'NeurIPS2021:{unit}:feature:{i}' for i in range(len(var))]
            feature['assay_name_key']=var.feature_types.astype(str).to_numpy()+':'+var.index.astype(str)
            feature['gene_id_missing_original']=var.gene_id.isna().to_numpy()
            feature['duplicate_name_across_native_axis']=var.index.duplicated(keep=False)
            feature['published_count_sum']=count_metrics['col_sum']
            feature['cells_with_nonzero_counts']=count_metrics['col_n']
            feature['usable_feature']=count_metrics['col_n']>0
            assert feature.feature_uid.is_unique and feature.assay_name_key.is_unique
            feature.to_csv(OUT/f'{unit}_feature_QC.tsv.gz',sep='\t')
            collisions[unit]={'duplicate_name_rows':int(var.index.duplicated(keep=False).sum()),'duplicate_name_groups':int(var.index[var.index.duplicated(keep=False)].nunique()),'unique_assay_name_keys':bool(feature.assay_name_key.is_unique),'names_never_merged':True}
            for mod in sorted(set(mods)):
                columns_mod=np.flatnonzero(mods==mod)
                assays.append({'unit':unit,'modality':mod,'rows':len(obs),'features':len(columns_mod),'usable_pairs':int(valid.sum()),
                    'usable_features':int(feature.loc[mods==mod,'usable_feature'].sum()),'source_file':str(path),'counts_slot':'layers/counts',
                    'row_index':f'{unit}_usable_pairs.tsv.gz','column_index':f'{unit}_feature_QC.tsv.gz','column_selector':'feature_types == '+mod+' and usable_feature',
                    'semantics':{'GEX':'published_unnormalized_RNA_counts','ADT':'published_unnormalized_antibody_counts','ATAC':'published_unnormalized_peak_counts'}[mod]})
            units.append({'unit':unit,'source_rows':len(obs),'source_features':len(var),'usable_paired_cells':int(valid.sum()),'quarantined_cells':int((~valid).sum()),'donors':int(obs.DonorID.nunique()),'batches':int(obs.batch.nunique()),'mt_genes_in_published_axis':int(mt.sum()),'source_path':str(path),'official_split_available':'is_train' in obs})
            donor=obs[['DonorID','DonorNumber','Site','batch','Samplename','DonorAge','DonorGender','VendorLot']].copy()
            donor.insert(0,'unit',unit)
            donor['official_split']=qc.official_split.to_numpy()
            donor['usable_paired_counts']=valid
            donor_tables.append(donor)
    pd.DataFrame(matrices).to_csv(OUT/'MATRIX_QC_SUMMARY.tsv',sep='\t',index=False)
    pd.DataFrame(axes).to_csv(OUT/'DERIVED_COORDINATE_QC.tsv.gz',sep='\t',index=False)
    pd.DataFrame(bad_aux,columns=['unit','slot','source_row','source_cell_id','nonfinite_values']).to_csv(OUT/'DERIVED_NONFINITE_ROWS.tsv.gz',sep='\t',index=False)
    pd.DataFrame(units).to_csv(OUT/'UNIT_QC_SUMMARY.tsv',sep='\t',index=False)
    pd.DataFrame(assays).to_csv(OUT/'ASSAY_INDEX.tsv',sep='\t',index=False)
    pd.DataFrame(missings).to_csv(OUT/'SOURCE_METADATA_MISSING.tsv',sep='\t',index=False)
    donors=pd.concat(donor_tables,ignore_index=True)
    grouped=donors.groupby(['unit','DonorID','DonorNumber','Site','batch','Samplename','official_split'],observed=True,dropna=False).agg(cells=('usable_paired_counts','size'),usable_pairs=('usable_paired_counts','sum')).reset_index()
    grouped.to_csv(OUT/'DONOR_SITE_SPLIT_SUMMARY.tsv',sep='\t',index=False)
    donor_audit=donors.groupby('DonorID',observed=True).agg(cells=('unit','size'),units=('unit','nunique'),sites=('Site','nunique'),age_values=('DonorAge','nunique'),gender_values=('DonorGender','nunique'),donor_number_values=('DonorNumber','nunique'),vendor_lot_values=('VendorLot','nunique')).reset_index()
    assert donor_audit[['age_values','gender_values','donor_number_values','vendor_lot_values']].eq(1).all().all()
    donor_audit.to_csv(OUT/'DONOR_AUDIT.tsv',sep='\t',index=False)
    save_json(OUT/'SOURCE_METADATA_DTYPES.json',metadata_types)
    save_json(OUT/'IDENTIFIER_AUDIT.json',collisions)
    policy={'dataset':'NeurIPS2021_GSE194122','source':'GSE194122 curated BMMC CITEseq and Multiome releases',
        'counts_location':'layers/counts; X and all obsm representations retain their published processed meaning',
        'originals_read_only':True,'primary_delivery':'original/expanded H5AD plus exact usable row/column indices, no duplicated large matrices',
        'count_requirements':['finite','nonnegative','integer valued','nonzero per paired modality','required source identity fields present'],
        'features':'Every source column retained; all-zero columns marked outside usable-feature indices; same names in GEX and ADT are distinct',
        'missing_gene_IDs':'Do not guess or fill; required Ensembl IDs checked for GEX only, ADT and peaks retain native identifiers',
        'derived_nonfinite':'Preserve original values; exclude any derived column with nonfinite entries from usable coordinate indices; never fill zeros or delete cells solely for missing embeddings',
        'missing_pseudotime':'Preserve missingness; optional author annotations are not mandatory cell quality criteria',
        'pairing':'RNA/ADT or RNA/ATAC share rows within a source object; CITE and Multiome objects are separate cells and not paired by matching barcode text',
        'split':'Preserve cite is_train train/test/iid_holdout exactly; no official row-level split supplied in this multiome object, do not invent one',
        'independence':'Same donor across sites/modalities remains one donor; official challenge split is not claimed donor-independent',
        'HCM_rules_applied':False,'HCM_labels_assigned':False,'new_normalization':False,'new_doublet_filter':False,'imputation':False,'model_training':False}
    save_json(OUT/'qc_policy.json',policy)
    for name in ['GEO_SOURCE.txt','SOURCE_LOOKUP.json','STAGED_INPUTS.json','INPUT_INSPECTION.json','TARGETED_INSPECTION.json']:
        if (WORK/name).exists():
            shutil.copyfile(WORK/name,OUT/name)
    report='# NeurIPS2021：分模态检查与可用索引\n\n'
    report+='来源为 GEO GSE194122 的健康成人骨髓 BMMC，作者已做质控和双细胞去除；不是 HCM 或心肌队列。GEO 整体设计描述12位供体，本次下载的两个发布对象实际覆盖 '+str(donors.DonorID.nunique())+' 位不同 DonorID，不补造缺失供体。\n\n'
    for item in units:
        report+=f"- {item['unit']}：{item['source_rows']} 细胞、{item['source_features']} 特征；可用双模态配对 {item['usable_paired_cells']}，隔离 {item['quarantined_cells']}；{item['donors']} 供体、{item['batches']} 捕获/站点批次。\n"
    report+='\n计数只取 layers/counts，完整扫描确认其数值域；X、ATAC gene activity、同型对照和嵌入不冒充未归一化计数。所有源特征、作者元数据和计数保持原样，大矩阵通过原件/展开件与行列索引交付。MT仅在发布基因轴计算并注明范围，不据此应用HCM阈值。\n\n'
    report+='CITE 的36个重复名称组各含RNA/ADT两个不同测量，共72列；以模态和原列号区分，不合并。98个ADT及所有ATAC峰未提供gene_id是原件事实，保留缺失，不把它当成RNA基因注释缺失。伪时间的缺失保留，不填值、不删细胞。\n\n'
    report+='Multiome 的 ATAC_lsi_red 存在非有限坐标；DERIVED_COORDINATE_QC和DERIVED_NONFINITE_ROWS逐项记录，仅有限且有信息的维度标为可用。原坐标不改，其他完整表示可通过索引使用，缺失坐标不构成计数细胞删除理由。\n\n'
    report+='CITE train/test/iid_holdout逐行保留；当前Multiome发布对象无is_train，不创造官方划分。跨站点同DonorID绑定，不能声称官方划分是供体独立测试。两个对象中的细胞不因条码文本一致就配对；配对只发生在同对象同一行的两种模态。未进行新的归一化、整合、插补、差异表达、模型训练或HCM标签赋值。\n\n'
    report+='本准备阶段产物须经独立脚本重读全矩阵、逐行/逐列复算和来源摘要核查后才能认证完成。\n'
    (OUT/'REPORT.md').write_text(report)
    manifest=[]
    for path in sorted(OUT.iterdir()):
        if path.is_file():
            manifest.append({'file':path.name,'server_path':str(path),'bytes':path.stat().st_size,'MD5':digest(path,'md5'),'SHA256':digest(path),'purpose':'full audit or usable index; large data remain in source files'})
    pd.DataFrame(manifest).to_csv(OUT/'OUTPUT_MANIFEST.tsv',sep='\t',index=False)
    save_json(OUT/'EXPORT_SUCCESS.json',{'status':'prepared_not_independently_verified','units':units,'matrices':len(matrices),'completed_utc':datetime.now(timezone.utc).isoformat()})
    print('Preparation completed; independent verification still required.',flush=True)


if __name__=='__main__':
    main()
