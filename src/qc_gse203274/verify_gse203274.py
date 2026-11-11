"""Independent source-to-output verification for all published pediatric nuclei."""
from datetime import datetime, timezone
from pathlib import Path
import hashlib
import json
import re
import tarfile
import traceback

import anndata as ad
import numpy as np
import pandas as pd
from scipy.io import mmread
from threadpoolctl import threadpool_limits

ROOT=Path.home()/'xiaoe-internship-data'
WORK=ROOT/'_processing/all_data_cleaning_20261111/GSE203274'
INPUT=WORK/'AllNuclei_snRNA_rawCount'
OUT=ROOT/'processed/pediatric_GSE203274_qc_20261111'
REF=ROOT/'processed/SCP1303_review_20261111/identity'
CHECKS=[]


def check(name,passed,detail=None):
    item={'name':name,'passed':bool(passed)}
    if detail is not None:item['detail']=detail
    CHECKS.append(item)
    if not passed:raise AssertionError(name)


def digest(path,algorithm='sha256'):
    with Path(path).open('rb') as handle:return hashlib.file_digest(handle,algorithm).hexdigest()


def same(x,y):
    if x.shape!=y.shape:return False
    difference=(x-y).tocsr();difference.eliminate_zeros()
    return difference.nnz==0


def age_parts(text):
    m=re.fullmatch(r'(\d+)y(?:[_/](\d+)m)?(?:[_/](\d+)d)?',str(text))
    if m is None:raise ValueError(text)
    return tuple(int(v or 0) for v in m.groups())


def run():
    baseline=json.loads((WORK/'ORIGINAL_DOWNLOAD_BASELINE.json').read_text(encoding='utf-8-sig'))
    check('two_original_files_only_cardiac_modality',len(baseline)==2)
    for item in baseline:
        path=Path(item['server_path'])
        check('original_bytes_'+path.name,path.stat().st_size==int(item['bytes']))
        check('original_MD5_'+path.name,digest(path,'md5')==item['md5'])
    print('Original download sizes and MD5 verified',flush=True)
    archive=ROOT/'validation/pediatric_GSE203274/GSE203274_AllNuclei_snRNA_rawCount.tar.gz'
    archived_files=[]
    with tarfile.open(archive,'r:gz') as tar:
        for member in tar:
            if not member.isfile():continue
            path=WORK/member.name
            archived_files.append(member.name)
            with tar.extractfile(member) as handle:archived_sha=hashlib.file_digest(handle,'sha256').hexdigest()
            check('extracted_member_'+Path(member.name).name,path.stat().st_size==member.size and digest(path)==archived_sha)
    check('all_three_count_members_verified',len(archived_files)==3)
    genes=pd.read_csv(INPUT/'genes.tsv.gz',sep='\t',header=None,keep_default_na=False)[0].astype(str)
    cells=pd.Index(pd.read_csv(INPUT/'barcodes.tsv.gz',sep='\t',header=None,keep_default_na=False)[0].astype(str))
    metadata=pd.read_csv(ROOT/'validation/pediatric_GSE203274/GSE203274_AllNuclei_snRNA_metadata.csv.gz',index_col=0,keep_default_na=False).loc[cells]
    with threadpool_limits(limits=4):source=mmread(INPUT/'expression_matrix.mtx.gz').T.tocsr()
    source.sum_duplicates();source.eliminate_zeros();source.sort_indices()
    check('real_integer_counts_not_scaled_expression',np.isfinite(source.data).all() and (source.data>0).all() and np.equal(source.data,np.rint(source.data)).all())
    raw=ad.read_h5ad(OUT/'all_source_counts.h5ad')
    check('source_shape_and_nnz',source.shape==(157293,29266) and source.nnz==448995451 and raw.shape==source.shape)
    check('full_source_cell_axis_unchanged',raw.obs_names.equals(cells) and cells.is_unique)
    check('full_source_gene_axis_unchanged',np.array_equal(raw.var.gene_symbol.to_numpy(),genes.to_numpy()) and np.array_equal(raw.var.source_row,np.arange(len(genes))))
    check('full_source_every_count_value_unchanged',same(raw.X,source) and same(raw.layers['counts'],source))
    print('Every original published matrix value matches full source output',flush=True)
    total=np.asarray(source.sum(axis=1)).ravel();detected=np.diff(source.indptr)
    mito=genes.str.startswith('MT-').to_numpy()
    mt_counts=np.asarray(source[:,mito].sum(axis=1)).ravel();pct=100*mt_counts/total
    strict=(detected>=500)&(detected<=5000)&(pct<5)
    check('MT_axis_and_true_percentage',mito.sum()==13 and np.allclose(pct,metadata['percent.mt'],rtol=0,atol=1e-8) and not np.allclose(pct/100,metadata['percent.mt'],rtol=0,atol=1e-8))
    check('author_count_and_detected_fields_exact',np.array_equal(total,pd.to_numeric(metadata.nCount_RNA)) and np.array_equal(detected,pd.to_numeric(metadata.nFeature_RNA)))
    check('recomputed_QC_metrics',np.array_equal(raw.obs.total_counts,total) and np.array_equal(raw.obs.n_genes_by_counts,detected) and np.array_equal(raw.obs.mt_counts,mt_counts) and np.allclose(raw.obs.pct_counts_mt,pct,rtol=0,atol=1e-10))
    check('strict_baseline_recomputed',np.array_equal(raw.obs.pass_strict_numeric,strict) and strict.sum()==141301)
    check('inherited_author_QC_is_visible',detected.min()>=700 and detected.max()<=10000 and pct.max()<=.2+1e-8)
    # Metadata preservation is checked field by field, retaining case-distinct diagnosis columns.
    for column in metadata:
        expected=metadata[column]
        actual=raw.obs['author_'+column]
        if expected.dtype.kind in 'fiu':valid=np.allclose(pd.to_numeric(actual),expected,rtol=0,atol=1e-10)
        else:valid=np.array_equal(actual.astype(str).to_numpy(),expected.astype(str).to_numpy())
        check('original_author_metadata_'+column,valid)
    check('no_case_folding_clinical_columns','author_diagnosis' in raw.obs and 'author_Diagnosis' in raw.obs)
    groups={'Donor':('Control',0),'HCM':('HCM',1),'DCM':('DCM',-1),'TOF':('TOF',-1),'Neo_HLHS':('HLHS',-1),'HF_HLHS':('HLHS',-1)}
    expected_disease=metadata.ClinicalRank.map({k:v[0] for k,v in groups.items()})
    expected_label=metadata.ClinicalRank.map({k:v[1] for k,v in groups.items()})
    check('clinical_disease_not_procedure_status',np.array_equal(raw.obs.disease,expected_disease) and np.array_equal(raw.obs.label,expected_label) and metadata.ClinicalRank.eq(metadata['Diagnosis']).all())
    check('other_pathologies_never_assigned_control',raw.obs.loc[raw.obs.disease.isin(['DCM','HLHS','TOF']),'label'].eq(-1).all())
    check('source_capture_prefix_matches_metadata',np.array_equal(raw.obs.sample_id.astype(str),cells.str.rsplit('_',n=1).str[0]))
    check('donor_binding_uses_author_labID',np.array_equal(raw.obs.donor_id.astype(str),'GSE203274_'+metadata.labID.astype(str)) and raw.obs.donor_id.nunique()==13)
    check('single_HCM_donor_two_captures',raw.obs.loc[raw.obs.disease.eq('HCM'),'donor_id'].nunique()==1 and raw.obs.loc[raw.obs.disease.eq('HCM'),'sample_id'].nunique()==2 and raw.obs.loc[raw.obs.disease.eq('HCM'),'donor_original'].eq('P86').all())
    check('four_control_donors',raw.obs.loc[raw.obs.disease.eq('Control'),'donor_id'].nunique()==4)
    regions=raw.obs.loc[raw.obs.donor_original.eq('13_198')]
    check('same_WU198_LV_RV_never_split_as_two_donors',set(regions.tissue_region)=={'LV','RV'} and regions.donor_id.nunique()==1 and regions.source_patientID.nunique()==2)
    check('echoEF_is_captureID_not_numeric_EF',metadata.echoEF.astype(str).eq(metadata['orig.ident'].astype(str)).all() and not raw.obs.numeric_ejection_fraction_available.any())
    ages=np.array([age_parts(v) for v in metadata.age])
    check('author_age_components_not_coarse_ageCont',np.array_equal(raw.obs.age_year_component,ages[:,0]) and np.array_equal(raw.obs.age_month_component,ages[:,1]) and np.array_equal(raw.obs.age_day_component,ages[:,2]))
    check('approximate_ages_explicit',np.allclose(raw.obs.age_years_approx_from_author_text,ages[:,0]+ages[:,1]/12+ages[:,2]/365.25,rtol=0,atol=1e-12))
    geo=json.loads((WORK/'GEO_SAMPLES_GSE203274.json').read_text(encoding='utf-8-sig'))
    geo={r['accession']:r for r in geo if r['source_name_ch1']==['Cardiac']}
    check('only_thirty_cardiac_GEO_records_used',len(geo)==30)
    samples=pd.read_csv(OUT/'sample_metadata.tsv',sep='\t')
    for row in samples.itertuples():
        obs=raw.obs.loc[raw.obs.sample_id.eq(row.sample_id)]
        expected_alias=(('WU198LV' if row.tissue_region=='LV' else 'WU198RV') if row.donor_original=='13_198' else 'WU13235' if row.donor_original=='13_235' else row.donor_original)
        candidate_set={key for key,r in geo.items() if r['title'][0].split(',')[0].strip()==expected_alias}
        check(row.sample_id+'_same_donor_region_candidate_GSMs',set(row.candidate_geo_accessions.split(';'))==candidate_set and len(candidate_set)>=2 and obs.capture_GSM_mapping_status.eq('candidate_set_within_same_donor_and_region_not_unique').all())
        characteristics=[dict(v.split(': ',1) for v in geo[gsm]['characteristics_ch1']) for gsm in candidate_set]
        expected={'Control':'Donor','HCM':'HCM','DCM':'DCM','TOF':'TOF'}
        valid=(all(c['diagnosis']==expected[row.disease] for c in characteristics) if row.disease in expected else all('HLH' in c['diagnosis'] for c in characteristics))
        check(row.sample_id+'_diagnosis_matches_GEO',valid)
        geo_ages={c['age'] for c in characteristics}
        check(row.sample_id+'_both_source_ages_retained',len(geo_ages)==1 and obs.geo_age_original.astype(str).isin(geo_ages).all() and obs.age_original.astype(str).eq(row.author_age).all())
        flag=age_parts(row.author_age)!=age_parts(next(iter(geo_ages)))
        check(row.sample_id+'_true_bool_age_conflict',raw.obs.age_source_conflict.dtype.kind=='b' and obs.age_source_conflict.eq(flag).all())
        check(row.sample_id+'_numeric_sample_loss',len(obs)==row.source_nuclei and int(obs.pass_strict_numeric.sum())==row.strict_nuclei)
    conflicts=set(raw.obs.loc[raw.obs.age_source_conflict,'donor_original'].astype(str))
    check('UK_age_discrepancies_preserved_only_where_real',conflicts=={'UK1','UK2'})
    repair=json.loads((WORK/'METADATA_REPAIR_VERIFICATION.json').read_text())
    check('metadata_dtype_repair_did_not_change_counts',repair['counts_unchanged'] and repair['before_X_SHA256']==repair['after_X_SHA256'])
    print('Original metadata, clinical groups, donor/region and age conflicts verified',flush=True)
    mapping=pd.read_csv(WORK/'reference_gene_mapping.tsv',sep='\t',keep_default_na=False)
    reference=pd.read_csv(REF/'reference_genes.tsv',sep='\t',keep_default_na=False)
    source_unique=~genes.duplicated(keep=False)
    ref_unique=reference.loc[~reference.gene_symbol.duplicated(keep=False)&reference.gene_symbol.ne('')]
    expected_mapping=ref_unique.loc[ref_unique.gene_symbol.isin(genes.loc[source_unique])]
    check('unique_exact_reference_symbols_only',mapping.gene_id.tolist()==expected_mapping.gene_id.tolist() and mapping.gene_symbol.tolist()==expected_mapping.gene_symbol.tolist())
    positions=mapping.source_column.to_numpy(dtype=int)
    check('reference_mapping_preserves_true_source_columns',np.array_equal(genes.iloc[positions].to_numpy(),mapping.gene_symbol.to_numpy()) and np.array_equal(raw.var.iloc[positions].reference_gene_id.to_numpy(),mapping.gene_id.to_numpy()))
    check('all_unmapped_source_genes_retained',raw.n_vars==len(genes) and raw.var.reference_gene_id.ne('').sum()==len(mapping) and len(mapping)==18393)
    rmeta=pd.read_csv(REF/'reference_metadata.tsv',sep='\t');cmeta=pd.read_csv(REF/'calibration_metadata.tsv',sep='\t')
    check('fixed_reference_donors_disjoint',not rmeta.donor_id.isin(cmeta.donor_id).any())
    strict_cells=cells[strict]
    check('reference_query_cell_order',pd.Index((WORK/'strict_query_cells.txt').read_text().splitlines()).equals(strict_cells))
    with threadpool_limits(limits=4):query=mmread(WORK/'strict_reference_query.mtx').T.tocsr()
    check('all_query_raw_values_correspond_exactly',same(query,source[strict][:,positions]))
    del query
    a=raw[strict].copy()
    identity=pd.read_csv(WORK/'identity.tsv',sep='\t',index_col='cell_id',keep_default_na=False)
    check('all_strict_identity_rows_unique_and_aligned',identity.index.is_unique and identity.index.equals(strict_cells))
    a.obs=a.obs.join(identity)
    audit=pd.read_csv(OUT/'cell_identity_audit.tsv.gz',sep='\t',index_col='cell_id',keep_default_na=False)
    check('identity_audit_exact_all_strict_rows',audit.index.equals(strict_cells))
    for column in identity:
        if identity[column].dtype.kind in 'fiu':valid=np.allclose(identity[column],audit[column],rtol=0,atol=1e-12)
        else:valid=identity[column].astype(str).eq(audit[column].astype(str)).all()
        check('R_prediction_preserved_'+column,valid)
    policy=json.loads((OUT/'qc_policy.json').read_text())
    check('frozen_strict_policy',policy['min_genes']==500 and policy['max_genes']==5000 and policy['mt_percent_exclusive_max']==5)
    cal=pd.read_csv(WORK/'reference_calibration.tsv',sep='\t',keep_default_na=False)
    ref_report=json.loads((OUT/'REFERENCE_CALIBRATION.json').read_text())
    check('calibration_cells_match_fixed_donor_holdout',set(cal.cell_id)==set(cmeta.reference_barcode))
    tail=np.zeros(len(a),dtype=bool)
    for label,g in cal.groupby('truth'):
        correct=g.loc[g.reference_label.eq(label)]
        if len(correct)<30:continue
        delta=float(correct.reference_delta.quantile(.05));next_delta=float(correct.reference_delta_next.quantile(.05))
        check('reference_tail_'+label,np.isclose(delta,ref_report['tails'][label]['delta_q05']) and np.isclose(next_delta,ref_report['tails'][label]['delta_next_q05']))
        tail|=(identity.reference_label.eq(label)&identity.reference_pruned.eq(label)&identity.reference_delta.ge(delta)&identity.reference_delta_next.ge(next_delta)).to_numpy()
    check('reference_tail_flags_independent',np.array_equal(audit.reference_tail_pass,tail))
    programs={};n1={};n2={};bounds=policy['mixed_program_bounds'];strict_total=total[strict]
    with threadpool_limits(limits=4):ref_x=mmread(REF/'calibration_counts.mtx.gz').T.tocsr()
    ref_total=np.asarray(ref_x.sum(axis=1)).ravel();truth=cmeta.broad_type.replace({'Pericyte':'Mural','Smooth_muscle':'Mural'})
    for label,panel in policy['panels'].items():
        values=a.X[:,a.var.gene_symbol.isin(panel).to_numpy()].toarray()
        n1[label]=(values>0).sum(axis=1);n2[label]=(values>=2).sum(axis=1)
        programs[label]=np.log1p(values*(10000/strict_total)[:,None]).mean(axis=1)
        check('every_cell_marker_evidence_'+label,np.array_equal(audit[label+'_markers'],n1[label]) and np.array_equal(audit[label+'_markers_ge2_UMI'],n2[label]) and np.allclose(audit[label+'_program'],programs[label],rtol=0,atol=1e-10))
        if label in bounds:
            v=ref_x[:,reference.gene_symbol.isin(panel).to_numpy()].toarray()
            score=np.log1p(v*(10000/np.maximum(ref_total,1))[:,None]).mean(axis=1)
            threshold=max(float(np.quantile(score[truth.eq(label)],.1)),float(np.quantile(score[truth.ne(label)],.99)))
            check('fixed_reference_bound_'+label,np.isclose(threshold,bounds[label]['threshold']))
    del ref_x
    strong=np.zeros(len(a),dtype=bool);weak=strong.copy()
    for left,right in policy['mixed_program_pairs']:
        ratio=np.minimum(programs[left],programs[right])/np.maximum(np.maximum(programs[left],programs[right]),1e-12)
        weak|=(n1[left]>=3)&(n1[right]>=3)&(ratio>=.25)
        strong|=(n2[left]>=3)&(n2[right]>=3)&(programs[left]>=bounds[left]['threshold'])&(programs[right]>=bounds[right]['threshold'])&(ratio>=.25)
    check('mixed_program_rules_recomputed',np.array_equal(audit.coherent_mixed_panel_flag,strong) and np.array_equal(audit.weak_coherent_mixed_panel_flag,weak))
    dominate=programs['Cardiomyocyte']>=np.max(np.column_stack([p for k,p in programs.items() if k!='Cardiomyocyte']),axis=1)
    check('CM_program_dominance_recomputed',np.array_equal(audit.cm_program_dominates,dominate))
    review=pd.read_csv(WORK/'RESIDUAL_REVIEW_SUMMARY.tsv',sep='\t').set_index('sample_id')
    caller_class=np.full(len(a),'',dtype=object);caller_score=np.full(len(a),np.nan)
    for sample in a.obs.sample_id.unique():
        mask=a.obs.sample_id.eq(sample).to_numpy();model=(raw.obs.sample_id.eq(sample)&raw.obs.n_genes_by_counts.ge(500)).to_numpy()
        count=int(strong[mask].sum());trigger=count>=10 and count/mask.sum()>=.01
        row=review.loc[sample]
        check(sample+'_residual_trigger_and_input_scope',bool(row.residual_model_triggered)==trigger and row.model_nuclei==model.sum() and row.strict_nuclei==mask.sum() and row.strong_mixed_nuclei==count)
        directory=WORK/'captures'/sample
        df=pd.read_csv(directory/'doublets.tsv',sep='\t',index_col='cell_id')
        check(sample+'_caller_ids_unique',df.index.is_unique and a.obs_names[mask].isin(df.index).all())
        if trigger:
            model_cells=pd.Index((directory/'residual_model_cells.txt').read_text().splitlines())
            with threadpool_limits(limits=4):model_x=mmread(directory/'residual_model_counts.mtx').T.tocsr()
            check(sample+'_residual_model_includes_high_complexity_and_all_types',model_cells.equals(raw.obs_names[model]) and same(model_x,source[model]))
            del model_x
        else:
            check(sample+'_inherit_author_QC_no_fabricated_scores',df.index.equals(a.obs_names[mask]) and df.scDblFinder_score.isna().all() and df.scDblFinder_class.eq('not_rerun_author_mixed_cluster_QC').all())
        caller_class[mask]=df.loc[a.obs_names[mask],'scDblFinder_class'].to_numpy()
        caller_score[mask]=df.loc[a.obs_names[mask],'scDblFinder_score'].to_numpy()
    check('caller_results_preserved',np.array_equal(audit.scDblFinder_class,caller_class) and np.allclose(pd.to_numeric(audit.scDblFinder_score,errors='coerce'),caller_score,equal_nan=True))
    excluded=(caller_class=='doublet')&strong
    check('only_joint_supported_doublets_excluded',np.array_equal(audit.exclude_suspected_doublet,excluded))
    print('Reference, markers, mixed programs and all30 residual decisions verified',flush=True)
    diagnostic=ad.read_h5ad(OUT/'diagnostic_normalized.h5ad')
    check('diagnostic_only_same_strict_axis',diagnostic.obs_names.equals(strict_cells) and diagnostic.var_names.equals(raw.var_names) and bool(diagnostic.uns['analysis_only']) and 'counts' not in diagnostic.layers)
    for start in range(0,len(a),5000):
        stop=min(start+5000,len(a))
        expected=a.X[start:stop].astype(np.float32).copy()
        expected=expected.multiply((10000/np.asarray(expected.sum(axis=1)).ravel())[:,None]).tocsr();expected.data=np.log1p(expected.data)
        difference=(diagnostic.X[start:stop]-expected).tocsr()
        check('diagnostic_log_CP10k_block_'+str(start),not difference.nnz or np.max(np.abs(difference.data))<3e-5)
    graph=diagnostic.obsp['connectivities'].tocsr()
    check('graph_is_finite_nonnegative_and_aligned',graph.shape==(len(a),len(a)) and np.isfinite(graph.data).all() and (graph.data>=0).all())
    anchor=identity.reference_pruned.eq('Cardiomyocyte').to_numpy()&(n1['Cardiomyocyte']>=3)&dominate&~strong
    support=np.asarray(graph@anchor.astype(float)).ravel()/np.maximum(np.asarray(graph.sum(axis=1)).ravel(),1e-12)
    check('CM_neighbor_support_recomputed',np.allclose(audit.cm_anchor_neighbor_fraction,support,rtol=0,atol=1e-10))
    cm=(identity.reference_label.eq('Cardiomyocyte')&identity.reference_pruned.eq('Cardiomyocyte')).to_numpy().copy()
    cm&=(n1['Cardiomyocyte']>=2)&dominate&~strong&~excluded
    cm&=tail|((n1['Cardiomyocyte']>=3)&(support>=.5))
    cm&=a.obs.author_MainCellType.eq('CM').to_numpy()
    candidate=a.obs.author_MainCellType.eq('CM').to_numpy()|identity.reference_label.eq('Cardiomyocyte').to_numpy()|((n1['Cardiomyocyte']>=3)&dominate)
    pending=candidate&~cm&~excluded;binary=a.obs.binary_disease_eligible.to_numpy()
    check('all_CM_core_and_pending_decisions_recomputed',np.array_equal(audit.qc_eligible_CM,cm) and np.array_equal(audit.qc_eligible_binary_CM,cm&binary))
    actions=np.select([excluded,cm,pending],['exclude_suspected_doublet','retain_CM','pending_CM'],default='retain_other_or_unknown')
    check('every_strict_review_action_recomputed',np.array_equal(audit.review_action,actions))
    source_actions=np.full(len(raw),'outside_strict_numeric',dtype=object);source_actions[strict]=actions
    check('every_source_final_action_recomputed',np.array_equal(raw.obs.final_action,source_actions))
    full_audit=pd.read_csv(OUT/'all_nuclei_audit.tsv.gz',sep='\t',index_col='cell_id',keep_default_na=False)
    check('full_nuclei_audit_order_and_dispositions',full_audit.index.equals(cells) and np.array_equal(full_audit.final_action,source_actions) and np.array_equal(full_audit.pass_strict_numeric,strict))
    expected_ids={'all_source_counts.h5ad':raw.obs_names,'all_cells_qc.h5ad':a.obs_names[~excluded],
                  'cm_qc_full.h5ad':a.obs_names[cm&binary],'cm_other_diseases.h5ad':a.obs_names[cm&~binary],
                  'cm_review_pending.h5ad':a.obs_names[pending],'suspected_doublets.h5ad':a.obs_names[excluded],
                  'outside_strict_numeric.h5ad':raw.obs_names[~strict]}
    export=json.loads((OUT/'EXPORT_SUCCESS.json').read_text())
    check('nine_h5ad_outputs_exist',len(export['hashes'])==9 and set(export['hashes'])==set(expected_ids)|{'cm_qc_light.h5ad','diagnostic_normalized.h5ad'})
    for filename,entry in export['hashes'].items():
        path=OUT/filename
        check(filename+'_digest',path.stat().st_size==entry['bytes'] and digest(path)==entry['sha256'])
        if filename=='diagnostic_normalized.h5ad':continue
        obj=ad.read_h5ad(path)
        if filename in expected_ids:check(filename+'_exact_selected_nuclei',obj.obs_names.equals(expected_ids[filename]))
        check(filename+'_full_native_gene_axis',obj.var_names.equals(raw.var_names) and np.array_equal(obj.var.gene_symbol,genes))
        check(filename+'_all_values_unchanged',same(obj.X,raw[obj.obs_names].X) and same(obj.layers['counts'],obj.X))
        check(filename+'_disease_and_donor_labels',obj.obs_names.is_unique and np.array_equal(obj.obs.disease.astype(str),raw.obs.loc[obj.obs_names,'disease'].astype(str)) and np.array_equal(obj.obs.label,raw.obs.loc[obj.obs_names,'label']))
        check(filename+'_donor_capture_region_and_clinical_metadata',all(np.array_equal(obj.obs[column].astype(str),raw.obs.loc[obj.obs_names,column].astype(str)) for column in ['donor_id','donor_original','sample_id','source_patientID','tissue_region','clinical_group','age_original','geo_age_original','age_source_conflict','candidate_geo_accessions','split_status']))
        check(filename+'_native_and_derived_gene_metadata',obj.var.equals(raw.var))
        if filename not in ['all_source_counts.h5ad','outside_strict_numeric.h5ad']:
            positions_in_strict=a.obs_names.get_indexer(obj.obs_names)
            check(filename+'_recomputed_review_actions',np.array_equal(obj.obs.review_action,actions[positions_in_strict]) and np.array_equal(obj.obs.qc_eligible_CM,cm[positions_in_strict]) and np.array_equal(obj.obs.qc_eligible_binary_CM,(cm&binary)[positions_in_strict]))
        if filename in ['cm_qc_full.h5ad','cm_qc_light.h5ad']:
            check(filename+'_only_true_HCM_Control',obj.obs.label.isin([0,1]).all() and obj.obs.disease.isin(['HCM','Control']).all() and obj.obs.qc_eligible_binary_CM.all())
        if filename=='cm_other_diseases.h5ad':check('other_CM_not_used_as_binary_controls',obj.obs.label.eq(-1).all() and obj.obs.disease.isin(['DCM','HLHS','TOF']).all())
        del obj
    light=ad.read_h5ad(OUT/'cm_qc_light.h5ad');core=ad.read_h5ad(OUT/'cm_qc_full.h5ad')
    check('light_unique_max20k_subset',light.obs_names.is_unique and light.n_obs==min(20000,core.n_obs) and light.obs_names.isin(core.obs_names).all())
    selected=light.obs.groupby('donor_id',observed=True).size();available=core.obs.groupby('donor_id',observed=True).size()
    not_saturated=selected.loc[selected<available]
    check('light_donor_balancing_without_replicate_inflation',len(not_saturated)<2 or not_saturated.max()-not_saturated.min()<=1)
    coverage=pd.read_csv(OUT/'sample_qc_summary.tsv',sep='\t')
    check('all_thirty_samples_reported',set(coverage.sample_id)==set(samples.sample_id) and coverage.source_nuclei.sum()==len(raw) and coverage.strict_nuclei.sum()==len(a))
    check('author_CM_loss_accounted_separately',coverage.author_CM.sum()==raw.obs.author_MainCellType.eq('CM').sum() and coverage.author_CM_strict.sum()==a.obs.author_MainCellType.eq('CM').sum() and coverage.author_CM_core.sum()==cm.sum())
    for row in coverage.itertuples():
        mask=a.obs.sample_id.eq(row.sample_id).to_numpy();raw_mask=raw.obs.sample_id.eq(row.sample_id).to_numpy()
        check(row.sample_id+'_core_pending_loss_counts',row.CM_core_all_diagnoses==cm[mask].sum() and row.CM_core_binary==(cm&binary)[mask].sum() and row.pending_CM==pending[mask].sum() and row.new_suspected_doublets==excluded[mask].sum() and row.author_CM==raw.obs.loc[raw_mask,'author_MainCellType'].eq('CM').sum())
    check('training_still_unassigned_single_pediatric_HCM_donor',raw.obs.split_status.eq('unassigned_bind_donor_across_replicates_and_LV_RV').all() and raw.obs.training_use.str.contains('requires_pediatric_design').all())
    code_hashes={}
    for line in (WORK/'CODE_EXPECTED_SHA256SUMS.txt').read_text().splitlines():
        expected,filename=line.split(maxsplit=1);filename=filename.strip();actual=digest(WORK/filename)
        check('executed_code_SHA256_'+filename,actual.lower()==expected.lower());code_hashes[filename]=actual
    for name in ['stage','repair_metadata','R','finalize']:
        check('job_exit_'+name,(WORK/(name+'_exit_code.txt')).read_text().strip()=='0')
    completed=datetime.now(timezone.utc).isoformat()
    summary={'dataset':'GSE203274','status':'complete_verified','checked_utc':completed,'checks_passed':len(CHECKS),
             'source_nuclei':len(raw),'genes':raw.n_vars,'strict_nuclei':len(a),'source_donors':13,'source_samples':30,
             'core_binary_CM':core.n_obs,'core_HCM':int(core.obs.disease.eq('HCM').sum()),'core_Control':int(core.obs.disease.eq('Control').sum()),
             'core_binary_donors':int(core.obs.donor_id.nunique()),'light_CM':light.n_obs,
             'core_other_disease_CM':int((cm&~binary).sum()),'pending_CM':int(pending.sum()),'strong_mixed':int(strong.sum()),
             'new_suspected_doublets':int(excluded.sum()),'author_CM_source':int(raw.obs.author_MainCellType.eq('CM').sum()),
             'author_CM_strict':int(a.obs.author_MainCellType.eq('CM').sum()),'author_CM_core_all_diagnoses':int(cm.sum()),
             'author_CM_pending':int((pending&a.obs.author_MainCellType.eq('CM').to_numpy()).sum()),
             'age_conflicting_donors':sorted(conflicts),'all_checks':CHECKS,'verified_code_sha256':code_hashes,
             'suitability':'Pediatric cohort with one HCM donor; retain donor/region/age and all other diagnoses. QC success does not establish independent multi-donor HCM generalization or authorize a training split.'}
    (OUT/'INDEPENDENT_VERIFICATION.json').write_text(json.dumps(summary,indent=2)+'\n')
    disease_loss=coverage.groupby('disease',observed=True)[['source_nuclei','strict_nuclei','author_CM','author_CM_strict','author_CM_core','author_CM_pending','CM_core_binary','pending_CM']].sum()
    report=f'''# GSE203274：儿童心脏单核数据清洗与独立验收

{completed}；{len(CHECKS)}项验收通过。输出：`{OUT}/`。

原发布157293核、29266个源基因；严格数值候选{len(a)}核。HCM/Control心肌核心{core.n_obs}核（HCM {summary['core_HCM']}，Control {summary['core_Control']}），覆盖{summary['core_binary_donors']}位来源供体；轻量版{light.n_obs}核。其他疾病心肌核心{summary['core_other_disease_CM']}核单独保存、label=-1。待定{summary['pending_CM']}核，强混合程序{summary['strong_mixed']}核，新增联合疑似双细胞移除{summary['new_suspected_doublets']}核。

## 原件与既有处理

Hill et al., Nature 2022, PMID35732239/PMC10405779，GEO GSE203274。原下载两个文件的大小和MD5复查一致；归档中3个矩阵/轴文件与暂存版本按流式SHA256核对，所有157293×29266发布计数与各个交付矩阵逐值一致。rawCount是作者文件名，内容是CellBender v0.1.0去背景和作者QC后的未归一化整数RNA计数，不能称为未经修正的CellRanger分子计数或完整液滴。

论文方法：CellRanger3.0.1以hg19 pre-mRNA计数，CellBender去背景，700–10000基因、MT<0.2%，Harmony及全局/细胞类型内多轮混合谱系簇删除。发布前被删核和基因无法恢复；本轮不重跑CellBender、EmptyDrops或机械重复作者双细胞过滤。项目严格版保持500≤基因≤5000、MT<5%；所有指标从完整发布计数重算。作者percent.mt为百分数，与重算逐值一致，不能再乘100；既有0.2%阈值比项目MT规则更严格。基因上限5000是本轮主要新增限制，不将高基因数直接判为双细胞。

所有29266个原始符号行和顺序保留。仅18393个在源基因轴与现有成人心室参考中均唯一、完全同名的符号，获得派生reference_gene_id并用于身份适配；原计数hg19与参考GRCh38来源区分，不能把派生ID冒充作者原生Ensembl ID，也未按映射成功与否删除基因。未归一化计数X和counts均保持原值；诊断log1p(CP10k)对象单独保存，不进行疾病整合、插补、交付基因HVG裁切、差异表达、tokenization或模型训练。

## 疾病、重复和年龄

共30个orig.ident批次，13位labID供体：仅P86为HCM，年龄11y3m12d，两次捕获不等于两位HCM患者。4位供体对照；其余为DCM、TOF与新生儿/心衰HLHS。二分类只对HCM=1、Control=0赋值；其他疾病label=-1另存。大小写不同的diagnosis与Diagnosis都保留：小写列混有LVAD/HTx手术状态，按ClinicalRank/Diagnosis并核对GEO临床诊断定疾病，绝不把其他心脏病当健康。

13_198的LV/RV在作者patientID中被拆成两个字符串，本轮按同一labID供体绑定，保留部位与原patientID。重复捕获不跨将来训练/验证划分；P33的9k/10k等批次不猜配rep1/rep2，保存同供体同部位候选GSM集合并标记未唯一对应，避免编造运行级溯源。全队列年龄范围约2月至17岁余；原作者age文本、分解年/月/日、明确标识的近似年数及GEO年龄均保留，UK1、UK2两位作者/GEO日期不一致，真实布尔冲突标记保留，不自行改成某一个值。作者ageCont仅作为原注释。echoEF整列实际等于批次ID，原值保留并标记numeric_ejection_fraction_available=False，没有伪造射血分数。

## 身份与残余混合

SingleR使用固定供体分隔的成人LV参考和独立留出供体，不把参考分数当后验概率；成人→儿童尤其婴幼儿的域差异仍可能影响身份判定。核心须作者CM、SingleR最佳/剪枝CM、多标记与心肌程序优势、无强混合，并满足参考尾部或至少3个CM标记加邻域支持。参考、标记和邻域是相关证据，不宣称对真实儿童身份进行了独立准确率验证；冲突候选另存，不由一个工具或单基因阴性强删。

弱单UMI混合组合仅作诊断；强混合须每谱系≥3基因各≥2UMI、双程序达到固定参考max(本类q10,其他类q99)、程序比例≥0.25。按样本只有强混合≥10且≥1%严格候选才补测scDblFinder，模型须包含该发布批次全部≥500基因的类型与高基因数核；加载量未知时率参数仅启发式。实际没有批次触发补测，沿用作者已有混合簇过滤，新增移除0；不能由此宣称所有核无双细胞或排除了同型双细胞。强混合者不进入心肌核心，但保留在复核审计中。

## 严格规则与作者CM损失

```text
{disease_loss.to_string()}
```

作者CM {summary['author_CM_source']}→严格 {summary['author_CM_strict']}→全部疾病核心 {summary['author_CM_core_all_diagnoses']}；作者CM中的待定{summary['author_CM_pending']}。总待定还可能含作者非CM的身份冲突候选，不能把总待定当作心肌细胞损失。各样本/供体/作者簇的筛除和保留分别列出，必须评估HCM与供体/部位/年龄的选择偏差。原发布计数全部留在all_source_counts.h5ad，不满足数值规则的15992核另存outside_strict_numeric.h5ad。

## 交付与使用边界

all_cells_qc.h5ad保留全部源疾病的通过数值QC核；cm_qc_full.h5ad和cm_qc_light.h5ad仅HCM/Control核心，cm_other_diseases.h5ad保存其他疾病的心肌核心，cm_review_pending.h5ad保留身份待定，suspected_doublets.h5ad反映新增联合排除的真实空结果。空排除结果不代表真实双细胞率为0。诊断对象只用于QC，不作为原始计数交付使用。

本套仅1位HCM供体，不能用大量同供体核宣称多患者HCM泛化能力。儿童与成人队列、LV/RV/RVOT部位、临床阶段、年龄和跨研究供体重叠须在训练前复核；训练/验证划分尚未分配，不能把P86两个批次分到两侧。清洗验收完成不等于已经验证模型或批准直接混入成人微调。

INDEPENDENT_VERIFICATION.json保存原件/归档/全矩阵逐值证据、临床与供体关系、参考/规则/图支持/残余模型检查和源码SHA256；FINAL_COMPLETION_GATE.json仅在全部通过后写入。逐核表及矩阵留服务器，Git只存源码与文本汇总。
'''
    (OUT/'REPORT.md').write_text(report)
    gate={'dataset':'GSE203274','status':'complete_verified','verified_utc':completed,'checks_passed':len(CHECKS),
          'verification_sha256':digest(OUT/'INDEPENDENT_VERIFICATION.json'),'report_sha256':digest(OUT/'REPORT.md'),
          'output_directory':str(OUT),'code_sha256':code_hashes,'training_split_assigned':False,
          'limitation':'one pediatric HCM donor; repeats and LV/RV donor binding retained; other diagnoses never mapped to Control'}
    for directory in [WORK,OUT]:(directory/'FINAL_COMPLETION_GATE.json').write_text(json.dumps(gate,indent=2)+'\n')
    print(json.dumps({k:v for k,v in summary.items() if k not in ['all_checks','verified_code_sha256']},indent=2),flush=True)


if __name__=='__main__':
    try:run()
    except BaseException as exc:
        (WORK/'VERIFICATION_FAILURE.json').write_text(json.dumps({'status':'verification_failed','error':str(exc),'traceback':traceback.format_exc(),'checks':CHECKS},indent=2)+'\n')
        raise
