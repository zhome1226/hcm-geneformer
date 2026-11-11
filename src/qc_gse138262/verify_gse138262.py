"""Independently reread all corrected tables and verify every delivered count value."""
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
from scipy.io import mmread
from threadpoolctl import threadpool_limits

ROOT=Path.home()/'xiaoe-internship-data'
WORK=ROOT/'_processing/all_data_cleaning_20261111/gse138262'
OUT=ROOT/'processed/wehrens_GSE138262_qc_20261111'
REF=ROOT/'processed/SCP1303_review_20261111/identity'
CHECKS=[]


def check(name,passed,detail=None):
    item={'name':name,'passed':bool(passed)}
    if detail is not None:item['detail']=detail
    CHECKS.append(item)
    if not passed:raise AssertionError(name)


def digest(path,algorithm='sha256'):
    with Path(path).open('rb') as handle:
        return hashlib.file_digest(handle,algorithm).hexdigest()


def same(x,y):
    if x.shape!=y.shape:return False
    d=(x-y).tocsr();d.eliminate_zeros()
    return d.nnz==0


def original_table(path):
    # Independent parser: ordinary whitespace split and NumPy numeric parsing.
    with Path(path).open() as handle:
        cells=next(handle).split()
        features=[];values=[]
        for line in handle:
            feature,numbers=line.strip().split(maxsplit=1)
            row=np.fromstring(numbers,sep=' ',dtype=np.int64)
            assert len(row)==len(cells)
            features.append(feature);values.append(row)
    return pd.Index(cells),pd.Index(features),sparse.csr_matrix(np.vstack(values).T)


def r_h5(path):
    with h5py.File(path) as f:
        m=f['matrix']
        x=sparse.csc_matrix((m['data'][:],m['indices'][:],m['indptr'][:]),shape=m['shape'][:]).T.tocsr()
        return pd.Index(m['barcodes'].asstr()[:]),pd.Index(m['features/id'].asstr()[:]),x


def run():
    record=json.loads((WORK/'ZENODO_RECORD_18206553.json').read_text())
    download=json.loads((WORK/'DOWNLOAD_VERIFICATION.json').read_text())
    source=Path(download['source_directory'])
    check('author_corrected_record_version',record['id']==18206554 and record['metadata']['doi']=='10.5281/zenodo.18206554')
    check('five_real_download_files',len(record['files'])==5 and len(download['files'])==5)
    for item in record['files']:
        path=source/item['key']
        check('download_bytes_'+item['key'],path.stat().st_size==item['size'])
        check('download_MD5_'+item['key'],digest(path,'md5')==item['checksum'].removeprefix('md5:'))
    baseline=json.loads((WORK/'ORIGINAL_DOWNLOAD_BASELINE.json').read_text(encoding='utf-8-sig'))
    check('all_53_original_download_records_retained',len(baseline)==53)
    for item in baseline:
        path=Path(item['server_path'])
        check('original_unchanged_'+item['file'],path.stat().st_size==int(item['bytes']) and digest(path,'md5')==item['md5'])
    original=ad.read_h5ad(WORK/'corrected_native_union_counts.h5ad')
    raw=ad.read_h5ad(OUT/'all_source_counts.h5ad')
    check('union_identity_unique',raw.obs_names.is_unique and raw.var_names.is_unique and raw.var.gene_id.is_unique)
    check('full_union_axes_preserved',raw.obs_names.equals(original.obs_names) and raw.var_names.equals(original.var_names))
    check('full_union_counts_preserved',same(raw.X,original.X) and same(raw.layers['counts'],original.X))
    check('integer_finite_nonnegative_raw',np.isfinite(raw.X.data).all() and (raw.X.data>=0).all() and np.equal(raw.X.data,np.rint(raw.X.data)).all())
    samples=pd.read_csv(WORK/'samples.tsv',sep='\t')
    check('sample_and_donor_scope',len(samples)==32 and samples.donor_id.nunique()==12 and samples.loc[samples.disease.eq('HCM'),'donor_id'].nunique()==5 and samples.loc[samples.disease.eq('Control'),'donor_id'].nunique()==7)
    geo={r['accession']:r for r in json.loads((WORK/'GEO_SAMPLES_GSE138262.json').read_text(encoding='utf-8-sig'))}
    author=[]
    for kind in ['ALL','ROOIJonly']:
        cells,features,x=r_h5(WORK/f'author_{kind}/RNA_counts.h5')
        author.append((kind,cells,features,x))
    author_meta=pd.read_csv(WORK/'author_ALL/metadata.tsv.gz',sep='\t',index_col='barcode',keep_default_na=False)
    checked_cells=set();checked_features=set();author_comparisons=0;mito_comparisons=[]
    reference_features=pd.read_csv(REF/'reference_genes.tsv',sep='\t')
    for row in samples.itertuples():
        cells,features,x=original_table(row.source_path)
        check(row.sample_id+'_source_SHA256',digest(row.source_path)==row.source_sha256)
        check(row.sample_id+'_native_unique_axes',cells.is_unique and features.is_unique)
        source_cells=pd.Index([row.sample_id+':'+c for c in cells])
        obs=raw.obs.loc[source_cells]
        native=ad.read_h5ad(WORK/row.sample_id/'native_counts.h5ad')
        check(row.sample_id+'_native_axes',native.obs_names.equals(source_cells) and native.var_names.equals(features))
        check(row.sample_id+'_native_all_values',same(native.X,x))
        locations=raw.var_names.get_indexer(features)
        coo=x.tocoo()
        expanded=sparse.csr_matrix((coo.data,(coo.row,locations[coo.col])),shape=(len(cells),raw.n_vars))
        check(row.sample_id+'_union_every_value_and_missing_zeros',same(raw[source_cells].X,expanded))
        checked_cells.update(source_cells)
        checked_features.update(features)
        sample_pos=list(raw.uns['feature_presence_sample_order']).index(row.sample_id)
        check(row.sample_id+'_native_feature_presence_mask',np.array_equal(raw.varm['native_feature_presence_by_sample'][:,sample_pos],raw.var_names.isin(features)))
        split=[feature.split('_',1) for feature in features]
        gene_ids=pd.Index([s[0] for s in split]);symbols=pd.Index([s[1] for s in split])
        check(row.sample_id+'_authoritative_gene_mapping',np.array_equal(raw.var.loc[features,'gene_id'].to_numpy(),gene_ids) and np.array_equal(raw.var.loc[features,'gene_symbol'].fillna('').to_numpy(),symbols))
        check(row.sample_id+'_MT_genes_are_endogenous',symbols.str.startswith('MT-').sum()==13)
        total=np.asarray(x.sum(axis=1)).ravel();detected=np.asarray((x>0).sum(axis=1)).ravel()
        mt_counts=np.asarray(x[:,symbols.str.startswith('MT-')].sum(axis=1)).ravel()
        pct=np.divide(100.*mt_counts,total,out=np.full(len(total),np.nan),where=total>0)
        numeric=(detected>=500)&(detected<=5000)&(pct<5)
        check(row.sample_id+'_QC_metrics',np.array_equal(obs.total_counts.to_numpy(),total) and np.array_equal(obs.n_genes_by_counts.to_numpy(),detected) and np.allclose(obs.pct_counts_mt,pct,equal_nan=True,rtol=0,atol=1e-10))
        check(row.sample_id+'_strict_baseline',np.array_equal(obs.pass_strict_numeric.to_numpy(),numeric))
        author_cells=pd.Index([row.author_prefix+'_'+c for c in cells])
        check(row.sample_id+'_author_cell_relation',np.array_equal(obs.author_cell_id.to_numpy(),author_cells) and np.array_equal(obs.author_analysis_included.to_numpy(),author_cells.isin(author_meta.index)))
        check(row.sample_id+'_source_disease_metadata',obs.donor_id.eq(row.donor_id).all() and obs.disease.eq(row.disease).all() and obs.label.eq(1 if row.disease=='HCM' else 0).all() and obs.study_origin.eq(row.study_origin).all())
        if row.disease=='HCM':
            check(row.sample_id+'_official_diagnosis',f'patient id: {row.donor_original}' in geo[row.sample_id]['characteristics_ch1'] and 'patient diagnosis: hypertrophic cardiomyopathy' in geo[row.sample_id]['characteristics_ch1'])
        else:
            expected_donor=Path(row.source_path).name.split('.')[1]
            check(row.sample_id+'_Wang_donor_binding',row.donor_id=='Wang_'+expected_donor)
        mapped=pd.Index([f.replace('_',':') for f in features])
        for kind,a_cells,a_genes,a_x in author:
            if kind=='ROOIJonly' and row.disease=='Control':continue
            shared_cells=np.flatnonzero(author_cells.isin(a_cells));shared_genes=np.flatnonzero(mapped.isin(a_genes))
            check(row.sample_id+'_'+kind+'_has_shared_counts',len(shared_cells)>0 and len(shared_genes)>0)
            check(row.sample_id+'_'+kind+'_every_shared_RNA_value',same(x[shared_cells][:,shared_genes],a_x[a_cells.get_indexer(author_cells[shared_cells])][:,a_genes.get_indexer(mapped[shared_genes])]))
            author_comparisons+=1
        included=np.flatnonzero(author_cells.isin(author_meta.index))
        if len(included):
            declared=author_meta.loc[author_cells[included]]
            expected_paper='vRooij' if row.disease=='HCM' else 'Hu'
            expected_condition='HCM' if row.disease=='HCM' else 'Ctrl1'
            expected_patient=('R.P'+str(row.donor_original)) if row.disease=='HCM' else ('H.'+str(row.donor_original))
            check(row.sample_id+'_author_original_study_and_condition',declared.annotation_paper_str.eq(expected_paper).all() and declared.annotation_paper_beatified.eq(expected_condition).all() and declared.annotation_patient_str.eq(expected_patient).all())
            published=pd.to_numeric(author_meta.loc[author_cells[included],'nCount_rnaMitoCounts'],errors='raise').to_numpy()
            delta=np.abs(mt_counts[included]-published)
            mito_comparisons.append({'sample_id':row.sample_id,'cells':len(included),'max_absolute_count_difference':float(delta.max())})
            check(row.sample_id+'_author_mitochondrial_counts',np.array_equal(mt_counts[included],published))
    check('every_input_cell_read_directly_from_corrected_sources',checked_cells==set(raw.obs_names))
    check('all_48_author_object_plate_comparisons_repeated',author_comparisons==48)
    check('all_native_genes_retained',set(raw.var_names)==checked_features)
    rmeta=pd.read_csv(REF/'reference_metadata.tsv',sep='\t');cmeta=pd.read_csv(REF/'calibration_metadata.tsv',sep='\t')
    check('reference_training_and_holdout_donors_disjoint',not rmeta.donor_id.isin(cmeta.donor_id).any())
    common=(WORK/'reference_common_genes.txt').read_text().splitlines()
    expected_common=reference_features.gene_id[reference_features.gene_id.isin(raw.var.gene_id)].tolist()
    check('reference_adapted_to_authoritative_Ensembl_genes',common==expected_common)
    query_cells=pd.Index((WORK/'strict_query_cells.txt').read_text().splitlines())
    strict_cells=raw.obs_names[raw.obs.pass_strict_numeric.to_numpy()]
    check('identity_only_unchanged_strict_query',query_cells.equals(strict_cells))
    with threadpool_limits(limits=2):query=mmread(WORK/'strict_reference_query.mtx').T.tocsr()
    query_positions=pd.Index(raw.var.gene_id).get_indexer(common)
    check('SingleR_query_is_unnormalized_corrected_counts',same(query,raw[strict_cells].X[:,query_positions]))
    a=ad.read_h5ad(OUT/'all_cells_qc.h5ad')
    check('all_strict_cells_retained_without_loading_model_exclusion',a.obs_names.equals(strict_cells) and a.obs.doublet_assessment.eq('not_modeled_plate_sorted_no_10x_loading_assumption').all())
    check('strict_counts_exact',a.var_names.equals(raw.var_names) and same(a.X,raw[strict_cells].X) and same(a.layers['counts'],a.X))
    policy=json.loads((OUT/'qc_policy.json').read_text())
    check('frozen_project_numeric_rules',policy['min_genes']==500 and policy['max_genes']==5000 and policy['mt_percent_exclusive_max']==5)
    identity=pd.read_csv(WORK/'identity.tsv',sep='\t',index_col='cell_id',keep_default_na=False)
    check('reference_predictions_preserved',identity.index.equals(a.obs_names) and all(np.array_equal(a.obs[c].to_numpy(),identity[c].to_numpy()) if identity[c].dtype.kind not in 'f' else np.allclose(a.obs[c],identity[c],rtol=0,atol=1e-14) for c in identity.columns))
    cal=pd.read_csv(WORK/'reference_calibration.tsv',sep='\t',keep_default_na=False)
    ref_report=json.loads((OUT/'REFERENCE_CALIBRATION.json').read_text())
    check('calibration_uses_all_fixed_holdout_cells',set(cal.cell_id)==set(cmeta.reference_barcode))
    tail=np.zeros(a.n_obs,dtype=bool)
    for label,g in cal.groupby('truth'):
        correct=g.loc[g.reference_label.eq(label)]
        if len(correct)<30:continue
        delta=float(correct.reference_delta.quantile(.05));next_delta=float(correct.reference_delta_next.quantile(.05))
        check('reference_tail_'+label,np.isclose(delta,ref_report['tails'][label]['delta_q05']) and np.isclose(next_delta,ref_report['tails'][label]['delta_next_q05']))
        tail|=(a.obs.reference_label.eq(label)&a.obs.reference_pruned.eq(label)&a.obs.reference_delta.ge(delta)&a.obs.reference_delta_next.ge(next_delta)).to_numpy()
    check('reference_tail_flags_recomputed',np.array_equal(a.obs.reference_tail_pass,tail))
    with threadpool_limits(limits=2):ref_x=mmread(REF/'calibration_counts.mtx.gz').T.tocsr()
    ref_total=np.asarray(ref_x.sum(axis=1)).ravel();truth=cmeta.broad_type.replace({'Pericyte':'Mural','Smooth_muscle':'Mural'})
    bounds=policy['mixed_program_bounds'];programs={};n2={};n1={}
    total=np.asarray(a.X.sum(axis=1)).ravel()
    for label,panel in policy['panels'].items():
        values=a.X[:,a.var.gene_symbol.isin(panel).to_numpy()].toarray()
        programs[label]=np.log1p(values*(10000/np.maximum(total,1))[:,None]).mean(axis=1)
        n1[label]=(values>0).sum(axis=1);n2[label]=(values>=2).sum(axis=1)
        check('independent_marker_evidence_'+label,np.array_equal(a.obs[label+'_markers'],n1[label]) and np.array_equal(a.obs[label+'_markers_ge2_UMI'],n2[label]) and np.allclose(a.obs[label+'_program'],programs[label],rtol=0,atol=1e-10))
        if label in bounds:
            ref_values=ref_x[:,reference_features.gene_symbol.isin(panel).to_numpy()].toarray()
            score=np.log1p(ref_values*(10000/np.maximum(ref_total,1))[:,None]).mean(axis=1)
            threshold=max(float(np.quantile(score[truth.eq(label)],.1)),float(np.quantile(score[truth.ne(label)],.99)))
            check('fixed_reference_mixed_bound_'+label,np.isclose(threshold,bounds[label]['threshold']))
    mixed=np.zeros(a.n_obs,dtype=bool)
    for left,right in policy['mixed_program_pairs']:
        ratio=np.minimum(programs[left],programs[right])/np.maximum(np.maximum(programs[left],programs[right]),1e-12)
        mixed|=(n2[left]>=3)&(n2[right]>=3)&(programs[left]>=bounds[left]['threshold'])&(programs[right]>=bounds[right]['threshold'])&(ratio>=.25)
    check('strong_mixed_flags_recomputed',np.array_equal(a.obs.coherent_mixed_panel_flag,mixed))
    dominate=programs['Cardiomyocyte']>=np.max(np.column_stack([s for k,s in programs.items() if k!='Cardiomyocyte']),axis=1)
    check('CM_program_dominance_recomputed',np.array_equal(a.obs.cm_program_dominates,dominate))
    check('no_small_cohort_graph_rescue',not a.obs.neighbor_evidence_enabled.any() and a.obs.cm_anchor_neighbor_fraction.eq(0).all() and len(a)<50)
    expected_core=(a.obs.reference_label.eq('Cardiomyocyte')&a.obs.reference_pruned.eq('Cardiomyocyte')).to_numpy().copy()
    expected_core &= (n1['Cardiomyocyte']>=2)&dominate&~mixed&tail
    core=ad.read_h5ad(OUT/'cm_qc_full.h5ad');pending=ad.read_h5ad(OUT/'cm_review_pending.h5ad')
    check('CM_core_identity_recomputed',set(core.obs_names)==set(a.obs_names[expected_core]) and np.array_equal(a.obs.qc_eligible_CM,expected_core))
    check('strict_candidates_partition_without_loss',set(core.obs_names).isdisjoint(pending.obs_names) and set(core.obs_names)|set(pending.obs_names)==set(a.obs_names))
    files=json.loads((OUT/'EXPORT_SUCCESS.json').read_text())['hashes']
    for filename,entry in files.items():
        path=OUT/filename
        check(filename+'_file_digest',path.stat().st_size==entry['bytes'] and digest(path)==entry['sha256'])
        obj=ad.read_h5ad(path)
        if filename.startswith('diagnostic_'):
            expected=raw[obj.obs_names].X.astype(np.float32).copy()
            expected=expected.multiply((10000/np.asarray(expected.sum(axis=1)).ravel())[:,None]).tocsr()
            expected.data=np.log1p(expected.data)
            delta=(obj.X-expected).tocsr()
            check(filename+'_explicit_log_CP10k_only',bool(obj.uns['analysis_only']) and 'counts' not in obj.layers and (not delta.nnz or np.max(np.abs(delta.data))<3e-5))
        else:
            check(filename+'_entire_count_matrix_preserved',obj.var_names.equals(raw.var_names) and same(obj.X,raw[obj.obs_names].X) and same(obj.layers['counts'],obj.X))
            check(filename+'_unique_cells_and_disease_labels',obj.obs_names.is_unique and obj.obs.label.eq(obj.obs.disease.map({'Control':0,'HCM':1})).all())
    outside=ad.read_h5ad(OUT/'outside_strict_numeric.h5ad')
    check('outside_strict_cells_retained_exactly',set(outside.obs_names)==set(raw.obs_names[~raw.obs.pass_strict_numeric.to_numpy()]))
    light=ad.read_h5ad(OUT/'cm_qc_light.h5ad')
    check('light_unique_and_not_inflated',light.obs_names.equals(core.obs_names) and light.n_obs<=20000)
    check('no_training_approval_from_severely_selected_data',raw.obs.training_use.str.startswith('not_approved_').all() and a.obs.training_use.str.startswith('not_approved_').all() and core.obs.disease.eq('HCM').all())
    coverage=pd.read_csv(OUT/'sample_qc_summary.tsv',sep='\t')
    check('all_32_samples_including_zero_survivors_reported',set(coverage.sample_id)==set(samples.sample_id) and coverage.source_cells.sum()==raw.n_obs and coverage.strict_cells.sum()==a.n_obs and coverage.CM_core.sum()==core.n_obs)
    external=pd.read_csv(WORK/'author_external_or_unmatched_cells.tsv.gz',sep='\t',index_col='barcode',keep_default_na=False)
    check('author_external_cells_inventoried_not_reused_as_independent_controls',set(external.index)==set(author_meta.index)-set(raw.obs.author_cell_id) and not raw.obs.study_origin.str.contains('Litvinukova|Teichmann').any())
    source_hashes={}
    for line in (WORK/'CODE_EXPECTED_SHA256SUMS.txt').read_text().splitlines():
        expected,filename=line.split(maxsplit=1);filename=filename.strip()
        actual=digest(WORK/filename)
        check('executed_source_SHA256_'+filename,actual.lower()==expected.lower())
        source_hashes[filename]=actual
    for name in ['download','inspect','stage','identity','finalize']:
        check('job_exit_'+name,(WORK/(name+'_exit_code.txt')).read_text().strip()=='0')
    completed=datetime.now(timezone.utc).isoformat()
    result={'dataset':'GSE138262','status':'complete_verified','checked_utc':completed,'checks_passed':len(CHECKS),
            'source_cells':raw.n_obs,'genes':raw.n_vars,'strict_cells':a.n_obs,'core_CM':core.n_obs,'core_HCM':int(core.obs.disease.eq('HCM').sum()),
            'core_Control':int(core.obs.disease.eq('Control').sum()),'core_donors':int(core.obs.donor_id.nunique()),'pending_identity':pending.n_obs,
            'author_count_comparisons':author_comparisons,'source_files_unchanged':len(baseline),'author_MT_count_comparisons':mito_comparisons,
            'training_use':'NOT APPROVED: severe strict whole-cell attrition; no Control; protocol/design and overlap review required',
            'all_checks':CHECKS,'verified_code_sha256':source_hashes}
    (OUT/'INDEPENDENT_VERIFICATION.json').write_text(json.dumps(result,indent=2)+'\n')
    by_disease=coverage.groupby('disease')[['source_cells','genes_500to5000','MT_lt5pct','strict_cells','CM_core','pending_identity']].sum()
    report=f'''# GSE138262：修正版计数清洗与独立验收

{completed}；{len(CHECKS)}项检查通过。输出：`{OUT}/`。

已清洗验收，不等于适合训练。严格基线仅留下{a.n_obs}个数值候选、{core.n_obs}个心肌核心、{pending.n_obs}个身份待定；核心均为HCM，只有{core.obs.donor_id.nunique()}位供体，无健康对照。本套不批准进入HCM/Control微调或独立疾病验证。没有为凑数量放宽规则、复制细胞或混入参考图谱。

## 输入更正与计数来源

作者仓库 https://github.com/vanrooij-lab/scRNAseq-HCM-human 明确更正：GEO旧TSV不是Wehrens et al., Cell Reports 2022实际使用的计数表。按作者链接取得概念DOI 10.5281/zenodo.18206553指向的版本DOI 10.5281/zenodo.18206554。5个文件全部符合Zenodo大小和MD5；32个完整修正UMI计数表保留。原下载清单的53个文件大小/MD5再次一致，旧GEO数据留作来源存档，没有覆盖。

16个CEL-seq2分选完整心肌细胞板来自5位HCM切肌患者，GSM及患者映射核对GEO和作者文件名；另外16个iCell8板/部分来自7位Wang健康供体，是作者按LV细胞群选取后重处理的对照伴随数据，保留其独立研究/技术来源。N2的多板、各板part1/2均绑定原供体；不把板当独立患者。Seurat ALL中另外的Litvinukova细胞单独列入来源审计，不重复加入当前队列，避免与现有参考图谱重用。

每个修正计数表均独立解析并与native、并集、最终各对象逐值核对；48组共享板/作者Seurat RNA计数全部一致。RNA assay已经删去线粒体基因；作者保存的线粒体计数元数据也与完整修正计数重算一致。MT百分比来自完整表，绝非从裁掉MT后的RNA层算出0%。原始命名、Ensembl ID及空gene symbol全部保留；并集{raw.n_vars}个基因，native轴与逐样本存在掩码保留，不以共有基因替代完整计数。

## 固定严格规则与真实损失

500≤检出基因≤5000、MT<5%是既定项目基线，不是所有HCM完整细胞的统一质量定义。逐疾病统计如下：

```text
{by_disease.to_string()}
```

完整心肌细胞的线粒体RNA占比可远高于单核数据；本套强烈筛除，不能据此称其余细胞均为坏细胞。完整7537个来源细胞及所有计数保存在all_source_counts.h5ad；不满足严格数值规则者单独保存outside_strict_numeric.h5ad，逐细胞及32个样本/12位供体损失有表可查。未来若制定完整细胞专用阈值，应另发版本并与本严格版比较，本轮未改阈值。

## 身份、混合信号与适用边界

SingleR使用27664个共同Ensembl基因和固定供体分隔的成人左心室单核参考；参考域与本完整细胞并不完全一致。联合多基因程序、程序优势及参考分数尾部判断心肌身份；评分不是后验概率。只有8个严格候选，不使用小样本邻域去挽救身份。2个身份未通过心肌核心门槛的细胞保留待定，不强制冒充心肌。

本数据为板式分选CEL-seq2和已选iCell8细胞，不运行EmptyDrops，不套用10x装载量双细胞率；高基因数不等于双细胞。强跨谱系程序按固定参考阈值、多基因/多UMI支持复核，结果仅为可疑混合信号，不能排除同型双细胞或证明单细胞纯度。未经校准的单一工具标签不用于强删。诊断log1p(CP10k)副本有明确标识，计数对象X与counts逐值不变，未插补、整合、按HVG裁切或开展差异表达/训练。

此套HCM与健康伴随来源具有研究、技术、组织取材及选择规则差异；跨研究供体重叠仍未解决，任何后续划分须按供体绑定。本严格版缺乏可用对照，6个核心不足以代表原5位HCM供体，更不支持疾病建模结论。

## 产物与验收

all_source_counts.h5ad：全部来源计数及审计；all_cells_qc.h5ad：8个数值候选；cm_qc_full.h5ad与cm_qc_light.h5ad：6个真实心肌核心；cm_review_pending.h5ad：2个待定。mixed_lineage_review.h5ad为真实规则过滤结果，本轮0个，不是伪造候选。outside_strict_numeric.h5ad保留7529个数值规则外的细胞。所有计数对象保留完整并集基因轴与原值，诊断对象单独标识。

INDEPENDENT_VERIFICATION.json记录逐值验收范围、原件与源码摘要、供体分隔、数值/身份/混合规则和小样本限制；FINAL_COMPLETION_GATE.json只在全部检查通过后写入。参考校准表、sample_qc_summary.tsv、donor_qc_summary.tsv、sample_metadata.tsv、all_cell_audit.tsv.gz与cell_identity_audit.tsv.gz可追溯。只将代码与文本汇总放Git，矩阵及逐细胞表留服务器。
'''
    (OUT/'REPORT.md').write_text(report)
    gate={'dataset':'GSE138262','status':'complete_verified','verified_utc':completed,'verification_sha256':digest(OUT/'INDEPENDENT_VERIFICATION.json'),
          'report_sha256':digest(OUT/'REPORT.md'),'checks_passed':len(CHECKS),'training_approved':False,
          'training_limitation':result['training_use'],'output_directory':str(OUT),'code_sha256':source_hashes}
    (WORK/'FINAL_COMPLETION_GATE.json').write_text(json.dumps(gate,indent=2)+'\n')
    (OUT/'FINAL_COMPLETION_GATE.json').write_text(json.dumps(gate,indent=2)+'\n')
    print(json.dumps({k:v for k,v in result.items() if k not in ['all_checks','verified_code_sha256','author_MT_count_comparisons']},indent=2),flush=True)


if __name__=='__main__':
    try:
        run()
    except BaseException as exc:
        result={'dataset':'GSE138262','status':'verification_failed','checked_utc':datetime.now(timezone.utc).isoformat(),
                'checks':CHECKS,'error':str(exc),'traceback':traceback.format_exc()}
        (WORK/'VERIFICATION_FAILURE.json').write_text(json.dumps(result,indent=2)+'\n')
        raise
