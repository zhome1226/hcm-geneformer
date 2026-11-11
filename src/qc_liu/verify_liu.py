"""Reopen final Liu objects and verify every delivered count against released MTX."""
from pathlib import Path
from datetime import datetime, timezone
import hashlib
import json

import anndata as ad
import numpy as np
import pandas as pd
from scipy.io import mmread
from threadpoolctl import threadpool_limits

ROOT = Path.home() / 'xiaoe-internship-data'
WORK = ROOT / '_processing/all_data_cleaning_20261111/liu2023'
OUT = ROOT / 'processed/liu2023_qc_20261111'
checks = {}


def check(name, value):
    checks[name] = bool(value)


def equal(a,b):
    return a.shape == b.shape and (a != b).nnz == 0


def file_hashes(path):
    md5, sha = hashlib.md5(), hashlib.sha256()
    with path.open('rb') as stream:
        while chunk := stream.read(32 * 1024**2):
            md5.update(chunk); sha.update(chunk)
    return {'bytes':path.stat().st_size,'md5':md5.hexdigest(),'sha256':sha.hexdigest()}


samples = pd.read_csv(WORK/'samples.tsv',sep='\t')
meta = pd.read_csv(WORK/'author_metadata.tsv.gz',sep='\t',index_col='barcode',keep_default_na=False)
preflight = json.loads((WORK/'PREFLIGHT.json').read_text())
policy = json.loads((OUT/'qc_policy.json').read_text())
export = json.loads((OUT/'EXPORT_SUCCESS.json').read_text())
names = ['all_cells_qc.h5ad','cm_qc_full.h5ad','cm_qc_light.h5ad','cm_review_pending.h5ad','suspected_doublets.h5ad']
objects = {name:ad.read_h5ad(OUT/name) for name in names}
allcell, core, light, pending, dbl = [objects[name] for name in names]
identity = pd.read_csv(OUT/'cell_identity_audit.tsv.gz',sep='\t',index_col='barcode',keep_default_na=False)
combined = ad.concat([allcell,dbl],join='inner',merge='same')[identity.index].copy()
obs = combined.obs
check('policy_strict_limits',policy['min_genes']==500 and policy['max_genes']==5000 and policy['mt_percent_exclusive_max']==5)
check('no_emptydrops_claim', 'No EmptyDrops' in policy['cell_calling'] and not any('emptydrops' in c for c in obs.columns))
check('strict_audit_partition',set(identity.index)==set(allcell.obs_names)|set(dbl.obs_names) and not dbl.obs_names.isin(allcell.obs_names).any())
check('audit_actions_match_objects',np.array_equal(obs.review_action.astype(str),identity.review_action))
check('core_matches_audit',set(core.obs_names)==set(identity.index[identity.review_action.eq('retain_CM')]))
check('pending_matches_audit',set(pending.obs_names)==set(identity.index[identity.review_action.eq('pending_CM')]))
check('light_unique_subset_target',light.obs_names.is_unique and light.obs_names.isin(core.obs_names).all() and light.n_obs==min(20000,core.n_obs))
check('input_sample_donor_count',samples.sample_id.nunique()==13 and samples.donor_id.nunique()==12)
check('author_groups_nuclei',len(meta)==55122 and meta.Group.eq('HCM').sum()==39183 and meta.Group.eq('HEALTHY').sum()==15939)
check('shared_control_donor',meta.loc[meta['orig.ident'].isin(['CL0412A','CL0412B']),'Subject'].eq('CL0412').all())
check('no_split_assigned',all('split' not in x.obs for x in objects.values()))
check('study_separate',all(x.obs.dataset.eq('liu2023').all() for x in objects.values()))
for name,x in objects.items():
    check(name+'_unique_axes',x.obs_names.is_unique and x.var_names.is_unique)
    check(name+'_strict_bounds',x.obs.n_genes_by_counts.between(500,5000).all() and x.obs.pct_counts_mt.lt(5).all())
    check(name+'_integer_nonnegative',np.isfinite(x.X.data).all() and (x.X.data>=0).all() and (x.X.data==np.round(x.X.data)).all())
    check(name+'_counts_layer_unchanged',equal(x.X,x.layers['counts']))
    check(name+'_all_released_genes',x.n_vars==31944 and x.var_names.equals(combined.var_names))
    check(name+'_author_prefilter_provenance',x.obs.source_pre_filtered.all())
for name in ['cm_qc_full.h5ad','cm_qc_light.h5ad','cm_review_pending.h5ad']:
    x=objects[name]
    check(name+'_subset_values',equal(x.X,allcell[x.obs_names].X))
check('allcell_no_exclusions',not allcell.obs.exclude_suspected_doublet.any())
check('doublet_two_evidence',dbl.obs.scDblFinder_class.eq('doublet').all() and dbl.obs.coherent_mixed_panel_flag.all())

# Independently recompute marker programs/mixed flags from delivered integer counts.
panels = {
 'Cardiomyocyte':['TNNT2','TTN','TNNI3','ACTC1','MYH7','MYBPC3','MYL2','TNNC1','RYR2','PLN','ACTN2','MYH6'],
 'Fibroblast':['COL1A1','COL1A2','DCN','LUM','PDGFRA'],
 'Endothelial':['PECAM1','VWF','KDR','EMCN','CLDN5'],
 'Immune':['PTPRC','LST1','TYROBP','CD3D','CD79A'],
 'Mural':['RGS5','PDGFRB','CSPG4','MYH11','ACTA2','TAGLN'],
 'Neural':['NRXN1','SNAP25','SYT1','SLC1A3','PLP1'],
 'Adipocyte':['ADIPOQ','FABP4','PLIN1','LPL','CIDEC'],
 'Mesothelial':['WT1','MSLN','UPK3B','KRT8','KRT18']}
pairs=[('Cardiomyocyte','Fibroblast'),('Cardiomyocyte','Endothelial'),('Cardiomyocyte','Immune'),('Cardiomyocyte','Mural'),('Fibroblast','Endothelial'),('Fibroblast','Immune'),('Endothelial','Immune'),('Immune','Mural')]
totals=np.asarray(combined.X.sum(axis=1)).ravel()
n, n2, programs={},{},{}
for label,genes in panels.items():
    values=combined.X[:,combined.var.gene_symbol.isin(genes).to_numpy()].toarray()
    n[label]=(values>0).sum(axis=1)
    n2[label]=(values>=2).sum(axis=1)
    programs[label]=np.log1p(values*(10000/np.maximum(totals,1))[:,None]).mean(axis=1)
    check(label+'_marker_counts',np.array_equal(n[label],obs[label+'_markers']))
    check(label+'_marker_ge2_counts',np.array_equal(n2[label],obs[label+'_markers_ge2_UMI']))
    check(label+'_program_values',np.allclose(programs[label],obs[label+'_program']))
dominates=programs['Cardiomyocyte']>=np.max(np.vstack([v for k,v in programs.items() if k!='Cardiomyocyte']),axis=0)
refdir=ROOT/'processed/SCP1303_review_20261111/identity'
mixed_genes=pd.read_csv(refdir/'reference_genes.tsv',sep='\t')
mixed_meta=pd.read_csv(refdir/'calibration_metadata.tsv',sep='\t')
with threadpool_limits(limits=2):
    mixed_ref=mmread(refdir/'calibration_counts.mtx.gz').T.tocsr()
mixed_size=np.asarray(mixed_ref.sum(axis=1)).ravel()
mixed_truth=mixed_meta.broad_type.replace({'Pericyte':'Mural','Smooth_muscle':'Mural'})
bounds={}
for label in ['Cardiomyocyte','Fibroblast','Endothelial','Immune','Mural']:
    values=mixed_ref[:,mixed_genes.gene_symbol.isin(panels[label]).to_numpy()].toarray()
    score=np.log1p(values*(10000/np.maximum(mixed_size,1))[:,None]).mean(axis=1)
    bounds[label]=max(np.quantile(score[mixed_truth.eq(label)],.1),np.quantile(score[mixed_truth.ne(label)],.99))
    check(label+'_mixed_threshold_from_reference',np.isclose(bounds[label],policy['mixed_program_bounds'][label]['threshold']))
mixed=np.zeros(combined.n_obs,dtype=bool)
weak=mixed.copy()
for left,right in pairs:
    ratio=np.minimum(programs[left],programs[right])/np.maximum(np.maximum(programs[left],programs[right]),1e-12)
    valid=(n[left]>=3)&(n[right]>=3)&(ratio>=.25)
    weak|=valid
    mixed|=valid&(n2[left]>=3)&(n2[right]>=3)&(programs[left]>=bounds[left])&(programs[right]>=bounds[right])
check('mixed_flags_recomputed',np.array_equal(mixed,obs.coherent_mixed_panel_flag))
check('weak_flags_diagnostic_only',np.array_equal(weak,obs.weak_coherent_mixed_panel_flag))
check('cm_dominance_recomputed',np.array_equal(dominates,obs.cm_program_dominates))
exclude=obs.scDblFinder_class.eq('doublet').to_numpy()&mixed
check('doublet_exclusion_membership',np.array_equal(exclude,obs.exclude_suspected_doublet))
expected_core=(obs.author_cellType.eq('CM')&obs.reference_label.eq('Cardiomyocyte')&obs.reference_pruned.eq('Cardiomyocyte')&
               (n['Cardiomyocyte']>=2)&dominates&~mixed&
               (obs.reference_tail_pass|((n['Cardiomyocyte']>=3)&obs.cm_anchor_neighbor_fraction.ge(.5))))
check('core_rule_recomputed',set(obs.index[expected_core])==set(core.obs_names))
candidate=obs.author_cellType.eq('CM')|obs.reference_label.eq('Cardiomyocyte')|((n['Cardiomyocyte']>=3)&dominates)
check('pending_rule_recomputed',set(obs.index[candidate&~expected_core&~exclude])==set(pending.obs_names))

# Verify adapted reference and original donor-separated training/holdout inputs.
refdir=ROOT/'processed/SCP1303_review_20261111/identity'
rmeta=pd.read_csv(refdir/'reference_metadata.tsv',sep='\t')
cmeta=pd.read_csv(refdir/'calibration_metadata.tsv',sep='\t')
check('reference_holdout_donors_disjoint',not rmeta.donor_id.isin(cmeta.donor_id).any())
refgenes=pd.read_csv(refdir/'reference_genes.tsv',sep='\t')
mapping=pd.read_csv(WORK/'reference_gene_mapping.tsv',sep='\t')
expected_map=refgenes[~refgenes.gene_symbol.duplicated(keep=False)&refgenes.gene_symbol.isin(combined.var_names)]
check('reference_unique_symbol_mapping',mapping.reset_index(drop=True).equals(expected_map.reset_index(drop=True)))
cal=pd.read_csv(WORK/'reference_calibration.tsv',sep='\t',keep_default_na=False)
check('calibration_truth_matches_heldout',cal.barcode.is_unique and set(cal.barcode)==set(cmeta.reference_barcode) and np.array_equal(cal.truth,cmeta.set_index('reference_barcode').loc[cal.barcode,'broad_type']))
tails={}
for label,frame in cal.groupby('truth'):
    correct=frame[frame.reference_label.eq(label)]
    if len(correct)>=30:
        tails[label]=(correct.reference_delta.quantile(.05),correct.reference_delta_next.quantile(.05))
tail=np.zeros(len(obs),dtype=bool)
for label,(delta,nxt) in tails.items():
    tail|=(obs.reference_label.eq(label)&obs.reference_pruned.eq(label)&obs.reference_delta.ge(delta)&obs.reference_delta_next.ge(nxt)).to_numpy()
check('reference_tail_flags_recomputed',np.array_equal(tail,obs.reference_tail_pass))

strict_expected=set()
source_summaries=[]
for row in samples.itertuples():
    sample=row.sample_id; directory=Path(row.raw_directory)
    var=pd.read_csv(directory/'features.tsv.gz',sep='\t',header=None,dtype=str,keep_default_na=False)
    barcodes=pd.read_csv(directory/'barcodes.tsv.gz',header=None,dtype=str)[0]
    with threadpool_limits(limits=2):
        raw=mmread(directory/'matrix.mtx.gz').T.tocsr()
    raw.sum_duplicates();raw.eliminate_zeros();raw.sort_indices()
    count=np.asarray(raw.sum(axis=1)).ravel();genes=np.diff(raw.indptr)
    mt=100*np.asarray(raw[:,var[1].str.startswith('MT-').to_numpy()].sum(axis=1)).ravel()/np.maximum(count,1)
    strict=(genes>=500)&(genes<=5000)&(mt<5)
    ids=pd.Index(sample+':'+barcodes)
    strict_expected.update(ids[strict])
    am=meta.loc[barcodes]
    check(sample+'_metadata_complete_source',am['orig.ident'].eq(sample).all() and am.Subject.nunique()==1)
    disease='Control' if am.Group.eq('HEALTHY').all() else 'HCM'
    label=0 if disease=='Control' else 1
    donor='Liu2023_'+am.Subject.iloc[0]
    check(sample+'_sample_table_matches_author',row.donor_id==donor and row.disease==disease)
    check(sample+'_author_metrics_and_MT_fraction',np.array_equal(count,am.nCount_RNA) and np.array_equal(genes,am.nFeature_RNA) and np.allclose(mt,100*am['percent.mito']))
    stage=json.loads((WORK/sample/'STAGED.json').read_text())
    for name,record in stage['sources'].items():
        check(sample+'_'+name+'_unchanged',file_hashes(Path(record['path']))['sha256']==record['sha256'])
    model=pd.read_csv(WORK/sample/'doublets.tsv',sep='\t',keep_default_na=False)
    check(sample+'_full_multitype_model_input_ge500',set(model.barcode)==set(ids[genes>=500]) and model.barcode.is_unique)
    trigger=str(row.trigger_residual_doublet_screen).lower()=='true'
    strict_obs=obs.loc[obs.sample_id.eq(sample)]
    nm=int(strict_obs.weak_coherent_mixed_panel_flag.sum())
    check(sample+'_residual_screen_trigger',trigger==(nm>=10 and nm/max(1,int(strict.sum()))>=.01) and nm==row.residual_mixed_strict)
    if trigger:
        check(sample+'_finite_model_scores',np.isfinite(pd.to_numeric(model.scDblFinder_score)).all() and model.scDblFinder_class.isin(['singlet','doublet']).all())
    else:
        check(sample+'_no_mechanical_doublet_rerun',model.scDblFinder_class.eq('not_rerun_author_QC').all() and model.scDblFinder_score.eq('').all())
    prediction=pd.read_csv(WORK/sample/'identity.tsv',sep='\t',index_col='barcode',keep_default_na=False)
    check(sample+'_identity_complete_strict_input',set(prediction.index)==set(ids[strict]))
    for column in ['reference_label','reference_pruned']:
        check(sample+'_'+column+'_traceable',np.array_equal(strict_obs[column].astype(str),prediction.loc[strict_obs.index,column]))
    audit=pd.read_csv(OUT/f'{sample}_released_nuclei_audit.tsv.gz',sep='\t',index_col='barcode')
    check(sample+'_all_released_nuclei_audited',audit.index.equals(ids) and np.array_equal(audit.n_genes_by_counts,genes) and np.array_equal(audit.total_counts,count) and np.allclose(audit.pct_counts_mt,mt) and np.array_equal(audit.pass_strict_numeric,strict))
    check(sample+'_audit_non_strict_quarantined',audit.loc[~strict,'final_action'].eq('outside_strict_numeric').all())
    for name,obj in objects.items():
        part=obj[obj.obs.sample_id.eq(sample)]
        positions=pd.Index(barcodes).get_indexer(part.obs.barcode_original.astype(str))
        stem=sample+'_'+name
        check(stem+'_original_gene_axis',np.array_equal(part.var_names.astype(str),var[0]) and np.array_equal(part.var.gene_symbol.astype(str),var[1]))
        check(stem+'_original_barcodes',(positions>=0).all() and np.array_equal(part.obs_names.astype(str),sample+':'+part.obs.barcode_original.astype(str)))
        check(stem+'_every_original_count',equal(part.X,raw[positions]))
        check(stem+'_recomputed_QC',np.array_equal(part.obs.total_counts,count[positions]) and np.array_equal(part.obs.n_genes_by_counts,genes[positions]) and np.allclose(part.obs.pct_counts_mt,mt[positions]))
        check(stem+'_donor_disease_source',part.obs.donor_id.eq(donor).all() and part.obs.disease.eq(disease).all() and part.obs.label.eq(label).all() and part.obs.capture_id.eq(sample).all())
        pm=am.loc[part.obs.barcode_original.astype(str)]
        check(stem+'_author_identity_preserved',np.array_equal(part.obs.author_cellType.astype(str),pm.cellType) and np.array_equal(part.obs.author_subcluster.astype(str),pm.subcluster))
    source_summaries.append({'sample_id':sample,'released_nuclei':len(barcodes),'strict_nuclei':int(strict.sum()),'donor_id':donor})
    print('Verified original counts and metadata: '+sample,flush=True)
check('every_strict_original_nucleus_accounted',strict_expected==set(combined.obs_names))
check('expected_released_input_count',sum(x['released_nuclei'] for x in source_summaries)==55122)
for name,record in preflight['originals'].items():
    hashes=file_hashes(ROOT/'validation/liu2023'/name)
    check('source_'+name+'_unchanged',hashes==record)

diagnostic=ad.read_h5ad(OUT/'diagnostic_normalized.h5ad')
check('diagnostic_explicit_analysis_copy',diagnostic.uns.get('analysis_only',False) and 'Diagnostic log1p CP10k' in diagnostic.uns.get('count_source',''))
check('diagnostic_full_strict_axes',set(diagnostic.obs_names)==set(combined.obs_names) and diagnostic.var_names.equals(combined.var_names))
rx=combined[diagnostic.obs_names].X.astype(float)
size=np.asarray(rx.sum(axis=1)).ravel()
expected=rx.multiply((10000/np.maximum(size,1))[:,None]).tocsr();expected.data=np.log1p(expected.data)
diff=expected-diagnostic.X
check('diagnostic_declared_values',np.max(np.abs(diff.data),initial=0)<1e-5)
hashes={}
for name,record in export['hashes'].items():
    hashes[name]=file_hashes(OUT/name)
    check(name+'_export_digest',hashes[name]['bytes']==record['bytes'] and hashes[name]['sha256']==record['sha256'])
summary=pd.read_csv(OUT/'sample_qc_summary.tsv',sep='\t')
check('summary_input_and_core',summary.released_nuclei.sum()==55122 and summary.core_CM.sum()==core.n_obs and summary.pending_CM.sum()==pending.n_obs and summary.suspected_doublets.sum()==dbl.n_obs)
check('summary_complete_samples',set(summary.sample_id)==set(samples.sample_id) and summary.sample_id.is_unique)
result={'status':'passed' if all(checks.values()) else 'failed','dataset':'liu2023','verified_utc':datetime.now(timezone.utc).isoformat(),'checks':checks,'hashes':hashes,'input_nuclei':55122,'input_samples':13,'input_donors':12,'strict_nuclei':combined.n_obs,'all_cells_qc':allcell.n_obs,'core_CM':core.n_obs,'light_CM':light.n_obs,'pending_CM':pending.n_obs,'suspected_doublets':dbl.n_obs,'genes':core.n_vars,'core_donors':int(core.obs.donor_id.nunique()),'core_disease_counts':core.obs.groupby('disease',observed=True).size().to_dict(),'pending_disease_counts':pending.obs.groupby('disease',observed=True).size().to_dict(),'residual_screens_triggered':samples.loc[samples.trigger_residual_doublet_screen.astype(str).str.lower().eq('true'),'sample_id'].tolist(),'verifier_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest()}
(OUT/'INDEPENDENT_VERIFICATION.json').write_text(json.dumps(result,indent=2)+'\n')
if result['status']!='passed':
    raise AssertionError([k for k,v in checks.items() if not v])
for algo,filename in [('md5','MD5SUMS'),('sha256','SHA256SUMS')]:
    (OUT/filename).write_text(''.join(f'{r[algo]}  {n}\n' for n,r in hashes.items()))
missing=sorted(set(samples.donor_id)-set(core.obs.donor_id))
report=f'''# Liu2023 清洗与验收完成

验收：{result['verified_utc']}，{len(checks)}项检查通过。
服务器：`{OUT}/`。

| 对象 | 核数 |
| --- | ---: |
| 作者已QC输入 | 55122 |
| 本轮严格候选 | {combined.n_obs} |
| 全细胞QC | {allcell.n_obs} |
| 心肌核心 / 轻量 | {core.n_obs} / {light.n_obs} |
| 心肌身份待定（单独保留） | {pending.n_obs} |
| 多证据疑似双细胞（单独保留） | {dbl.n_obs} |

原发布31944基因及未归一化整数计数完整保留。心肌核心按疾病：{result['core_disease_counts']}；覆盖{result['core_donors']}/12位供体，核心缺失供体：{missing or '无'}。少于20k不复制细胞。

来源：Liu et al., Cell Discovery 2023, DOI 10.1038/s41421-022-00490-3；10 HCM、2健康供体，室间隔snRNA-seq。13个发布样本中CL0412A/B的作者Subject均为CL0412，按同一生物供体绑定；不能声称13位供体或把A/B跨训练验证划分。两样本究竟属于技术还是不同组织重复尚未确定。

作者已按UMI、基因、MT和核糖体比例过滤，并用Scrublet分数>0.35和混合谱系过滤；本轮输入是已QC核，不能恢复原作者删除的核/基因，不能重新运行EmptyDrops或宣称重新识别原始液滴。13个发布MTX矩阵与Seurat RNA counts层逐值一致，未使用integrated assay。作者percent.mito单位为0–1比例，本轮从计数重算0–100百分数，使用项目严格基线500≤基因≤5000、MT<5%。这不是HCM统一行业硬阈值。

先以宽松跨谱系信号触发补充筛查（≥10且≥1%严格候选），实际补测样本：{result['residual_screens_triggered']}。宽松单UMI标记会标记23815/53538个输入候选及307/2603个正常参考核，不能直接当作双细胞证据。最终采用强程序：每谱系至少3基因各≥2 UMI，程序强度达到成人参考的max(本类q10,其他类q99)，双程序比例≥.25；仅14个严格候选通过强程序条件，作者CM为0。阈值来自同一固定参考留出集，不能把该集的表观低假阳性率当作独立性能证明。

仅补测工具调用及强混合程序共同支持时另存疑似双细胞。原始加载率和真实残余双细胞率未知，默认加载假设仅用于筛查，不能实验确证同型/异型双细胞；工具单独调用或弱跨谱系表达不导致删除。只有强混合程序而无调用支持的CM候选进入待定。第一次宽松程序导出已单独存档，最终版本为Liu2023_strict_v2。

心肌核心要求作者CM、适配成人左心室单核参考SingleR、多个心肌标记及程序一致，并有参考尾部或邻域支持。身份冲突、混合程序候选另存待定。参考训练/留出供体固定分开；阈值和表观参考性能共享同一留出集，不能当独立HCM验证或概率。原作者CM1/CM2及其他注释均保留，损失见author_subcluster_retention.tsv与sample_qc_summary.tsv。

独立验收重新读取全部13个原MTX，逐值核对5个计数对象的基因、条码、计数、重算QC与疾病/供体/作者注释；检查全部严格候选均在全细胞或疑似双细胞对象中，审计覆盖55122个发布核，原文件MD5/SHA256不变。诊断副本明确log1p(CP10k)，核对数值且不当原始计数输入。MD5SUMS/SHA256SUMS覆盖6个H5AD。

诊断副本之外不插补、归一化、整合或裁切HVG；不把疾病标签作为细胞身份输入。疾病为来源供体诊断，不是逐核病理确证。仅2位健康供体，逐供体数量与待定损失需在后续评估中报告。跨研究供体重叠及最终训练划分尚未解决；本轮不训练虚拟细胞模型。
'''
(OUT/'REPORT.md').write_text(report)
(WORK/'VERIFICATION_SUMMARY.json').write_text(json.dumps({k:v for k,v in result.items() if k not in ['checks','hashes']},indent=2)+'\n')
print(f'VERIFIED Liu2023: {len(checks)} checks, CM={core.n_obs}',flush=True)
