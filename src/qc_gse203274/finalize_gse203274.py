"""Join pediatric source labels, independent identity evidence and immutable counts."""
from datetime import datetime, timezone
from importlib.metadata import version
from pathlib import Path
import hashlib
import json
import shutil

import anndata as ad
import numpy as np
import pandas as pd
import scanpy as sc
from qc_common import PANELS,PAIRS

ROOT=Path.home()/'xiaoe-internship-data'
WORK=ROOT/'_processing/all_data_cleaning_20261111/GSE203274'
OUT=ROOT/'processed/pediatric_GSE203274_qc_20261111'
SEED=20261111


def digest(path):
    with path.open('rb') as handle:return hashlib.file_digest(handle,'sha256').hexdigest()


def light_ids(obs):
    target=min(20000,len(obs));groups={str(k):v.index.to_numpy() for k,v in obs.groupby('donor_id',observed=True)}
    allocation={k:0 for k in groups}
    for _ in range(target):
        available=[k for k in groups if allocation[k]<len(groups[k])]
        chosen=min(available,key=lambda k:(allocation[k],k));allocation[chosen]+=1
    rng=np.random.default_rng(SEED);selected=[]
    for donor in sorted(groups):selected.extend(rng.choice(groups[donor],allocation[donor],replace=False))
    return sorted(selected)


def main():
    assert (WORK/'R_exit_code.txt').read_text().strip()=='0'
    assert (WORK/'RESIDUAL_REVIEW_SUCCESS.txt').exists()
    OUT.mkdir(exist_ok=True,parents=True)
    assert not (OUT/'EXPORT_SUCCESS.json').exists(),'Inspect successful export before rerunning'
    raw=ad.read_h5ad(WORK/'author_full_counts.h5ad')
    a=raw[raw.obs.pass_strict_numeric.to_numpy()].copy()
    identities=pd.read_csv(WORK/'identity.tsv',sep='\t',index_col='cell_id',keep_default_na=False)
    markers=pd.read_csv(WORK/'strict_marker_evidence.tsv.gz',sep='\t',index_col='cell_id',keep_default_na=False)
    assert identities.index.equals(a.obs_names) and markers.index.equals(a.obs_names)
    a.obs=a.obs.join(identities).join(markers[[c for c in markers.columns if c not in a.obs]])
    doublets=[]
    for sample in a.obs.sample_id.unique():
        table=pd.read_csv(WORK/'captures'/sample/'doublets.tsv',sep='\t',index_col='cell_id')
        doublets.append(table)
    doubles=pd.concat(doublets)
    assert doubles.index.is_unique and a.obs_names.isin(doubles.index).all()
    a.obs=a.obs.join(doubles)
    a.obs['exclude_suspected_doublet']=a.obs.scDblFinder_class.eq('doublet')&a.obs.coherent_mixed_panel_flag
    cal=pd.read_csv(WORK/'reference_calibration.tsv',sep='\t',keep_default_na=False)
    tails={};performance=[]
    for label,g in cal.groupby('truth'):
        correct=g.loc[g.reference_label.eq(label)]
        if len(correct)>=30:
            tails[label]={'n_correct':len(correct),'delta_q05':float(correct.reference_delta.quantile(.05)),
                          'delta_next_q05':float(correct.reference_delta_next.quantile(.05))}
        truth=cal.truth.eq(label);prediction=cal.reference_label.eq(label)
        performance.append({'label':label,'reference_truth':int(truth.sum()),
                            'precision':float((truth&prediction).sum()/max(1,prediction.sum())),
                            'recall':float((truth&prediction).sum()/max(1,truth.sum()))})
    ref_report={'common_genes':len((WORK/'reference_common_genes.txt').read_text().splitlines()),
                'accuracy':float(cal.truth.eq(cal.reference_label).mean()),'tails':tails,'per_label':performance,
                'reference':'Existing adult LV nuclei; fixed donor-separated training/holdout; unique exact source symbol matches only',
                'limitation':'Pediatric domain shift unresolved; empirical reference tails are not posterior probabilities or independent disease validation; threshold and apparent tail metrics share the fixed holdout.'}
    a.obs['reference_tail_pass']=False
    for label,b in tails.items():
        passed=a.obs.reference_label.eq(label)&a.obs.reference_pruned.eq(label)
        passed &= a.obs.reference_delta.ge(b['delta_q05'])&a.obs.reference_delta_next.ge(b['delta_next_q05'])
        a.obs.loc[passed,'reference_tail_pass']=True
    print('Computing pediatric diagnostic PCA and exact neighbors on strict candidates',flush=True)
    d=a.copy();d.X=d.X.astype(np.float32)
    sc.settings.n_jobs=4
    sc.pp.normalize_total(d,target_sum=1e4);sc.pp.log1p(d)
    sc.pp.highly_variable_genes(d,n_top_genes=3000,flavor='seurat')
    d.var.loc[d.var.mt,'highly_variable']=False
    pcs=min(40,d.n_obs-1,int(d.var.highly_variable.sum())-1)
    sc.pp.pca(d,n_comps=pcs,mask_var='highly_variable',random_state=SEED)
    sc.pp.neighbors(d,n_neighbors=15,n_pcs=pcs,transformer='sklearn',random_state=SEED)
    sc.tl.leiden(d,resolution=.6,flavor='igraph',n_iterations=2,random_state=SEED)
    anchor=a.obs.reference_pruned.eq('Cardiomyocyte')&a.obs.Cardiomyocyte_markers.ge(3)&a.obs.cm_program_dominates&~a.obs.coherent_mixed_panel_flag
    graph=d.obsp['connectivities'].tocsr()
    a.obs['cm_anchor_neighbor_fraction']=np.asarray(graph@anchor.to_numpy(dtype=float)).ravel()/np.maximum(np.asarray(graph.sum(axis=1)).ravel(),1e-12)
    a.obs['diagnostic_leiden']=d.obs.leiden.astype(str)
    a.obsm['X_pca_qc']=d.obsm['X_pca']
    cm=a.obs.reference_label.eq('Cardiomyocyte')&a.obs.reference_pruned.eq('Cardiomyocyte')
    cm &= a.obs.Cardiomyocyte_markers.ge(2)&a.obs.cm_program_dominates&~a.obs.coherent_mixed_panel_flag&~a.obs.exclude_suspected_doublet
    cm &= a.obs.reference_tail_pass|(a.obs.Cardiomyocyte_markers.ge(3)&a.obs.cm_anchor_neighbor_fraction.ge(.5))
    cm &= a.obs.author_MainCellType.eq('CM')
    candidate=a.obs.author_MainCellType.eq('CM')|a.obs.reference_label.eq('Cardiomyocyte')|(a.obs.Cardiomyocyte_markers.ge(3)&a.obs.cm_program_dominates)
    pending=candidate&~cm&~a.obs.exclude_suspected_doublet
    a.obs['cell_type_qc']='Unknown'
    for label in a.obs.reference_label.unique():
        panel='Mural' if label in ['Pericyte','Smooth_muscle'] else label
        if label!='Cardiomyocyte' and panel in PANELS:
            passed=a.obs.reference_pruned.eq(label)&a.obs[panel+'_markers'].ge(2)
            a.obs.loc[passed,'cell_type_qc']=label
    a.obs.loc[cm,'cell_type_qc']='Cardiomyocyte'
    a.obs['review_action']=np.select([a.obs.exclude_suspected_doublet,cm,pending],
                                   ['exclude_suspected_doublet','retain_CM','pending_CM'],default='retain_other_or_unknown')
    a.obs['qc_eligible_CM']=cm
    a.obs['qc_eligible_binary_CM']=cm&a.obs.binary_disease_eligible
    a.obs['training_use']='requires_pediatric_design_and_donor_overlap_split_review_single_HCM_donor'
    assert a.obs.author_echoEF.astype(str).eq(a.obs.sample_id.astype(str)).all()
    a.obs['numeric_ejection_fraction_available']=False
    raw.obs['numeric_ejection_fraction_available']=False
    mixed=json.loads((WORK/'MIXED_PROGRAM_REVIEW.json').read_text())
    residual=json.loads((WORK/'RESIDUAL_REVIEW_POLICY.json').read_text())
    policy={'version':'GSE203274_pediatric_strict_v1','seed':SEED,'min_genes':500,'max_genes':5000,'mt_percent_exclusive_max':5,
            'count_source':'GEO AllNuclei rawCount file: released unnormalized integer RNA counts after author CellBender v0.1.0 background correction and author QC. Not uncorrected CellRanger molecular counts, scaled data, Harmony coordinates or full raw droplets.',
            'author_QC':'Hill et al., Nature 2022, PMID35732239/PMC10405779: hg19 pre-mRNA CellRanger3.0.1, CellBender0.1.0, 700-10000 genes, MT<0.2%, Harmony and repeated global/cell-type mixed-marker cluster removal; published set157293 nuclei. Source restrictions inherited and excluded original nuclei/genes cannot be recovered.',
            'numeric':'Recompute from released full counts: 500<=detected genes<=5000 and MT%<5. Published percent.mt is percentage, not fraction; exactly matches recomputation. Inherited MT<0.2% is already stricter; this project upper-gene threshold adds selection.',
            'gene_axis':'All29266 native source symbols retained including zeros/duplicates. Derived reference_gene_id assigned only for18393 unique exact symbol matches to existing GRCh38 adult LV reference; source assembly is hg19, no native Ensembl IDs claimed or guessed; unmapped/ambiguous genes remain.',
            'disease_map':{'Control':0,'HCM':1,'DCM':-1,'HLHS':-1,'TOF':-1},
            'disease_source':'Use author ClinicalRank and Diagnosis with GEO clinical diagnoses. Lowercase diagnosis mixes LVAD/HTx/procedure status with disease and cannot define HCM/Control alone. Non-HCM pathology is never relabeled Control.',
            'donor_binding':'Use13 author labID donors,30 orig.ident batches. Same13_198 LV/RV is one donor even though author patientID splits sites; repeats/plate counts are not independent patients. P33 9k/10k and other batches retain candidate GEO sets, no invented unique GSM association.',
            'age':'Pediatric source text retained and split into year/month/day components; approximate years explicitly derived. UK1/UK2 author age and GEO age disagree, both retained with true conflict flags; ageCont is source metadata and not preferred exact age.',
            'clinical_field_issue':'author echoEF equals orig.ident/capture ID for every nucleus; retained unchanged, numeric_ejection_fraction_available=False, no invented EF values.',
            'reference_identity':'SingleR adult LV reference with fixed donor-separated holdout on unique symbol/Ensembl matches. Source author labels, multiple markers and neighborhood support combined; neonatal/pediatric reference domain shift remains. Reference scores/tails and correlated graph support are not independent confidence probabilities.',
            'CM_core':'Author MainCellType=CM, SingleR best and pruned Cardiomyocyte, >=2 CM-panel genes, CM-program dominance, no strong mixed or joint suspected-doublet exclusion; and empirical reference-tail pass OR >=3 CM genes with >=.5 CM-anchor neighbor fraction.',
            'CM_pending':'Author CM OR best reference CM OR >=3 CM genes with CM program dominance, but not core and not jointly excluded; all diagnoses retained with source label.',
            'doublet_model':residual['trigger_rule'],'doublet_exclusion':residual['exclusion_rule'],
            'doublet_rate_limitation':'If triggered, dbr=NULL,dbr.per1k=.008 is a supplemental residual-screening heuristic; published post-QC batch size is not known original loading or a measured residual rate. No EmptyDrops or repeated CellBender.',
            'panels':PANELS,'mixed_program_pairs':PAIRS,'mixed_program_bounds':mixed['bounds'],
            'strong_mixed_rule':'Each distinct lineage >=3 genes with >=2 UMIs; both mean log1p(CP10k) programs>=max(fixed reference positive q10,negative q99), min/max program ratio>=.25. Weak one-UMI combinations are diagnostic only, not deletion criteria.',
            'graph':'Pediatric cohort alone; diagnostic log1p(CP10k),3000HVGs excluding MT,40PCs,15 exact sklearn neighbors,Leiden.6; no disease integration or disease labels in neighbor construction. Original-count objects retain all genes.',
            'output_scope':'all_cells_qc retains all source diagnoses; cm_qc_full/light contain ONLY Control/HCM core; cm_other_diseases retains other verified CM with label=-1. Pending, suspected doublets and numeric exclusions kept separately; no duplicate nuclei to fill20k.',
            'training':'QC completed, train/test unassigned. Only one HCM donor P86 at11y3m12d; two captures are not independent HCM patients. Pediatric/external-reference/region and cross-study donor-overlap review needed before training or validation; not an independent multi-donor HCM benchmark.',
            'scope':'No tokenization, virtual-cell training, differential expression, target ranking, count normalization/imputation/integration or HVG truncation in delivered counts.'}
    a.layers['counts']=a.X.copy();a.uns['qc_policy_json']=json.dumps(policy);a.uns['reference_calibration_json']=json.dumps(ref_report)
    raw.obs['final_action']='outside_strict_numeric'
    raw.obs.loc[a.obs_names,'final_action']=a.obs.review_action
    raw.obs['qc_eligible_CM']=False;raw.obs.loc[a.obs_names,'qc_eligible_CM']=cm.to_numpy()
    raw.obs['qc_eligible_binary_CM']=False;raw.obs.loc[a.obs_names,'qc_eligible_binary_CM']=a.obs.qc_eligible_binary_CM.to_numpy()
    raw.obs['training_use']='requires_pediatric_design_and_donor_overlap_split_review_single_HCM_donor'
    raw.layers['counts']=raw.X.copy();raw.uns['qc_policy_json']=json.dumps(policy)
    raw.obs.to_csv(OUT/'all_nuclei_audit.tsv.gz',sep='\t',index_label='cell_id')
    a.obs.to_csv(OUT/'cell_identity_audit.tsv.gz',sep='\t',index_label='cell_id')
    binary_core=cm&a.obs.binary_disease_eligible
    masks={'all_cells_qc.h5ad':~a.obs.exclude_suspected_doublet,'cm_qc_full.h5ad':binary_core,
           'cm_other_diseases.h5ad':cm&~a.obs.binary_disease_eligible,'cm_review_pending.h5ad':pending,
           'suspected_doublets.h5ad':a.obs.exclude_suspected_doublet}
    for filename,mask in masks.items():a[mask.to_numpy()].copy().write_h5ad(OUT/filename,compression='gzip')
    light=a[light_ids(a.obs.loc[binary_core])].copy();light.write_h5ad(OUT/'cm_qc_light.h5ad',compression='gzip')
    raw.write_h5ad(OUT/'all_source_counts.h5ad',compression='gzip')
    raw[~raw.obs.pass_strict_numeric.to_numpy()].copy().write_h5ad(OUT/'outside_strict_numeric.h5ad',compression='gzip')
    d.obs['review_action']=a.obs.review_action
    d.uns['analysis_only']=True;d.uns['count_source']='Diagnostic X=log1p(CP10k); unchanged full-gene author counts are delivered in separate count objects'
    d.write_h5ad(OUT/'diagnostic_normalized.h5ad',compression='gzip')
    samples=pd.read_csv(WORK/'samples.tsv',sep='\t')
    residual_summary=pd.read_csv(WORK/'RESIDUAL_REVIEW_SUMMARY.tsv',sep='\t').set_index('sample_id')
    rows=[]
    for row in samples.itertuples():
        original=raw.obs.loc[raw.obs.sample_id.eq(row.sample_id)];current=a.obs.loc[a.obs.sample_id.eq(row.sample_id)]
        rows.append({'sample_id':row.sample_id,'donor_id':row.donor_id,'disease':row.disease,'clinical_group':row.clinical_group,
                     'region':row.tissue_region,'source_nuclei':len(original),'genes_500to5000':int(original.n_genes_by_counts.between(500,5000).sum()),
                     'MT_lt5pct':int(original.pct_counts_mt.lt(5).sum()),'strict_nuclei':len(current),
                     'author_CM':int(original.author_MainCellType.eq('CM').sum()),'author_CM_strict':int(current.author_MainCellType.eq('CM').sum()),
                     'author_CM_core':int((current.author_MainCellType.eq('CM')&current.qc_eligible_CM).sum()),
                     'author_CM_pending':int((current.author_MainCellType.eq('CM')&current.review_action.eq('pending_CM')).sum()),
                     'CM_core_all_diagnoses':int(current.qc_eligible_CM.sum()),'CM_core_binary':int(current.qc_eligible_binary_CM.sum()),
                     'pending_CM':int(current.review_action.eq('pending_CM').sum()),'strong_mixed':int(current.coherent_mixed_panel_flag.sum()),
                     'new_suspected_doublets':int(current.exclude_suspected_doublet.sum()),
                     'residual_model_triggered':bool(residual_summary.loc[row.sample_id,'residual_model_triggered']),
                     'age_source_conflict':row.age_source_conflict})
    coverage=pd.DataFrame(rows);coverage.to_csv(OUT/'sample_qc_summary.tsv',sep='\t',index=False)
    counts=[c for c in coverage.columns if c not in ['sample_id','donor_id','disease','clinical_group','region','residual_model_triggered','age_source_conflict']]
    coverage.groupby(['donor_id','disease','clinical_group'],observed=True)[counts].sum().to_csv(OUT/'donor_qc_summary.tsv',sep='\t')
    samples.to_csv(OUT/'sample_metadata.tsv',sep='\t',index=False)
    a.obs.groupby(['donor_id','disease','author_MainCellType','review_action'],observed=True).size().to_csv(OUT/'author_identity_retention.tsv',sep='\t',header=['nuclei'])
    a.obs.groupby(['disease','author_Cluster','review_action'],observed=True).size().to_csv(OUT/'author_cluster_retention.tsv',sep='\t',header=['nuclei'])
    raw.obs.groupby(['disease','numeric_exclusion_reason'],observed=True).size().to_csv(OUT/'numeric_loss_by_disease.tsv',sep='\t',header=['nuclei'])
    for filename in ['STAGE_SUCCESS.json','METADATA_REPAIR_VERIFICATION.json','RESIDUAL_REVIEW_SUMMARY.tsv','RESIDUAL_REVIEW_POLICY.json','reference_gene_mapping.tsv','METHODS_EXCERPTS.txt']:
        shutil.copy2(WORK/filename,OUT/filename)
    (OUT/'qc_policy.json').write_text(json.dumps(policy,indent=2)+'\n')
    (OUT/'REFERENCE_CALIBRATION.json').write_text(json.dumps(ref_report,indent=2)+'\n')
    (OUT/'MIXED_PROGRAM_REVIEW.json').write_text(json.dumps(mixed,indent=2)+'\n')
    core_obs=a.obs.loc[binary_core]
    result={'dataset':'GSE203274','status':'exported_pending_independent_verification','completed_utc':datetime.now(timezone.utc).isoformat(),
            'source_nuclei':raw.n_obs,'genes':raw.n_vars,'strict_nuclei':a.n_obs,'source_donors':13,'source_samples':30,
            'core_binary_CM':int(binary_core.sum()),'core_HCM':int((binary_core&a.obs.disease.eq('HCM')).sum()),
            'core_Control':int((binary_core&a.obs.disease.eq('Control')).sum()),'core_binary_donors':int(core_obs.donor_id.nunique()),
            'core_other_disease_CM':int((cm&~a.obs.binary_disease_eligible).sum()),'pending_CM':int(pending.sum()),
            'strong_mixed':int(a.obs.coherent_mixed_panel_flag.sum()),'new_suspected_doublets':int(a.obs.exclude_suspected_doublet.sum()),
            'light_CM':light.n_obs,'training_split':'unassigned_single_HCM_donor_pediatric_design_review_required',
            'hashes':{p.name:{'bytes':p.stat().st_size,'sha256':digest(p)} for p in sorted(OUT.glob('*.h5ad'))},
            'python_versions':{p:version(p) for p in ['scanpy','anndata','numpy','scipy','pandas']}}
    (OUT/'EXPORT_SUCCESS.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps({k:v for k,v in result.items() if k!='hashes'},indent=2),flush=True)


if __name__=='__main__':main()
