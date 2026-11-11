"""Apply the frozen strict baseline without discarding excluded whole-cell counts."""
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

ROOT = Path.home()/'xiaoe-internship-data'
WORK = ROOT/'_processing/all_data_cleaning_20261111/gse138262'
OUT = ROOT/'processed/wehrens_GSE138262_qc_20261111'
SEED = 20261111
PANELS = {
 'Cardiomyocyte':['TNNT2','TTN','TNNI3','ACTC1','MYH7','MYBPC3','MYL2','TNNC1','RYR2','PLN','ACTN2','MYH6'],
 'Fibroblast':['COL1A1','COL1A2','DCN','LUM','PDGFRA'],
 'Endothelial':['PECAM1','VWF','KDR','EMCN','CLDN5'],
 'Immune':['PTPRC','LST1','TYROBP','CD3D','CD79A'],
 'Mural':['RGS5','PDGFRB','CSPG4','MYH11','ACTA2','TAGLN'],
 'Neural':['NRXN1','SNAP25','SYT1','SLC1A3','PLP1'],
 'Adipocyte':['ADIPOQ','FABP4','PLIN1','LPL','CIDEC'],
 'Mesothelial':['WT1','MSLN','UPK3B','KRT8','KRT18']}
PAIRS = [('Cardiomyocyte','Fibroblast'),('Cardiomyocyte','Endothelial'),('Cardiomyocyte','Immune'),
 ('Cardiomyocyte','Mural'),('Fibroblast','Endothelial'),('Fibroblast','Immune'),('Endothelial','Immune'),('Immune','Mural')]


def digest(path):
    with path.open('rb') as handle:
        return hashlib.file_digest(handle,'sha256').hexdigest()


def main():
    assert (WORK/'IDENTITY_SUCCESS.txt').exists()
    OUT.mkdir(exist_ok=True,parents=True)
    assert not (OUT/'EXPORT_SUCCESS.json').exists(), 'Existing final export requires inspection before rerun'
    raw = ad.read_h5ad(WORK/'corrected_native_union_counts.h5ad')
    a = raw[raw.obs.pass_strict_numeric.to_numpy()].copy()
    predictions = pd.read_csv(WORK/'identity.tsv',sep='\t',index_col='cell_id',keep_default_na=False)
    assert set(predictions.index)==set(a.obs_names)
    a.obs = a.obs.join(predictions)
    calibration = pd.read_csv(WORK/'reference_calibration.tsv',sep='\t',keep_default_na=False)
    mixed = json.loads((WORK/'MIXED_PROGRAM_REVIEW.json').read_text())
    tails = {}
    for label, group in calibration.groupby('truth'):
        correct = group.loc[group.reference_label.eq(label)]
        if len(correct)>=30:
            tails[label]={'n_correct':len(correct),'delta_q05':float(correct.reference_delta.quantile(.05)),
                          'delta_next_q05':float(correct.reference_delta_next.quantile(.05))}
    reference_report = {'common_genes':len((WORK/'reference_common_genes.txt').read_text().splitlines()),
                        'accuracy':float(calibration.truth.eq(calibration.reference_label).mean()),'tails':tails,
                        'limitation':'Adult LV nuclear reference applied to whole cells; empirical tails are correlated heuristics, not posterior probabilities or disease validation. Same fixed donor holdout supplies thresholds and their apparent metrics.'}
    a.obs['reference_tail_pass']=False
    for label,b in tails.items():
        flag=a.obs.reference_label.eq(label)&a.obs.reference_pruned.eq(label)
        flag &= a.obs.reference_delta.ge(b['delta_q05'])&a.obs.reference_delta_next.ge(b['delta_next_q05'])
        a.obs.loc[flag,'reference_tail_pass']=True
    total=np.asarray(a.X.sum(axis=1)).ravel()
    for label,panel in PANELS.items():
        loc=np.flatnonzero(a.var.gene_symbol.isin(panel).to_numpy())
        assert len(loc)>=3,(label,len(loc))
        values=a.X[:,loc].toarray()
        a.obs[label+'_markers']=(values>0).sum(axis=1)
        a.obs[label+'_markers_ge2_UMI']=(values>=2).sum(axis=1)
        a.obs[label+'_program']=np.log1p(values*(10000/np.maximum(total,1))[:,None]).mean(axis=1)
    a.obs['cm_program_dominates']=a.obs.Cardiomyocyte_program.ge(a.obs[[p+'_program' for p in PANELS if p!='Cardiomyocyte']].max(axis=1))
    strong=np.zeros(a.n_obs,dtype=bool);reason=np.full(a.n_obs,'',dtype=object)
    for left,right in PAIRS:
        l=a.obs[left+'_program'].to_numpy();r=a.obs[right+'_program'].to_numpy()
        flag=(a.obs[left+'_markers_ge2_UMI'].to_numpy()>=3)&(a.obs[right+'_markers_ge2_UMI'].to_numpy()>=3)
        flag &= (l>=mixed['bounds'][left]['threshold'])&(r>=mixed['bounds'][right]['threshold'])
        flag &= np.minimum(l,r)/np.maximum(np.maximum(l,r),1e-12)>=.25
        reason[flag&~strong]=left+'+'+right;strong|=flag
    a.obs['coherent_mixed_panel_flag']=strong;a.obs['mixed_panel_pair']=reason
    a.obs['cm_anchor_neighbor_fraction']=0.
    a.obs['neighbor_evidence_enabled']=False
    diagnostic_files=[]
    for study in a.obs.study_origin.unique():
        ids=a.obs_names[a.obs.study_origin.eq(study)]
        d=a[ids].copy();d.layers.clear()
        d.X=d.X.astype(np.float32)
        sc.pp.normalize_total(d,target_sum=1e4);sc.pp.log1p(d)
        if d.n_obs>=5:
            sc.pp.highly_variable_genes(d,n_top_genes=min(3000,d.n_vars),flavor='seurat')
            d.var.loc[d.var.mt,'highly_variable']=False
            pcs=min(40,d.n_obs-1,int(d.var.highly_variable.sum())-1)
            sc.pp.pca(d,n_comps=pcs,mask_var='highly_variable',random_state=SEED)
            sc.pp.neighbors(d,n_neighbors=min(15,d.n_obs-1),n_pcs=pcs,transformer='sklearn',random_state=SEED)
            # Fewer than 50 cells do not supply an adequate cohort graph for identity rescue.
            if d.n_obs>=50:
                anchor=(a.obs.loc[ids,'reference_pruned'].eq('Cardiomyocyte')&a.obs.loc[ids,'Cardiomyocyte_markers'].ge(3)
                        &a.obs.loc[ids,'cm_program_dominates']&~a.obs.loc[ids,'coherent_mixed_panel_flag'])
                graph=d.obsp['connectivities'].tocsr()
                support=np.asarray(graph@anchor.to_numpy(dtype=float)).ravel()/np.maximum(np.asarray(graph.sum(axis=1)).ravel(),1e-12)
                a.obs.loc[ids,'cm_anchor_neighbor_fraction']=support
                a.obs.loc[ids,'neighbor_evidence_enabled']=True
        d.uns['analysis_only']=True
        d.uns['count_source']='Diagnostic X=log1p(CP10k), no counts layer; unchanged counts are delivered separately'
        d.uns['study_graph_policy']='Each source study separately; no disease integration or cross-study neighbors; no identity rescue below 50 strict cells'
        path=OUT/('diagnostic_'+str(study)+'.h5ad')
        d.write_h5ad(path,compression='gzip');diagnostic_files.append(path.name)
    cm=a.obs.reference_label.eq('Cardiomyocyte')&a.obs.reference_pruned.eq('Cardiomyocyte')
    cm &= a.obs.Cardiomyocyte_markers.ge(2)&a.obs.cm_program_dominates&~a.obs.coherent_mixed_panel_flag
    cm &= a.obs.reference_tail_pass|(a.obs.neighbor_evidence_enabled&a.obs.Cardiomyocyte_markers.ge(3)&a.obs.cm_anchor_neighbor_fraction.ge(.5))
    a.obs['cell_type_qc']='Unknown'
    for label in a.obs.reference_label.unique():
        panel='Mural' if label in ['Pericyte','Smooth_muscle'] else label
        if panel in PANELS and label!='Cardiomyocyte':
            flag=a.obs.reference_pruned.eq(label)&a.obs[panel+'_markers'].ge(2)
            a.obs.loc[flag,'cell_type_qc']=label
    a.obs.loc[cm,'cell_type_qc']='Cardiomyocyte'
    # Both source collections preselected CM candidates. Unsupported identities remain reviewable.
    a.obs['review_action']=np.where(cm,'retain_CM','pending_identity_or_mixture')
    a.obs['qc_eligible_CM']=cm
    a.obs['training_use']='not_approved_strict_low_yield_no_control_requires_protocol_and_design_review'
    policy={'version':'GSE138262_corrected_strict_v1','seed':SEED,'min_genes':500,'max_genes':5000,'mt_percent_exclusive_max':5,
            'count_source':'32 full corrected UMI count tables at Zenodo version DOI 10.5281/zenodo.18206554; retain native Ensembl/symbol IDs, missing symbols, all cells and genes. RNA assays in released Seurat objects have already removed mitochondrial genes; never derive MT% from those assays.',
            'numeric':'500 <= detected genes <= 5000, MT percentage <5 from full corrected count tables; frozen project baseline, not a universal whole-cell heart threshold. No silent relaxation even if yield is near zero.',
            'source_scope':'5 HCM donors/16 CEL-seq2 plates from GSE138262; 7 separately identified Wang healthy donors/16 iCell8 plate-parts are author-selected LV control companions. Author ALL also contains Litvinukova cells, inventoried separately; not appended or declared independent validation.',
            'disease_map':{'Control':0,'HCM':1},'disease_definition':'Source donor diagnosis; not independently measured per-cell pathology.',
            'cell_calling_and_doublets':'FACS-sorted CEL-seq2 and selected iCell8 whole cells: no EmptyDrops, no 10x-loading scDblFinder/Scrublet model, no claim that high gene counts or mixed programs prove doublets. Strong mixed programs stay in review; homotypic doublets remain unresolved.',
            'identity':'Adapted adult LV nuclear SingleR reference with fixed donor-separated calibration, >=2 CM panel genes, CM-program dominance, no coherent mixed panels, and empirical reference-tail pass or >=3 CM genes plus >=.5 anchor-neighbor fraction. Neighbor rescue enabled only for a source study with >=50 strict cells; this dataset has too few cells. Domain shift and small sample limit confidence.',
            'mixed_program_bounds':mixed['bounds'],'mixed_program_pairs':PAIRS,'panels':PANELS,
            'mixed_rule':'Each lineage >=3 genes with >=2 UMIs; both mean log1p(CP10k) programs >= max(fixed reference positive q10,negative q99); min/max score ratio >=.25. Heuristic correlated evidence, not independent doublet confidence.',
            'graph':'Per source study only, diagnostic log1p(CP10k), <=3000 HVGs excluding MT, <=40 PCs, <=15 exact sklearn neighbors; graphs do not replace raw counts or integrate disease.',
            'retention':'All 7537 input cells retained in all_source_counts.h5ad, all 28663 union features retained in every count object; per-native-axis presence mask and native plate H5AD files preserved.',
            'training':'Not approved under the frozen strict baseline: only 8 numeric candidates, no healthy controls; no train/test split, tokenization, training, differential expression or target ranking.',
            'light':'At most 20000 distinct core cells, no duplication; this strict result is already below the limit.',
            'author_comparison':'All 48 shared plate/object RNA comparisons match corrected source counts exactly; only shared genes/cells are comparable because published RNA assays are filtered.',
            'unmatched_features':'Empty author gene symbols retained with authoritative Ensembl identifiers, not guessed or discarded.',
            'author_external_reference_overlap':'Author ALL incorporates Litvinukova2020 reference cells also used elsewhere; no independent disease-validation claim or automatic train inclusion.'}
    a.layers['counts']=a.X.copy()
    a.uns['qc_policy_json']=json.dumps(policy)
    a.uns['reference_calibration_json']=json.dumps(reference_report)
    raw.obs['final_action']='outside_strict_numeric'
    raw.obs.loc[a.obs_names,'final_action']=a.obs.review_action
    raw.obs['qc_eligible_CM']=False;raw.obs.loc[a.obs_names,'qc_eligible_CM']=cm.to_numpy()
    raw.obs['training_use']='not_approved_strict_low_yield_no_control_requires_protocol_and_design_review'
    raw.layers['counts']=raw.X.copy();raw.uns['qc_policy_json']=json.dumps(policy)
    raw.obs.to_csv(OUT/'all_cell_audit.tsv.gz',sep='\t',index_label='cell_id')
    a.obs.to_csv(OUT/'cell_identity_audit.tsv.gz',sep='\t',index_label='cell_id')
    objects={'all_source_counts.h5ad':raw,'all_cells_qc.h5ad':a,'cm_qc_full.h5ad':a[cm.to_numpy()].copy(),
             'cm_review_pending.h5ad':a[~cm.to_numpy()].copy(),'mixed_lineage_review.h5ad':a[strong].copy(),
             'outside_strict_numeric.h5ad':raw[~raw.obs.pass_strict_numeric.to_numpy()].copy()}
    core=objects['cm_qc_full.h5ad']
    assert core.n_obs<=20000
    objects['cm_qc_light.h5ad']=core.copy()
    for filename,obj in objects.items():
        obj.write_h5ad(OUT/filename,compression='gzip')
    samples=pd.read_csv(WORK/'samples.tsv',sep='\t')
    summaries=[]
    for row in samples.itertuples():
        obs=raw.obs.loc[raw.obs.sample_id.eq(row.sample_id)]
        strict_obs=a.obs.loc[a.obs.sample_id.eq(row.sample_id)]
        summaries.append({'sample_id':row.sample_id,'donor_id':row.donor_id,'disease':row.disease,'technology':row.technology,
                          'source_cells':len(obs),'genes_500to5000':int(obs.n_genes_by_counts.between(500,5000).sum()),
                          'MT_lt5pct':int(obs.pct_counts_mt.lt(5).sum()),'strict_cells':len(strict_obs),
                          'CM_core':int(strict_obs.qc_eligible_CM.sum()),'pending_identity':int((~strict_obs.qc_eligible_CM).sum()),
                          'strong_mixed':int(strict_obs.coherent_mixed_panel_flag.sum()),
                          'median_MT_all':float(obs.pct_counts_mt.median()),'median_genes_all':float(obs.n_genes_by_counts.median()),
                          'numeric_retention_fraction':float(obs.pass_strict_numeric.mean()),
                          'author_analysis_cells':int(obs.author_analysis_included.sum()),
                          'author_analysis_strict_cells':int((obs.author_analysis_included&obs.pass_strict_numeric).sum())})
    sample_summary=pd.DataFrame(summaries)
    sample_summary.to_csv(OUT/'sample_qc_summary.tsv',sep='\t',index=False)
    aggregation={'source_cells':'sum','genes_500to5000':'sum','MT_lt5pct':'sum','strict_cells':'sum','CM_core':'sum','pending_identity':'sum','strong_mixed':'sum','author_analysis_cells':'sum','author_analysis_strict_cells':'sum'}
    sample_summary.groupby(['donor_id','disease']).agg(aggregation).to_csv(OUT/'donor_qc_summary.tsv',sep='\t')
    samples.to_csv(OUT/'sample_metadata.tsv',sep='\t',index=False)
    for filename in ['AUTHOR_COUNT_COMPARISONS.tsv','DOWNLOAD_VERIFICATION.json','STAGE_SUCCESS.json']:
        shutil.copy2(WORK/filename,OUT/filename)
    (OUT/'qc_policy.json').write_text(json.dumps(policy,indent=2)+'\n')
    (OUT/'REFERENCE_CALIBRATION.json').write_text(json.dumps(reference_report,indent=2)+'\n')
    (OUT/'MIXED_PROGRAM_REVIEW.json').write_text(json.dumps(mixed,indent=2)+'\n')
    raw.obs.groupby(['disease','numeric_exclusion_reason'],observed=True).size().to_csv(OUT/'numeric_loss_by_disease.tsv',sep='\t',header=['cells'])
    result={'dataset':'GSE138262','status':'exported_pending_independent_verification','completed_utc':datetime.now(timezone.utc).isoformat(),
            'source_cells':raw.n_obs,'genes':raw.n_vars,'strict_cells':a.n_obs,'core_CM':core.n_obs,
            'pending_identity':int((~cm).sum()),'strong_mixed':int(strong.sum()),'donors_source':int(raw.obs.donor_id.nunique()),
            'core_donors':int(core.obs.donor_id.nunique()),'core_HCM':int(core.obs.disease.eq('HCM').sum()),
            'core_Control':int(core.obs.disease.eq('Control').sum()),'training_use':'not_approved_strict_low_yield_no_Control',
            'diagnostic_files':diagnostic_files,'python_versions':{p:version(p) for p in ['scanpy','anndata','numpy','scipy','pandas']},
            'hashes':{p.name:{'bytes':p.stat().st_size,'sha256':digest(p)} for p in sorted(OUT.glob('*.h5ad'))}}
    (OUT/'EXPORT_SUCCESS.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps({k:v for k,v in result.items() if k!='hashes'},indent=2),flush=True)


if __name__=='__main__':
    main()
