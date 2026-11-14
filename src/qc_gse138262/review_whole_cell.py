"""Versioned whole-cell review using author depth QC and conservative CM identity."""
import argparse
from datetime import datetime, timezone
import gzip
import hashlib
import json
import os
from pathlib import Path

os.environ.setdefault('OPENBLAS_NUM_THREADS', '2')
os.environ.setdefault('OMP_NUM_THREADS', '2')
import anndata as ad
import numpy as np
import pandas as pd
from scipy import sparse
from scipy.io import mmread, mmwrite
from threadpoolctl import threadpool_limits

SOURCE = Path('/home/zhome/xiaoe-internship-data/processed/wehrens_GSE138262_qc_20261111/all_source_counts.h5ad')
OLD_CORE = SOURCE.parent/'cm_qc_full.h5ad'
REF = Path('/home/zhome/xiaoe-internship-data/processed/SCP1303_review_20261111/identity')
PANELS = {
    'Cardiomyocyte': ['TNNT2','TTN','TNNI3','ACTC1','MYH7','MYBPC3','MYL2','TNNC1','RYR2','PLN','ACTN2','MYH6'],
    'Fibroblast': ['COL1A1','COL1A2','DCN','LUM','PDGFRA'],
    'Endothelial': ['PECAM1','VWF','KDR','EMCN','CLDN5'],
    'Immune': ['PTPRC','LST1','TYROBP','CD3D','CD79A'],
    'Mural': ['RGS5','PDGFRB','CSPG4','MYH11','ACTA2','TAGLN'],
    'Neural': ['NRXN1','SNAP25','SYT1','SLC1A3','PLP1'],
    'Adipocyte': ['ADIPOQ','FABP4','PLIN1','LPL','CIDEC'],
    'Mesothelial': ['WT1','MSLN','UPK3B','KRT8','KRT18']}


def digest(path, algorithm='sha256'):
    with path.open('rb') as stream:
        return hashlib.file_digest(stream, algorithm).hexdigest()


def save_json(path, value):
    path.write_text(json.dumps(value, ensure_ascii=False, indent=2)+'\n', encoding='utf-8')


def write_tsv(frame, path, index=False):
    frame.astype(object).where(frame.notna(), 'NA').replace('', 'NA').to_csv(path, sep='\t', index=index)


def programs(x, symbols, nonmt_total):
    result = pd.DataFrame(index=np.arange(x.shape[0]))
    for label, panel in PANELS.items():
        positions = np.flatnonzero(pd.Index(symbols).isin(panel))
        assert len(positions) >= 3, (label, len(positions))
        values = x[:, positions].toarray()
        result[label+'_genes_ge2'] = (values >= 2).sum(axis=1)
        result[label+'_program'] = np.log1p(values*10000/np.maximum(nonmt_total, 1)[:, None]).mean(axis=1)
    return result


def stage(work):
    assert not work.exists(), 'Use a fresh work directory'
    work.mkdir(parents=True)
    policy = {
        'version': 'GSE138262_author_depth_whole_cell_v2_20261114', 'created_utc': datetime.now(timezone.utc).isoformat(),
        'author_repo_commit': '621c874b6e1f98f6eee83084c6424a52e24d4f97',
        'author_script': 'Projects/HCM_SCS_2021_06_SeuratRevisedAnalysis_v2_UmiTools.R',
        'author_lines': '626-635, 738-777; executed filter line 769, not commented line 766',
        'numeric_gate': 'author_analysis_included AND total non-MT UMIs > 1000',
        'MT_cutoff': None, 'upper_gene_cutoff': None,
        'numeric_rationale': 'Source-specific author protocol for sorted whole cells; MT-heavy whole cells differ from nuclei. MT retained in stored counts.',
        'outlier_review': '3 raw MAD within each native table/plate-part (>=20 candidates; else study fallback). High MT AND low log1p non-MT detected genes goes to review; high depth/complexity alone is a flag.',
        'identity': 'SingleR broad and pruned Cardiomyocyte labels AND >=3 CM panel genes each >=2 UMIs AND CM-program dominance; no graph rescue.',
        'identity_normalization': 'Panel scores: diagnostic log1p(CP10k) using full non-MT library size. SingleR: common non-MT gene subset normalized monotonically within each cell. Stored X/counts remain untouched.',
        'mixed_review': 'For supported reference lineages: CM and another lineage each >=3 genes with >=2 UMIs, both programs >=max(reference truth positive q10,negative q99), weaker/stronger score >=0.25. For rare lineages without >=20 reference positives: >=3 genes each >=2 UMIs in both panels and other/CM program >=0.5 is an uncalibrated review flag. Neither proves a doublet.',
        'doublet_scope': 'No 10x loading model, EmptyDrops, Scrublet or automatic high-gene doublet call for these plate-sorted/selected whole cells; homotypic doublets unresolved.',
        'reference': 'Same fixed donor-separated adult LV nuclear reference and calibration; retrain SingleR on common non-MT genes; reference domain shift remains.',
        'use': 'Supplementary whole-cell resource; disease, study and technology confounded. Not approved as independent binary validation or automatically pooled adult snRNA training.',
        'split': 'unassigned', 'originals': {str(SOURCE): digest(SOURCE), str(OLD_CORE): digest(OLD_CORE)},
        'new_DGE_or_training': False, 'new_tokenization': False, 'preserve_all_source_features': True}
    save_json(work/'POLICY.json', policy)  # Freeze before identity predictions or outcome inspection.
    a = ad.read_h5ad(SOURCE)
    assert sparse.issparse(a.X) and 'counts' in a.layers
    assert (a.X != a.layers['counts']).nnz == 0
    x = a.X.tocsr()
    mt = a.var['mt'].to_numpy(dtype=bool)
    nonmt_total = np.asarray(x[:, ~mt].sum(axis=1)).ravel()
    nonmt_genes = np.asarray((x[:, ~mt] > 0).sum(axis=1)).ravel()
    a.obs['whole_cell_nonmt_counts'] = nonmt_total
    a.obs['whole_cell_nonmt_genes'] = nonmt_genes
    included = a.obs['author_analysis_included'].to_numpy(dtype=bool)
    known = pd.to_numeric(a.obs['author_nCount.nMT'], errors='coerce').to_numpy()
    assert np.equal(nonmt_total[included], known[included]).all(), 'Recomputed non-MT counts differ from author metadata'
    eligible = included & (nonmt_total > 1000)
    a.obs['pass_author_depth_whole_cell'] = eligible
    a.obs['qc_joint_outlier_review'] = False
    a.obs['qc_high_complexity_flag'] = False
    bounds = []
    for sample, ids in a.obs.loc[eligible].groupby('sample_id', observed=True).groups.items():
        fit = ids
        fallback = False
        if len(fit) < 20:
            origin = a.obs.loc[ids[0], 'study_origin']
            fit = a.obs_names[eligible & a.obs.study_origin.eq(origin).to_numpy()]
            fallback = True
        metrics = {'mt': a.obs.loc[fit,'pct_counts_mt'].to_numpy(),
                   'genes': np.log1p(a.obs.loc[fit,'whole_cell_nonmt_genes'].to_numpy())}
        threshold = {}
        for key, values in metrics.items():
            center = float(np.median(values)); mad = float(np.median(np.abs(values-center)))
            if mad == 0:
                mad = float(np.subtract(*np.quantile(values, [.75, .25]))/1.349)
            threshold[key] = {'median': center, 'raw_MAD': mad, 'low': center-3*mad, 'high': center+3*mad}
        log_genes = np.log1p(a.obs.loc[ids,'whole_cell_nonmt_genes'].to_numpy())
        a.obs.loc[ids,'qc_joint_outlier_review'] = (a.obs.loc[ids,'pct_counts_mt'].to_numpy() > threshold['mt']['high']) & (log_genes < threshold['genes']['low'])
        a.obs.loc[ids,'qc_high_complexity_flag'] = log_genes > threshold['genes']['high']
        bounds.append({'sample': str(sample), 'fit_cells': len(fit), 'study_fallback': fallback, 'bounds': threshold})
    save_json(work/'OUTLIER_BOUNDS.json', bounds)
    query = a[eligible].copy()
    scores = programs(query.X, query.var.gene_symbol, nonmt_total[eligible])
    scores.index = query.obs_names
    query.obs = query.obs.join(scores)
    other = [p+'_program' for p in PANELS if p != 'Cardiomyocyte']
    query.obs['cm_program_dominates'] = query.obs.Cardiomyocyte_program.ge(query.obs[other].max(axis=1))
    refgenes = pd.read_csv(REF/'reference_genes.tsv', sep='\t')
    common = refgenes.loc[refgenes.gene_id.isin(query.var.gene_id) & ~refgenes.gene_symbol.fillna('').str.startswith('MT-'), 'gene_id'].tolist()
    positions = pd.Index(query.var.gene_id).get_indexer(common)
    assert len(common) > 5000 and (positions >= 0).all()
    with threadpool_limits(limits=2):
        mmwrite(work/'identity_query.mtx', query.X[:, positions].T)
    (work/'identity_genes.txt').write_text('\n'.join(common)+'\n')
    (work/'identity_cells.txt').write_text('\n'.join(query.obs_names)+'\n')
    with gzip.open(REF/'calibration_counts.mtx.gz', 'rb') as handle:
        calibration = mmread(handle).tocsr().T.tocsr()
    metadata = pd.read_csv(REF/'calibration_metadata.tsv', sep='\t')
    assert calibration.shape == (len(metadata), len(refgenes))
    total = np.asarray(calibration[:, ~refgenes.gene_symbol.fillna('').str.startswith('MT-').to_numpy()].sum(axis=1)).ravel()
    cal_scores = programs(calibration, refgenes.gene_symbol, total)
    mixed = {}
    for label in PANELS:
        values = cal_scores[label+'_program']
        positive = values[metadata.broad_type.eq(label).to_numpy()]
        negative = values[metadata.broad_type.ne(label).to_numpy()]
        supported = len(positive) >= 20 and len(negative) >= 20
        mixed[label] = {'positive_n': len(positive), 'negative_n': len(negative), 'calibration_supported': supported,
                        'positive_q10': float(positive.quantile(.1)) if supported else None,
                        'negative_q99': float(negative.quantile(.99)) if supported else None,
                        'threshold': float(max(positive.quantile(.1), negative.quantile(.99))) if supported else None}
    assert mixed['Cardiomyocyte']['calibration_supported']
    save_json(work/'MIXED_BOUNDS.json', mixed)
    query.write_h5ad(work/'author_depth_candidates.h5ad', compression='gzip')
    a.obs.to_csv(work/'all_cell_audit_before_identity.tsv.gz', sep='\t', index_label='cell_id')
    summary = a.obs.groupby(['disease','donor_id'], observed=True).agg(
        source_cells=('disease','size'), author_included=('author_analysis_included','sum'),
        author_depth_pass=('pass_author_depth_whole_cell','sum'), strict_baseline_pass=('pass_strict_numeric','sum'))
    write_tsv(summary.reset_index(), work/'STAGE_DONOR_COUNTS.tsv')
    save_json(work/'STAGE_SUCCESS.json', {'source_cells': a.n_obs, 'source_genes': a.n_vars,
              'candidate_cells': query.n_obs, 'common_nonmt_genes': len(common),
              'candidate_disease_counts': {str(k): int(v) for k,v in query.obs.disease.value_counts().items()}})
    print(json.dumps(json.loads((work/'STAGE_SUCCESS.json').read_text())), flush=True)


def finalize(work, output):
    assert not output.exists(), 'Use a fresh final output directory'
    assert (work/'IDENTITY_SUCCESS.json').exists()
    policy = json.loads((work/'POLICY.json').read_text())
    for name, value in policy['originals'].items():
        assert digest(Path(name)) == value, 'Source changed during review: '+name
    a = ad.read_h5ad(work/'author_depth_candidates.h5ad')
    predictions = pd.read_csv(work/'identity.tsv', sep='\t', index_col='cell_id', keep_default_na=False)
    assert set(predictions.index) == set(a.obs_names) and predictions.index.is_unique
    a.obs = a.obs.join(predictions)
    bounds = json.loads((work/'MIXED_BOUNDS.json').read_text())
    mixed = np.zeros(a.n_obs, dtype=bool)
    reason = np.full(a.n_obs, '', dtype=object)
    cm_score = a.obs.Cardiomyocyte_program.to_numpy()
    for label in PANELS:
        if label == 'Cardiomyocyte':
            continue
        other = a.obs[label+'_program'].to_numpy()
        flag = a.obs.Cardiomyocyte_genes_ge2.ge(3).to_numpy() & a.obs[label+'_genes_ge2'].ge(3).to_numpy()
        if bounds[label]['calibration_supported']:
            flag &= (cm_score >= bounds['Cardiomyocyte']['threshold']) & (other >= bounds[label]['threshold'])
            flag &= np.minimum(cm_score, other)/np.maximum(np.maximum(cm_score, other), 1e-12) >= .25
            prefix = ''
        else:
            flag &= other/np.maximum(cm_score, 1e-12) >= .5
            prefix = 'uncalibrated_review:'
        reason[flag & ~mixed] = prefix+'Cardiomyocyte+'+label
        mixed |= flag
    a.obs['coherent_mixed_lineage_review'] = mixed
    a.obs['mixed_lineage_reason'] = reason
    a.obs['qc_eligible_CM'] = (a.obs.reference_label.eq('Cardiomyocyte') & a.obs.reference_pruned.eq('Cardiomyocyte')
                              & a.obs.Cardiomyocyte_genes_ge2.ge(3) & a.obs.cm_program_dominates
                              & ~a.obs.coherent_mixed_lineage_review & ~a.obs.qc_joint_outlier_review)
    a.obs['review_action'] = np.where(a.obs.qc_eligible_CM, 'whole_cell_CM_core', 'identity_or_mixed_QC_review')
    a.obs['label'] = a.obs.disease.map({'Control': 0, 'HCM': 1}).astype(int)
    a.obs['split_status'] = 'unassigned_requires_donor_and_model_history_review'
    a.obs['training_use'] = 'supplementary_whole_cell_only_study_technology_disease_confounded'
    a.uns['qc_policy_json'] = json.dumps(policy)
    a.uns['whole_cell_mixed_bounds_json'] = json.dumps(bounds)
    a.uns['whole_cell_identity_method'] = 'SingleR on common non-MT genes; nuclear reference, domain shift unresolved; heuristic consensus, not calibrated posterior.'
    output.mkdir(parents=True)
    core = a[a.obs.qc_eligible_CM.to_numpy()].copy()
    pending = a[~a.obs.qc_eligible_CM.to_numpy()].copy()
    objects = {'author_depth_candidates.h5ad': a, 'cm_qc_full.h5ad': core,
               'cm_review_pending.h5ad': pending}
    manifests = []
    for filename, obj in objects.items():
        obj.write_h5ad(output/filename, compression='gzip')
        manifests.append({'file': filename, 'server_path': str(output/filename), 'cells': obj.n_obs,
                          'genes': obj.n_vars, 'bytes': (output/filename).stat().st_size,
                          'MD5': digest(output/filename, 'md5'), 'SHA256': digest(output/filename)})
    audit = pd.read_csv(work/'all_cell_audit_before_identity.tsv.gz', sep='\t', index_col='cell_id')
    audit['whole_cell_final_action'] = 'outside_author_depth_gate'
    audit.loc[a.obs_names, 'whole_cell_final_action'] = a.obs.review_action
    audit.to_csv(output/'all_cell_audit.tsv.gz', sep='\t', index_label='cell_id')
    a.obs.to_csv(output/'cell_identity_audit.tsv.gz', sep='\t', index_label='cell_id')
    donor = audit.groupby(['disease','donor_id'], observed=True).agg(source_cells=('disease','size'),
             author_included=('author_analysis_included','sum'), author_depth_candidates=('pass_author_depth_whole_cell','sum'),
             old_strict_numeric=('pass_strict_numeric','sum'))
    donor['CM_core'] = core.obs.groupby(['disease','donor_id'], observed=True).size().reindex(donor.index, fill_value=0)
    donor['identity_QC_review'] = pending.obs.groupby(['disease','donor_id'], observed=True).size().reindex(donor.index, fill_value=0)
    donor['CM_fraction_of_source'] = donor.CM_core/donor.source_cells
    write_tsv(donor.reset_index(), output/'DONOR_RETENTION.tsv')
    sample = audit.groupby(['disease','donor_id','sample_id'], observed=True).agg(source_cells=('disease','size'), author_depth_candidates=('pass_author_depth_whole_cell','sum'))
    sample['CM_core'] = core.obs.groupby(['disease','donor_id','sample_id'], observed=True).size().reindex(sample.index, fill_value=0)
    write_tsv(sample.reset_index(), output/'SAMPLE_RETENTION.tsv')
    write_tsv(pd.DataFrame(manifests), output/'OUTPUT_MANIFEST.tsv')
    (output/'MD5SUMS').write_text(''.join(r['MD5']+'  '+r['file']+'\n' for r in manifests))
    for filename in ['POLICY.json','OUTLIER_BOUNDS.json','MIXED_BOUNDS.json','STAGE_SUCCESS.json','IDENTITY_SUCCESS.json']:
        (output/filename).write_bytes((work/filename).read_bytes())
    summary = {'exported_utc': datetime.now(timezone.utc).isoformat(), 'source_cells': len(audit),
               'author_depth_candidates': a.n_obs, 'core_cells': core.n_obs, 'pending_cells': pending.n_obs,
               'core_disease_cells': {str(k): int(v) for k,v in core.obs.disease.value_counts().items()},
               'core_donors': {str(k): int(v) for k,v in core.obs.groupby('disease', observed=True).donor_id.nunique().items()},
               'source_donors': 12, 'genes': core.n_vars, 'mixed_review_candidates': int(mixed.sum()),
               'joint_numeric_outliers': int(a.obs.qc_joint_outlier_review.sum()),
               'high_complexity_flag_not_automatic_doublet': int(a.obs.qc_high_complexity_flag.sum()),
               'old_strict_core_cells': 6, 'source_counts_preserved': True, 'split': 'unassigned',
               'training_or_DGE_run': False, 'tokenization_run': False,
               'limitation': policy['use'], 'verification_status': 'pending_independent_verifier'}
    save_json(output/'EXPORT_SUMMARY.json', summary)
    print(json.dumps(summary), flush=True)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('mode', choices=['stage','finalize'])
    parser.add_argument('--work', type=Path, required=True)
    parser.add_argument('--output', type=Path)
    args = parser.parse_args()
    if args.mode == 'stage':
        stage(args.work)
    else:
        assert args.output is not None
        finalize(args.work, args.output)
