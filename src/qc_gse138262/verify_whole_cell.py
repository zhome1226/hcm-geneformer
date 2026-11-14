"""Independently compare the new whole-cell exports to the complete source counts."""
import argparse
from datetime import datetime, timezone
import hashlib
import json
import os
from pathlib import Path

os.environ.setdefault('OPENBLAS_NUM_THREADS', '2')
import anndata as ad
import numpy as np
import pandas as pd


def sha(path):
    with path.open('rb') as f:
        return hashlib.file_digest(f, 'sha256').hexdigest()


def check(condition, name, checks):
    assert bool(condition), name
    checks.append(name)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    root = args.output
    assert not (root/'FINAL_COMPLETION_GATE.json').exists()
    policy = json.loads((root/'POLICY.json').read_text())
    checks = []
    for path, digest in policy['originals'].items():
        check(sha(Path(path)) == digest, 'source_and_old_strict_file_unchanged_'+Path(path).name, checks)
    source_path = next(Path(p) for p in policy['originals'] if Path(p).name == 'all_source_counts.h5ad')
    source = ad.read_h5ad(source_path)
    x = source.X.tocsr()
    check(source.obs_names.is_unique and source.var_names.is_unique, 'source_axes_unique', checks)
    mt = source.var.gene_symbol.str.startswith('MT-').to_numpy()
    check(np.array_equal(mt, source.var.mt.to_numpy()), 'source_MT_annotation_matches_symbols', checks)
    nonmt_counts = np.asarray(x[:, ~mt].sum(axis=1)).ravel()
    expected = source.obs.author_analysis_included.to_numpy(dtype=bool) & (nonmt_counts > 1000)
    expected_ids = source.obs_names[expected]
    objects = {}
    for filename in ['author_depth_candidates.h5ad','cm_qc_full.h5ad','cm_review_pending.h5ad']:
        a = ad.read_h5ad(root/filename)
        objects[filename] = a
        prefix = filename+':'
        check(a.obs_names.is_unique and a.var_names.equals(source.var_names), prefix+'axes_unique_full_source_genes', checks)
        rows = source.obs_names.get_indexer(a.obs_names)
        check((rows >= 0).all(), prefix+'all_cells_real_source_cells', checks)
        difference = (a.X - x[rows]).tocsr(); difference.eliminate_zeros()
        check(difference.nnz == 0, prefix+'every_count_matches_source', checks)
        difference = (a.layers['counts'] - x[rows]).tocsr(); difference.eliminate_zeros()
        check(difference.nnz == 0, prefix+'every_counts_layer_value_matches_source', checks)
        check(np.isfinite(a.X.data).all() and (a.X.data >= 0).all() and np.equal(a.X.data, np.rint(a.X.data)).all(), prefix+'all_stored_counts_finite_integer_nonnegative', checks)
        for col in ['donor_id','sample_id','disease','author_cell_id','technology','study_origin']:
            check(np.array_equal(a.obs[col].astype(str),source.obs.iloc[rows][col].astype(str)), prefix+'source_metadata_'+col, checks)
        check(a.obs.label.eq(a.obs.disease.map({'Control':0,'HCM':1})).all(), prefix+'label_diagnosis_match', checks)
        check((a.obs.groupby('donor_id', observed=True).disease.nunique() == 1).all(), prefix+'one_diagnosis_per_donor', checks)
        mask = 'native_feature_presence_by_sample'
        check(np.array_equal(a.varm[mask], source.varm[mask]), prefix+'full_native_feature_missingness_mask_preserved', checks)
        check(np.array_equal(a.uns['feature_presence_sample_order'],source.uns['feature_presence_sample_order']), prefix+'native_mask_sample_order_preserved', checks)
        totals = np.asarray(a.X.sum(axis=1)).ravel()
        mt_counts = np.asarray(a.X[:,mt].sum(axis=1)).ravel()
        check(np.allclose(a.obs.pct_counts_mt, 100*mt_counts/totals, atol=1e-9), prefix+'MT_from_complete_counts', checks)
        check(a.obs.whole_cell_nonmt_counts.gt(1000).all() and a.obs.author_analysis_included.all(), prefix+'author_depth_gate_exact', checks)
        check(a.obs.split_status.eq('unassigned_requires_donor_and_model_history_review').all(), prefix+'no_split_assigned', checks)
    candidates = objects['author_depth_candidates.h5ad']
    core = objects['cm_qc_full.h5ad']
    pending = objects['cm_review_pending.h5ad']
    check(candidates.obs_names.equals(expected_ids), 'candidate_ids_exact_author_depth_rule', checks)
    check(len(set(core.obs_names) & set(pending.obs_names)) == 0, 'core_review_disjoint', checks)
    check(set(core.obs_names) | set(pending.obs_names) == set(expected_ids), 'core_review_partition_all_candidates', checks)
    independent_core = (candidates.obs.reference_label.eq('Cardiomyocyte') & candidates.obs.reference_pruned.eq('Cardiomyocyte')
                        & candidates.obs.Cardiomyocyte_genes_ge2.ge(3) & candidates.obs.cm_program_dominates
                        & ~candidates.obs.coherent_mixed_lineage_review & ~candidates.obs.qc_joint_outlier_review)
    check(core.obs_names.equals(candidates.obs_names[independent_core]), 'core_exact_frozen_identity_rule', checks)
    check(core.obs.qc_eligible_CM.all() and not pending.obs.qc_eligible_CM.any(), 'core_review_labels_match', checks)
    check(core.n_obs > 0 and core.obs.disease.nunique() == 2, 'core_has_real_HCM_and_control_cells', checks)
    check(core.obs.groupby('disease',observed=True).donor_id.nunique().min() >= 2, 'at_least_two_donors_per_group', checks)
    check(core.n_vars == source.n_vars == 28663 and core.var.gene_id.is_unique, 'full_gene_axis_not_HVG_or_normalized_subset', checks)
    table = pd.read_csv(root/'DONOR_RETENTION.tsv', sep='\t')
    check(len(table) == 12 and table.source_cells.sum() == source.n_obs, 'all_source_donors_in_retention_table', checks)
    check(table.CM_core.sum() == core.n_obs, 'retention_table_matches_export', checks)
    audit = pd.read_csv(root/'all_cell_audit.tsv.gz', sep='\t', index_col='cell_id')
    check(audit.index.equals(source.obs_names), 'all_7537_source_cells_auditable', checks)
    for path, digest in policy['originals'].items():
        check(sha(Path(path)) == digest, 'source_unchanged_after_verification_'+Path(path).name, checks)
    report = {'checked_utc': datetime.now(timezone.utc).isoformat(), 'checks_passed': len(checks), 'checks': checks,
              'source_cells': source.n_obs, 'candidate_cells': candidates.n_obs, 'core_cells': core.n_obs,
              'pending_cells': pending.n_obs, 'core_disease_cells': {str(k):int(v) for k,v in core.obs.disease.value_counts().items()},
              'core_donors': {str(k):int(v) for k,v in core.obs.groupby('disease',observed=True).donor_id.nunique().items()},
              'source_and_old_strict_preserved': True, 'full_output_count_comparison_to_source': True,
              'independence_or_training_approval': False, 'QC_domain_shift': 'nuclear reference applied to whole cells; heuristic consensus, no validated sensitivity/specificity',
              'training_DGE_tokenization_run': False}
    (root/'INDEPENDENT_VERIFICATION.json').write_text(json.dumps(report, indent=2)+'\n')
    summary = json.loads((root/'EXPORT_SUMMARY.json').read_text())
    summary.update(verification_status='passed', verification_checks_passed=len(checks),
                   verification_report='INDEPENDENT_VERIFICATION.json')
    (root/'EXPORT_SUMMARY.json').write_text(json.dumps(summary, indent=2)+'\n')
    (root/'FINAL_COMPLETION_GATE.json').write_text(json.dumps({'complete': True, 'verified_UTC': report['checked_utc'],
                 'checks_passed': len(checks), 'scope': 'source-conserving supplementary whole-cell QC export; not independent disease validation'},indent=2)+'\n')
    print(json.dumps({k:v for k,v in report.items() if k != 'checks'}), flush=True)


if __name__ == '__main__':
    main()
