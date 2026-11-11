"""Freeze evidence-based review decisions and export a separate strict QC version."""
from datetime import datetime, timezone
from pathlib import Path
import argparse
import gc
import hashlib
import json

import anndata as ad
import h5py
import numpy as np
import pandas as pd

root = Path.home() / 'xiaoe-internship-data'
parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument('--output-name', default='SCP1303_qc_reviewed_20261111')
parser.add_argument('--use-group-evidence', action='store_true')
args = parser.parse_args()
if '/' in args.output_name or '\\' in args.output_name:
    raise ValueError('Output name must be a single directory component')
old = root / 'processed/SCP1303_qc_20261111'
review = root / 'processed/SCP1303_review_20261111/identity'
out = root / 'processed' / args.output_name
out.mkdir(exist_ok=True)
if (out / 'EXPORT_SUCCESS.json').exists():
    raise RuntimeError('Completed export exists; verify it instead of overwriting')
audit = pd.read_csv(old / 'reports/cell_qc_audit.tsv.gz', sep='\t', index_col='barcode', low_memory=False)
e = pd.read_csv(review / 'joined_identity_evidence.tsv.gz', sep='\t', index_col='barcode', low_memory=False)
doublets = pd.read_csv(review / 'doublets/secondary_doublet_predictions.tsv.gz', sep='\t', index_col='barcode')
assert doublets.index.is_unique
e = e.join(doublets[['scDblFinder_score', 'scDblFinder_class']])
is_cm = e['cell_type_original'].astype(str).str.startswith('Cardiomyocyte')
controls = e['review_context_anchor'] & is_cm & e['reference_labels'].eq('Cardiomyocyte') & e['cardiomyocyte_markers_detected'].ge(3)
neighbor_floor = float(e.loc[controls, 'author_cm_neighbor_fraction'].quantile(.05))
noncm = ['fibroblast', 'endothelial', 'immune', 'pericyte', 'smooth_muscle', 'neuronal', 'adipocyte', 'mesothelial']
nc_fraction = e[[x + '_cellbender_fraction' for x in noncm]].max(axis=1)
cm_dominant = e['cardiomyocyte_cellbender_fraction'].ge(nc_fraction) & e['cardiomyocyte_cellbender_fraction'].gt(0)
neighborhood = e['author_cm_neighbor_fraction'].ge(neighbor_floor)
reference_weak = ~e['reference_empirical_reference_tail_pass']
supported = is_cm & (e['independent_cm_support'] |
    (e['reference_labels'].eq('Cardiomyocyte') & e['cardiomyocyte_markers_detected'].ge(3) & cm_dominant & neighborhood) |
    (reference_weak & e['cardiomyocyte_markers_detected'].ge(5) & cm_dominant & neighborhood))
group_supported = pd.Series(False, index=e.index)
group_calibration = None
if args.use_group_evidence:
    groupdir = review / 'group_review'
    gp = pd.read_csv(groupdir / 'query_group_predictions.tsv', sep='\t').set_index('profile_id')
    gm = pd.read_csv(groupdir / 'query_group_metadata.tsv', sep='\t', dtype={'cluster':str}).set_index('profile_id')
    cp = pd.read_csv(groupdir / 'calibration_group_predictions.tsv', sep='\t').set_index('profile_id')
    cm = pd.read_csv(groupdir / 'calibration_group_metadata.tsv', sep='\t').set_index('profile_id')
    cp = cm.join(cp)
    true_cm, pred_cm = cp['truth'].eq('Cardiomyocyte'), cp['label'].eq('Cardiomyocyte')
    group_calibration = {'reference_holdout_profiles': len(cp),
        'accuracy': float(cp['truth'].eq(cp['label']).mean()),
        'CM_truth_profiles': int(true_cm.sum()),
        'CM_precision': float((true_cm&pred_cm).sum()/max(1,pred_cm.sum())),
        'CM_recall': float((true_cm&pred_cm).sum()/max(1,true_cm.sum()))}
    if group_calibration['CM_precision']<.95 or group_calibration['CM_recall']<.95:
        raise ValueError('Group-level cardiomyocyte reference applicability gate failed')
    groups = gm.join(gp)
    cluster_records = []
    for cluster, part in groups.groupby('cluster', sort=True):
        aggregate = part.loc[part['donor_id'].eq('ALL')]
        donors = part.loc[~part['donor_id'].eq('ALL')]
        fraction = float(donors.loc[donors['label'].eq('Cardiomyocyte'),'nuclei'].sum()/max(1,donors['nuclei'].sum()))
        cluster_records.append({'cluster':str(cluster), 'group_label':str(aggregate.iloc[0]['label']),
            'donor_profile_count':len(donors), 'CM_donor_profile_nucleus_fraction':fraction})
    cluster_table = pd.DataFrame(cluster_records).set_index('cluster')
    cluster_table.to_csv(out / 'group_identity_consistency.tsv', sep='\t')
    program = pd.read_csv(groupdir / 'query_CM_program.tsv', sep='\t', index_col='barcode')
    e = e.join(program)
    e['group_label'] = e['leiden_1.0'].astype(str).map(cluster_table['group_label'])
    e['group_CM_donor_fraction'] = e['leiden_1.0'].astype(str).map(cluster_table['CM_donor_profile_nucleus_fraction'])
    e['group_donor_profile_count'] = e['leiden_1.0'].astype(str).map(cluster_table['donor_profile_count'])
    group_supported = is_cm & e['reference_labels'].eq('Cardiomyocyte') & ~e['independent_noncm_support']
    group_supported &= e['reference_CM_program_genes_detected'].ge(3) & e['group_label'].eq('Cardiomyocyte')
    group_supported &= e['group_CM_donor_fraction'].ge(.8) & e['group_donor_profile_count'].ge(3)
    supported |= group_supported
    (out / 'GROUP_REFERENCE_APPLICABILITY.json').write_text(json.dumps(group_calibration, indent=2)+'\n')
matching_panel = {'Fibroblast': 'fibroblast', 'Endothelial': 'endothelial', 'Immune': 'immune',
    'Pericyte': 'pericyte', 'Smooth_muscle': 'smooth_muscle', 'Neural': 'neuronal',
    'Adipocyte': 'adipocyte', 'Mesothelial': 'mesothelial'}
matching_detected = np.zeros(len(e), dtype=int)
matching_fraction = np.zeros(len(e), dtype=float)
for label, panel in matching_panel.items():
    mask = e['reference_labels'].eq(label).to_numpy()
    matching_detected[mask] = e.loc[mask, panel + '_markers_detected'].to_numpy()
    matching_fraction[mask] = e.loc[mask, panel + '_cellbender_fraction'].to_numpy()
noncm_supported = is_cm & e['independent_noncm_support'] & (matching_detected >= 2)
noncm_supported &= e['cardiomyocyte_markers_detected'].le(1) & (matching_fraction > e['cardiomyocyte_cellbender_fraction'])
residual_doublet = is_cm & e['scDblFinder_class'].eq('doublet') & e['coherent_mixed_panel_flag']
residual_doublet &= e['author_doublet_capture_percentile'].ge(.95) & e['marker_mixed_warning']
state = pd.Series('context_only', index=e.index)
reason = pd.Series('non-target context for review', index=e.index)
state.loc[is_cm] = 'pending'
reason.loc[is_cm] = 'identity evidence insufficient or conflicting; kept separately'
state.loc[supported] = 'retain_supported_cm'
reason.loc[supported] = 'author identity plus applicable reference or extended cardiomyocyte program and calibrated neighbor context'
reason.loc[group_supported] = 'individual CM best reference label plus >=3 independent-reference-derived CM program genes and consistent multi-donor CM group; sparse per-cell panel not required'
state.loc[noncm_supported] = 'exclude_noncm'
reason.loc[noncm_supported] = 'reference-tail-supported non-CM type plus matching lineage program, minimal cardiomyocyte panel, and non-CM count dominance'
state.loc[residual_doublet] = 'exclude_suspected_doublet'
reason.loc[residual_doublet] = 'scDblFinder call plus original high capture-relative score and coherent mixed programs; not experimental confirmation'
e['review_identity_state'] = state
e['review_reason'] = reason
e.to_csv(review / 'final_joined_evidence.tsv.gz', sep='\t', index_label='barcode', compression='gzip')

eligible_high = audit['pass_target_disease'] & audit['pass_cm_annotation'] & audit['pass_author_subcluster']
eligible_high &= audit['pass_min_genes'] & audit['pass_mito'] & audit['n_genes_by_counts'].gt(5000)
target = audit['kept_full'] | eligible_high
decisions = audit.loc[target, ['donor_id', 'sample_id', 'capture_id', 'disease', 'label',
    'cell_type_original', 'SubCluster', 'n_genes_by_counts', 'pct_counts_mt', 'doublet_score',
    'marker_mixed_warning', 'n_cm_markers_detected', 'kept_full', 'selected_light', 'split']].copy()
decisions['review_identity_state'] = 'retain_author_qc_no_new_signal'
decisions['review_reason'] = 'author QC reused; no registered new identity concern'
reviewed = decisions.index.intersection(e.index[is_cm])
decisions.loc[reviewed, 'review_identity_state'] = e.loc[reviewed, 'review_identity_state']
decisions.loc[reviewed, 'review_reason'] = e.loc[reviewed, 'review_reason']
decisions['eligible_strict_5000'] = decisions['kept_full']
decisions['eligible_upper_6000'] = decisions['n_genes_by_counts'].le(6000)
decisions['kept_reviewed_core'] = decisions['eligible_strict_5000'] & decisions['review_identity_state'].str.startswith('retain_')
decisions['action'] = decisions['review_identity_state']
high = ~decisions['eligible_strict_5000']
decisions.loc[high, 'action'] = 'exclude_project_upper_limit'
decisions.loc[high, 'review_reason'] = 'above unchanged project upper limit; identity state separately recorded: ' + decisions.loc[high, 'review_identity_state']
decisions['scDblFinder_score'] = doublets['scDblFinder_score'].reindex(decisions.index)
decisions['scDblFinder_class'] = doublets['scDblFinder_class'].reindex(decisions.index).fillna('not_rerun_author_qc_reused')
expected_warning = audit.index[audit['kept_full'] & (audit['marker_mixed_warning'] | audit['n_cm_markers_detected'].eq(0))]
assert len(expected_warning) == 987 and expected_warning.isin(decisions.index).all()
assert int(eligible_high.sum()) == 5043 and audit.index[eligible_high].isin(decisions.index).all()
policy = {'version': 'SCP1303_reviewed_strict_v3' if args.use_group_evidence else 'SCP1303_reviewed_strict_v2', 'created_utc': datetime.now(timezone.utc).isoformat(),
    'seed': 20261111, 'min_genes': 500, 'max_genes': 5000, 'mt_percent_exclusive_max': 5,
    'count_layer': 'author CellBender X and counts; original CellRanger retained',
    'disease_map': {'NF': 0, 'HCM': 1}, 'reference': 'Litvinukova2020 adult left ventricular nuclei; SingleR with reference-donor holdout calibration',
    'neighbor_support_floor': neighbor_floor,
    'neighbor_floor_derivation': 'Fifth percentile of author-CM unflagged context anchors with >=3 expanded CM genes and CM as best reference label; contextual support, not independent confidence.',
    'retain_rule': 'Applicable reference-tail CM support OR CM best reference label + >=3 expanded CM genes + CM program count dominance + neighbor support OR uncertain reference + >=5 expanded CM genes + same dominance/neighbors; author CM required.',
    'noncm_rule': 'Reference-tail non-CM support + >=2 matching lineage genes + <=1 CM panel gene + matching non-CM counts dominate CM counts.',
    'doublet_rule': 'scDblFinder doublet + original sample-relative score >=95th percentile + >=2 CM and >=3 non-CM genes + original mixed marker warning; no single-tool-only removal.',
    'uncertain_rule': 'Other reviewed author cardiomyocytes remain pending and outside reviewed core, retained in a separate H5AD.',
    'unchallenged_cells': 'Reuse author QC with explicit status; not claimed to have independent per-cell reference confirmation.',
    'upper_limit_decision': 'Maintain user-confirmed 5000 primary limit. 6000 is an identity/coverage sensitivity report only, not silently adopted.',
    'split_status': 'Original donor split preserved; provisional pending final training protocol',
    'secondary_doublet_rate_assumption': 'Residual dbr=.01, dbr.sd=.015 on author-retained complete library cell-type composition; original droplets/doublet rate not reconstructed.'}
if args.use_group_evidence:
    policy['group_reference_applicability'] = group_calibration
    policy['additional_sparse_CM_rule'] = 'Individual best reference label CM + >=3 genes from a CM-specific program derived only from training-reference counts + CM global cluster and >=80% weighted donor-profile CM consistency across >=3 donor profiles; non-CM/doublet exclusions override. Group profiles are correlated evidence, not independent votes.'
    policy['program_source_sha256'] = hashlib.sha256((review / 'group_review/reference_CM_program.tsv').read_bytes()).hexdigest()
(out / 'qc_policy.json').write_text(json.dumps(policy, indent=2) + '\n')
decisions.to_csv(out / 'review_decisions.tsv.gz', sep='\t', index_label='barcode', compression='gzip')
decisions.groupby(['disease', 'donor_id', 'cell_type_original', 'action'], observed=True).agg(nuclei=('label', 'size')).to_csv(out / 'review_summary.tsv', sep='\t')
impact = []
for upper in [5000, 6000]:
    for (disease, donor), group in decisions.groupby(['disease', 'donor_id'], observed=True):
        keep = group['n_genes_by_counts'].le(upper) & group['review_identity_state'].str.startswith('retain_')
        impact.append({'upper_limit': upper, 'disease': disease, 'donor_id': donor, 'eligible_identity_supported_nuclei': int(keep.sum())})
pd.DataFrame(impact).to_csv(out / 'upper_limit_sensitivity.tsv', sep='\t', index=False)

print('Loading original strict object for count-preserving export', flush=True)
full = ad.read_h5ad(old / 'cm_qc_full.h5ad')
barcodes = full.obs['barcode_original'].astype(str)
assert set(barcodes) == set(decisions.index[decisions['eligible_strict_5000']])
for column in ['review_identity_state', 'review_reason', 'action', 'kept_reviewed_core', 'scDblFinder_score', 'scDblFinder_class']:
    full.obs[column] = decisions[column].reindex(barcodes).to_numpy()
full.uns['qc_review_policy_json'] = json.dumps(policy)
full.uns['author_qc_reused'] = 'Original author processing reused; selective independent identity and residual-doublet review recorded per cell.'
core = full[full.obs['kept_reviewed_core'].to_numpy()].copy()
pending = full[full.obs['review_identity_state'].eq('pending').to_numpy()].copy()
with h5py.File(old / 'cm_qc_light.h5ad') as handle:
    old_light = ad.io.read_elem(handle['obs'])
rng = np.random.default_rng(20261111)
chosen = [name for name in old_light.index if name in core.obs_names]
chosen_set = set(chosen)
replacements = []
for removed in old_light.index:
    if removed in chosen_set:
        continue
    original = old_light.loc[removed]
    available = core.obs.loc[~core.obs.index.isin(chosen_set)]
    pool = available.loc[available['donor_id'].astype(str).eq(str(original['donor_id'])) &
        available['cell_type_original'].astype(str).eq(str(original['cell_type_original']))]
    match_level = 'same_donor_and_author_subtype'
    if len(pool) == 0:
        pool = available.loc[available['donor_id'].astype(str).eq(str(original['donor_id']))]
        match_level = 'same_donor'
    if len(pool) == 0:
        pool = available.loc[available['disease'].astype(str).eq(str(original['disease']))]
        match_level = 'same_disease_donor_loss_reported'
    if len(pool) == 0:
        raise ValueError('Cannot replace removed light nucleus without duplication')
    replacement = str(rng.choice(pool.index.to_numpy()))
    chosen.append(replacement); chosen_set.add(replacement)
    replacements.append({'removed': removed, 'replacement': replacement, 'match_level': match_level,
        'donor_id_original': str(original['donor_id']), 'donor_id_replacement': str(core.obs.loc[replacement, 'donor_id'])})
assert len(chosen) == 20000 and len(chosen_set) == 20000
core.obs['selected_light'] = core.obs.index.isin(chosen_set)
light = core[chosen].copy()
light.uns['sampling'] = {'seed': 20261111, 'target': 20000, 'actual': 20000,
    'method': 'Retain eligible original light nuclei; replace exclusions within donor and author subtype where possible; no duplication'}
pd.DataFrame(replacements, columns=['removed', 'replacement', 'match_level', 'donor_id_original', 'donor_id_replacement']).to_csv(out / 'light_replacements.tsv', sep='\t', index=False)
for name, data in [('cm_qc_full.h5ad', core), ('cm_qc_light.h5ad', light), ('cm_review_pending.h5ad', pending)]:
    print(f'Writing {name}: {data.shape}', flush=True)
    temporary = out / (name + '.tmp')
    data.write_h5ad(temporary, compression='gzip', compression_opts=2)
    temporary.replace(out / name)
core.obs.groupby(['disease', 'donor_id', 'cell_type_original'], observed=True).size().to_csv(out / 'core_coverage.tsv', sep='\t', header=['nuclei'])
summary = {'status': 'exported_pending_independent_verification', 'completed_utc': datetime.now(timezone.utc).isoformat(),
    'full_nuclei': core.n_obs, 'light_nuclei': light.n_obs, 'genes': core.n_vars,
    'full_by_disease': core.obs['disease'].value_counts().to_dict(),
    'pending_nuclei': pending.n_obs, 'core_donors': int(core.obs['donor_id'].nunique()),
    'light_replacements': len(replacements), 'actions_strict_baseline': decisions.loc[decisions['eligible_strict_5000'], 'action'].value_counts().to_dict(),
    'actions_high_gene_identity': decisions.loc[~decisions['eligible_strict_5000'], 'review_identity_state'].value_counts().to_dict(),
    'candidate_decision_rows': len(decisions), 'hashes': {}}
for path in out.glob('*.h5ad'):
    with path.open('rb') as handle:
        digest = hashlib.file_digest(handle, 'sha256').hexdigest()
    summary['hashes'][path.name] = {'bytes': path.stat().st_size, 'sha256': digest}
(out / 'EXPORT_SUCCESS.json').write_text(json.dumps(summary, indent=2) + '\n')
print(json.dumps(summary, indent=2), flush=True)
