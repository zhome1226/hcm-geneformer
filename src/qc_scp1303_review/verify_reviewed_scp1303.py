"""Independently reopen exports; verify every retained row, count layer and QC invariant."""
from datetime import datetime, timezone
from pathlib import Path
import argparse
import hashlib
import json

import anndata as ad
import numpy as np
import pandas as pd

root = Path.home() / 'xiaoe-internship-data'
parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument('--output-name', default='SCP1303_qc_reviewed_20261111')
args = parser.parse_args()
old = root / 'processed/SCP1303_qc_20261111'
out = root / 'processed' / args.output_name
summary = json.loads((out / 'EXPORT_SUCCESS.json').read_text())
policy = json.loads((out / 'qc_policy.json').read_text())
decisions = pd.read_csv(out / 'review_decisions.tsv.gz', sep='\t', index_col='barcode', low_memory=False)
baseline = ad.read_h5ad(old / 'cm_qc_full.h5ad')
core = ad.read_h5ad(out / 'cm_qc_full.h5ad')
light = ad.read_h5ad(out / 'cm_qc_light.h5ad')
pending = ad.read_h5ad(out / 'cm_review_pending.h5ad')
checks = {}
def check(name, condition):
    checks[name] = bool(condition)
    if not checks[name]:
        raise AssertionError(name)

check('declared_strict_limits', policy['min_genes']==500 and policy['max_genes']==5000 and policy['mt_percent_exclusive_max']==5)
check('decision_barcodes_unique', decisions.index.is_unique)
check('all_baseline_and_high_candidates_covered', len(decisions)==114927+5043 and int(decisions['eligible_strict_5000'].sum())==114927)
expected_core = set(decisions.index[decisions['kept_reviewed_core']])
expected_pending = set(decisions.index[decisions['eligible_strict_5000'] & decisions['review_identity_state'].eq('pending')])
check('core_matches_decision_set', set(core.obs['barcode_original'].astype(str)) == expected_core)
check('pending_matches_decision_set', set(pending.obs['barcode_original'].astype(str)) == expected_pending)
check('pending_disjoint_core', not set(pending.obs_names).intersection(core.obs_names))
check('all_original_genes_preserved', core.var_names.equals(baseline.var_names) and light.var_names.equals(baseline.var_names) and pending.var_names.equals(baseline.var_names))
for name, obj in [('full', core), ('light', light), ('pending', pending)]:
    loc = baseline.obs_names.get_indexer(obj.obs_names)
    check(name + '_unique_ids', obj.obs_names.is_unique and obj.obs['barcode_original'].is_unique)
    check(name + '_source_members', bool((loc >= 0).all()))
    original = baseline[loc]
    for field in ['barcode_original', 'donor_id', 'sample_id', 'disease', 'label', 'split']:
        check(name + '_source_metadata_' + field, np.array_equal(obj.obs[field].astype(str).to_numpy(), original.obs[field].astype(str).to_numpy()))
    for key in ['X', 'counts', 'cellranger_raw']:
        matrix = obj.X if key=='X' else obj.layers[key]
        source = original.X if key=='X' else original.layers[key]
        check(name + '_' + key + '_all_values_match_original', (matrix != source).nnz==0)
        check(name + '_' + key + '_valid_values', np.isfinite(matrix.data).all() and (matrix.data>=0).all())
    check(name + '_x_counts_identical', (obj.X != obj.layers['counts']).nnz==0)
    genes = np.asarray((obj.X>0).sum(axis=1)).ravel()
    raw = obj.layers['cellranger_raw']
    totals = np.asarray(raw.sum(axis=1)).ravel()
    mito = np.asarray(raw[:, obj.var_names.str.startswith('MT-')].sum(axis=1)).ravel()*100/np.maximum(totals, 1)
    check(name + '_strict_gene_limits', ((genes>=500)&(genes<=5000)).all())
    check(name + '_strict_mt_limit', (mito<5).all())
    check(name + '_qc_metadata', np.array_equal(genes, obj.obs['n_genes_by_counts']) and np.allclose(mito, obj.obs['pct_counts_mt'], atol=1e-6))
    check(name + '_disease_labels', set(obj.obs['disease']).issubset({'NF','HCM'}) and np.array_equal(obj.obs['disease'].astype(str).map({'NF':0,'HCM':1}), obj.obs['label']))
    check(name + '_donor_split_consistent', obj.obs.groupby('donor_id', observed=True)['split'].nunique().max()<=1)
check('full_31_donors', core.obs['donor_id'].nunique()==31)
check('light_20000_unique_31_donors', light.n_obs==20000 and light.obs['donor_id'].nunique()==31)
check('light_balanced_labels', light.obs['label'].value_counts().to_dict()=={0:10000,1:10000})
for key in ['X', 'counts', 'cellranger_raw']:
    part = core[light.obs_names]
    a = part.X if key=='X' else part.layers[key]
    b = light.X if key=='X' else light.layers[key]
    check('light_full_exact_' + key, (a!=b).nnz==0)
check('core_contains_only_retained_states', core.obs['review_identity_state'].astype(str).str.startswith('retain_').all())
check('pending_contains_only_pending_states', pending.obs['review_identity_state'].eq('pending').all())
doublet_exclusions = decisions.loc[decisions['action'].eq('exclude_suspected_doublet')]
check('doublet_exclusions_have_secondary_call_and_original_warning', doublet_exclusions['scDblFinder_class'].eq('doublet').all() and doublet_exclusions['marker_mixed_warning'].all())
check('high_gene_candidates_not_silently_admitted', decisions.loc[~decisions['eligible_strict_5000'], 'action'].eq('exclude_project_upper_limit').all())
hashes = {}
for name, entry in summary['hashes'].items():
    path = out / name
    md5, sha = hashlib.md5(), hashlib.sha256()
    with path.open('rb') as handle:
        while chunk := handle.read(32*1024**2):
            md5.update(chunk); sha.update(chunk)
    hashes[name] = {'bytes': path.stat().st_size, 'md5': md5.hexdigest(), 'sha256': sha.hexdigest()}
    check(name + '_size_and_export_sha256', hashes[name]['bytes']==entry['bytes'] and hashes[name]['sha256']==entry['sha256'])
for name, expected in {'cm_qc_full.h5ad':'16456fea621116003b7f6e35039e2d7e', 'cm_qc_light.h5ad':'5517c21c2595bc42394bb1af3f09a949'}.items():
    with (old / name).open('rb') as handle:
        actual = hashlib.file_digest(handle, 'md5').hexdigest()
    check('original_baseline_unchanged_' + name, actual==expected)
for algorithm in ['md5', 'sha256']:
    (out / (algorithm.upper()+'SUMS')).write_text(''.join(f'{entry[algorithm]}  {name}\n' for name, entry in hashes.items()))
report = {'status': 'passed', 'verified_utc': datetime.now(timezone.utc).isoformat(),
    'checks': checks, 'hashes': hashes, 'core_by_disease': core.obs['disease'].value_counts().to_dict(),
    'core_nuclei': core.n_obs, 'light_nuclei': light.n_obs, 'pending_nuclei': pending.n_obs,
    'verifier_sha256': hashlib.sha256(Path(__file__).read_bytes()).hexdigest()}
(out / 'INDEPENDENT_VERIFICATION.json').write_text(json.dumps(report, indent=2) + '\n')
print(json.dumps({key:report[key] for key in ['status','core_nuclei','light_nuclei','pending_nuclei','core_by_disease']}, indent=2), flush=True)
