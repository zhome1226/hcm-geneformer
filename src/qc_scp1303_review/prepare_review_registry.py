"""Build review membership, targeted marker evidence and an unchanged count subset."""
from datetime import datetime, timezone
import gc
import gzip
import hashlib
import json
from pathlib import Path

import anndata as ad
import h5py
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.io import mmwrite

root = Path.home()/'xiaoe-internship-data'
work = root/'_processing/SCP1303_review_20261111'
out = root/'processed/SCP1303_review_20261111'
out.mkdir(parents=True, exist_ok=True)
seed = 20261111
print('Reading existing audit and output membership', flush=True)
o = pd.read_csv(root/'processed/SCP1303_qc_20261111/reports/cell_qc_audit.tsv.gz',
                sep='\t', index_col='barcode', low_memory=False)
for name in ['full', 'light']:
    with h5py.File(root/f'processed/SCP1303_qc_20261111/cm_qc_{name}.h5ad') as f:
        saved = ad.io.read_elem(f['obs'])
    expected = o.index[o['kept_full' if name == 'full' else 'selected_light']]
    if set(saved['barcode_original'].astype(str)) != set(expected):
        raise AssertionError('Existing audit membership differs from ' + name)
o['original_row'] = np.arange(len(o))
o['review_mixed'] = o['kept_full'] & o['marker_mixed_warning']
o['review_no_five_markers'] = o['kept_full'] & o['n_cm_markers_detected'].eq(0)
o['review_baseline_warning'] = o['review_mixed'] | o['review_no_five_markers']
eligible_identity = o['pass_target_disease'] & o['pass_cm_annotation'] & o['pass_author_subcluster']
eligible_lower = eligible_identity & o['pass_min_genes'] & o['pass_mito']
o['review_high_gene'] = eligible_lower & o['n_genes_by_counts'].gt(5000)
o['review_p1735'] = eligible_lower & o['donor_id'].eq('P1735')
o['review_cm_iii_context'] = eligible_lower & o['cell_type_original'].eq('Cardiomyocyte_III')
expected_counts = {'baseline_warnings_full': 987, 'baseline_warnings_light': 185,
    'above_5000': 5043, '5000_to_6000': 3438, 'above_6000': 1605,
    'p1735_above_5000': 327}
actual = {'baseline_warnings_full': int(o['review_baseline_warning'].sum()),
    'baseline_warnings_light': int((o['review_baseline_warning'] & o['selected_light']).sum()),
    'above_5000': int(o['review_high_gene'].sum()),
    '5000_to_6000': int((o['review_high_gene'] & o['n_genes_by_counts'].le(6000)).sum()),
    'above_6000': int((o['review_high_gene'] & o['n_genes_by_counts'].gt(6000)).sum()),
    'p1735_above_5000': int((o['review_high_gene'] & o['review_p1735']).sum())}
if actual != expected_counts:
    raise AssertionError(f'Review registry differs from known audit: {actual}')

candidate = o['review_baseline_warning'] | o['review_high_gene']
# Same-capture context retains multiple cell types; unflagged cells are not gold labels.
anchor_pool = o.loc[o['pass_target_disease'] & o['pass_baseline'] &
                    o['pass_author_subcluster'] & ~candidate]
rng = np.random.default_rng(seed)
anchor_ids = []
for _, group in anchor_pool.groupby(['sample_id', 'cell_type_original'], sort=True, observed=True):
    anchor_ids.extend(rng.choice(group.index.to_numpy(), size=min(20, len(group)), replace=False))
o['review_context_anchor'] = o.index.isin(anchor_ids)
selected = candidate | o['review_context_anchor'] | o['review_cm_iii_context'] | o['review_p1735']
selected_rows = np.flatnonzero(selected.to_numpy())
query = o.loc[selected].copy()
query['review_action'] = np.where(candidate.loc[query.index], 'pending_review', 'context_only')
o.loc[candidate].to_csv(out/'review_candidates.tsv.gz', sep='\t', index_label='barcode', compression='gzip')
query.to_csv(out/'query_metadata.tsv', sep='\t', index_label='barcode')
for column in ['donor_id', 'sample_id', 'cell_type_original']:
    o.loc[selected].groupby([column, 'disease'], observed=True).agg(
        selected=('original_row', 'size'), baseline_warning=('review_baseline_warning', 'sum'),
        high_gene=('review_high_gene', 'sum'), anchors=('review_context_anchor', 'sum')).to_csv(
            out/f'context_by_{column}.tsv', sep='\t')

print(f'Creating count-preserving review context: {len(query)} nuclei', flush=True)
with h5py.File(root/'_processing/SCP1303_20261111/author_input.h5ad') as f:
    source_obs = ad.io.read_elem(f['obs'])
    var = ad.io.read_elem(f['var'])
    if not source_obs.index.equals(o.index):
        raise AssertionError('Audit row ordering differs from source')
    matrix = ad.io.read_elem(f['X'])
    counts = matrix[selected_rows].copy()
    del matrix
    gc.collect()
    panels = {'cardiomyocyte': ['TNNT2', 'TTN', 'TNNI3', 'ACTC1', 'MYH7', 'MYBPC3',
                'MYL2', 'TNNC1', 'RYR2', 'PLN', 'ACTN2', 'MYH6'],
        'fibroblast': ['COL1A1', 'COL1A2', 'DCN', 'LUM', 'PDGFRA'],
        'endothelial': ['PECAM1', 'VWF', 'KDR', 'EMCN', 'CLDN5'],
        'immune': ['PTPRC', 'LST1', 'TYROBP', 'CD3D', 'CD79A'],
        'pericyte': ['RGS5', 'PDGFRB', 'CSPG4'],
        'smooth_muscle': ['MYH11', 'ACTA2', 'TAGLN'],
        'neuronal': ['NRXN1', 'SNAP25', 'SYT1'],
        'adipocyte': ['ADIPOQ', 'FABP4', 'PLIN1'],
        'mesothelial': ['WT1', 'MSLN', 'UPK3B']}
    genes = list(dict.fromkeys(g for panel in panels.values() for g in panel if g in var.index))
    positions = var.index.get_indexer(genes)
    values = counts[:, positions].toarray()
    raw = ad.io.read_elem(f['layers']['cellranger_raw'])
    raw_markers = raw[selected_rows][:, positions].toarray()
    del raw
    gc.collect()

gene_ids = var['gene_ids'].astype(str)
if not gene_ids.is_unique:
    raise ValueError('Source Ensembl IDs duplicated; resolve explicit mapping before export')
pd.DataFrame({'gene_id': gene_ids.to_numpy(), 'gene_symbol': var.index}).to_csv(
    out/'query_genes.tsv', sep='\t', index=False)
with gzip.open(out/'query_counts.mtx.gz', 'wb') as handle:
    mmwrite(handle, counts.T, comment='Unnormalized author CellBender counts; unchanged selected source rows')
logged = np.log1p(values * (1e4/np.maximum(query['total_counts'].to_numpy(), 1))[:, None])
evidence = query[['donor_id', 'sample_id', 'disease', 'cell_type_original', 'selected_light',
    'review_baseline_warning', 'review_high_gene', 'review_context_anchor', 'review_p1735',
    'n_genes_by_counts', 'total_counts', 'pct_counts_mt', 'doublet_score']].copy()
for panel, symbols in panels.items():
    ids = [genes.index(g) for g in symbols if g in genes]
    evidence[panel+'_markers_detected'] = (values[:, ids] > 0).sum(axis=1)
    evidence[panel+'_mean_log1p_cp10k'] = logged[:, ids].mean(axis=1)
extra = ['MYBPC3', 'RYR2', 'PLN', 'TNNC1', 'ACTN2', 'MYH6']
five = ['TNNT2', 'TTN', 'ACTC1', 'TNNI3', 'MYH7']
no_five = ~(values[:, [genes.index(g) for g in five]] > 0).any(axis=1)
no_extra = ~(values[:, [genes.index(g) for g in extra]] > 0).any(axis=1)
priority = evidence['selected_light'].to_numpy() & no_five & no_extra
if int(priority.sum()) != 26:
    raise AssertionError('Priority two-panel-negative light nuclei differ from previous review')
evidence['priority_two_panels_absent_light'] = priority
evidence.to_csv(out/'marker_program_evidence.tsv.gz', sep='\t', index_label='barcode', compression='gzip')
evidence.loc[priority].to_csv(out/'priority_26.tsv', sep='\t', index_label='barcode')
paired = pd.DataFrame({'gene': genes, 'mean_cellbender': values.mean(axis=0),
    'mean_cellranger': raw_markers.mean(axis=0),
    'fraction_detected_cellbender': (values > 0).mean(axis=0),
    'fraction_detected_cellranger': (raw_markers > 0).mean(axis=0)})
paired.to_csv(out/'paired_marker_summary.tsv', sep='\t', index=False)

group_names = query['disease'].astype(str) + ':' + query['cell_type_original'].astype(str)
means = pd.DataFrame(logged, columns=genes).groupby(group_names.to_numpy(), observed=True).mean()
means.to_csv(out/'marker_means_by_author_type.tsv', sep='\t')
fig, ax = plt.subplots(figsize=(16, max(6, .25*len(means))))
im = ax.imshow(means.to_numpy(), aspect='auto', cmap='viridis')
ax.set_xticks(np.arange(len(genes)), genes, rotation=75, ha='right', fontsize=8)
ax.set_yticks(np.arange(len(means)), means.index, fontsize=8)
ax.set_title('Review context: multi-lineage marker evidence (diagnostic only)')
fig.colorbar(im, ax=ax, label='Mean log1p(counts / total * 10000)')
fig.tight_layout(); fig.savefig(out/'marker_program_context.png', dpi=150); plt.close(fig)
part = o.loc[o['review_p1735']]
fig, axes = plt.subplots(1, 2, figsize=(11, 4))
for mask, title, color in [(part['n_genes_by_counts'].le(5000), '<=5000', '#247BA0'),
                            (part['n_genes_by_counts'].gt(5000), '>5000', '#D95D39')]:
    axes[0].scatter(part.loc[mask, 'total_counts'], part.loc[mask, 'n_genes_by_counts'],
                    s=9, alpha=.6, color=color, label=title)
    axes[1].scatter(part.loc[mask, 'n_genes_by_counts'], part.loc[mask, 'doublet_score'],
                    s=9, alpha=.6, color=color)
axes[0].set_xscale('log'); axes[0].set_xlabel('CellBender total counts'); axes[0].set_ylabel('Detected genes')
axes[1].set_xlabel('Detected genes'); axes[1].set_ylabel('Author Scrublet score')
axes[0].legend(); fig.suptitle('P1735: upper-threshold review; no new doublet call')
fig.tight_layout(); fig.savefig(out/'P1735_evidence.png', dpi=160); plt.close(fig)

summary = {'status': 'registry_and_initial_evidence_complete', 'created_utc': datetime.now(timezone.utc).isoformat(),
    'seed': seed, 'membership_checks': actual, 'context_nuclei': len(query),
    'context_cell_types': int(query['cell_type_original'].nunique()),
    'context_donors': int(query['donor_id'].nunique()), 'priority_light_nuclei': int(priority.sum()),
    'source_gene_ids': len(gene_ids), 'panels': panels,
    'decision_status': 'All candidates remain pending; no new filters applied',
    'files': []}
for name in ['query_counts.mtx.gz', 'query_genes.tsv', 'query_metadata.tsv', 'review_candidates.tsv.gz']:
    path = out/name
    with path.open('rb') as handle:
        digest = hashlib.file_digest(handle, 'sha256').hexdigest()
    summary['files'].append({'file': name, 'bytes': path.stat().st_size, 'sha256': digest})
(out/'REGISTRY_SUMMARY.json').write_text(json.dumps(summary, indent=2)+'\n')
print(json.dumps({k:v for k,v in summary.items() if k not in ['panels', 'files']}, indent=2), flush=True)
