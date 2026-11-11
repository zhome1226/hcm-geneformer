"""Calibrate reference-tail flags and join independent, marker and graph evidence."""
from datetime import datetime, timezone
from pathlib import Path
import json

import anndata as ad
import h5py
import numpy as np
import pandas as pd

root = Path.home() / 'xiaoe-internship-data'
base = root / 'processed/SCP1303_review_20261111'
out = base / 'identity'
cal = pd.read_csv(out / 'calibration_predictions.tsv', sep='\t')
query = pd.read_csv(out / 'query_predictions.tsv', sep='\t', index_col='barcode')
meta = pd.read_csv(base / 'query_metadata.tsv', sep='\t', index_col='barcode', low_memory=False)
panels = pd.read_csv(base / 'marker_program_evidence.tsv.gz', sep='\t', index_col='barcode')
graph = pd.read_csv(out / 'query_graph_evidence.tsv', sep='\t', index_col='barcode')
background = pd.read_csv(out / 'paired_background_evidence.tsv.gz', sep='\t', index_col='barcode')
assert query.index.equals(meta.index) and meta.index.equals(panels.index) and meta.index.equals(graph.index)
thresholds = {}
metrics = []
for label in sorted(cal['truth'].unique()):
    correct = cal.loc[cal['truth'].eq(label) & cal['labels'].eq(label)]
    true = cal['truth'].eq(label)
    assigned = cal['labels'].eq(label)
    metrics.append({'label': label, 'truth_nuclei': int(true.sum()),
        'predicted_nuclei': int(assigned.sum()),
        'recall': float((true & assigned).sum() / max(1, true.sum())),
        'precision': float((true & assigned).sum() / max(1, assigned.sum()))})
    if len(correct) >= 30:
        thresholds[label] = {'correct_calibration_nuclei': len(correct),
            'delta_q05': float(correct['delta'].quantile(.05)),
            'delta_next_q05': float(correct['delta_next'].quantile(.05))}
calibration = {'completed_utc': datetime.now(timezone.utc).isoformat(),
    'accuracy': float(cal['truth'].eq(cal['labels']).mean()), 'per_label': metrics,
    'empirical_tail_thresholds': thresholds,
    'interpretation': 'Thresholds flag lower tails of correctly annotated reference-donor holdout nuclei. They do not estimate posterior confidence or validate disease classification. No disease labels used.'}
cal_pass = pd.Series(False, index=cal.index)
for label, values in thresholds.items():
    cal_pass |= (cal['labels'].eq(label) & cal['pruned_labels'].fillna('').eq(label)
        & cal['delta'].ge(values['delta_q05']) & cal['delta_next'].ge(values['delta_next_q05']))
calibration['calibration_set_tail_flag_metrics'] = []
for label in sorted(thresholds):
    accepted = cal['labels'].eq(label) & cal_pass
    true = cal['truth'].eq(label)
    calibration['calibration_set_tail_flag_metrics'].append({'label': label,
        'accepted_nuclei': int(accepted.sum()),
        'accepted_precision_on_calibration_set': float((accepted & true).sum()/max(1, accepted.sum())),
        'accepted_recall_on_calibration_set': float((accepted & true).sum()/max(1, true.sum()))})
calibration['tail_metrics_limitation'] = 'Thresholds and these apparent metrics use the same reference calibration set; not an independent validation of the thresholds.'
(out / 'REFERENCE_CALIBRATION.json').write_text(json.dumps(calibration, indent=2) + '\n')
pd.crosstab(cal['truth'], cal['labels']).to_csv(out / 'reference_holdout_confusion.tsv', sep='\t')
query['empirical_reference_tail_pass'] = False
for label, values in thresholds.items():
    mask = query['labels'].eq(label) & query['pruned_labels'].fillna('').eq(label)
    mask &= query['delta'].ge(values['delta_q05']) & query['delta_next'].ge(values['delta_next_q05'])
    query.loc[mask, 'empirical_reference_tail_pass'] = True
joined = meta.join(query.add_prefix('reference_')).join(panels.filter(regex='markers_detected|mean_log1p_cp10k|priority_'))
joined = joined.join(graph).join(background)
is_cm = joined['cell_type_original'].astype(str).str.startswith('Cardiomyocyte')
joined['independent_cm_support'] = joined['reference_labels'].eq('Cardiomyocyte') & joined['reference_empirical_reference_tail_pass']
joined['independent_noncm_support'] = joined['reference_labels'].ne('Cardiomyocyte') & joined['reference_empirical_reference_tail_pass']
noncm = ['fibroblast', 'endothelial', 'immune', 'pericyte', 'smooth_muscle', 'neuronal', 'adipocyte', 'mesothelial']
detected = joined[[name + '_markers_detected' for name in noncm]]
programs = joined[[name + '_mean_log1p_cp10k' for name in noncm]]
joined['max_noncm_markers_detected'] = detected.max(axis=1)
joined['max_noncm_program'] = programs.max(axis=1)
joined['dominant_noncm_program'] = programs.idxmax(axis=1).str.replace('_mean_log1p_cp10k', '', regex=False)
joined['coherent_mixed_panel_flag'] = joined['cardiomyocyte_markers_detected'].ge(2) & joined['max_noncm_markers_detected'].ge(3)
audit = pd.read_csv(root / 'processed/SCP1303_qc_20261111/reports/cell_qc_audit.tsv.gz',
    sep='\t', usecols=['barcode', 'capture_id', 'doublet_score'], index_col='barcode')
rank = audit.groupby('capture_id')['doublet_score'].rank(pct=True)
joined['author_doublet_capture_percentile'] = rank.reindex(joined.index)
joined['secondary_doublet_review_trigger'] = is_cm & (
    (joined['coherent_mixed_panel_flag'] & joined['author_doublet_capture_percentile'].ge(.95)) |
    (joined['independent_noncm_support'] & joined['cardiomyocyte_markers_detected'].ge(2) & joined['max_noncm_markers_detected'].ge(3)))
joined.to_csv(out / 'joined_identity_evidence.tsv.gz', sep='\t', index_label='barcode', compression='gzip')
groups = {
    'baseline_warnings_full': joined['review_baseline_warning'],
    'baseline_warnings_light': joined['review_baseline_warning'] & joined['selected_light'],
    'priority_26': joined['priority_two_panels_absent_light'],
    'HCM_CM_III': joined['disease'].eq('HCM') & joined['cell_type_original'].eq('Cardiomyocyte_III'),
    'high_5000_to_6000': joined['review_high_gene'] & joined['n_genes_by_counts'].le(6000),
    'high_above_6000': joined['review_high_gene'] & joined['n_genes_by_counts'].gt(6000),
    'P1735_at_most_5000': joined['review_p1735'] & joined['n_genes_by_counts'].le(5000),
    'P1735_above_5000': joined['review_p1735'] & joined['n_genes_by_counts'].gt(5000),
    'author_cm_context_controls': joined['review_context_anchor'] & is_cm,
}
summary = {}
for name, mask in groups.items():
    part = joined.loc[mask]
    summary[name] = {'nuclei': len(part),
        'reference_labels': part['reference_labels'].value_counts().to_dict(),
        'empirical_cm_support': int(part['independent_cm_support'].sum()),
        'empirical_noncm_support': int(part['independent_noncm_support'].sum()),
        'reference_tail_flagged_or_uncalibrated': int((~part['reference_empirical_reference_tail_pass']).sum()),
        'secondary_doublet_review_trigger': int(part['secondary_doublet_review_trigger'].sum()),
        'author_cm_neighbors_median': float(part['author_cm_neighbor_fraction'].median())}
joined.loc[joined['priority_two_panels_absent_light']].to_csv(out / 'priority_26_identity.tsv', sep='\t', index_label='barcode')
joined.loc[joined['secondary_doublet_review_trigger']].groupby('capture_id').size().to_csv(out / 'doublet_trigger_captures.tsv', sep='\t', header=['trigger_nuclei'])
report = {'reference_calibration': calibration, 'groups': summary,
    'doublet_trigger_nuclei': int(joined['secondary_doublet_review_trigger'].sum()),
    'doublet_trigger_captures': int(joined.loc[joined['secondary_doublet_review_trigger'], 'capture_id'].nunique()),
    'new_filters_applied': False}
(out / 'IDENTITY_SUMMARY.json').write_text(json.dumps(report, indent=2) + '\n')
print(json.dumps({'calibration_accuracy': calibration['accuracy'], 'groups': summary,
    'doublet_trigger_nuclei': report['doublet_trigger_nuclei'],
    'doublet_trigger_captures': report['doublet_trigger_captures']}, indent=2), flush=True)
