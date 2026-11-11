"""Independently stream all official counts and verify every read index and split role."""
from datetime import datetime, timezone
from pathlib import Path
import csv
import hashlib
import json
import traceback
import anndata as ad
import h5py
import numpy as np
import pandas as pd
from scipy.sparse import csr_matrix

ROOT = Path.home() / 'xiaoe-internship-data'
SOURCE = ROOT / 'benchmarks/competitions/vcc2025'
WORK = ROOT / '_processing/all_data_cleaning_20261111/VCC2025'
OUT = ROOT / 'processed/VCC2025_qc_20261111'
FILES = {'train': ('train', 'Training'), 'validation': ('validation', 'Validation'), 'test': ('test', 'Test')}
CHECKS = []


def check(name, passed):
    CHECKS.append({'name': name, 'passed': bool(passed)})
    if not passed:
        raise AssertionError(name)


def digest(path, algorithm='sha256'):
    with path.open('rb') as handle:
        return hashlib.file_digest(handle, algorithm).hexdigest()


def read_index(name):
    return pd.read_csv(OUT / name, sep='\t', index_col=0, keep_default_na=False)


def run():
    baseline = json.loads((WORK / 'ORIGINAL_DOWNLOAD_BASELINE.json').read_text())
    export = json.loads((OUT / 'EXPORT_SUCCESS.json').read_text())
    original_records = pd.read_csv(OUT / 'ORIGINAL_FILES.tsv', sep='\t', keep_default_na=False).set_index('server_path')
    check('all_seven_original_files_accounted_for', len(baseline) == len(original_records) == 7)
    for row in baseline:
        path = Path(row['server_path'])
        actual_md5 = digest(path, 'md5')
        actual_sha256 = digest(path)
        recorded = original_records.loc[str(path)]
        check('original_bytes_' + str(path.relative_to(SOURCE)), path.stat().st_size == int(row['bytes']) == recorded['bytes'])
        check('original_MD5_' + path.name, actual_md5 == row['md5'] == recorded.MD5)
        check('original_new_SHA256_record_' + path.name, actual_sha256 == recorded.SHA256)
        if row.get('sha256'):
            check('original_historical_SHA256_' + path.name, actual_sha256 == row['sha256'])
        check('honest_SHA256_baseline_availability_' + path.name, bool(recorded.historical_SHA256_available) == bool(row.get('sha256')))
    print('All7 original files independently rehashed', flush=True)
    with (SOURCE / 'gene_names.csv').open() as handle:
        genes = pd.Index([row[0] for row in csv.reader(handle)], name='gene_symbol')
    feature = read_index('feature_index_QC.tsv.gz')
    check('fixed18080_gene_panel_unique_and_unmodified', len(genes) == 18080 and genes.is_unique and feature.index.equals(genes))
    check('native_gene_positions_and_MT_flags', np.array_equal(feature.source_column, np.arange(len(genes))) and np.array_equal(feature.mt_symbol_within_Flex_panel, genes.str.startswith('MT-')))
    mt_mask = genes.str.startswith('MT-')
    cell_tables = {}
    source_obs = {}
    target_sets = {}
    full_gene_totals = {}
    full_gene_cells = {}
    for split, (directory, suffix) in FILES.items():
        path = SOURCE / directory / ('adata_' + suffix + '.h5ad')
        audit = read_index(split + '_cell_QC.tsv.gz')
        cell_tables[split] = audit
        with h5py.File(path) as handle:
            obs = ad.io.read_elem(handle['obs'])
            var = ad.io.read_elem(handle['var'])
            source_obs[split] = obs
            x = handle['X']
            pointer = x['indptr'][:]
            n_cells, n_genes = map(int, x.attrs['shape'])
            check(split + '_real_CSR_dimensions_and_indptr', x.attrs['encoding-type'] == 'csr_matrix' and n_cells == len(obs) and n_genes == len(var) == len(genes) and len(pointer) == n_cells + 1 and pointer[0] == 0 and pointer[-1] == len(x['data']) == len(x['indices']) and (np.diff(pointer) >= 0).all())
            check(split + '_no_fabricated_raw_or_counts_layer', 'raw' not in handle and not list(handle['layers']))
            check(split + '_original_cells_and_gene_axis_unique_and_aligned', obs.index.is_unique and var.index.equals(genes) and audit.index.equals(obs.index))
            check(split + '_native_metadata_fields_exactly_preserved', all(audit[column].astype(str).tolist() == obs[column].astype(str).tolist() for column in ['target_gene', 'guide_id', 'batch']))
            check(split + '_exact_source_row_and_path_mapping', np.array_equal(audit.source_row, np.arange(n_cells)) and audit.count_source_path.eq(str(path)).all())
            if split == 'train':
                check('native_train_geneIDs_used_only_on_exact_panel', np.array_equal(feature.gene_id_from_train_native_var, var.gene_id.astype(str)) and feature.gene_id_source.eq('author_train_var_same_exact_gene_axis').all())
            total = np.zeros(n_cells, dtype=np.int64)
            detected = np.zeros(n_cells, dtype=np.int32)
            mt_counts = np.zeros(n_cells, dtype=np.int64)
            gene_sum = np.zeros(n_genes, dtype=np.int64)
            gene_n = np.zeros(n_genes, dtype=np.int64)
            stored_zeros = 0
            noncanonical_blocks = 0
            maximum_count = 0
            for start in range(0, n_cells, 1536):
                stop = min(start + 1536, n_cells)
                left, right = int(pointer[start]), int(pointer[stop])
                values = x['data'][left:right]
                columns = x['indices'][left:right]
                check(split + '_full_numeric_block_' + str(start), np.isfinite(values).all() and (values >= 0).all() and np.equal(values, np.rint(values)).all())
                check(split + '_full_column_bounds_block_' + str(start), not len(columns) or (columns.min() >= 0 and columns.max() < n_genes))
                stored_zeros += int(np.count_nonzero(values == 0))
                if len(values):
                    maximum_count = max(maximum_count, int(values.max()))
                matrix = csr_matrix((values.astype(np.int64), columns, pointer[start:stop + 1] - left), shape=(stop - start, n_genes))
                noncanonical_blocks += int(not matrix.has_canonical_format)
                matrix.sum_duplicates()
                matrix.sort_indices()
                matrix.eliminate_zeros()
                # Use direct reduceat/bincount rather than the preparation row/column sums.
                occupied = np.diff(matrix.indptr) > 0
                if occupied.any():
                    result = np.zeros(stop - start, dtype=np.int64)
                    result[occupied] = np.add.reduceat(matrix.data, matrix.indptr[:-1][occupied])
                    total[start:stop] = result
                detected[start:stop] = np.diff(matrix.indptr)
                gene_sum += np.asarray(matrix.sum(axis=0)).ravel()
                gene_n += np.bincount(matrix.indices, minlength=n_genes)
                mt_counts[start:stop] = np.asarray(matrix[:, mt_mask].sum(axis=1)).ravel()
                controls = obs.target_gene.iloc[start:stop].eq('non-targeting').to_numpy()
                for local in np.flatnonzero(controls):
                    low, high = matrix.indptr[local:local + 2]
                    fingerprint = hashlib.sha256()
                    fingerprint.update(np.asarray(matrix.indices[low:high], dtype='<i4').tobytes())
                    fingerprint.update(np.asarray(matrix.data[low:high], dtype='<i8').tobytes())
                    if fingerprint.hexdigest() != audit.control_expression_SHA256.iloc[start + local]:
                        raise AssertionError(split + '_control_value_mismatch_row_' + str(start + local))
                if start % (1536 * 32) == 0:
                    print(split, 'independent full scan', stop, '/', n_cells, flush=True)
            check(split + '_every_row_total_and_detected_count', np.array_equal(audit.total_counts, total) and np.array_equal(audit.n_genes_by_counts, detected))
            check(split + '_within_panel_MT_recomputed', np.array_equal(audit.mt_counts_within_published_panel, mt_counts) and np.allclose(audit.pct_MT_within_published_panel, np.divide(100 * mt_counts, total, out=np.zeros(n_cells), where=total > 0), rtol=0, atol=1e-10))
            check(split + '_every_gene_total_and_detection_count', np.array_equal(feature[split + '_total_counts'], gene_sum) and np.array_equal(feature[split + '_n_cells_detected'], gene_n) and np.array_equal(feature[split + '_all_zero'], gene_n == 0))
            full_gene_totals[split] = gene_sum
            full_gene_cells[split] = gene_n
            reported = next(row for row in export['split_summaries'] if row['split'] == split)
            check(split + '_matrix_summary_independently_recomputed', reported['stored_zeros'] == stored_zeros and reported['maximum_count'] == maximum_count and reported['nnz'] == int(pointer[-1]))
            check(split + '_canonical_storage_status_recomputed', (reported['noncanonical_blocks'] > 0) == (noncanonical_blocks > 0))
        control = obs.target_gene.eq('non-targeting').to_numpy()
        valid_metadata = obs[['target_gene', 'guide_id', 'batch']].notna().all(axis=1).to_numpy() & obs[['target_gene', 'guide_id', 'batch']].astype(str).ne('').all(axis=1).to_numpy()
        usable = (total > 0) & (detected > 0) & valid_metadata
        check(split + '_only_fixed_technical_QC_applied', np.array_equal(audit.qc_usable, usable) and np.array_equal(audit.is_control, control) and np.array_equal(audit.shared_control_reference, control))
        expected_reason = np.select([~valid_metadata, total == 0], ['missing_source_annotation', 'zero_information_cell'], default='retained_publisher_release')
        check(split + '_quarantine_reasons_recomputed', np.array_equal(audit.qc_reason, expected_reason))
        check(split + '_roles_preserve_official_splits', audit.official_split.eq(split).all() and np.array_equal(audit.training_eligible_under_official_split, usable & (split == 'train')) and np.array_equal(audit.held_out_perturbation_evaluation, usable & ~control & (split != 'train')))
        expected_role = np.where(control, 'shared_non_targeting_reference', 'training_perturbation' if split == 'train' else 'held_out_' + split + '_perturbation')
        check(split + '_shared_reference_never_claimed_independent_outcome', np.array_equal(audit.use_role, expected_role) and np.array_equal(audit.source_occurrences_across_official_splits, np.where(control, 3, 1)))
        check(split + '_cellUID_stable_and_no_HCM_labels', audit.cell_uid.tolist() == ('VCC2025:' + obs.index.astype(str)).tolist() and not audit.disease_label_assigned.any() and not audit.patient_donor_split_applicable.any() and audit.cell_context.eq('H1_hESC').all())
        usable_table = read_index(split + '_usable_index.tsv.gz')
        check(split + '_usable_index_is_exact_rows_and_metadata', usable_table.index.equals(obs.index[usable]) and usable_table.astype(str).equals(audit.loc[usable].astype(str)))
        excluded = read_index(split + '_quarantine_index.tsv')
        check(split + '_quarantine_rows_exact_and_not_placeholder', excluded.index.tolist() == obs.index[~usable].tolist())
        if split != 'train':
            heldout = read_index(split + '_held_out_perturbation_index.tsv.gz')
            check(split + '_heldout_index_exactly_excludes_shared_controls', heldout.index.equals(obs.index[usable & ~control]) and heldout.astype(str).equals(audit.loc[usable & ~control].astype(str)))
        official = pd.read_csv(SOURCE / directory / ('pert_counts_' + suffix + '.csv')).set_index('target_gene')
        target_sets[split] = set(official.index)
        check(split + '_official_target_set_has_no_controls', official.index.is_unique and set(official.index) == set(obs.target_gene) - {'non-targeting'})
        for target, row in official.iterrows():
            positions = obs.target_gene.eq(target).to_numpy()
            check(split + '_official_count_and_median_UMI_' + target, positions.sum() == row.n_cells and float(np.median(total[positions])) == row.median_umi_per_cell)
        summary = pd.read_csv(OUT / (split + '_perturbation_summary.tsv'), sep='\t').set_index('target_gene')
        actual_summary = audit.groupby('target_gene', observed=True).agg(n_cells=('cell_uid', 'size'), median_umi_per_cell=('total_counts', 'median'))
        check(split + '_full_perturbation_summary_recomputed', summary.index.equals(actual_summary.index) and np.array_equal(summary.n_cells, actual_summary.n_cells) and np.array_equal(summary.median_umi_per_cell, actual_summary.median_umi_per_cell))
        print('Independently verified split:', split, flush=True)
    finish(baseline, export, original_records, genes, feature, cell_tables, source_obs, target_sets)


def finish(baseline, export, original_records, genes, feature, cell_tables, source_obs, target_sets):
    check('official_target_counts_150_50_100', {key: len(value) for key, value in target_sets.items()} == {'train': 150, 'validation': 50, 'test': 100})
    check('native_gene_mapping_present_and_nonempty', feature.gene_id_from_train_native_var.ne('').all() and feature.gene_id_from_train_native_var.ne('nan').all())
    expected_rows = []
    for left, right in [('train', 'validation'), ('train', 'test'), ('validation', 'test')]:
        a, b = source_obs[left], source_obs[right]
        shared = a.index.intersection(b.index).sort_values()
        check(left + '_' + right + '_target_genes_disjoint', not (target_sets[left] & target_sets[right]))
        check(left + '_' + right + '_only_shared38176_controls', len(shared) == 38176 and a.loc[shared, 'target_gene'].eq('non-targeting').all() and b.loc[shared, 'target_gene'].eq('non-targeting').all())
        check(left + '_' + right + '_shared_source_metadata_identical', a.loc[shared].astype(str).equals(b.loc[shared].astype(str)))
        ca, cb = cell_tables[left].loc[shared], cell_tables[right].loc[shared]
        check(left + '_' + right + '_all_shared_count_rows_fingerprinted_identically', ca.control_expression_SHA256.str.len().eq(64).all() and ca.control_expression_SHA256.equals(cb.control_expression_SHA256))
    combined = pd.concat(cell_tables.values())
    noncontrol = combined.loc[~combined.is_control]
    check('noncontrol_cell_identifiers_unique_across_splits', noncontrol.index.is_unique)
    check('each_guide_maps_to_exactly_one_target', (combined.groupby('guide_id', observed=True).target_gene.nunique() == 1).all())
    check('noncontrol_hash_field_is_empty', noncontrol.control_expression_SHA256.eq('').all())
    check('cellUID_and_sourceID_identify_same_unique_cells', combined.cell_uid.nunique() == combined.index.nunique())
    unique_controls = read_index('shared_controls_unique_index.tsv.gz')
    expected_controls = cell_tables['train'].loc[cell_tables['train'].is_control]
    check('unique_control_index_exact_train_reference', len(unique_controls) == 38176 and unique_controls[expected_controls.columns].astype(str).equals(expected_controls.astype(str)))
    check('unique_control_copy_choice_explicit', unique_controls.reference_count_source_choice.eq('train_copy_same_expression_verified_in_all_three_releases').all())
    check('total_export_counts_recomputed', export['file_observations'] == len(combined) and export['unique_cells'] == combined.index.nunique() and export['unique_perturbed_cells'] == len(noncontrol) and export['shared_controls'] == len(unique_controls) and export['genes'] == len(genes))
    dataset_index = pd.read_csv(OUT / 'DATASET_INDEX.tsv', sep='\t').set_index('official_split')
    summary = pd.read_csv(OUT / 'split_QC_summary.tsv', sep='\t').set_index('split')
    perturbations = pd.read_csv(OUT / 'PERTURBATION_SPLITS.tsv', sep='\t', keep_default_na=False)
    check('exact_three_splits_in_dataset_and_summary', dataset_index.index.tolist() == summary.index.tolist() == list(FILES))
    for split, audit in cell_tables.items():
        entry = dataset_index.loc[split]
        directory, suffix = FILES[split]
        expected_source = str(SOURCE / directory / ('adata_' + suffix + '.h5ad'))
        check(split + '_dataset_index_paths_dimensions_and_no_redundant_matrix', entry.dataset == 'VCC2025' and entry.source_count_file == expected_source and entry.usable_index == str(OUT / (split + '_usable_index.tsv.gz')) and entry.source_cells == len(audit) and entry.usable_cells == int(audit.qc_usable.sum()) and entry.genes == len(genes) and entry.matrix_handling == 'read_validated_original_counts_no_rewrite')
        actual = {
            'source_cells': len(audit), 'qc_usable': int(audit.qc_usable.sum()),
            'genes_preserved': len(genes), 'target_perturbations': len(target_sets[split]),
            'control_cells': int(audit.is_control.sum()),
            'unique_noncontrol_cells': int((~audit.is_control).sum()),
            'zero_information_cells': int(audit.total_counts.eq(0).sum()),
            'zero_features_preserved': int(feature[split + '_all_zero'].sum()),
            'median_total_counts': float(audit.total_counts.median()),
            'median_detected_genes': float(audit.n_genes_by_counts.median()),
        }
        exported = next(row for row in export['split_summaries'] if row['split'] == split)
        check(split + '_summary_all_cell_and_feature_counts_recomputed', all(summary.loc[split, key] == value == exported[key] for key, value in actual.items()))
        target_summary = pd.read_csv(OUT / (split + '_perturbation_summary.tsv'), sep='\t', keep_default_na=False)
        expected_rows.append(target_summary)
    check('combined_perturbation_split_table_exact', perturbations.astype(str).equals(pd.concat(expected_rows, ignore_index=True).astype(str)))
    policy = json.loads((OUT / 'qc_policy.json').read_text())
    check('policy_preserves_modality_and_count_semantics', policy['dataset'] == 'VCC2025' and 'H1 human embryonic' in policy['context'] and 'not HCM' in policy['modality_boundary'] and 'read-only' in policy['count_semantics'])
    check('HCM_filters_and_training_explicitly_not_applied', {'HCM500to5000_gene_filter', 'MT_lt5pct', 'CM_identity_filter', 'training', 'tokenization'}.issubset(policy['not_applied']))
    inspection = json.loads((OUT / 'INPUT_INSPECTION.json').read_text())
    check('metadata_inspection_unique_cells_matches_full_audit', inspection['unique_cell_identifiers'] == combined.index.nunique())
    check('baseline_copy_exact', json.loads((OUT / 'ORIGINAL_DOWNLOAD_BASELINE.json').read_text()) == baseline)
    for name, expected in export['hashes'].items():
        path = OUT / name
        check('exported_payload_bytes_' + name, path.stat().st_size == expected['bytes'])
        check('exported_payload_SHA256_' + name, digest(path) == expected['SHA256'])
        check('exported_payload_MD5_' + name, digest(path, 'md5') == expected['MD5'])
    prepare_exit_code = (WORK / 'prepare_exit_code.txt').read_text().strip()
    check('prepare_job_exited_successfully', prepare_exit_code == '0')
    completed = datetime.now(timezone.utc).isoformat()
    verification = {
        'dataset': 'VCC2025', 'status': 'complete_verified', 'verified_utc': completed,
        'checks_passed': len(CHECKS), 'checks_failed': 0, 'checks': CHECKS,
        'scope': 'Independent full scan of every stored count value, source hashes, every exported cell/feature metric, official split and source metadata alignment; no biological efficacy or singlet certification.',
        'unique_cells': int(combined.index.nunique()), 'unique_perturbed_cells': len(noncontrol),
        'shared_controls': len(unique_controls), 'published_file_observations': len(combined),
        'genes_preserved': len(genes), 'retained_observations': int(combined.qc_usable.sum()),
        'quarantined_observations': int((~combined.qc_usable).sum()),
    }
    (OUT / 'INDEPENDENT_VERIFICATION.json').write_text(json.dumps(verification, indent=2) + '\n')
    artifacts = []
    for original in original_records.reset_index().to_dict('records'):
        artifacts.append({'dataset': 'VCC2025', 'artifact_kind': 'unchanged_original',
            'file': original['file'], 'server_path': original['server_path'], 'bytes': original['bytes'],
            'MD5': original['MD5'], 'SHA256': original['SHA256']})
    for path in sorted(OUT.iterdir()):
        if path.is_file() and path.name != 'OUTPUT_MANIFEST.tsv':
            artifacts.append({'dataset': 'VCC2025', 'artifact_kind': 'validated_index_or_audit',
                'file': path.name, 'server_path': str(path), 'bytes': path.stat().st_size,
                'MD5': digest(path, 'md5'), 'SHA256': digest(path)})
    pd.DataFrame(artifacts).to_csv(OUT / 'OUTPUT_MANIFEST.tsv', sep='\t', index=False)
    gate = {'dataset': 'VCC2025', 'status': 'complete_verified', 'verified_utc': completed,
        'independent_checks_passed': len(CHECKS), 'independent_checks_failed': 0,
        'all7_original_files_rehashed_against_download_baseline': True,
        'all_stored_values_independently_scanned': True, 'official_splits_preserved': True,
        'shared_control_expression_and_metadata_identical': True,
        'original_source_files_used_read_only_via_explicit_indices': True,
        'disease_label_assigned': False, 'biological_knockdown_or_singlet_certified': False,
        'output_directory': str(OUT), 'verification_file_SHA256': digest(OUT / 'INDEPENDENT_VERIFICATION.json'),
        'output_manifest_SHA256': digest(OUT / 'OUTPUT_MANIFEST.tsv'),
        'source_prepare_exit_code': prepare_exit_code}
    (WORK / 'FINAL_COMPLETION_GATE.json').write_text(json.dumps(gate, indent=2) + '\n')
    print(json.dumps({key: value for key, value in verification.items() if key != 'checks'}, indent=2), flush=True)


if __name__ == '__main__':
    try:
        run()
    except BaseException as error:
        (WORK / 'VERIFICATION_FAILURE.json').write_text(json.dumps({'status': 'failed', 'error': str(error),
            'traceback': traceback.format_exc(), 'checks': CHECKS}, indent=2) + '\n')
        raise
