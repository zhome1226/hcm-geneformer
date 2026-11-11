"""Independently re-read all assays and prove the published indices and boundaries."""
from datetime import datetime, timezone
from pathlib import Path
import hashlib
import json
import traceback
import anndata as ad
import h5py
import numpy as np
import pandas as pd
from scipy import sparse

ROOT = Path.home() / 'xiaoe-internship-data'
SOURCE = ROOT / 'benchmarks/competitions/neurips2023'
WORK = ROOT / '_processing/all_data_cleaning_20261111/NeurIPS2023'
OUT = ROOT / 'processed/NeurIPS2023_qc_20261111'
EXPECTED = {
    'sc_original': ('sc_counts.h5ad', '', 'layers/counts'),
    'sc_reannotated': ('sc_counts_reannotated_with_counts.h5ad', '', 'raw/X'),
    'baseline_rna': ('multiome_counts.h5mu', 'mod/rna', 'raw/X'),
    'baseline_atac': ('multiome_counts.h5mu', 'mod/atac', 'layers/binarized'),
    'de_train': ('2023-09-14_kaggle_upload/2023-09-12_de_by_cell_type_train.h5ad', '', 'X'),
    'de_test': ('2023-09-14_kaggle_upload/2023-09-12_de_by_cell_type_test.h5ad', '', 'X'),
}
CHECKS = []


def check(name, result):
    CHECKS.append({'name': name, 'passed': bool(result)})
    if not result:
        raise AssertionError(name)


def digest(path, algorithm='sha256'):
    with path.open('rb') as handle:
        return hashlib.file_digest(handle, algorithm).hexdigest()


def read_table(name, indexed=False):
    result = pd.read_csv(OUT / name, sep='\t', index_col=0 if indexed else None,
        keep_default_na=False, float_precision='round_trip')
    if indexed:
        result.index = result.index.astype(str)
    return result


def same_table(a, b):
    # Compare by column to avoid allocating a huge fixed-width Unicode matrix.
    return a.columns.tolist() == b.columns.tolist() and a.index.astype(str).tolist() == b.index.astype(str).tolist() and all(np.array_equal(a[name].astype(str).to_numpy(), b[name].astype(str).to_numpy()) for name in a)


def native_columns_preserved(native, saved):
    if native.index.astype(str).tolist() != saved.index.astype(str).tolist():
        return False
    for name in native:
        if name not in saved:
            return False
        column = native[name]
        if pd.api.types.is_numeric_dtype(column) and not pd.api.types.is_bool_dtype(column):
            # CSV float32 tokens round-trip to float32, not to the expanded
            # float64 representation of the original value. Require exactness
            # in the publisher's actual dtype rather than relaxing tolerance.
            original = column.to_numpy()
            restored = pd.to_numeric(saved[name], errors='coerce').to_numpy().astype(original.dtype)
            if not np.array_equal(restored, original, equal_nan=True):
                return False
        else:
            expected = column.astype(object).where(column.notna(), '').astype(str).to_numpy()
            if not np.array_equal(saved[name].astype(str), expected):
                return False
    return True


def numeric_audit(node, semantics, mt_mask=None, fingerprint=False):
    if isinstance(node, h5py.Group) and node.attrs.get('encoding-type') == 'csr_matrix':
        pointer = node['indptr'][:]
        shape = tuple(node.attrs['shape'])
        assert len(pointer) == shape[0] + 1 and pointer[0] == 0 and pointer[-1] == len(node['data']) == len(node['indices']) and (np.diff(pointer) >= 0).all()
        reader = ad.io.sparse_dataset(node)
        is_sparse_storage = True
    elif isinstance(node, h5py.Group):
        reader = ad.io.read_elem(node)
        assert isinstance(reader, pd.DataFrame)
        shape = reader.shape
        is_sparse_storage = False
    else:
        reader, shape = node, node.shape
        assert len(shape) == 2
        is_sparse_storage = False
    n, p = shape
    row_sum, row_n, row_mt = np.zeros(n), np.zeros(n, dtype=np.int32), np.zeros(n)
    gene_sum, gene_n = np.zeros(p), np.zeros(p, dtype=np.int64)
    hashes = np.full(n, '', dtype=object)
    summary = {'nonzero_values': 0, 'negative_values': 0, 'fractional_values': 0, 'stored_sparse_zeros': 0}
    minimum, maximum, noncanonical = np.inf, -np.inf, False
    for start in range(0, n, 1024):
        stop = min(start + 1024, n)
        if isinstance(reader, pd.DataFrame):
            matrix = sparse.csr_matrix(reader.iloc[start:stop].to_numpy(dtype=np.float64))
        elif is_sparse_storage:
            # SparseDataset/astype may sort indices on read. Audit the native
            # on-disk order separately, excluding transitions between rows.
            left, right = int(pointer[start]), int(pointer[stop])
            native_columns = node['indices'][left:right]
            if len(native_columns) > 1:
                not_increasing = native_columns[1:] <= native_columns[:-1]
                boundaries = pointer[start + 1:stop] - left
                boundaries = boundaries[(boundaries > 0) & (boundaries < len(native_columns))]
                not_increasing[boundaries - 1] = False
                noncanonical |= bool(not_increasing.any())
            matrix = reader[start:stop, :].astype(np.float64)
        else:
            matrix = sparse.csr_matrix(reader[start:stop, :].astype(np.float64))
        values = matrix.data
        assert np.isfinite(values).all() and (not len(matrix.indices) or (matrix.indices.min() >= 0 and matrix.indices.max() < p))
        if semantics in {'counts', 'binary', 'processed_nonnegative', 'annotation_score', 'pvalue'}:
            assert (values >= 0).all()
        if semantics in {'counts', 'binary'}:
            assert (values == np.floor(values)).all() and (values < np.iinfo(np.int64).max).all()
        if semantics in {'binary', 'pvalue'}:
            assert (values <= 1).all()
        summary['negative_values'] += int((values < 0).sum())
        summary['fractional_values'] += int((values != np.floor(values)).sum())
        if is_sparse_storage:
            summary['stored_sparse_zeros'] += int((values == 0).sum())
        if len(values):
            minimum, maximum = min(minimum, float(values.min())), max(maximum, float(values.max()))
        noncanonical |= not matrix.has_canonical_format
        matrix.sum_duplicates()
        matrix.eliminate_zeros()
        matrix.sort_indices()
        # A different block size, source reader, and direct row reductions are used.
        lengths = np.diff(matrix.indptr)
        nonempty = lengths > 0
        result = np.zeros(stop - start)
        if nonempty.any():
            result[nonempty] = np.add.reduceat(matrix.data, matrix.indptr[:-1][nonempty])
        row_sum[start:stop], row_n[start:stop] = result, lengths
        gene_sum += np.asarray(matrix.sum(axis=0)).ravel()
        gene_n += np.bincount(matrix.indices, minlength=p)
        if mt_mask is not None:
            row_mt[start:stop] = np.asarray(matrix[:, mt_mask].sum(axis=1)).ravel()
        if fingerprint:
            for local in range(stop - start):
                low, high = matrix.indptr[local:local + 2]
                hasher = hashlib.sha256()
                hasher.update(matrix.indices[low:high].astype('<i4').tobytes())
                hasher.update(matrix.data[low:high].astype('<i8').tobytes())
                hashes[start + local] = hasher.hexdigest()
        summary['nonzero_values'] += matrix.nnz
    summary.update(rows=n, features=p, zero_entries=int(n * p - summary['nonzero_values']),
        minimum_stored_nonzero_or_sparse_value=float(minimum) if np.isfinite(minimum) else 0.0,
        maximum_stored_nonzero_or_sparse_value=float(maximum) if np.isfinite(maximum) else 0.0,
        zero_information_rows=int((row_n == 0).sum()), all_zero_features=int((gene_n == 0).sum()),
        row_sum_total=float(row_sum.sum()), has_noncanonical_storage=noncanonical)
    return summary, (row_sum, row_n, row_mt, hashes, gene_sum, gene_n)


def expected_semantics(unit, slot):
    if slot == 'X':
        return 'signed_DE_score' if unit.startswith('de_') else 'processed_nonnegative'
    if slot.startswith('raw/') or slot in {'layers/counts', 'obsm/protein_counts'}:
        return 'counts'
    if slot in {'layers/P.Value', 'layers/adj.P.Val', 'obsm/ora_pvals'}:
        return 'pvalue'
    if slot in {'layers/is_de', 'layers/is_de_adj', 'layers/binarized'}:
        return 'binary'
    return {'layers/logFC': 'signed_logFC', 'obsm/HTO_clr': 'signed_hashing_CLR', 'obsm/ora_estimate': 'annotation_score'}.get(slot, 'author_embedding')


def main():
    baseline = json.loads((WORK / 'ORIGINAL_DOWNLOAD_BASELINE.json').read_text())
    originals = read_table('ORIGINAL_FILES.tsv').set_index('server_path')
    export = json.loads((OUT / 'EXPORT_SUCCESS.json').read_text())
    check('all_eight_originals_accounted_for', len(baseline) == len(originals) == 8)
    for item in baseline:
        path = Path(item['server_path'])
        check('source_bytes_' + path.name, path.stat().st_size == int(item['bytes']) == originals.loc[str(path), 'bytes'])
        check('source_MD5_' + path.name, digest(path, 'md5') == item['md5'] == originals.loc[str(path), 'MD5'])
        check('source_SHA256_' + path.name, digest(path) == item['sha256'] == originals.loc[str(path), 'SHA256'])
    print('Eight originals independently rehashed', flush=True)
    dataset = read_table('DATASET_INDEX.tsv').set_index('unit')
    units = read_table('UNIT_QC_SUMMARY.tsv').set_index('unit')
    matrices = read_table('MATRIX_QC_SUMMARY.tsv').set_index(['unit', 'relative_slot'])
    check('exact_six_assay_units_and_unique_matrix_slots', dataset.index.tolist() == units.index.tolist() == list(EXPECTED) and matrices.index.is_unique)
    audits, native_obs, native_var, cell_hashes = {}, {}, {}, {}
    recomputed_matrix_slots = []
    meta_differences = {}
    for unit, (relative, group_name, primary) in EXPECTED.items():
        source_file = SOURCE / relative
        with h5py.File(source_file, 'r') as handle:
            group = handle[group_name] if group_name else handle
            obs, var = ad.io.read_elem(group['obs']), ad.io.read_elem(group['var'])
            native_obs[unit], native_var[unit] = obs, var
            audit = read_table(unit + '_observation_QC.tsv.gz', indexed=True)
            feature = read_table(unit + '_feature_QC.tsv.gz', indexed=True)
            audits[unit] = audit
            check(unit + '_unique_unmodified_source_axes', obs.index.is_unique and var.index.is_unique and audit.index.tolist() == obs.index.astype(str).tolist() and feature.index.tolist() == var.index.astype(str).tolist())
            check(unit + '_every_author_observation_field_preserved', native_columns_preserved(obs, audit))
            check(unit + '_every_author_feature_field_preserved', native_columns_preserved(var, feature))
            check(unit + '_exact_read_only_count_or_score_locations', audit.source_file.eq(str(source_file)).all() and audit.source_group.eq(group_name or '/').all() and audit.primary_matrix_slot.eq(primary).all() and np.array_equal(audit.source_row, np.arange(len(obs))) and np.array_equal(feature.source_column, np.arange(len(var))))
            if primary.startswith('raw/'):
                check(unit + '_raw_count_gene_axis_matches_native_var', ad.io.read_elem(group['raw']['var']).index.equals(var.index))
            slots = ['X'] + ['layers/' + name for name in group.get('layers', {})]
            if 'raw' in group:
                slots.append('raw/X')
            slots += ['obsm/' + name for name in group.get('obsm', {})]
            primary_metrics = None
            for slot in slots:
                semantics = expected_semantics(unit, slot)
                row = matrices.loc[(unit, slot)]
                node = group[slot]
                if isinstance(node, h5py.Group) and node.attrs.get('encoding-type') == 'dataframe':
                    table = ad.io.read_elem(node)
                    check(unit + '_' + slot + '_native_table_alignment', table.index.equals(obs.index) and table.columns.is_unique)
                    if slot == 'obsm/protein_counts':
                        check(unit + '_CMO_tags_are_hashing_not_protein_ADT', len(table.columns) == 12 and table.columns.str.startswith('CMO_').all() and row.assay_interpretation == 'cell_hashing_CMO_counts')
                        cell_hashes[unit] = table
                mt = var.index.astype(str).str.startswith('MT-') if slot == primary and unit not in {'baseline_atac', 'de_train', 'de_test'} else None
                actual, metrics = numeric_audit(node, semantics, mt_mask=mt, fingerprint=mt is not None)
                check(unit + '_' + slot + '_all_values_and_sparse_structure_checked', actual['rows'] == len(obs) and row.semantics == semantics and row.source_file == str(source_file))
                exact = [name for name in actual if name not in {'has_noncanonical_storage', 'row_sum_total'}]
                check(unit + '_' + slot + '_full_numeric_summary_recomputed', all(actual[name] == row[name] for name in exact) and np.isclose(actual['row_sum_total'], row.row_sum_total, rtol=1e-11, atol=1e-6) and actual['has_noncanonical_storage'] == (row.noncanonical_blocks > 0))
                if slot in {'X', primary} or slot.startswith('layers/'):
                    check(unit + '_' + slot + '_full_gene_axis_shape', actual['features'] == len(var))
                if slot == primary:
                    primary_metrics = metrics
                recomputed_matrix_slots.append((unit, slot))
                print('Independent matrix verified:', unit, slot, flush=True)
            total, detected, mt, fingerprints, gene_sum, gene_n = primary_metrics
            is_de = unit.startswith('de_')
            compare = (lambda a, b: np.allclose(a, b, rtol=1e-11, atol=1e-8)) if is_de else np.array_equal
            check(unit + '_every_primary_row_sum_and_detection', compare(audit.primary_row_sum, total) and np.array_equal(audit.primary_nonzero_features, detected) and np.array_equal(audit.primary_zero_information, detected == 0))
            check(unit + '_every_primary_feature_sum_and_detection', compare(feature.primary_feature_sum, gene_sum) and np.array_equal(feature.primary_nonzero_observations, gene_n) and np.array_equal(feature.all_zero_in_primary_matrix, gene_n == 0))
            if unit.startswith('sc_'):
                columns = ['donor_id', 'cell_type', 'split', 'sm_name', 'dose_uM', 'timepoint_hr', 'plate_name', 'well']
                valid = obs[columns].notna().all(axis=1).to_numpy(copy=True)
                valid &= obs[columns].astype(str).apply(lambda column: column.str.strip().ne('')).all(axis=1).to_numpy()
                valid &= np.asarray(np.isfinite(obs.dose_uM) & (obs.dose_uM >= 0) & np.isfinite(obs.timepoint_hr))
                check(unit + '_original_splits_preserved', audit.official_split.tolist() == obs['split'].astype(str).tolist() and set(obs['split']) == {'train', 'control', 'public_test', 'private_test'})
                check(unit + '_stable_shared_release_cellUID', audit.cell_uid.tolist() == ('NeurIPS2023:perturbation:' + obs.index.astype(str)).tolist() and audit.canonical_unique_cell_source.eq(unit == 'sc_original').all())
                tags = cell_hashes[unit]
                check(unit + '_every_CMO_total_and_nonzero_count', np.array_equal(audit.cell_hashing_CMO_total, tags.to_numpy().sum(axis=1)) and np.array_equal(audit.cell_hashing_CMO_nonzero_features, np.count_nonzero(tags.to_numpy(), axis=1)))
            elif unit.startswith('baseline_'):
                root = ad.io.read_elem(handle['obs'])
                valid = obs.cell_type.notna().to_numpy() & obs.cell_type.astype(str).str.strip().ne('').to_numpy() & root.donor_id.notna().to_numpy() & root.donor_id.astype(str).str.strip().ne('').to_numpy()
                check(unit + '_paired_root_donors_exact', root.index.equals(obs.index) and audit.source_donor_id.tolist() == root.donor_id.astype(str).tolist() and set(root.donor_id) == {'Donor 1', 'Donor 2', 'Donor 3'})
                check(unit + '_baseline_role_not_invented_holdout_or_patient', audit.official_split.eq('baseline_reference_no_cell_level_challenge_split').all() and audit.cell_uid.tolist() == ('NeurIPS2023:baseline_multiome:' + obs.index.astype(str)).tolist() and audit.canonical_unique_cell_source.eq(unit == 'baseline_rna').all())
                if unit == 'baseline_atac':
                    locations = var.index.to_series().str.extract(r'^([^:]+):(\d+)-(\d+)$')
                    check('every_ATAC_peak_coordinate_preserved_and_valid', locations.notna().all().all() and (locations[2].astype(int) > locations[1].astype(int)).all() and feature.peak_coordinate_valid.all())
                    meta_differences['ATAC_author_feature_count_vs_published_binary_differences'] = int(np.count_nonzero(np.asarray(obs.n_features_per_cell) != detected))
                else:
                    meta_differences['RNA_author_total_counts_vs_published_gene_total_differences'] = int(np.count_nonzero(np.asarray(obs.total_counts) != total))
                    meta_differences['RNA_author_n_counts_vs_published_gene_total_differences'] = int(np.count_nonzero(np.asarray(obs.n_counts) != total))
                    meta_differences['RNA_author_n_genes_vs_published_detection_differences'] = int(np.count_nonzero(np.asarray(obs.n_genes) != detected))
            else:
                columns = ['cell_type', 'sm_lincs_id', 'sm_name', 'SMILES', 'dose_uM', 'split']
                valid = obs[columns].notna().all(axis=1).to_numpy() & obs[columns].astype(str).apply(lambda column: column.str.strip().ne('')).all(axis=1).to_numpy()
                expected_uid = 'NeurIPS2023:condition:' + obs[['cell_type', 'sm_lincs_id', 'dose_uM']].astype(str).agg('||'.join, axis=1)
                check(unit + '_condition_units_and_official_splits_exact', audit.condition_uid.tolist() == expected_uid.tolist() and audit.condition_uid.is_unique and audit.official_split.tolist() == obs['split'].astype(str).tolist() and set(obs['split']) == ({'train', 'control'} if unit == 'de_train' else {'private_test', 'public_test'}))
                score, p, fc = group['X'][:], group['layers/P.Value'][:], group['layers/logFC'][:]
                mask = p > 0
                expected_score = np.sign(fc[mask]) * np.negative(np.log10(p[mask]))
                check(unit + '_published_signed_P_statistic_on_every_positive_P', np.allclose(score[mask], expected_score, rtol=1e-7, atol=1e-7))
                recorded = read_table('published_DE_consistency.tsv').set_index('unit').loc[unit]
                check(unit + '_reported_DE_score_formula_error_exact', np.isclose(recorded.maximum_score_formula_error_positive_P, np.max(np.abs(score[mask] - expected_score)), rtol=0, atol=1e-12))
                check(unit + '_no_guessed_P_zero_floor_or_DE_flag_rule', recorded.positive_P_entries == int(mask.sum()) and recorded.zero_P_entries == int((~mask).sum()) and recorded.is_de_vs_P_lt_0_05_mismatches == int(np.count_nonzero(group['layers/is_de'][:] != (p < .05))) and recorded.is_de_adj_vs_adj_P_lt_0_05_mismatches == int(np.count_nonzero(group['layers/is_de_adj'][:] != (group['layers/adj.P.Val'][:] < .05))))
            usable = valid & ((detected > 0) if not is_de else True)
            reason = np.select([~valid, (detected == 0) & (not is_de)], ['incomplete_required_source_metadata', 'zero_information_cell'], default='retained_publisher_curated_release')
            check(unit + '_fixed_metadata_and_information_QC_only', np.array_equal(audit.technical_qc_usable, usable) and np.array_equal(audit.qc_reason, reason) and not audit.disease_label_assigned.any())
            expected_role = np.select([audit.official_split.eq('train'), audit.official_split.eq('control'), audit.official_split.eq('public_test'), audit.official_split.eq('private_test')], ['official_training', 'publisher_control_reference', 'held_out_public_test', 'held_out_private_test'], default='paired_baseline_reference')
            check(unit + '_training_test_reference_roles_exact', np.array_equal(audit.training_role, expected_role))
            check(unit + '_usable_index_every_row_and_metadata_exact', same_table(read_table(unit + '_usable_index.tsv.gz', indexed=True), audit.loc[usable]))
            check(unit + '_quarantine_exact_not_placeholder', same_table(read_table(unit + '_quarantine_index.tsv', indexed=True), audit.loc[~usable]))
            if not is_de and unit != 'baseline_atac':
                check(unit + '_every_count_row_canonical_fingerprint_verified', np.array_equal(audit.count_row_SHA256, fingerprints))
                check(unit + '_every_MT_value_recomputed_on_released_genes', np.array_equal(feature.MT_symbol_flag, var.index.astype(str).str.startswith('MT-')) and np.array_equal(audit.MT_counts_in_published_genes, mt) and np.allclose(audit.MT_percentage_in_published_genes, np.divide(100 * mt, total, out=np.zeros(len(total)), where=total > 0), rtol=0, atol=1e-10))
            entry = dataset.loc[unit]
            check(unit + '_dataset_index_points_to_real_originals_and_indices', entry.source_file == str(source_file) and entry.source_group == (group_name or '/') and entry.primary_matrix_slot == primary and entry.observations == len(obs) and entry.features == len(var) and entry.usable_index == str(OUT / (unit + '_usable_index.tsv.gz')) and entry.feature_index == str(OUT / (unit + '_feature_QC.tsv.gz')) and entry.matrix_handling == 'validated_original_read_only_with_explicit_indices')
            expected_assay = 'signed_DE_score' if is_de else 'ATAC_binary_accessibility' if unit == 'baseline_atac' else 'RNA_counts'
            expected_sum_semantics = 'signed_score_sum_not_UMI' if is_de else 'binary_accessible_peak_count' if unit == 'baseline_atac' else 'unnormalized_RNA_UMI_counts'
            check(unit + '_unit_summary_recomputed', units.loc[unit, 'observations'] == len(obs) and units.loc[unit, 'features'] == len(var) and units.loc[unit, 'usable_observations'] == int(usable.sum()) and units.loc[unit, 'quarantined_observations'] == int((~usable).sum()) and units.loc[unit, 'zero_information_rows'] == int((detected == 0).sum()) and units.loc[unit, 'all_zero_features_preserved'] == int((gene_n == 0).sum()) and units.loc[unit, 'assay'] == expected_assay and entry.assay == expected_assay and audit.assay_semantics.eq(expected_assay).all() and units.loc[unit, 'primary_sum_semantics'] == expected_sum_semantics)
    check('all_source_assay_and_obsm_matrices_accounted_for', set(recomputed_matrix_slots) == set(matrices.index) and len(recomputed_matrix_slots) == export['matrix_components_scanned'])
    finish(baseline, export, originals, audits, native_obs, native_var, cell_hashes, meta_differences)


def finish(baseline, export, originals, audits, native_obs, native_var, cell_hashes, differences):
    original, updated = audits['sc_original'], audits['sc_reannotated']
    check('release_size_subset_and_unmodified_gene_axes', len(original) == 301785 and len(updated) == 298087 and updated.index.isin(original.index).all() and native_var['sc_original'].index.equals(native_var['sc_reannotated'].index))
    check('every_shared_count_row_identical_across_release_versions', original.loc[updated.index, 'count_row_SHA256'].equals(updated.count_row_SHA256))
    check('every_shared_CMO_count_identical_across_release_versions', cell_hashes['sc_original'].loc[updated.index].equals(cell_hashes['sc_reannotated']))
    for field in ['split', 'donor_id', 'dose_uM', 'timepoint_hr', 'plate_name', 'well', 'sm_name']:
        check('shared_source_field_preserved_' + field, original.loc[updated.index, field].astype(str).tolist() == updated[field].astype(str).tolist())
    correspondence = read_table('sc_release_correspondence.tsv.gz', indexed=True)
    shared = original.index.isin(updated.index)
    expected_updated_rows = updated.source_row.reindex(original.index)
    saved_updated_rows = pd.to_numeric(correspondence.reannotated_source_row, errors='coerce')
    check('release_correspondence_every_row_exact', correspondence.index.equals(original.index) and np.array_equal(correspondence.original_source_row, original.source_row) and np.array_equal(saved_updated_rows, expected_updated_rows, equal_nan=True) and np.array_equal(correspondence.present_in_author_reannotated_release, shared) and np.array_equal(correspondence.count_identity_verified_when_shared, shared))
    check('release_correspondence_preserves_both_annotations', correspondence.original_cell_type.tolist() == original.cell_type.astype(str).tolist() and correspondence.loc[updated.index, 'author_reannotated_cell_type'].tolist() == updated.cell_type.astype(str).tolist() and correspondence.original_split.tolist() == original['split'].astype(str).tolist())
    excluded = read_table('original_only_after_author_reannotation.tsv.gz', indexed=True)
    check('author_removed3698_are_traceable_without_deleting_originals', len(excluded) == 3698 and same_table(excluded, correspondence.loc[~shared]))
    transitions = pd.crosstab(original.loc[updated.index, 'cell_type'].rename('original_cell_type'), updated.cell_type.rename('author_reannotated_cell_type')).stack()
    expected_transitions = transitions.loc[transitions > 0].reset_index()
    expected_transitions.columns = ['original_cell_type', 'author_reannotated_cell_type', 'cells']
    actual_transitions = read_table('annotation_transition_summary.tsv')
    check('every_annotation_transition_count_recomputed', same_table(actual_transitions, expected_transitions))
    differences['shared_cells_with_changed_author_cell_type'] = int((original.loc[updated.index, 'cell_type'].astype(str).to_numpy() != updated.cell_type.astype(str).to_numpy()).sum())
    if 'cell_type_orig' in updated:
        differences['reannotated_cell_type_orig_vs_original_release_differences'] = int((original.loc[updated.index, 'cell_type'].astype(str).to_numpy() != updated.cell_type_orig.astype(str).to_numpy()).sum())
    paired = read_table('multiome_paired_index.tsv.gz', indexed=True)
    check('baseline_and_perturbation_share_same_three_publisher_donorIDs', set(original.donor_id) == set(paired.donor_id) == {'Donor 1', 'Donor 2', 'Donor 3'})
    with h5py.File(SOURCE / 'multiome_counts.h5mu', 'r') as handle:
        root_obs, root_var = ad.io.read_elem(handle['obs']), ad.io.read_elem(handle['var'])
        check('baseline_pair26335_native_root_metadata_exact', len(paired) == 26335 and native_columns_preserved(root_obs, paired) and np.array_equal(paired.source_root_row, np.arange(len(paired))))
        for modality in ['rna', 'atac']:
            obsmap = ad.io.read_elem(handle['obsmap/' + modality]).ravel()
            varmap = ad.io.read_elem(handle['varmap/' + modality]).ravel()
            native = native_var['baseline_' + modality]
            check(modality + '_every_native_cell_pair_map_exact', np.array_equal(obsmap, np.arange(1, len(paired) + 1)) and np.array_equal(paired[modality + '_source_row'], obsmap - 1) and paired.index.equals(audits['baseline_' + modality].index))
            present = varmap > 0
            check(modality + '_every_native_feature_map_exact', np.array_equal(np.sort(varmap[present]), np.arange(1, len(native) + 1)) and root_var.index[present].equals(native.index[varmap[present] - 1]))
            check(modality + '_paired_QC_exact', np.array_equal(paired[modality + '_technical_qc_usable'], audits['baseline_' + modality].technical_qc_usable))
        check('both_modalities_same_cell_identity_and_annotation', audits['baseline_rna'].cell_uid.equals(audits['baseline_atac'].cell_uid) and audits['baseline_rna'].cell_type.equals(audits['baseline_atac'].cell_type))
        check('paired_usable_exact_and_reference_role_explicit', np.array_equal(paired.paired_usable, paired.rna_technical_qc_usable & paired.atac_technical_qc_usable) and paired.reference_role.eq('same_three_donors_baseline_not_independent_patient_test').all())
    combined_de = pd.concat([audits['de_train'], audits['de_test']])
    check('869_disjoint_DE_condition_units_and_fixed18211_genes', len(combined_de) == 869 and combined_de.condition_uid.is_unique and native_var['de_train'].index.equals(native_var['de_test'].index) and len(native_var['de_train']) == 18211 and native_var['de_train'].index.isin(native_var['sc_original'].index).all())
    idmap = pd.read_csv(SOURCE / 'workflow-resources/neurips-2023-kaggle/id_map.csv')
    mapped = read_table('official_test_id_mapping.tsv')
    check('official_test_idmap_all255_rows_and_fields_exact', len(idmap) == 255 and np.array_equal(idmap.id, np.arange(255)) and all(np.array_equal(mapped[name], idmap[name]) for name in idmap.columns) and mapped.cell_type.tolist() == audits['de_test'].cell_type.tolist() and mapped.sm_name.tolist() == audits['de_test'].sm_name.tolist() and mapped.condition_uid.tolist() == audits['de_test'].condition_uid.tolist() and mapped.official_split.tolist() == audits['de_test'].official_split.tolist() and np.array_equal(mapped.source_de_test_row, np.arange(255)))
    keys = ['cell_type', 'sm_name', 'dose_uM', 'split']
    counts = native_obs['sc_original'].groupby(keys, observed=True).agg(source_cells=('donor_id', 'size'), source_donors=('donor_id', 'nunique')).reset_index()
    expected = combined_de[keys + ['condition_uid']].merge(counts, how='left', on=keys, validate='one_to_one')
    conditions = read_table('DE_condition_source_summary.tsv')
    check('every_DE_condition_link_uses_original_six_type_labels', all(np.array_equal(conditions[name], expected[name]) for name in keys + ['condition_uid']) and np.array_equal(pd.to_numeric(conditions.source_cells, errors='coerce'), expected.source_cells, equal_nan=True) and np.array_equal(pd.to_numeric(conditions.source_donors, errors='coerce'), expected.source_donors, equal_nan=True) and np.array_equal(conditions.source_condition_available, expected.source_cells.notna()) and conditions.source_annotation_for_linking.eq('sc_counts_original_six_cell_types').all())
    moa = pd.read_csv(SOURCE / 'moa_annotations.csv', keep_default_na=False)
    compound = read_table('compound_annotation_QC.tsv')
    check('all147_source_compound_annotations_preserved', len(moa) == len(compound) == 147 and all(np.array_equal(moa[name].astype(str), compound[name].astype(str)) for name in moa.columns) and np.array_equal(compound.source_row, np.arange(147)))
    check('missing_and_ambiguous_drug_annotations_flagged_not_invented', np.array_equal(compound.sm_lincs_id_duplicate, moa.sm_lincs_id.duplicated(keep=False)) and np.array_equal(compound.sm_name_duplicate, moa.sm_name.duplicated(keep=False)) and np.array_equal(compound.moa_annotation_missing, moa.moa.eq('')) and np.array_equal(compound.moa_citation_missing, moa.moa_citation.eq('')))
    controls = read_table('publisher_control_definition.tsv')
    expected_controls = native_obs['sc_original'].loc[native_obs['sc_original']['split'].eq('control')].groupby(['sm_name', 'dose_uM'], observed=True).size().rename('cells').reset_index()
    check('publisher_control_split_keeps_all_drug_dose_definitions', same_table(controls, expected_controls))
    policy = json.loads((OUT / 'qc_policy.json').read_text())
    check('modality_count_and_disease_boundaries_documented', policy['count_sources'] == {'sc_original': 'layers/counts', 'sc_reannotated': 'raw/X', 'baseline_rna': 'mod/rna/raw/X'} and 'not HCM' in policy['source_context'] and 'not raw fragment counts' in policy['ATAC_boundary'] and 'CMO' in policy['hashing_boundary'])
    check('unrequested_processing_not_applied', {'HCM_gene500to5000', 'MT_lt5percent', 'doublet_recalling', 'new_DE', 'model_training'}.issubset(policy['not_applied']))
    check('source_baseline_copy_exact', json.loads((OUT / 'ORIGINAL_DOWNLOAD_BASELINE.json').read_text()) == baseline)
    check('export_counts_are_honest_and_not_double_counted', export['primary_unique_perturbation_cells'] == len(original) and export['paired_baseline_cells'] == len(paired) and export['author_reannotated_subset_cells'] == len(updated) and export['DE_condition_rows'] == len(combined_de) and export['DE_conditions_missing_source_match'] == int(expected.source_cells.isna().sum()))
    for name, record in export['hashes'].items():
        path = OUT / name
        check('payload_bytes_' + name, path.stat().st_size == record['bytes'])
        check('payload_MD5_' + name, digest(path, 'md5') == record['MD5'])
        check('payload_SHA256_' + name, digest(path) == record['SHA256'])
    prepare_exit = (WORK / 'prepare_exit_code.txt').read_text().strip()
    check('prepare_exited_successfully', prepare_exit == '0')
    schema = {unit: {'obs': {name: str(frame[name].dtype) for name in frame},
        'var': {name: str(native_var[unit][name].dtype) for name in native_var[unit]}}
        for unit, frame in native_obs.items()}
    (OUT / 'SOURCE_METADATA_DTYPES.json').write_text(json.dumps({'comparison': 'Numeric CSV fields restore the native source dtype before exact comparison; float32 is not compared to an expanded float64 representation. Original HDF5 metadata retains native dtypes.', 'units': schema}, indent=2) + '\n')
    result = {'dataset': 'NeurIPS2023', 'status': 'complete_verified', 'verified_utc': datetime.now(timezone.utc).isoformat(),
        'checks_passed': len(CHECKS), 'checks_failed': 0, 'checks': CHECKS,
        'scope': 'All8 source hashes, all native assay/raw/layer/obsm matrix values, full row/feature metrics and source metadata, official splits, author release correspondence, every multiome pair/feature map and source DE layer semantics. Publisher biological QC is retained, not recertified.',
        'primary_unique_perturbation_cells': len(original), 'paired_baseline_cells': len(paired),
        'author_reannotated_subset_cells_not_additional_cells': len(updated), 'DE_condition_rows_not_cells': len(combined_de),
        'published_matrix_components': export['matrix_components_scanned'], 'source_metadata_differences_preserved': differences,
        'metadata_numeric_comparison': 'exact_after_restoring_source_dtype_not_tolerance_relaxation',
        'quarantined_by_unit': {key: int((~value.technical_qc_usable).sum()) for key, value in audits.items()}}
    (OUT / 'INDEPENDENT_VERIFICATION.json').write_text(json.dumps(result, indent=2) + '\n')
    manifest = []
    for row in originals.reset_index().to_dict('records'):
        manifest.append({'dataset': 'NeurIPS2023', 'artifact_kind': 'unchanged_original', 'file': row['file'],
            'server_path': row['server_path'], 'bytes': row['bytes'], 'MD5': row['MD5'], 'SHA256': row['SHA256']})
    for path in sorted(OUT.iterdir()):
        if path.is_file() and path.name != 'OUTPUT_MANIFEST.tsv':
            manifest.append({'dataset': 'NeurIPS2023', 'artifact_kind': 'validated_index_or_audit', 'file': path.name,
                'server_path': str(path), 'bytes': path.stat().st_size, 'MD5': digest(path, 'md5'), 'SHA256': digest(path)})
    pd.DataFrame(manifest).to_csv(OUT / 'OUTPUT_MANIFEST.tsv', sep='\t', index=False)
    gate = {'dataset': 'NeurIPS2023', 'status': 'complete_verified', 'verified_utc': result['verified_utc'],
        'checks_passed': len(CHECKS), 'checks_failed': 0, 'all8_originals_rehashed': True,
        'all_native_assay_and_obsm_matrices_fully_scanned': True, 'official_splits_and_all_pair_maps_preserved': True,
        'source_originals_read_only_with_real_usable_indices': True, 'reannotated_subset_not_double_counted': True,
        'disease_labels_assigned': False, 'biological_response_or_singlet_certified': False, 'prepare_exit_code': prepare_exit,
        'output_directory': str(OUT), 'verification_SHA256': digest(OUT / 'INDEPENDENT_VERIFICATION.json'),
        'output_manifest_SHA256': digest(OUT / 'OUTPUT_MANIFEST.tsv')}
    (WORK / 'FINAL_COMPLETION_GATE.json').write_text(json.dumps(gate, indent=2) + '\n')
    print(json.dumps({key: value for key, value in result.items() if key != 'checks'}, indent=2), flush=True)


if __name__ == '__main__':
    try:
        main()
    except BaseException as error:
        (WORK / 'VERIFICATION_FAILURE.json').write_text(json.dumps({'status': 'failed', 'error': str(error),
            'traceback': traceback.format_exc(), 'checks': CHECKS}, indent=2) + '\n')
        raise
