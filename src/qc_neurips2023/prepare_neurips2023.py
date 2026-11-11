"""Validate the publisher's curated assays; publish faithful row/feature indices."""
from datetime import datetime, timezone
from pathlib import Path
import hashlib
import json
import shutil
import anndata as ad
import h5py
import numpy as np
import pandas as pd
from scipy.sparse import csr_matrix

ROOT = Path.home() / 'xiaoe-internship-data'
SOURCE = ROOT / 'benchmarks/competitions/neurips2023'
WORK = ROOT / '_processing/all_data_cleaning_20261111/NeurIPS2023'
OUT = ROOT / 'processed/NeurIPS2023_qc_20261111'
UNITS = [
    ('sc_original', 'sc_counts.h5ad', '', 'layers/counts', 'RNA_counts'),
    ('sc_reannotated', 'sc_counts_reannotated_with_counts.h5ad', '', 'raw/X', 'RNA_counts'),
    ('baseline_rna', 'multiome_counts.h5mu', 'mod/rna', 'raw/X', 'RNA_counts'),
    ('baseline_atac', 'multiome_counts.h5mu', 'mod/atac', 'layers/binarized', 'ATAC_binary_accessibility'),
    ('de_train', '2023-09-14_kaggle_upload/2023-09-12_de_by_cell_type_train.h5ad', '', 'X', 'signed_DE_score'),
    ('de_test', '2023-09-14_kaggle_upload/2023-09-12_de_by_cell_type_test.h5ad', '', 'X', 'signed_DE_score'),
]
BLOCK = 2048


def digest(path, algorithm='sha256'):
    with path.open('rb') as handle:
        return hashlib.file_digest(handle, algorithm).hexdigest()


def read_group(handle, name):
    return handle[name] if name else handle


def frame_ok(frame, columns):
    return (frame[columns].notna().all(axis=1) & frame[columns].astype(str).apply(lambda col: col.str.strip().ne('')).all(axis=1)).to_numpy(dtype=bool, copy=True)


def specs_for(group, primary, assay):
    specs = {'X': 'signed_DE_score' if assay == 'signed_DE_score' else 'processed_nonnegative'}
    for name in group.get('layers', {}):
        specs['layers/' + name] = {'counts': 'counts', 'binarized': 'binary', 'P.Value': 'pvalue',
            'adj.P.Val': 'pvalue', 'is_de': 'binary', 'is_de_adj': 'binary', 'logFC': 'signed_logFC'}[name]
    if 'raw' in group:
        specs['raw/X'] = 'counts'
    for name in group.get('obsm', {}):
        specs['obsm/' + name] = ('counts' if name == 'protein_counts' else 'pvalue' if name == 'ora_pvals'
            else 'annotation_score' if name == 'ora_estimate' else 'signed_hashing_CLR' if name == 'HTO_clr' else 'author_embedding')
    assert primary in specs
    return specs


def blocks(matrix):
    if isinstance(matrix, h5py.Group) and matrix.attrs.get('encoding-type') == 'dataframe':
        table = ad.io.read_elem(matrix)
        for start in range(0, len(table), BLOCK):
            yield start, csr_matrix(table.iloc[start:start + BLOCK].to_numpy(dtype=np.float64))
    elif isinstance(matrix, h5py.Group):
        assert matrix.attrs['encoding-type'] == 'csr_matrix'
        shape = tuple(matrix.attrs['shape'])
        pointer = matrix['indptr'][:]
        assert len(pointer) == shape[0] + 1 and pointer[0] == 0
        assert pointer[-1] == len(matrix['data']) == len(matrix['indices']) and (np.diff(pointer) >= 0).all()
        for start in range(0, shape[0], BLOCK):
            stop = min(start + BLOCK, shape[0])
            low, high = int(pointer[start]), int(pointer[stop])
            values, columns = matrix['data'][low:high], matrix['indices'][low:high]
            assert not len(columns) or (columns.min() >= 0 and columns.max() < shape[1])
            yield start, csr_matrix((values.astype(np.float64), columns, pointer[start:stop + 1] - low), shape=(stop - start, shape[1]))
    else:
        assert matrix.ndim == 2
        for start in range(0, matrix.shape[0], BLOCK):
            yield start, csr_matrix(matrix[start:start + BLOCK].astype(np.float64))


def scan(matrix, semantics, mt_mask=None, fingerprint=False):
    if isinstance(matrix, h5py.Group):
        if matrix.attrs.get('encoding-type') == 'dataframe':
            table = ad.io.read_elem(matrix)
            shape = table.shape
        else:
            shape = tuple(matrix.attrs['shape'])
    else:
        shape = matrix.shape
    rows, columns = shape
    row_sum, detected, mt = np.zeros(rows), np.zeros(rows, dtype=np.int32), np.zeros(rows)
    feature_sum, feature_n = np.zeros(columns), np.zeros(columns, dtype=np.int64)
    fingerprints = np.full(rows, '', dtype=object)
    minimum, maximum = np.inf, -np.inf
    nonzero, negative, fractional, stored_zeros, noncanonical = 0, 0, 0, 0, 0
    for start, chunk in blocks(matrix):
        data = chunk.data
        assert np.isfinite(data).all(), (matrix.name, start, 'nonfinite')
        if semantics in {'counts', 'binary', 'processed_nonnegative', 'annotation_score', 'pvalue'}:
            assert (data >= 0).all(), (matrix.name, start, 'negative')
        if semantics in {'counts', 'binary'}:
            assert np.equal(data, np.rint(data)).all() and (data < np.iinfo(np.int64).max).all()
        if semantics in {'binary', 'pvalue'}:
            assert (data <= 1).all(), (matrix.name, start, 'outside_unit_interval')
        if len(data):
            minimum, maximum = min(minimum, float(data.min())), max(maximum, float(data.max()))
        fractional += int(np.count_nonzero(data != np.rint(data)))
        negative += int(np.count_nonzero(data < 0))
        stored_zeros += int(np.count_nonzero(data == 0))
        noncanonical += int(not chunk.has_canonical_format)
        chunk.sum_duplicates()
        chunk.sort_indices()
        chunk.eliminate_zeros()
        stop = start + chunk.shape[0]
        row_sum[start:stop] = np.asarray(chunk.sum(axis=1)).ravel()
        detected[start:stop] = chunk.getnnz(axis=1)
        feature_sum += np.bincount(chunk.indices, weights=chunk.data, minlength=columns)
        feature_n += np.bincount(chunk.indices, minlength=columns)
        nonzero += chunk.nnz
        if mt_mask is not None:
            mt[start:stop] = np.asarray(chunk[:, mt_mask].sum(axis=1)).ravel()
        if fingerprint:
            for local in range(chunk.shape[0]):
                low, high = chunk.indptr[local:local + 2]
                hasher = hashlib.sha256()
                hasher.update(np.asarray(chunk.indices[low:high], dtype='<i4').tobytes())
                hasher.update(np.asarray(chunk.data[low:high], dtype='<i8').tobytes())
                fingerprints[start + local] = hasher.hexdigest()
    # Dense/dataframe zeros were dropped during CSR construction; zero counts below
    # are semantic zeros, while stored_zeros specifically describes sparse storage.
    info = {'matrix_slot': matrix.name, 'rows': rows, 'features': columns, 'semantics': semantics,
        'nonzero_values': nonzero, 'zero_entries': int(rows * columns - nonzero),
        'negative_values': negative, 'fractional_values': fractional,
        'stored_sparse_zeros': stored_zeros, 'noncanonical_blocks': noncanonical,
        'minimum_stored_nonzero_or_sparse_value': float(minimum) if np.isfinite(minimum) else 0.0,
        'maximum_stored_nonzero_or_sparse_value': float(maximum) if np.isfinite(maximum) else 0.0,
        'zero_information_rows': int((detected == 0).sum()), 'all_zero_features': int((feature_n == 0).sum()),
        'row_sum_total': float(row_sum.sum())}
    return info, (row_sum, detected, mt, fingerprints, feature_sum, feature_n)


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    assert not (OUT / 'EXPORT_SUCCESS.json').exists(), 'Inspect successful outputs before rerunning'
    baseline = json.loads((WORK / 'ORIGINAL_DOWNLOAD_BASELINE.json').read_text())
    assert len(baseline) == 8
    originals = []
    for item in baseline:
        path = Path(item['server_path'])
        md5, sha = digest(path, 'md5'), digest(path)
        assert path.stat().st_size == int(item['bytes']) and md5 == item['md5'] and sha == item['sha256']
        originals.append({'file': path.name, 'server_path': str(path), 'bytes': path.stat().st_size,
            'MD5': md5, 'SHA256': sha, 'source_url': item['source_url']})
        print('Original file verified:', path.name, flush=True)
    pd.DataFrame(originals).to_csv(OUT / 'ORIGINAL_FILES.tsv', sep='\t', index=False)
    indices, variables, author_obs, hashing = {}, {}, {}, {}
    summaries, matrix_summaries, dataset_index, de_rules = [], [], [], []
    for unit, relative, group_name, primary, assay in UNITS:
        path = SOURCE / relative
        with h5py.File(path, 'r') as handle:
            group = read_group(handle, group_name)
            obs, var = ad.io.read_elem(group['obs']), ad.io.read_elem(group['var'])
            assert obs.index.is_unique and var.index.is_unique
            assert (obs.index.astype(str).str.strip() != '').all() and (var.index.astype(str).str.strip() != '').all()
            if primary.startswith('raw/'):
                assert ad.io.read_elem(group['raw']['var']).index.equals(var.index)
            is_de = assay == 'signed_DE_score'
            mt_mask = var.index.astype(str).str.startswith('MT-') if assay == 'RNA_counts' else None
            primary_metrics = None
            for slot, semantics in specs_for(group, primary, assay).items():
                matrix = group[slot]
                if isinstance(matrix, h5py.Group) and matrix.attrs.get('encoding-type') == 'dataframe':
                    table = ad.io.read_elem(matrix)
                    assert table.index.equals(obs.index) and table.columns.is_unique
                    if slot == 'obsm/protein_counts':
                        assert table.columns.str.startswith('CMO_').all()
                        hashing[unit] = table
                result, metrics = scan(matrix, semantics, mt_mask if slot == primary else None,
                    fingerprint=(slot == primary and assay == 'RNA_counts'))
                assert result['rows'] == len(obs)
                if slot in {'X', primary} or slot.startswith('layers/'):
                    assert result['features'] == len(var)
                result.update(unit=unit, source_file=str(path), relative_slot=slot,
                    assay_interpretation='cell_hashing_CMO_counts' if slot == 'obsm/protein_counts' else semantics)
                matrix_summaries.append(result)
                if slot == primary:
                    primary_metrics = metrics
                print('Full matrix scan:', unit, slot, result['rows'], result['features'], flush=True)
            total, detected, mt, fingerprints, feature_sum, feature_n = primary_metrics
            index = obs.copy()
            index.index = index.index.astype(str)
            index.index.name = 'source_observation_id'
            index['source_row'] = np.arange(len(index))
            index['source_file'] = str(path)
            index['source_group'] = group_name or '/'
            index['primary_matrix_slot'] = primary
            index['assay_semantics'] = assay
            if unit.startswith('sc_'):
                required = ['donor_id', 'cell_type', 'split', 'sm_name', 'dose_uM', 'timepoint_hr', 'plate_name', 'well']
                valid = frame_ok(obs, required)
                valid &= (np.isfinite(obs.dose_uM) & (obs.dose_uM >= 0) & np.isfinite(obs.timepoint_hr)).to_numpy()
                assert set(obs.split) == {'train', 'control', 'public_test', 'private_test'}
                index['cell_uid'] = 'NeurIPS2023:perturbation:' + index.index
                index['official_split'] = obs.split.astype(str).to_numpy()
                index['canonical_unique_cell_source'] = unit == 'sc_original'
                index['release_version'] = unit
                index['control_definition'] = 'publisher_control_split_includes_DMSO_and_low_dose_positive_controls'
            elif unit.startswith('baseline_'):
                root_obs = ad.io.read_elem(handle['obs'])
                assert root_obs.index.equals(obs.index)
                index['source_donor_id'] = root_obs.donor_id.astype(str).to_numpy()
                valid = frame_ok(obs, ['cell_type']) & frame_ok(root_obs, ['donor_id'])
                index['cell_uid'] = 'NeurIPS2023:baseline_multiome:' + index.index
                index['official_split'] = 'baseline_reference_no_cell_level_challenge_split'
                index['canonical_unique_cell_source'] = unit == 'baseline_rna'
                index['paired_modality'] = 'baseline_atac' if unit == 'baseline_rna' else 'baseline_rna'
            else:
                valid = frame_ok(obs, ['cell_type', 'sm_lincs_id', 'sm_name', 'SMILES', 'dose_uM', 'split'])
                index['condition_uid'] = 'NeurIPS2023:condition:' + obs[['cell_type', 'sm_lincs_id', 'dose_uM']].astype(str).agg('||'.join, axis=1).to_numpy()
                assert index.condition_uid.is_unique
                index['official_split'] = obs.split.astype(str).to_numpy()
                assert set(obs.split) == ({'train', 'control'} if unit == 'de_train' else {'public_test', 'private_test'})
                # Validate the published statistic, without estimating new DE results.
                score, p, fc = group['X'][:], group['layers/P.Value'][:], group['layers/logFC'][:]
                positive_p = p > 0
                predicted = np.zeros_like(p)
                predicted[positive_p] = np.sign(fc[positive_p]) * -np.log10(p[positive_p])
                error = np.abs(score[positive_p] - predicted[positive_p])
                de_rules.append({'unit': unit, 'positive_P_entries': int(positive_p.sum()),
                    'zero_P_entries': int((~positive_p).sum()), 'maximum_score_formula_error_positive_P': float(error.max()) if len(error) else 0.0,
                    'is_de_vs_P_lt_0_05_mismatches': int(np.count_nonzero(group['layers/is_de'][:] != (p < 0.05))),
                    'is_de_adj_vs_adj_P_lt_0_05_mismatches': int(np.count_nonzero(group['layers/is_de_adj'][:] != (group['layers/adj.P.Val'][:] < 0.05)))})
                assert np.allclose(score[positive_p], predicted[positive_p], rtol=1e-7, atol=1e-7)
            index['primary_row_sum'] = total
            index['primary_nonzero_features'] = detected
            index['primary_zero_information'] = detected == 0
            index['technical_qc_usable'] = valid & ((detected > 0) if not is_de else True)
            index['qc_reason'] = np.select([~valid, (detected == 0) & (not is_de)],
                ['incomplete_required_source_metadata', 'zero_information_cell'], default='retained_publisher_curated_release')
            index['training_role'] = np.where(index.official_split.eq('train'), 'official_training',
                np.where(index.official_split.eq('control'), 'publisher_control_reference',
                np.where(index.official_split.eq('public_test'), 'held_out_public_test',
                np.where(index.official_split.eq('private_test'), 'held_out_private_test', 'paired_baseline_reference'))))
            index['disease_label_assigned'] = False
            if assay == 'RNA_counts':
                index['count_row_SHA256'] = fingerprints
                index['MT_counts_in_published_genes'] = mt
                index['MT_percentage_in_published_genes'] = np.divide(100 * mt, total, out=np.zeros(len(index)), where=total > 0)
            if unit.startswith('sc_'):
                tags = hashing[unit]
                index['cell_hashing_CMO_total'] = tags.sum(axis=1).to_numpy()
                index['cell_hashing_CMO_nonzero_features'] = (tags > 0).sum(axis=1).to_numpy()
            feature = var.copy()
            feature.index = feature.index.astype(str)
            feature.index.name = 'source_feature_id'
            feature['source_column'] = np.arange(len(var))
            feature['primary_feature_sum'] = feature_sum
            feature['primary_nonzero_observations'] = feature_n
            feature['all_zero_in_primary_matrix'] = feature_n == 0
            if mt_mask is not None:
                feature['MT_symbol_flag'] = mt_mask
            if unit == 'baseline_atac':
                positions = feature.index.to_series().str.extract(r'^([^:]+):(\d+)-(\d+)$')
                assert positions.notna().all().all() and (positions[2].astype(int) > positions[1].astype(int)).all()
                feature['peak_coordinate_valid'] = True
            indices[unit], variables[unit], author_obs[unit] = index, feature, obs
            index.to_csv(OUT / (unit + '_observation_QC.tsv.gz'), sep='\t')
            index.loc[index.technical_qc_usable].to_csv(OUT / (unit + '_usable_index.tsv.gz'), sep='\t')
            index.loc[~index.technical_qc_usable].to_csv(OUT / (unit + '_quarantine_index.tsv'), sep='\t')
            feature.to_csv(OUT / (unit + '_feature_QC.tsv.gz'), sep='\t')
            summary = {'unit': unit, 'assay': assay, 'observations': len(obs), 'features': len(var),
                'usable_observations': int(index.technical_qc_usable.sum()), 'quarantined_observations': int((~index.technical_qc_usable).sum()),
                'zero_information_rows': int((detected == 0).sum()), 'all_zero_features_preserved': int((feature_n == 0).sum()),
                'donors': int(obs.donor_id.nunique()) if 'donor_id' in obs else int(index.source_donor_id.nunique()) if 'source_donor_id' in index else None,
                'primary_matrix_slot': primary, 'primary_sum_semantics': 'signed_score_sum_not_UMI' if is_de else 'binary_accessible_peak_count' if unit == 'baseline_atac' else 'unnormalized_RNA_UMI_counts'}
            summaries.append(summary)
            dataset_index.append({'unit': unit, 'source_file': str(path), 'source_group': group_name or '/',
                'primary_matrix_slot': primary, 'assay': assay, 'usable_index': str(OUT / (unit + '_usable_index.tsv.gz')),
                'feature_index': str(OUT / (unit + '_feature_QC.tsv.gz')), 'observations': len(obs), 'features': len(var),
                'matrix_handling': 'validated_original_read_only_with_explicit_indices'})
            print('Unit prepared:', json.dumps(summary), flush=True)
    finish(indices, variables, author_obs, hashing, summaries, matrix_summaries, dataset_index, de_rules)


def finish(indices, variables, author_obs, hashing, summaries, matrix_summaries, dataset_index, de_rules):
    original, updated = indices['sc_original'], indices['sc_reannotated']
    assert variables['sc_original'].index.equals(variables['sc_reannotated'].index)
    assert updated.index.isin(original.index).all() and len(original) - len(updated) == 3698
    assert original.loc[updated.index, 'count_row_SHA256'].equals(updated.count_row_SHA256)
    assert (original.loc[updated.index, 'split'].astype(str).to_numpy() == updated['split'].astype(str).to_numpy()).all()
    metadata_columns = ['donor_id', 'dose_uM', 'timepoint_hr', 'plate_name', 'well', 'sm_name']
    assert all(original.loc[updated.index, name].astype(str).tolist() == updated[name].astype(str).tolist() for name in metadata_columns)
    assert hashing['sc_original'].columns.equals(hashing['sc_reannotated'].columns)
    assert np.array_equal(hashing['sc_original'].loc[updated.index].to_numpy(), hashing['sc_reannotated'].to_numpy())
    release = pd.DataFrame({'original_source_row': original.source_row,
        'reannotated_source_row': updated.source_row.reindex(original.index).astype('Int64'),
        'present_in_author_reannotated_release': original.index.isin(updated.index),
        'original_cell_type': original.cell_type.astype(str), 'author_reannotated_cell_type': updated.cell_type.astype(str).reindex(original.index),
        'original_split': original['split'].astype(str), 'count_source_for_unique_cells': 'sc_counts.h5ad:layers/counts'})
    release['count_identity_verified_when_shared'] = release.present_in_author_reannotated_release
    release.to_csv(OUT / 'sc_release_correspondence.tsv.gz', sep='\t', index_label='source_cell_id')
    transition = release.loc[release.present_in_author_reannotated_release].groupby(['original_cell_type', 'author_reannotated_cell_type'], observed=True).size().rename('cells').reset_index()
    transition.to_csv(OUT / 'annotation_transition_summary.tsv', sep='\t', index=False)
    release.loc[~release.present_in_author_reannotated_release].to_csv(OUT / 'original_only_after_author_reannotation.tsv.gz', sep='\t', index_label='source_cell_id')
    with h5py.File(SOURCE / 'multiome_counts.h5mu', 'r') as handle:
        root_obs, root_var = ad.io.read_elem(handle['obs']), ad.io.read_elem(handle['var'])
        obsmap, varmap = ad.io.read_elem(handle['obsmap']), ad.io.read_elem(handle['varmap'])
        paired = root_obs.copy()
        paired.index.name = 'source_cell_id'
        paired['source_root_row'] = np.arange(len(paired))
        for modality in ['rna', 'atac']:
            unit = 'baseline_' + modality
            mapping = np.asarray(obsmap[modality]).ravel()
            feature_mapping = np.asarray(varmap[modality]).ravel()
            assert np.array_equal(mapping, np.arange(1, len(paired) + 1))
            assert indices[unit].index.equals(root_obs.index)
            selected = feature_mapping > 0
            assert np.array_equal(np.sort(feature_mapping[selected]), np.arange(1, len(variables[unit]) + 1))
            assert root_var.index[selected].equals(variables[unit].index[feature_mapping[selected] - 1])
            paired[modality + '_source_row'] = mapping - 1
            paired[modality + '_technical_qc_usable'] = indices[unit].technical_qc_usable.to_numpy()
        assert indices['baseline_rna'].cell_uid.equals(indices['baseline_atac'].cell_uid)
        assert indices['baseline_rna'].cell_type.astype(str).equals(indices['baseline_atac'].cell_type.astype(str))
        paired['paired_usable'] = paired.rna_technical_qc_usable & paired.atac_technical_qc_usable
        paired['reference_role'] = 'same_three_donors_baseline_not_independent_patient_test'
        paired.to_csv(OUT / 'multiome_paired_index.tsv.gz', sep='\t')
    de_combined = pd.concat([indices['de_train'], indices['de_test']])
    assert de_combined.condition_uid.is_unique
    assert variables['de_train'].index.equals(variables['de_test'].index)
    assert variables['de_train'].index.isin(variables['sc_original'].index).all()
    idmap = pd.read_csv(SOURCE / 'workflow-resources/neurips-2023-kaggle/id_map.csv')
    test = indices['de_test']
    assert np.array_equal(idmap.id, np.arange(len(test)))
    assert idmap.sm_name.tolist() == test.sm_name.astype(str).tolist() and idmap.cell_type.tolist() == test.cell_type.astype(str).tolist()
    condition_keys = ['cell_type', 'sm_name', 'dose_uM', 'split']
    counts = author_obs['sc_original'].groupby(condition_keys, observed=True).agg(source_cells=('donor_id', 'size'), source_donors=('donor_id', 'nunique')).reset_index()
    conditions = de_combined[condition_keys + ['condition_uid']].merge(counts, on=condition_keys, how='left', validate='one_to_one')
    conditions['source_annotation_for_linking'] = 'sc_counts_original_six_cell_types'
    conditions['source_condition_available'] = conditions.source_cells.notna()
    conditions.to_csv(OUT / 'DE_condition_source_summary.tsv', sep='\t', index=False)
    idmap['official_split'] = test.official_split.to_numpy()
    idmap['source_de_test_row'] = test.source_row.to_numpy()
    idmap['condition_uid'] = test.condition_uid.to_numpy()
    idmap.to_csv(OUT / 'official_test_id_mapping.tsv', sep='\t', index=False)
    moa = pd.read_csv(SOURCE / 'moa_annotations.csv', keep_default_na=False)
    moa['source_row'] = np.arange(len(moa))
    moa['sm_lincs_id_duplicate'] = moa.sm_lincs_id.duplicated(keep=False)
    moa['sm_name_duplicate'] = moa.sm_name.duplicated(keep=False)
    moa['moa_annotation_missing'] = moa.moa.eq('')
    moa['moa_citation_missing'] = moa.moa_citation.eq('')
    moa.to_csv(OUT / 'compound_annotation_QC.tsv', sep='\t', index=False)
    pd.DataFrame(summaries).to_csv(OUT / 'UNIT_QC_SUMMARY.tsv', sep='\t', index=False)
    pd.DataFrame(matrix_summaries).to_csv(OUT / 'MATRIX_QC_SUMMARY.tsv', sep='\t', index=False)
    pd.DataFrame(dataset_index).to_csv(OUT / 'DATASET_INDEX.tsv', sep='\t', index=False)
    pd.DataFrame(de_rules).to_csv(OUT / 'published_DE_consistency.tsv', sep='\t', index=False)
    controls = author_obs['sc_original'].loc[author_obs['sc_original']['split'].eq('control')].groupby(['sm_name', 'dose_uM'], observed=True).size().rename('cells').reset_index()
    controls.to_csv(OUT / 'publisher_control_definition.tsv', sep='\t', index=False)
    policy = {'dataset': 'NeurIPS2023', 'version': 'curated_release_validation_v1',
        'source_context': 'PBMC small-molecule treatment for24hours and paired baseline RNA/ATAC from the same3 healthy donors; not HCM or cardiomyocytes.',
        'primary_sources': ['https://github.com/openproblems-bio/task_perturbation_prediction', 'https://www.kaggle.com/competitions/open-problems-single-cell-perturbations'],
        'count_sources': {'sc_original': 'layers/counts', 'sc_reannotated': 'raw/X', 'baseline_rna': 'mod/rna/raw/X'},
        'processed_expression': 'All X expression arrays kept at their native paths and fully scanned; fractional published X is not relabelled as counts or overwritten.',
        'ATAC_boundary': 'X is processed accessibility; layers/binarized is binary accessible peaks, not raw fragment counts. Source fragment-count/TSS/nucleosome metadata preserved; no invented raw matrix or fragment files.',
        'hashing_boundary': 'obsm/protein_counts contains12 CMO cell-hashing tags, not disease labels or a biological ADT protein panel. HTO_clr is signed transformed hashing signal.',
        'DE_boundary': 'X=sign(logFC)*-log10(P.Value) on positive-P entries, signed values are legitimate; rows are cell-type/compound/dose conditions, not individual cells. Keep all5 source layers and official public/private test split. Zero P entries, if any, are explicitly counted without inventing the publisher floor.',
        'annotation_versions': 'sc_reannotated is a298087-cell subset of301785 original cells. Shared canonical counts and CMO counts checked exactly, release metadata differences kept; original six-type annotations link official DE. Do not concatenate the versions as independent cells.',
        'technical_QC': 'Retain curated release; nonzero information and required source metadata for assay cells. Zero DE response remains valid. All genes/peaks retained, including zero features; no effect-dependent selection.',
        'splits_and_donors': 'Preserve original train/control/public_test/private_test. Official split is not a held-out-donor split; baseline shares the3 donors. Positive low-dose controls must not be collapsed into DMSO. No disease0/1 assignment.',
        'not_applied': ['HCM_gene500to5000', 'MT_lt5percent', 'EmptyDrops', 'doublet_recalling', 'CM_identity', 'normalization', 'imputation', 'integration', 'HVG_truncation', 'new_DE', 'tokenization', 'model_training'],
        'delivery': 'Validated original H5AD/H5MU read-only plus real usable/feature/paired/release indices; no redundant large expression copies. Biological response efficacy, original droplet calling and donor independence across other studies are not certified.'}
    (OUT / 'qc_policy.json').write_text(json.dumps(policy, indent=2) + '\n')
    for name in ['ORIGINAL_DOWNLOAD_BASELINE.json', 'INPUT_INSPECTION.json']:
        shutil.copy2(WORK / name, OUT / name)
    files = sorted(path for path in OUT.iterdir() if path.is_file())
    result = {'dataset': 'NeurIPS2023', 'status': 'exported_pending_independent_verification',
        'completed_utc': datetime.now(timezone.utc).isoformat(), 'unit_summaries': summaries,
        'matrix_components_scanned': len(matrix_summaries), 'primary_unique_perturbation_cells': len(original),
        'paired_baseline_cells': len(paired), 'author_reannotated_subset_cells': len(updated),
        'DE_condition_rows': len(de_combined), 'DE_conditions_missing_source_match': int((~conditions.source_condition_available).sum()),
        'hashes': {path.name: {'bytes': path.stat().st_size, 'MD5': digest(path, 'md5'), 'SHA256': digest(path)} for path in files}}
    (OUT / 'EXPORT_SUCCESS.json').write_text(json.dumps(result, indent=2) + '\n')
    print(json.dumps({key: value for key, value in result.items() if key != 'hashes'}, indent=2), flush=True)


if __name__ == '__main__':
    main()
