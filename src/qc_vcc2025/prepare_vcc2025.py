"""Stream all official counts; preserve original files and publish audited read indices."""
from collections import Counter
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
SOURCE = ROOT / 'benchmarks/competitions/vcc2025'
WORK = ROOT / '_processing/all_data_cleaning_20261111/VCC2025'
OUT = ROOT / 'processed/VCC2025_qc_20261111'
FILES = {'train': ('train', 'Training'), 'validation': ('validation', 'Validation'), 'test': ('test', 'Test')}
BLOCK = 1024


def checksum(path, algorithm='sha256'):
    with path.open('rb') as handle:
        return hashlib.file_digest(handle, algorithm).hexdigest()


def row_fingerprint(indices, counts):
    digest = hashlib.sha256()
    digest.update(np.asarray(indices, dtype='<i4').tobytes())
    digest.update(np.asarray(counts, dtype='<i8').tobytes())
    return digest.hexdigest()


def main():
    OUT.mkdir(exist_ok=True, parents=True)
    assert not (OUT / 'EXPORT_SUCCESS.json').exists(), 'Inspect successful outputs before rerunning'
    baseline = json.loads((WORK / 'ORIGINAL_DOWNLOAD_BASELINE.json').read_text())
    assert len(baseline) == 7
    input_records = []
    for row in baseline:
        path = Path(row['server_path'])
        md5 = checksum(path, 'md5')
        sha256 = checksum(path)
        assert path.stat().st_size == int(row['bytes']) and md5 == row['md5'], path
        if row.get('sha256'):
            assert sha256 == row['sha256'], path
        input_records.append({'file': row['file'], 'server_path': str(path), 'bytes': path.stat().st_size,
            'MD5': md5, 'SHA256': sha256, 'historical_SHA256_available': bool(row.get('sha256')),
            'source_url': row['source_url']})
        print('Original checksum verified:', path.name, flush=True)
    genes = pd.Index(pd.read_csv(SOURCE / 'gene_names.csv', header=None)[0].astype(str), name='gene_symbol')
    assert genes.is_unique and len(genes) == 18080
    metadata = {}
    var_tables = {}
    source_paths = {}
    for split, (directory, suffix) in FILES.items():
        source_paths[split] = SOURCE / directory / ('adata_' + suffix + '.h5ad')
        with h5py.File(source_paths[split]) as handle:
            metadata[split] = ad.io.read_elem(handle['obs'])
            var_tables[split] = ad.io.read_elem(handle['var'])
            assert var_tables[split].index.equals(genes)
    shared = set(metadata['train'].index) & set(metadata['validation'].index) & set(metadata['test'].index)
    assert len(shared) == 38176
    for split, obs in metadata.items():
        assert obs.index.is_unique and set(obs.columns) == {'target_gene', 'guide_id', 'batch'}
        assert set(obs.loc[obs.target_gene.eq('non-targeting')].index) == shared
    features = pd.DataFrame({'source_column': np.arange(len(genes)), 'gene_symbol': genes,
        'gene_id_from_train_native_var': var_tables['train'].gene_id.astype(str).to_numpy(),
        'gene_id_source': 'author_train_var_same_exact_gene_axis',
        'mt_symbol_within_Flex_panel': genes.str.startswith('MT-')}, index=genes)
    assert features.gene_id_from_train_native_var.ne('').all()
    summaries = []
    all_indices = []
    manifest = []
    target_tables = []
    shared_hashes = {}
    for split, obs in metadata.items():
        path = source_paths[split]
        directory, suffix = FILES[split]
        n_cells, n_genes = len(obs), len(genes)
        totals = np.zeros(n_cells, dtype=np.int64)
        detected = np.zeros(n_cells, dtype=np.int32)
        mt = np.zeros(n_cells, dtype=np.int64)
        gene_totals = np.zeros(n_genes, dtype=np.int64)
        gene_cells = np.zeros(n_genes, dtype=np.int64)
        control_hash = np.full(n_cells, '', dtype=object)
        noncanonical_blocks = 0
        stored_zeros = 0
        observed_max = 0
        with h5py.File(path) as handle:
            x = handle['X']
            assert x.attrs['encoding-type'] == 'csr_matrix' and tuple(x.attrs['shape']) == (n_cells, n_genes)
            pointer = x['indptr'][:]
            assert len(pointer) == n_cells + 1 and pointer[0] == 0 and pointer[-1] == len(x['data']) == len(x['indices'])
            assert (np.diff(pointer) >= 0).all() and not list(handle['layers']) and 'raw' not in handle
            for start in range(0, n_cells, BLOCK):
                stop = min(start + BLOCK, n_cells)
                left, right = int(pointer[start]), int(pointer[stop])
                data = x['data'][left:right]
                indices = x['indices'][left:right]
                assert np.isfinite(data).all() and (data >= 0).all() and np.equal(data, np.rint(data)).all()
                assert not len(indices) or (indices.min() >= 0 and indices.max() < n_genes)
                stored_zeros += int((data == 0).sum())
                if len(data): observed_max = max(observed_max, int(data.max()))
                chunk = csr_matrix((data.astype(np.int64), indices, pointer[start:stop + 1] - left), shape=(stop - start, n_genes))
                if not chunk.has_canonical_format:
                    noncanonical_blocks += 1
                    chunk.sum_duplicates()
                    chunk.sort_indices()
                chunk.eliminate_zeros()
                totals[start:stop] = np.asarray(chunk.sum(axis=1)).ravel()
                detected[start:stop] = chunk.getnnz(axis=1)
                mt[start:stop] = np.asarray(chunk[:, features.mt_symbol_within_Flex_panel.to_numpy()].sum(axis=1)).ravel()
                gene_totals += np.bincount(chunk.indices, weights=chunk.data, minlength=n_genes).astype(np.int64)
                gene_cells += np.bincount(chunk.indices, minlength=n_genes)
                for row in np.flatnonzero(obs.index[start:stop].isin(shared)):
                    low, high = chunk.indptr[row:row + 2]
                    control_hash[start + row] = row_fingerprint(chunk.indices[low:high], chunk.data[low:high])
                if start % (BLOCK * 32) == 0:
                    print(split, 'full count scan', stop, '/', n_cells, flush=True)
        index = obs.astype(str).copy()
        index.index.name = 'source_cell_id'
        index['source_row'] = np.arange(n_cells)
        index['cell_uid'] = 'VCC2025:' + index.index.astype(str)
        index['official_split'] = split
        index['count_source_path'] = str(path)
        index['is_control'] = index.target_gene.eq('non-targeting')
        index['shared_control_reference'] = index.index.isin(shared)
        index['source_occurrences_across_official_splits'] = np.where(index.shared_control_reference, 3, 1)
        index['total_counts'] = totals
        index['n_genes_by_counts'] = detected
        index['mt_counts_within_published_panel'] = mt
        index['pct_MT_within_published_panel'] = np.divide(100 * mt, totals, out=np.zeros(n_cells), where=totals > 0)
        index['control_expression_SHA256'] = control_hash
        # Pandas may expose a read-only view; this mask is updated below.
        valid_meta = index[['target_gene', 'guide_id', 'batch']].ne('').all(axis=1).to_numpy(dtype=bool, copy=True)
        valid_meta &= ~obs[['target_gene', 'guide_id', 'batch']].isna().any(axis=1).to_numpy()
        index['qc_usable'] = (totals > 0) & (detected > 0) & valid_meta
        index['qc_reason'] = np.select([~valid_meta, totals == 0], ['missing_source_annotation', 'zero_information_cell'], default='retained_publisher_release')
        index['use_role'] = np.where(index.is_control, 'shared_non_targeting_reference',
                                     'training_perturbation' if split == 'train' else 'held_out_' + split + '_perturbation')
        index['training_eligible_under_official_split'] = index.qc_usable & (split == 'train')
        index['held_out_perturbation_evaluation'] = index.qc_usable & ~index.is_control & (split != 'train')
        index['cell_context'] = 'H1_hESC'
        index['perturbation_modality'] = 'dual_guide_CRISPRi'
        index['disease_label_assigned'] = False
        index['patient_donor_split_applicable'] = False
        index.to_csv(OUT / (split + '_cell_QC.tsv.gz'), sep='\t')
        index.loc[index.qc_usable].to_csv(OUT / (split + '_usable_index.tsv.gz'), sep='\t')
        index.loc[index.held_out_perturbation_evaluation].to_csv(OUT / (split + '_held_out_perturbation_index.tsv.gz'), sep='\t') if split != 'train' else None
        index.loc[~index.qc_usable].to_csv(OUT / (split + '_quarantine_index.tsv'), sep='\t')
        all_indices.append(index)
        shared_hashes[split] = index.loc[index.shared_control_reference, 'control_expression_SHA256'].sort_index()
        official = pd.read_csv(SOURCE / directory / ('pert_counts_' + suffix + '.csv')).set_index('target_gene')
        actual = index.groupby('target_gene', observed=True).agg(n_cells=('cell_uid', 'size'), median_umi_per_cell=('total_counts', 'median'))
        assert set(official.index) == set(actual.index) - {'non-targeting'}
        assert np.array_equal(official.n_cells, actual.loc[official.index, 'n_cells'])
        assert np.array_equal(official.median_umi_per_cell, actual.loc[official.index, 'median_umi_per_cell'])
        actual['official_split'] = split
        actual['is_control'] = actual.index == 'non-targeting'
        actual.reset_index().to_csv(OUT / (split + '_perturbation_summary.tsv'), sep='\t', index=False)
        target_tables.append(actual.reset_index())
        features[split + '_total_counts'] = gene_totals
        features[split + '_n_cells_detected'] = gene_cells
        features[split + '_all_zero'] = gene_cells == 0
        summaries.append({'split': split, 'source_cells': n_cells, 'qc_usable': int(index.qc_usable.sum()),
            'genes_preserved': n_genes, 'target_perturbations': len(official), 'control_cells': int(index.is_control.sum()),
            'unique_noncontrol_cells': int((~index.is_control).sum()), 'zero_information_cells': int((totals == 0).sum()),
            'zero_features_preserved': int((gene_cells == 0).sum()), 'stored_zeros': stored_zeros,
            'noncanonical_blocks': noncanonical_blocks, 'nnz': int(pointer[-1]), 'maximum_count': observed_max,
            'median_total_counts': float(np.median(totals)), 'median_detected_genes': float(np.median(detected))})
        manifest.append({'dataset': 'VCC2025', 'official_split': split, 'source_count_file': str(path),
            'usable_index': str(OUT / (split + '_usable_index.tsv.gz')), 'source_cells': n_cells,
            'usable_cells': int(index.qc_usable.sum()), 'genes': n_genes, 'matrix_handling': 'read_validated_original_counts_no_rewrite'})
        print('Split completed:', split, json.dumps(summaries[-1]), flush=True)
    assert shared_hashes['train'].equals(shared_hashes['validation']) and shared_hashes['train'].equals(shared_hashes['test'])
    combined = pd.concat(all_indices)
    noncontrol = combined.loc[~combined.is_control]
    assert noncontrol.index.is_unique
    assert (combined.groupby('guide_id', observed=True).target_gene.nunique() == 1).all()
    controls = all_indices[0].loc[all_indices[0].is_control].copy()
    controls['reference_count_source_choice'] = 'train_copy_same_expression_verified_in_all_three_releases'
    controls.to_csv(OUT / 'shared_controls_unique_index.tsv.gz', sep='\t')
    features.to_csv(OUT / 'feature_index_QC.tsv.gz', sep='\t')
    pd.DataFrame(summaries).to_csv(OUT / 'split_QC_summary.tsv', sep='\t', index=False)
    pd.DataFrame(manifest).to_csv(OUT / 'DATASET_INDEX.tsv', sep='\t', index=False)
    pd.concat(target_tables).to_csv(OUT / 'PERTURBATION_SPLITS.tsv', sep='\t', index=False)
    pd.DataFrame(input_records).to_csv(OUT / 'ORIGINAL_FILES.tsv', sep='\t', index=False)
    policy = {'dataset': 'VCC2025', 'version': 'official_release_QC_v1', 'context': 'H1 human embryonic stem cell, dual-guide CRISPRi, 10x Flex',
        'count_semantics': 'Official AnnData X count matrices, float32 storage but every stored value checked finite, nonnegative and integer-valued. No layers/raw; matrices retained read-only at original source paths.',
        'source_gene_space': 'Official gene_names.csv18080 symbols, exact same order in all three releases. Native train var gene_id supplies traceable mapping only; source val/test QC annotations preserved, not overwritten.',
        'primary_sources': ['https://github.com/ArcInstitute/arc-virtual-cell-atlas/blob/main/virtual-cell-challenge/README.md', 'https://arcinstitute.org/news/behind-the-data-virtual-cell-challenge'],
        'official_splits': '150 train/50 validation/100 test target genes. Perturbed cells and target sets disjoint. Preserve release split; validation/test labels used only for integrity checks, never threshold selection/training.',
        'shared_controls': 'Same38176 non-targeting cells repeated in all3 official files; identity, source metadata and canonical count-rowSHA256 must all agree. Count once for unique-cell statistics; shared reference is not independent held-out outcomes. Single canonical read source is train copy.',
        'cell_QC': 'Structural validity, finite/nonnegative integer counts, nonzero information and complete target_gene/guide_id/batch. Index/quarantine only; no disease-dependent thresholds, no deletion for weak knockdown or large genuine response.',
        'feature_QC': 'All18080 release features retained, including any zero features; duplicate/missing identifiers audited. No HVG truncation or gene aggregation.',
        'gene_panel_MT': 'MT metrics apply only to published Flex probe panel and are descriptive; no MT threshold.',
        'modality_boundary': 'H1 hESC CRISPRi perturbation benchmark, not HCM, not cardiomyocytes; non-targeting is experimental negative control, not healthy patient. No disease0/1 labels or invented donors.',
        'not_applied': ['HCM500to5000_gene_filter', 'MT_lt5pct', 'EmptyDrops', 'doublet_caller', 'CM_identity_filter', 'renormalization', 'imputation', 'batch_integration', 'DE_analysis', 'tokenization', 'training'],
        'publisher_QC_limit': 'Retain curated publisher release; original FASTQ/droplet screening and per-cell knockdown QC are not re-established. Clean guide annotation does not independently prove single-cell singlets or biological knockdown efficacy.',
        'delivery': 'Validated original count files plus explicit row/feature/role indices; no34GB redundant copies, no synthetic or placeholder matrices. Per-cell QC, unique controls and official held-out perturbation indices remain on server.'}
    (OUT / 'qc_policy.json').write_text(json.dumps(policy, indent=2) + '\n')
    for name in ['ORIGINAL_DOWNLOAD_BASELINE.json', 'INPUT_INSPECTION.json']:
        shutil.copy2(WORK / name, OUT / name)
    payload = sorted(path for path in OUT.iterdir() if path.is_file())
    result = {'dataset': 'VCC2025', 'status': 'exported_pending_independent_verification',
        'completed_utc': datetime.now(timezone.utc).isoformat(), 'split_summaries': summaries,
        'file_observations': len(combined), 'unique_cells': combined.cell_uid.nunique(), 'shared_controls': len(controls),
        'unique_perturbed_cells': len(noncontrol), 'genes': len(genes),
        'original_files': input_records,
        'hashes': {path.name: {'bytes': path.stat().st_size, 'SHA256': checksum(path), 'MD5': checksum(path, 'md5')} for path in payload}}
    (OUT / 'EXPORT_SUCCESS.json').write_text(json.dumps(result, indent=2) + '\n')
    print(json.dumps({key: value for key, value in result.items() if key not in ['hashes', 'original_files']}, indent=2), flush=True)


if __name__ == '__main__':
    main()
