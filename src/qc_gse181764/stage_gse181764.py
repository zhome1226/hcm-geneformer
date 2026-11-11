"""Stage GSE181764 without changing counts; preserve unequal feature spaces."""
from concurrent.futures import ThreadPoolExecutor
from datetime import datetime, timezone
from pathlib import Path
import gc
import gzip
import hashlib
import json
import re

import anndata as ad
import h5py
import numpy as np
import pandas as pd
from scipy.io import mmread
from threadpoolctl import threadpool_limits

ROOT = Path.home() / 'xiaoe-internship-data'
WORK = ROOT / '_processing/all_data_cleaning_20261111/gse181764'
WORK.mkdir(parents=True, exist_ok=True)
COHORTS = {'GSE181764': ('mixed', -1, 9)}


def digest(path):
    with path.open('rb') as handle:
        return hashlib.file_digest(handle, 'sha256').hexdigest()


def parse_samples(path):
    records = []
    current = None
    for line in path.read_text().splitlines():
        if line.startswith('^SAMPLE = '):
            current = {'accession': line.split(' = ', 1)[1], 'source_path': str(path)}
            records.append(current)
        elif current is not None and line.startswith('!Sample_') and ' = ' in line:
            key, value = line.removeprefix('!Sample_').split(' = ', 1)
            current.setdefault(key, []).append(value)
    return records


metadata = parse_samples(ROOT / 'validation/chin_GSE181764/sample_metadata.soft')
other_records = []
prior_samples = pd.read_csv(ROOT / '_processing/all_data_cleaning_20261111/chin_pair/samples.tsv', sep='\t', dtype=str)
for row in prior_samples.itertuples():
    other_records.append({'accession': row.sample_id, 'title': [row.donor_original]})
samples = []
for cohort, (disease, label, expected) in COHORTS.items():
    source_dir = ROOT / f'validation/chin_{cohort}/raw'
    for matrix in sorted(source_dir.glob('*_matrix.mtx.gz')):
        prefix = matrix.name.removesuffix('_matrix.mtx.gz')
        accession = prefix.split('_')[0]
        records = [r for r in metadata if r['accession'] == accession]
        assert len(records) == 1, accession
        record = records[0]
        characteristics = dict(x.split(': ', 1) for x in record['characteristics_ch1'])
        donor = record['title'][0].split(':')[0].strip()
        assert re.fullmatch(r'\d+', donor)
        original_disease = characteristics['disease status']
        disease = 'Control' if original_disease == 'Normal' else 'HCM'
        assert original_disease in ['Normal', 'Obstructive HCM', 'Non-Obstructive HCM']
        label = 0 if disease == 'Control' else 1
        overlap = [r['accession'] for r in other_records if r['title'][0].split(':')[0].strip() == donor]
        space_features = pd.read_csv(str(source_dir / prefix) + '_features.tsv.gz', sep='\t', header=None, dtype=str, keep_default_na=False)
        feature_space = hashlib.sha256(space_features.to_csv(sep='\t', header=False, index=False).encode()).hexdigest()[:12]
        samples.append({'dataset': f'chin_{cohort}', 'series': cohort, 'sample_id': accession,
                        'capture_id': accession, 'donor_id': 'Chin_' + donor,
                        'donor_original': donor, 'disease': disease, 'label': label,
                        'disease_subtype': original_disease,
                        'clinical_stage': 'end_stage_heart_failure_explant' if original_disease == 'Non-Obstructive HCM' else ('unused_donor_heart' if disease == 'Control' else 'obstructive_myectomy'),
                        'clinical_stage_source': 'https://pmc.ncbi.nlm.nih.gov/articles/PMC9555339/',
                        'feature_space_id': feature_space,
                        'source_gene_count': len(space_features),
                        'platform_id': record.get('platform_id', ['unspecified'])[0],
                        'sequencing_run_id': 'GSE181764_run16' if record['description'][0] == '16' else 'GSE181764_run14_15_17_22',
                        'known_other_Chin_accessions_same_donor_id': ','.join(overlap),
                        'tissue': 'interventricular septum', 'sex': characteristics['gender'],
                        'age_years': int(characteristics.get('age', characteristics.get('age (years)'))),
                        'technology': '10x 3prime snRNA-seq', 'source_prefix': str(source_dir / prefix),
                        'metadata_path': record['source_path'],
                        'disease_source_url': f'https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc={cohort}',
                        'capture_evidence': 'GEO sample-specific 10x matrix; shared sequencing run is not a shared capture',
                        'donor_independence': 'Bind same Chin donor identifier across series; absence of matching identifier is not proof of biological independence'})
    assert sum(s['series'] == cohort for s in samples) == expected
assert len({s['donor_id'] for s in samples}) == len(samples) == 9
assert sum(s['disease'] == 'Control' for s in samples) == 2
assert sum(s['disease_subtype'] == 'Obstructive HCM' for s in samples) == 1
assert sum(s['disease_subtype'] == 'Non-Obstructive HCM' for s in samples) == 6
pd.DataFrame(samples).to_csv(WORK / 'samples.tsv', sep='\t', index=False)
(WORK / 'SAMPLE_PROVENANCE.json').write_text(json.dumps(samples, indent=2) + '\n')


def stage(sample):
    directory = WORK / sample['sample_id']
    directory.mkdir(exist_ok=True)
    paths = {key: Path(sample['source_prefix'] + suffix) for key, suffix in
             [('matrix', '_matrix.mtx.gz'), ('features', '_features.tsv.gz'), ('barcodes', '_barcodes.tsv.gz')]}
    sources = {key: {'path': str(path), 'bytes': path.stat().st_size, 'sha256': digest(path)}
               for key, path in paths.items()}
    success = directory / 'STAGED.json'
    if success.exists():
        old = json.loads(success.read_text())
        if old['sources'] == sources and all((directory / p).exists() for p in
                                            ['raw_10x.h5', 'prefilter.h5ad', 'droplet_qc.tsv.gz']):
            print('Stage already complete: ' + sample['sample_id'], flush=True)
            return old
        raise RuntimeError('Existing stage does not match source: ' + sample['sample_id'])
    features = pd.read_csv(paths['features'], sep='\t', header=None, dtype=str, keep_default_na=False)
    barcodes = pd.read_csv(paths['barcodes'], sep='\t', header=None, dtype=str)[0]
    assert features[0].is_unique and barcodes.is_unique
    assert features.shape[1] == 3 and features[2].eq('Gene Expression').all()
    with gzip.open(paths['matrix'], 'rb') as handle, threadpool_limits(limits=2):
        original = mmread(handle)
        shape = original.shape
        original_nnz = original.nnz
        assert shape == (len(features), len(barcodes))
        assert np.isfinite(original.data).all() and (original.data >= 0).all()
        assert (original.data == np.round(original.data)).all()
        counts = original.T.tocsr()
    counts.sum_duplicates()
    counts.eliminate_zeros()
    assert counts.data.max(initial=0) <= np.iinfo(np.int32).max
    counts = counts.astype(np.int32)
    del original
    total = np.asarray(counts.sum(axis=1)).ravel()
    detected = np.diff(counts.indptr)
    mito = features[1].str.startswith('MT-').to_numpy()
    mito_counts = np.asarray(counts[:, mito].sum(axis=1)).ravel()
    assert mito.sum() >= 10
    pct_mito = np.divide(100 * mito_counts, total, out=np.zeros(len(total), dtype=float), where=total > 0)
    audit = pd.DataFrame({'barcode_original': barcodes.to_numpy(), 'source_row': np.arange(len(barcodes)),
                          'total_counts': total, 'n_genes_by_counts': detected,
                          'mt_counts': mito_counts, 'pct_counts_mt': pct_mito})
    audit.index = pd.Index(sample['sample_id'] + ':' + barcodes, name='barcode')
    audit['pass_strict_numeric'] = (detected >= 500) & (detected <= 5000) & (pct_mito < 5)
    audit.to_csv(directory / 'droplet_qc.tsv.gz', sep='\t', compression='gzip')
    # A CSR cell x gene matrix has exactly the CSC layout expected by 10x gene x barcode files.
    with h5py.File(directory / 'raw_10x.h5', 'w') as handle:
        group = handle.create_group('matrix')
        for name, values in [('data', counts.data), ('indices', counts.indices),
                             ('indptr', counts.indptr), ('shape', np.array(counts.shape[::-1]))]:
            group.create_dataset(name, data=values, compression='gzip')
        group.create_dataset('barcodes', data=np.asarray(barcodes, dtype='S'))
        f = group.create_group('features')
        for name, values in [('id', features[0]), ('name', features[1]), ('feature_type', features[2]),
                             ('genome', pd.Series(['GRCh38'] * len(features)))]:
            f.create_dataset(name, data=np.asarray(values, dtype='S'))
        f.create_dataset('_all_tag_keys', data=np.asarray(['genome'], dtype='S'))
    keep = detected >= 200
    obs = audit.loc[keep].copy()
    for key, value in sample.items():
        if key != 'source_prefix':
            obs[key] = value
    var = pd.DataFrame({'gene_symbol': features[1].to_numpy(), 'feature_type': features[2].to_numpy(), 'mt': mito},
                       index=pd.Index(features[0], name='gene_id'))
    query = ad.AnnData(counts[keep].copy(), obs=obs, var=var)
    query.uns['count_source'] = 'Unmodified GEO Cell Ranger 3.0 raw droplet gene-expression counts'
    query.write_h5ad(directory / 'prefilter.h5ad', compression='gzip')
    summary = {'sample_id': sample['sample_id'], 'dataset': sample['dataset'], 'sources': sources,
               'raw_shape_genes_droplets': list(shape), 'input_nnz': int(original_nnz),
               'canonical_nnz': int(counts.nnz), 'prefilter_nuclei_ge200_genes': int(keep.sum()),
               'strict_numeric_candidates_before_cell_calling': int(audit['pass_strict_numeric'].sum()),
               'mt_gene_count': int(mito.sum()), 'count_sum': int(total.sum()),
               'feature_space_id': sample['feature_space_id'],
               'ambient_positive_droplets_UMI_le100': int(((total > 0) & (total <= 100)).sum()),
               'total_UMI_quantiles': np.quantile(total, [0,.1,.5,.9,.99,1]).tolist(),
               'detected_gene_quantiles': np.quantile(detected, [0,.1,.5,.9,.99,1]).tolist(),
               'completed_utc': datetime.now(timezone.utc).isoformat()}
    success.write_text(json.dumps(summary, indent=2) + '\n')
    print(f"Staged {sample['sample_id']}: {int(keep.sum())} >=200-gene droplets; {summary['strict_numeric_candidates_before_cell_calling']} strict candidates", flush=True)
    del counts, query, audit
    gc.collect()
    return summary


with ThreadPoolExecutor(max_workers=3) as pool:
    summaries = list(pool.map(stage, samples))
(WORK / 'STAGING_SUMMARY.json').write_text(json.dumps(summaries, indent=2) + '\n')
feature_sets = {}
feature_frames = {}
for sample in samples:
    frame = pd.read_csv(sample['source_prefix'] + '_features.tsv.gz', sep='\t', header=None, dtype=str, keep_default_na=False)
    key = sample['feature_space_id']
    feature_sets[key] = set(frame[0])
    feature_frames[key] = frame
    directory = WORK / 'feature_axes'
    directory.mkdir(exist_ok=True)
    frame.to_csv(directory / (key + '.tsv'), sep='\t', header=False, index=False)
union = sorted(set.union(*feature_sets.values()))
common = sorted(set.intersection(*feature_sets.values()))
preferred = sorted(feature_frames, key=lambda key: (-len(feature_frames[key]), key))
symbols = {}
for key in reversed(preferred):
    symbols.update(dict(zip(feature_frames[key][0], feature_frames[key][1])))
index = pd.DataFrame({'gene_id': union, 'gene_symbol': [symbols[g] for g in union]})
for key in preferred:
    index['available_' + key] = index.gene_id.isin(feature_sets[key])
index.to_csv(WORK / 'UNION_FEATURES.tsv', sep='\t', index=False)
(WORK / 'COMMON_GENE_IDS.txt').write_text('\n'.join(common) + '\n')
originals = {p.name: {'bytes': p.stat().st_size, 'sha256': digest(p)} for p in (ROOT / 'validation/chin_GSE181764').iterdir() if p.is_file()}
result = {'status': 'staged_not_cleaned', 'samples': len(samples), 'donors': len({s['donor_id'] for s in samples}),
          'feature_spaces': {key: len(value) for key, value in feature_sets.items()}, 'union_genes': len(union), 'common_genes': len(common),
          'union_semantics': 'Sparse zero padding for unavailable features is structural missingness, not observed biological zero; feature availability mask is mandatory.',
          'originals': originals, 'known_donor_id_matches': {s['sample_id']: s['known_other_Chin_accessions_same_donor_id'] for s in samples if s['known_other_Chin_accessions_same_donor_id']}}
(WORK / 'PREFLIGHT.json').write_text(json.dumps(result, indent=2) + '\n')
print(json.dumps(result, indent=2), flush=True)
print('GSE181764_STAGING_COMPLETE', flush=True)
