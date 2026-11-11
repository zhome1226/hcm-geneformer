"""Preserve and audit the real GSE227734 counts; do not invent gene annotations."""
from datetime import datetime, timezone
from pathlib import Path
import gzip
import hashlib
import json

import anndata as ad
import h5py
import numpy as np
import pandas as pd
from scipy.sparse import csr_matrix

ROOT = Path.home() / 'xiaoe-internship-data'
SOURCE = ROOT / 'validation/cox7b_GSE227734'
WORK = ROOT / '_processing/all_data_cleaning_20261111/GSE227734'
OUT = ROOT / 'processed/cox7b_GSE227734_qc_20261111_pending_annotations'
OUT.mkdir(parents=True, exist_ok=True)


def sha256(path):
    with path.open('rb') as stream:
        return hashlib.file_digest(stream, 'sha256').hexdigest()


def parse_soft(path):
    records = []
    record = None
    for line in path.read_text().splitlines():
        if line.startswith('^SAMPLE = '):
            record = {'accession': line.split(' = ', 1)[1]}
            records.append(record)
        elif record is not None and line.startswith('!Sample_') and ' = ' in line:
            key, value = line.removeprefix('!Sample_').split(' = ', 1)
            record.setdefault(key, []).append(value)
    return records


source_records = {}
for path in sorted(SOURCE.iterdir()):
    if path.is_file():
        source_records[path.name] = {'path': str(path), 'bytes': path.stat().st_size, 'sha256': sha256(path)}
compressed = SOURCE / 'GSE227734_hcm_raw.h5ad.gz'
native = WORK / 'GSE227734_hcm_raw.h5ad'
with gzip.open(compressed, 'rb') as stream:
    assert hashlib.file_digest(stream, 'sha256').hexdigest() == sha256(native), 'Decompression differs from original'

# Inspect physical HDF5, rather than treating this installed AnnData API's
# None-keyed alias of X as a separate source layer. There is no layers group.
with h5py.File(native) as handle:
    assert set(handle.keys()) == {'X', 'obs', 'var'}
    counts = np.asarray(handle['X'])
    original_obs = ad.io.read_elem(handle['obs'])
    original_var = ad.io.read_elem(handle['var'])
assert counts.shape == (338, 9760)
assert np.isfinite(counts).all() and (counts >= 0).all() and (counts == np.rint(counts)).all()
assert original_obs.index.is_unique and original_var.index.is_unique
assert np.array_equal(original_var.index.astype(str), np.arange(9760).astype(str))
assert original_var.shape[1] == 0
totals = counts.sum(axis=1, dtype=np.float64)
detected = (counts > 0).sum(axis=1)
assert np.array_equal(totals, original_obs.nCount_RNA)
assert np.array_equal(detected, original_obs.nFeature_RNA)
assert (totals > 0).all() and (counts.sum(axis=0) > 0).all()

records = parse_soft(SOURCE / 'sample_metadata.soft')
mapping = {}
for record in records:
    alias = record['description'][0]
    characteristics = dict(value.split(': ', 1) for value in record['characteristics_ch1'])
    assert alias not in mapping and characteristics['cell type'] == 'CM cells'
    mapping[alias] = {'source_sample_alias': alias, 'sample_id': record['accession'],
                      'donor_id': 'GSE227734:' + alias, 'source_title': record['title'][0],
                      'age_years': int(characteristics['age']), 'tissue': characteristics['tissue'],
                      'disease': 'HCM', 'label': 1,
                      'disease_evidence': 'GEO series overall_design and PMID37979444 abstract: primary CM from HCM myectomy patients',
                      'technology': 'STRT-seq isolated primary cardiomyocytes; GEO reports Drop-Seq processing',
                      'donor_independence': 'Five GEO source patients; cross-study donor overlap unestablished',
                      'source_url': 'https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=' + record['accession']}
assert len(mapping) == 5 and set(original_obs['sample'].astype(str)) == set(mapping)
series_text = (SOURCE / 'series_metadata.soft').read_text()
assert 'Primary cardiocytes were isolated from human heart of HCM patients' in series_text
obs = original_obs.copy()
obs['source_cell_id'] = original_obs.index.astype(str)
obs['source_row'] = np.arange(len(obs))
for key in next(iter(mapping.values())):
    obs[key] = [mapping[str(alias)][key] for alias in obs['sample']]
assert np.array_equal(original_obs.age, obs.age_years)
assert np.array_equal(original_obs['X'].astype(str), original_obs.index.astype(str))
obs.index = pd.Index('GSE227734:' + original_obs.index.astype(str), name='cell_id')
obs['dataset'] = 'cox7b_GSE227734'
obs['sex'] = original_obs.gender.astype(str).to_numpy()
obs['total_counts'] = totals
obs['n_anonymous_features_by_counts'] = detected
obs['pass_gene_count_bounds'] = (detected >= 500) & (detected <= 5000)
obs['pct_counts_mt'] = np.full(len(obs), np.nan)
obs['gene_annotation_status'] = 'missing_in_official_published_matrix'
obs['mt_qc_status'] = 'not_computable_without_gene_mapping'
obs['cell_identity_status'] = 'author_CM_only_not_independently_verified'
obs['doublet_status'] = 'not_assessed_no_gene_mapping_or_full_capture_population'
obs['cleaning_status'] = 'blocked_missing_gene_annotation'
obs['eligible_for_training'] = False
obs['author_area_units'] = 'not_confirmed_from_available_metadata'
obs['author_plate_or_library_assignment'] = 'unavailable_do_not_infer_from_repeated_8base_barcode'
var = pd.DataFrame({'source_feature_position': np.arange(counts.shape[1]),
                    'feature_identity_status': 'anonymous_unmapped_column_not_gene_ID'},
                   index=pd.Index(original_var.index.astype(str), name='anonymous_feature_id'))
obj = ad.AnnData(csr_matrix(counts), obs=obs, var=var)
obj.layers['counts'] = obj.X.copy()
policy = {
    'version': 'GSE227734_partial_audit_v1', 'dataset_cleaning_complete': False,
    'dataset_status': 'blocked_missing_gene_annotation',
    'project_numeric_rule': '500<=detected genes<=5000 and MT<5%; only anonymous detected-feature count can currently be checked',
    'gene_mapping': 'Official HDF5 has only X/obs/var; var has 9760 positional identifiers 0..9759 and no columns. No raw/uns/layers gene annotation is present. Do not map columns by guessing alphabetical or reference-genome order.',
    'MT': 'Uncomputable; NaN means missing measurement, never 0% or passing MT QC',
    'identity': 'CM is an author/source annotation; marker/reference verification unavailable without gene names',
    'doublets': 'No 10x loading-rate assumption, EmptyDrops, or gene-count-only doublet exclusion. Reused 8-base sample barcodes do not uniquely identify cells; preserve original patient-prefixed cell IDs.',
    'disease': 'Five HCM patients, no published healthy cell group; area_class Large/Small is morphology and is never recoded as disease or perturbation',
    'counts': 'All 338 rows/9760 anonymous features, exact published integer-valued RNA counts; no normalization, imputation, gene deletion or integration',
    'scope': 'Validated preparation/quarantine only, not a clean training/validation object. No core/light CM placeholder files. Full molecular QC remains incomplete.',
    'required_recovery': 'Obtain an authoritative ordered 9760-column gene mapping or corrected matching matrix and validate exact axes/count values before MT, identity and appropriate doublet review',
}
obj.uns['count_source'] = 'Unmodified published GSE227734 X; source contains no gene-name mapping'
obj.uns['qc_policy_json'] = json.dumps(policy)
obj.uns['author_obs_columns_json'] = json.dumps(list(original_obs.columns))
obj.uns['source_files_json'] = json.dumps(source_records)
obj.uns['source_decompressed_sha256'] = sha256(native)
obj.uns['dataset_cleaning_complete'] = False
obj.uns['feature_identity_warning'] = policy['gene_mapping']
obj.uns['no_train_test_split_assigned'] = True
obj.write_h5ad(OUT / 'author_counts_annotation_pending.h5ad', compression='gzip')
obs.to_csv(OUT / 'cell_preparation_audit.tsv.gz', sep='\t', compression='gzip')
var.to_csv(OUT / 'feature_annotation_audit.tsv.gz', sep='\t', compression='gzip')
rows = []
for alias, record in mapping.items():
    part = obs.loc[obs['sample'].astype(str).eq(alias)]
    rows.append({**record, 'published_cells': len(part), 'sex_from_published_obs': ','.join(sorted(part.sex.unique())),
                 'gene_bound_candidates_only': int(part.pass_gene_count_bounds.sum()),
                 'above_5000_features': int(part.n_anonymous_features_by_counts.gt(5000).sum()),
                 'MT_assessable_cells': 0, 'fully_QC_passed_cells': 0, 'training_eligible_cells': 0,
                 'published_Large': int(part.area_class.eq('Large').sum()),
                 'published_Small': int(part.area_class.eq('Small').sum()),
                 'status': 'blocked_missing_gene_annotation'})
pd.DataFrame(rows).to_csv(OUT / 'sample_preparation_summary.tsv', sep='\t', index=False)
(OUT / 'qc_policy.json').write_text(json.dumps(policy, indent=2) + '\n')
result = {'dataset': 'cox7b_GSE227734', 'status': 'prepared_pending_independent_audit',
          'cleaning_complete': False, 'prepared_utc': datetime.now(timezone.utc).isoformat(),
          'published_cells': len(obs), 'anonymous_features': len(var), 'donors': len(mapping),
          'gene_bound_candidates_only': int(obs.pass_gene_count_bounds.sum()),
          'MT_assessable_cells': 0, 'independently_identity_verified_cells': 0,
          'counts_file_bytes': (OUT / 'author_counts_annotation_pending.h5ad').stat().st_size,
          'counts_file_sha256': sha256(OUT / 'author_counts_annotation_pending.h5ad'),
          'source_files': source_records, 'decompressed_sha256': sha256(native)}
(OUT / 'PREPARATION.json').write_text(json.dumps(result, indent=2) + '\n')
print(json.dumps({k: v for k, v in result.items() if k != 'source_files'}, indent=2), flush=True)
