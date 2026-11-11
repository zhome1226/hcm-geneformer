"""Prepare bulk expression and traceable metadata without altering published values."""
from collections import Counter
from datetime import datetime, timezone
from importlib.metadata import version
from pathlib import Path
import gzip
import hashlib
import json
import re
import shutil

import anndata as ad
import numpy as np
import pandas as pd
from scipy.sparse import csr_matrix

ROOT = Path.home() / 'xiaoe-internship-data'
SOURCE = ROOT / 'validation/bulk_GSE249925'
WORK = ROOT / '_processing/all_data_cleaning_20261111/GSE249925'
OUT = ROOT / 'processed/bulk_GSE249925_qc_20261111'


def digest(path, algorithm='sha256'):
    with path.open('rb') as handle:
        return hashlib.file_digest(handle, algorithm).hexdigest()


def soft_samples(path):
    records = []
    current = None
    for line in path.read_text().splitlines():
        if line.startswith('^SAMPLE = '):
            current = {'accession': line.split(' = ', 1)[1]}
            records.append(current)
        elif current is not None and line.startswith('!Sample_') and ' = ' in line:
            field, value = line[1:].split(' = ', 1)
            current.setdefault(field, []).append(value)
    return records


def mad_flags(values):
    logged = np.log1p(values)
    median = np.median(logged)
    deviation = np.median(np.abs(logged - median))
    return np.abs(logged - median) > 5 * 1.4826 * deviation if deviation > 0 else np.zeros(len(values), dtype=bool)


def main():
    OUT.mkdir(exist_ok=True, parents=True)
    assert not (OUT / 'EXPORT_SUCCESS.json').exists(), 'Inspect successful export before rerunning'
    baseline = json.loads((WORK / 'ORIGINAL_DOWNLOAD_BASELINE.json').read_text())
    assert len(baseline) == 3
    for row in baseline:
        path = Path(row['server_path'])
        assert path.stat().st_size == int(row['bytes']) and digest(path, 'md5') == row['md5']
    path = SOURCE / 'GSE249925_20231208-For-GEO-Raw-Counts.tsv.gz'
    table = pd.read_csv(path, sep='\t', dtype=str, keep_default_na=False)
    assert table.columns[:2].tolist() == ['Gene_Index', 'Length'] and table.Gene_Index.is_unique
    values = table.iloc[:, 2:].apply(pd.to_numeric, errors='raise').to_numpy()
    assert np.isfinite(values).all() and (values >= 0).all() and np.equal(values, np.rint(values)).all()
    values = values.astype(np.int64)
    native = pd.Index(table.Gene_Index, name='source_feature_id')
    components = table.Gene_Index.str.split('@', expand=True)
    assert components.shape[1] == 5
    var = pd.DataFrame({'source_row': np.arange(len(table)), 'source_Gene_Index': native,
                        'source_Length': table.Length.to_numpy(), 'chromosome': components[0].to_numpy(),
                        'gene_id': components[1].to_numpy(), 'gene_symbol': components[2].to_numpy(),
                        'start_original': components[3].to_numpy(), 'end_original': components[4].to_numpy()}, index=native)
    var['length_bases'] = pd.to_numeric(var.source_Length, errors='coerce')
    var['start'] = pd.to_numeric(var.start_original, errors='coerce')
    var['end'] = pd.to_numeric(var.end_original, errors='coerce')
    var['gene_id_is_Ensembl_format'] = var.gene_id.str.fullmatch(r'ENSG\d{11}(?:\.\d+)?')
    var['gene_id_nonunique'] = var.gene_id.duplicated(keep=False)
    var['gene_symbol_nonunique'] = var.gene_symbol.duplicated(keep=False)
    var['identity_valid'] = var.gene_id_is_Ensembl_format & var.chromosome.ne('') & var.gene_symbol.ne('')
    var['identity_valid'] &= var.start.ge(1) & var.end.ge(var.start) & var.length_bases.gt(0)
    var['published_expression_total'] = values.sum(axis=1)
    var['n_samples_detected'] = (values > 0).sum(axis=1)
    var['all_zero'] = var.n_samples_detected.eq(0)
    var['mt_by_source_chromosome'] = var.chromosome.isin(['chrM', 'chrMT', 'MT', 'M'])
    var['qc_usable'] = var.identity_valid & ~var.all_zero
    var['qc_reason'] = np.select([~var.identity_valid, var.all_zero], ['invalid_feature_identity', 'all_samples_zero'], default='retained')
    records = soft_samples(SOURCE / 'sample_metadata.soft')
    by_title = {record['Sample_title'][0]: record for record in records}
    sample_ids = pd.Index(table.columns[2:], name='sample_id')
    assert len(by_title) == len(records) == len(sample_ids) == 120 and set(sample_ids) == set(by_title)
    rows = []
    all_processing = set()
    for sample in sample_ids:
        record = by_title[sample]
        characteristics = dict(text.split(': ', 1) for text in record['Sample_characteristics_ch1'])
        relations = record.get('Sample_relation', [])
        biosamples = [text.rsplit('/', 1)[1] for text in relations if text.startswith('BioSample: ')]
        sra = [text.split('term=', 1)[1] for text in relations if text.startswith('SRA: ') and 'term=' in text]
        assert len(biosamples) == len(sra) == 1 and characteristics['group'] in ['HCM', 'Control']
        assert record['Sample_organism_ch1'] == ['Homo sapiens'] and record['Sample_library_strategy'] == ['RNA-Seq']
        all_processing.update(record.get('Sample_data_processing', []))
        rows.append({'sample_id': sample, 'GSM': record['accession'], 'source_title': sample,
                     'BioSample': biosamples[0], 'SRA_experiment': sra[0],
                     'disease': characteristics['group'], 'label': {'Control': 0, 'HCM': 1}[characteristics['group']],
                     'tissue_original': characteristics['tissue'], 'source_name': '|'.join(record['Sample_source_name_ch1']),
                     'organism': 'Homo sapiens', 'modality': 'bulk_RNA',
                     'instrument': '|'.join(record.get('Sample_instrument_model', [])),
                     'donor_id': '', 'donor_mapping_status': 'not_provided_per_sample',
                     'age_available': False, 'sex_available': False, 'batch_available': False,
                     'split_status': 'unassigned_subject_and_cross_study_overlap_review_required',
                     'training_eligible_single_cell': False, 'certified_unnormalized_counts': False,
                     'GEO_url': 'https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=' + record['accession']})
    obs = pd.DataFrame(rows).set_index('sample_id')
    assert obs.GSM.is_unique and dict(Counter(obs.disease)) == {'HCM': 97, 'Control': 23}
    obs['published_expression_total'] = values.sum(axis=0)
    obs['n_features_detected'] = (values > 0).sum(axis=0)
    obs['pct_published_expression_mt'] = 100 * values[var.mt_by_source_chromosome.to_numpy()].sum(axis=0) / np.maximum(values.sum(axis=0), 1)
    obs['log_total_5MAD_flag_only'] = mad_flags(obs.published_expression_total.to_numpy())
    obs['log_detected_5MAD_flag_only'] = mad_flags(obs.n_features_detected.to_numpy())
    hashes = [hashlib.sha256(np.ascontiguousarray(values[:, i]).tobytes()).hexdigest() for i in range(len(obs))]
    frequencies = Counter(hashes)
    obs['expression_vector_sha256'] = hashes
    obs['exact_expression_duplicate'] = [frequencies[value] > 1 for value in hashes]
    obs['biosample_nonunique'] = obs.BioSample.duplicated(keep=False)
    obs['sra_experiment_nonunique'] = obs.SRA_experiment.duplicated(keep=False)
    obs['qc_usable'] = obs.published_expression_total.gt(0) & obs.n_features_detected.gt(0)
    obs['qc_reason'] = np.where(obs.qc_usable, 'retained', 'zero_information_sample')
    assert any('Timmed means of M-values' in text for text in all_processing)
    assert any('Raw gene-count matrix' in text for text in all_processing)
    policy = {'dataset': 'GSE249925', 'version': 'bulk_QC_v1',
              'matrix_semantics': 'Published nonnegative integer expression; GEO raw-count label conflicts with statement that provided reads were TMM-normalized. Integer values do not resolve provenance. Not certified unnormalized counts.',
              'certified_unnormalized_counts': False,
              'published_processing': sorted(all_processing),
              'author_reference': 'GEO: hg38, Ensembl v78, STAR, HT-Seq, edgeR/TMM; PMID39255898, DOI10.1016/j.yjmcc.2024.09.005',
              'original_preservation': 'Original three files read-only; all native composite feature identifiers, lengths, sample order and values retained in all_published_expression.h5ad. X and published_values layer are identical integers, never called uncorrected raw counts.',
              'QC_features': 'Keep identifiable native features with any nonzero published value. All-zero and structurally invalid features retained in full source and separately indexed in feature_exclusions.tsv.gz; no low-abundance, disease-dependent or MT threshold filtering.',
              'feature_axis': 'Split native chromosome@Ensembl@symbol@start@end without renaming or aggregation; retain duplicate native Ensembl IDs and symbols as distinct source rows, mark nonunique identifiers. Length is source annotation, not inferred genomic span.',
              'QC_samples': 'Exact title-to-GSM mapping, nonnegative finite integer values, known source group and nonzero information. Expression/registry duplicates and 5MAD log totals/detected-feature flags are review flags only, never guessed donor identities or automatic deletions.',
              'disease': {'HCM': 1, 'Control': 0},
              'metadata': 'Preserve exact sample titles including numeric strings and R suffixes. GSM/BioSample/SRA identify source records, not proven unique patients. No individual donor, sex, age or batch map provided; no imputation from names or study sex totals.',
              'study': 'Garmany et al., J Mol Cell Cardiol 2024: reported97 symptomatic obstructive HCM myectomy patients and23 donor cardiac controls. Per-sample group/tissue from GEO; no fabricated specimen-level septal location or sex.',
              'training': 'Bulk tissue mixtures are not individual cardiomyocytes. No single-cell training eligibility or automatic substitution for snRNA validation; donor/sex/batch and raw-count provenance require review before count-model inference or splitting.',
              'not_applied': ['single-cell500to5000', 'MT_lt5pct', 'EmptyDrops', 'doublet_calling', 'CM_extraction', 'normalization', 'imputation', 'batch_correction', 'differential_expression', 'tokenization', 'training']}
    x = csr_matrix(values.T)
    full = ad.AnnData(X=x, obs=obs, var=var)
    full.layers['published_values'] = x.copy()
    full.uns['qc_policy_json'] = json.dumps(policy)
    full.uns['X_semantics'] = policy['matrix_semantics']
    full.uns['source_input_path'] = str(path)
    full.uns['count_certification'] = 'unresolved_raw_vs_TMM_source_statement_conflict'
    full.write_h5ad(OUT / 'all_published_expression.h5ad', compression='gzip')
    core = full[obs.qc_usable.to_numpy(), var.qc_usable.to_numpy()].copy()
    core.write_h5ad(OUT / 'bulk_qc_expression.h5ad', compression='gzip')
    selected = table.loc[var.qc_usable.to_numpy(), ['Gene_Index', 'Length'] + obs.index[obs.qc_usable].tolist()]
    selected.to_csv(OUT / 'bulk_qc_matrix.tsv.gz', sep='\t', index=False, compression='gzip')
    obs.to_csv(OUT / 'sample_metadata_QC.tsv', sep='\t')
    var.to_csv(OUT / 'feature_audit.tsv.gz', sep='\t')
    var.loc[~var.qc_usable].to_csv(OUT / 'feature_exclusions.tsv.gz', sep='\t')
    obs.loc[~obs.qc_usable].to_csv(OUT / 'sample_exclusions.tsv', sep='\t')
    obs.groupby('disease', observed=True).agg(samples=('GSM', 'size'), usable=('qc_usable', 'sum'),
        total_expression_min=('published_expression_total', 'min'), total_expression_max=('published_expression_total', 'max'),
        detected_min=('n_features_detected', 'min'), detected_max=('n_features_detected', 'max'),
        MAD_total_flags=('log_total_5MAD_flag_only', 'sum'), MAD_detected_flags=('log_detected_5MAD_flag_only', 'sum'),
        exact_expression_duplicates=('exact_expression_duplicate', 'sum')).to_csv(OUT / 'group_QC_summary.tsv', sep='\t')
    (OUT / 'qc_policy.json').write_text(json.dumps(policy, indent=2) + '\n')
    for name in ['ORIGINAL_DOWNLOAD_BASELINE.json', 'INPUT_INSPECTION.json']:
        shutil.copy2(WORK / name, OUT / name)
    outputs = ['all_published_expression.h5ad', 'bulk_qc_expression.h5ad', 'bulk_qc_matrix.tsv.gz',
               'sample_metadata_QC.tsv', 'feature_audit.tsv.gz', 'feature_exclusions.tsv.gz', 'sample_exclusions.tsv',
               'group_QC_summary.tsv', 'qc_policy.json']
    result = {'dataset': 'GSE249925', 'status': 'exported_pending_independent_verification',
              'completed_utc': datetime.now(timezone.utc).isoformat(),
              'source_features': full.n_vars, 'QC_features': core.n_vars, 'source_samples': full.n_obs, 'QC_samples': core.n_obs,
              'source_groups': dict(Counter(obs.disease)), 'all_zero_features': int(var.all_zero.sum()),
              'invalid_feature_identity': int((~var.identity_valid).sum()),
              'duplicate_Ensembl_rows': int(var.gene_id_nonunique.sum()), 'duplicate_symbol_rows': int(var.gene_symbol_nonunique.sum()),
              'registry_duplicates': int(obs.biosample_nonunique.sum() + obs.sra_experiment_nonunique.sum()),
              'expression_duplicate_samples': int(obs.exact_expression_duplicate.sum()),
              'certified_unnormalized_counts': False,
              'hashes': {name: {'bytes': (OUT / name).stat().st_size, 'MD5': digest(OUT / name, 'md5'),
                                 'SHA256': digest(OUT / name)} for name in outputs},
              'python_versions': {package: version(package) for package in ['anndata', 'numpy', 'pandas', 'scipy']}}
    (OUT / 'EXPORT_SUCCESS.json').write_text(json.dumps(result, indent=2) + '\n')
    print(json.dumps({key: value for key, value in result.items() if key != 'hashes'}, indent=2), flush=True)


if __name__ == '__main__':
    main()
