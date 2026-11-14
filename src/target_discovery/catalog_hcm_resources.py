"""Catalog existing QC objects without changing counts, QC thresholds or splits."""
import argparse
import hashlib
import json
import os
import platform
from datetime import datetime, timezone
from pathlib import Path

os.environ.setdefault('OPENBLAS_NUM_THREADS', '2')
os.environ.setdefault('OMP_NUM_THREADS', '2')
import h5py
import numpy as np
import pandas as pd
from preflight_inputs import plain, read_elem, sample_matrix, unique_values


def write_tsv(rows, path):
    # Explicit missing values preserve the schema and avoid trailing whitespace.
    pd.DataFrame(rows).fillna('NA').replace('', 'NA').to_csv(path, sep='\t', index=False)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--config', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    config = json.loads(args.config.read_text(encoding='utf-8'))
    assert not args.output.exists(), 'Choose a new output directory; do not overwrite a catalog'
    args.output.mkdir(parents=True)
    audit = {'checked_utc': datetime.now(timezone.utc).isoformat(),
             'hostname': platform.node(), 'scope': 'metadata plus sampled count values; no full revalidation',
             'new_qc_run': False, 'matrices_concatenated': False, 'dge_run': False,
             'training_run': False, 'split': 'unassigned', 'model_exposure': 'unresolved',
             'config_sha256': hashlib.sha256(args.config.read_bytes()).hexdigest(), 'datasets': {}}
    donors, samples, catalog = [], [], []
    for key, record in config['datasets'].items():
        path = Path(record['path'])
        before = path.stat()
        with h5py.File(path, 'r') as f:
            obs, var = read_elem(f['obs']), read_elem(f['var'])
            matrices = {'X': sample_matrix(f['X'])}
            matrices.update({'layers/' + k: sample_matrix(v) for k, v in f.get('layers', {}).items()})
            provenance = {k: plain(read_elem(f['uns'][k])) for k in
                          ['qc_parameters', 'qc_policy', 'count_source', 'count_layer_provenance',
                           'source', 'review_policy', 'author_qc_reused', 'feature_availability_captures']
                          if k in f.get('uns', {})}
            varm_shapes = {k: list(v.shape) for k, v in f.get('varm', {}).items()
                           if isinstance(v, h5py.Dataset)}
            required = ['donor_id', 'sample_id', 'disease', 'label']
            assert all(c in obs for c in required), key + ': missing metadata fields'
            assert not obs[required].isna().any().any(), key + ': missing metadata values'
            labels = pd.to_numeric(obs['label'], errors='raise')
            assert labels.eq(obs['disease'].astype(str).map({'HCM': 1, 'NF': 0, 'Control': 0})).all()
            assert obs.index.is_unique and var.index.is_unique, key + ': nonunique axis'
            assert (obs.groupby('donor_id', observed=True)['disease'].nunique() == 1).all()
            assert (obs.groupby('sample_id', observed=True)['donor_id'].nunique() == 1).all()
            assert len(obs) == record['expected_cells'], key + ': count changed from prior report'
            assert before.st_size == record['bytes'], key + ': size changed from prior manifest'
            assert 'layers/counts' in matrices, key + ': no counts layer'
            for name in ['X', 'layers/counts']:
                m = matrices[name]
                assert m['shape'] == [len(obs), len(var)], key + ': matrix shape mismatch'
                assert m.get('sample_finite') and m.get('sample_nonnegative') and m.get('sample_integer_valued'), key + ': count sample failed'
            disease_cells = {str(k): int(v) for k, v in obs['disease'].value_counts().items()}
            disease_donors = {str(k): int(v) for k, v in obs.groupby('disease', observed=True)['donor_id'].nunique().items()}
            for donor, rows in obs.groupby('donor_id', observed=True):
                donors.append({'dataset': key, 'donor_id': str(donor), 'disease': unique_values(rows, 'disease')[0],
                               'cells_or_nuclei': len(rows), 'sample_ids': ';'.join(unique_values(rows, 'sample_id')),
                               'resource_group': record['group'], 'split': 'unassigned',
                               'cross_study_person_id': 'unresolved'})
            for sample, rows in obs.groupby('sample_id', observed=True):
                samples.append({'dataset': key, 'sample_id': str(sample),
                                'donor_id': unique_values(rows, 'donor_id')[0],
                                'disease': unique_values(rows, 'disease')[0], 'cells_or_nuclei': len(rows)})
            audit['datasets'][key] = {'path': str(path), 'bytes': before.st_size, 'shape': [len(obs), len(var)],
                                      'resource_group': record['group'], 'disease_cells': disease_cells,
                                      'donors_per_disease': disease_donors, 'sample_count': int(obs['sample_id'].nunique()),
                                      'matrices': matrices, 'provenance': provenance, 'varm_shapes': varm_shapes,
                                      'uns_keys': list(f.get('uns', {}).keys()),
                                      'obs_columns': obs.columns.tolist()}
        after = path.stat()
        assert before.st_size == after.st_size and before.st_mtime_ns == after.st_mtime_ns
        audit['datasets'][key]['source_size_and_mtime_unchanged'] = True
        catalog.append({'dataset': key, 'resource_group': record['group'], 'cells_or_nuclei': len(obs),
                        'HCM': disease_cells.get('HCM', 0),
                        'Control_or_NF': disease_cells.get('Control', 0) + disease_cells.get('NF', 0),
                        'HCM_donors': disease_donors.get('HCM', 0),
                        'Control_or_NF_donors': disease_donors.get('Control', 0) + disease_donors.get('NF', 0),
                        'genes': len(var), 'server_path': str(path), 'bytes': before.st_size,
                        'qc_status': record['qc_status'], 'count_semantics': record['count_semantics'],
                        'limitation': record['limitation'], 'new_v1_tokenization_needed': record['new_v1_tokenization_needed'],
                        'split': 'unassigned', 'model_exposure': 'unresolved'})
        print(json.dumps({'dataset': key, 'cells': len(obs), 'donors': disease_donors}), flush=True)

    # Read existing full-cell QC metrics only; no filtering or diagnostic expression transform.
    source = Path(config['full_cell_reqc_source'])
    before = source.stat()
    with h5py.File(source, 'r') as f:
        obs = read_elem(f['obs'])
    metrics = [c for c in obs if c in ['n_genes_by_counts', 'pct_counts_mt', 'total_counts', 'pct_mito', 'n_genes']]
    source_summary = {'path': str(source), 'source_cells': len(obs), 'obs_columns': obs.columns.tolist(),
                      'disease_cells': plain(obs['disease'].value_counts().to_dict()),
                      'donors_per_disease': plain(obs.groupby('disease', observed=True)['donor_id'].nunique().to_dict()),
                      'existing_metric_quantiles': {}}
    for disease, rows in obs.groupby('disease', observed=True):
        source_summary['existing_metric_quantiles'][str(disease)] = {
            c: plain(pd.to_numeric(rows[c], errors='raise').quantile([0, .1, .25, .5, .75, .9, 1]).to_dict())
            for c in metrics}
    after = source.stat()
    assert before.st_size == after.st_size and before.st_mtime_ns == after.st_mtime_ns
    source_summary['source_size_and_mtime_unchanged'] = True
    audit['GSE138262_full_source'] = source_summary

    adult = [r for r in catalog if r['resource_group'] == 'adult_non_scp_existing_qc']
    audit['adult_non_scp_existing_qc'] = {
        'datasets': [r['dataset'] for r in adult], 'nuclei': sum(r['cells_or_nuclei'] for r in adult),
        'HCM': sum(r['HCM'] for r in adult), 'Control': sum(r['Control_or_NF'] for r in adult),
        'dataset_scoped_donor_entries': sum(r['HCM_donors'] + r['Control_or_NF_donors'] for r in adult),
        'unique_people_across_studies': 'unresolved', 'independent_of_model_training': 'unresolved',
        'approved_as_pooled_dge_or_training_input': False}
    links = args.output/'adult_non_scp_existing_qc'
    links.mkdir()
    for r in adult:
        destination = links/(r['dataset'] + '.h5ad')
        destination.symlink_to(r['server_path'])
        assert destination.is_file() and destination.resolve() == Path(r['server_path']).resolve()
    (links/'README.txt').write_text(
        'Links to existing per-study QC objects; this is a candidate resource collection.\n'
        'No new QC, concatenation, donor split, training or independent-validation claim.\n'
        'Preserve each source file, feature availability masks and study/donor metadata.\n'
        'GSE181764 has five control nuclei across two donors; not a balanced cohort.\n', encoding='utf-8')
    write_tsv(catalog, args.output/'QC_RESOURCE_CATALOG.tsv')
    write_tsv(donors, args.output/'DONOR_CATALOG.tsv')
    write_tsv(samples, args.output/'SAMPLE_CATALOG.tsv')
    (args.output/'CATALOG_AUDIT.json').write_text(json.dumps(plain(audit), ensure_ascii=False, indent=2)+'\n', encoding='utf-8')
    print(json.dumps(audit['adult_non_scp_existing_qc']), flush=True)


if __name__ == '__main__':
    main()
