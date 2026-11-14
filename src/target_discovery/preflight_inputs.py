"""Read-only metadata/count-layer preflight; no DGE, normalization or data rewrite."""
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
try:
    from anndata.io import read_elem
except ImportError:  # Existing, pinned project environment.
    from anndata._io.specs import read_elem


def plain(value):
    if isinstance(value, dict):
        return {str(k): plain(v) for k, v in value.items()}
    if isinstance(value, (list, tuple, np.ndarray)):
        return [plain(v) for v in value]
    if isinstance(value, np.generic):
        return value.item()
    if isinstance(value, bytes):
        return value.decode('utf-8')
    return value


def sample_matrix(node):
    """Check spread-out stored values, explicitly not a full matrix audit."""
    if isinstance(node, h5py.Group):
        values = node['data']
        encoding = plain(node.attrs.get('encoding-type', 'unknown'))
        shape = plain(node.attrs['shape'])
    else:
        shape = list(node.shape)
        return {'encoding': 'dense', 'shape': shape, 'sample_status': 'not_sampled_dense'}
    n = len(values)
    width = min(n, 4096)
    starts = np.unique(np.linspace(0, max(0, n-width), 24, dtype=np.int64))
    samples = np.concatenate([values[int(s):int(s)+width] for s in starts]) if n else np.array([])
    finite = bool(np.isfinite(samples).all())
    return {'encoding': encoding, 'shape': shape, 'dtype': str(values.dtype),
            'stored_values': n, 'sampled_values': len(samples),
            'sampling': 'up to 24 evenly spread slices of at most 4096 stored entries; slices may overlap',
            'sample_finite': finite,
            'sample_nonnegative': bool((samples >= 0).all()) if finite else False,
            'sample_integer_valued': bool(np.equal(samples, np.round(samples)).all()) if finite else False,
            'sample_min': float(samples.min()) if len(samples) else None,
            'sample_max': float(samples.max()) if len(samples) else None,
            'sample_status': 'sampled_only_not_full_audit'}


def unique_values(frame, column):
    if column not in frame:
        return []
    return sorted(frame[column].dropna().astype(str).unique().tolist())


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--config', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    assert args.config.is_file()
    args.output.mkdir(parents=True, exist_ok=True)
    assert not (args.output/'PREFLIGHT.json').exists(), 'Completed result exists; choose a new output directory'
    config = json.loads(args.config.read_text())
    evidence = {'checked_utc': datetime.now(timezone.utc).isoformat(),
                'hostname': platform.node(), 'datasets': {},
                'scope': 'metadata and spread-out value samples; no hashes or full matrix comparison',
                'input_files_opened_read_only': True, 'dge_run': False, 'model_training_run': False}
    donor_rows, sample_rows, subtype_rows, alias_rows = [], [], [], []
    for dataset, record in config.items():
        path = Path(record['path'])
        before = path.stat()
        with h5py.File(path, 'r') as file:
            obs = read_elem(file['obs'])
            var = read_elem(file['var'])
            x = sample_matrix(file['X'])
            layers = {name: sample_matrix(node) for name, node in file.get('layers', {}).items()}
            provenance = {}
            for key in ['count_layer_provenance', 'qc_parameters', 'qc_policy', 'author_qc_reused',
                        'tokenizer_preparation', 'source', 'count_source', 'review_policy']:
                if key in file.get('uns', {}):
                    provenance[key] = plain(read_elem(file['uns'][key]))
            item = {'path': str(path), 'bytes': before.st_size,
                    'bytes_match_prior_manifest': before.st_size == record['bytes'],
                    'shape': [len(obs), len(var)], 'obs_index_unique': bool(obs.index.is_unique),
                    'var_index_unique': bool(var.index.is_unique),
                    'obs_columns': obs.columns.tolist(), 'var_columns': var.columns.tolist(),
                    'uns_keys': list(file.get('uns', {}).keys()), 'provenance': provenance,
                    'X': x, 'layers': layers,
                    'fields': {}, 'warnings': []}
        if 'donor_id' not in obs or 'sample_id' not in obs:
            raise ValueError(dataset + ' missing donor_id/sample_id')
        for column in ['donor_id', 'sample_id', 'disease', 'label', 'sex', 'age', 'age_years',
                       'tissue', 'tissue_site', 'SubCluster', 'author_subcluster',
                       'author_SubCluster', 'author_subject', 'split', 'split_status']:
            if column in obs:
                item['fields'][column] = {
                    'nunique': int(obs[column].nunique()), 'missing': int(obs[column].isna().sum()),
                    'values': unique_values(obs, column) if obs[column].nunique() <= 80 else []}
        if obs[['donor_id', 'sample_id', 'disease', 'label']].isna().any().any():
            item['warnings'].append('missing_required_metadata')
        labels = pd.to_numeric(obs['label'], errors='coerce')
        expected = obs['disease'].astype(str).map({'NF': 0, 'Control': 0, 'HCM': 1})
        item['labels_match_disease'] = bool(labels.eq(expected).all())
        item['disease_nuclei'] = {str(k): int(v) for k,v in obs['disease'].value_counts().items()}
        item['donors_per_disease'] = {str(k): int(v) for k,v in obs.groupby('disease', observed=True)['donor_id'].nunique().items()}
        item['donors_with_multiple_diseases'] = int((obs.groupby('donor_id', observed=True)['disease'].nunique() > 1).sum())
        item['samples_with_multiple_donors'] = int((obs.groupby('sample_id', observed=True)['donor_id'].nunique() > 1).sum())
        for donor, rows in obs.groupby('donor_id', observed=True):
            row = {'dataset': dataset, 'donor_id': str(donor), 'disease': ';'.join(unique_values(rows, 'disease')),
                   'nuclei': len(rows), 'samples': ';'.join(unique_values(rows, 'sample_id')),
                   'n_samples': rows['sample_id'].nunique()}
            for col in ['sex', 'age', 'age_years', 'tissue', 'tissue_site', 'donor_original',
                        'donor_id_original', 'author_subject', 'subject', 'Subject']:
                if col in rows:
                    row[col] = ';'.join(unique_values(rows, col))
            donor_rows.append(row)
            for col in ['donor_id', 'donor_original', 'donor_id_original', 'author_subject', 'subject', 'Subject']:
                for value in unique_values(rows, col):
                    alias_rows.append({'dataset': dataset, 'donor_id': str(donor), 'field': col, 'alias': value})
        for sample, rows in obs.groupby('sample_id', observed=True):
            sample_rows.append({'dataset': dataset, 'sample_id': str(sample),
                                'donor_id': ';'.join(unique_values(rows, 'donor_id')),
                                'disease': ';'.join(unique_values(rows, 'disease')), 'nuclei': len(rows)})
        subtype_fields = [c for c in ['SubCluster', 'author_subcluster', 'author_SubCluster',
                                     'cell_type_original', 'author_celltype', 'author_cell_type'] if c in obs]
        for column in subtype_fields:
            counts = obs.groupby(['donor_id', 'disease', column], observed=True, dropna=False).size()
            for (donor, disease, subtype), count in counts.items():
                subtype_rows.append({'dataset': dataset, 'donor_id': str(donor), 'disease': str(disease),
                                     'annotation_field': column, 'subtype': str(subtype), 'nuclei': int(count)})
        gene_ids = var['ensembl_id'].astype(str) if 'ensembl_id' in var else var.index.astype(str)
        item['gene_axis'] = {'ensembl_field': 'ensembl_id' if 'ensembl_id' in var else 'index',
                             'ensembl_format_count': int(gene_ids.str.match(r'^ENSG\d+(?:\.\d+)?$').sum()),
                             'gene_ids_unique': bool(gene_ids.is_unique), 'source_genes': len(var)}
        if dataset == 'liu2023':
            item['warnings'].append('Only two independent control donors; CL0412A/B must remain together')
        if dataset.startswith('chin_'):
            item['warnings'].append('Accession completely confounded with disease in Chin pair')
        if dataset == 'SCP1303':
            item['warnings'].append('X/counts are author CellBender-adjusted; original CellRanger layer is separate')
        after = path.stat()
        item['size_and_mtime_unchanged'] = before.st_size == after.st_size and before.st_mtime_ns == after.st_mtime_ns
        evidence['datasets'][dataset] = item
        print(json.dumps({'checked': dataset, 'shape': item['shape'], 'donors_per_disease': item['donors_per_disease']}), flush=True)
    donors = pd.DataFrame(donor_rows)
    donors.to_csv(args.output/'DONOR_INVENTORY.tsv', sep='\t', index=False)
    pd.DataFrame(sample_rows).to_csv(args.output/'SAMPLE_INVENTORY.tsv', sep='\t', index=False)
    pd.DataFrame(subtype_rows).to_csv(args.output/'DONOR_SUBTYPES.tsv', sep='\t', index=False)
    aliases = pd.DataFrame(alias_rows).drop_duplicates()
    aliases.to_csv(args.output/'OBSERVED_DONOR_ALIASES.tsv', sep='\t', index=False)
    shared = []
    for alias, rows in aliases.groupby('alias'):
        if rows['dataset'].nunique() > 1:
            shared.append({'alias': alias, 'matches': rows.to_dict('records')})
    evidence['cross_dataset_observed_alias_matches'] = shared
    evidence['cross_study_donor_independence'] = 'unresolved; no common alias is not proof of independent people'
    evidence['total_nuclei'] = sum(d['shape'][0] for d in evidence['datasets'].values())
    evidence['total_dataset_scoped_donors'] = len(donors)
    evidence['source_config_sha256'] = hashlib.sha256(args.config.read_bytes()).hexdigest()
    (args.output/'PREFLIGHT.json').write_text(json.dumps(plain(evidence), ensure_ascii=False, indent=2)+'\n')
    print(json.dumps({'status': 'completed', 'output': str(args.output),
                      'total_nuclei': evidence['total_nuclei'], 'dataset_scoped_donors': len(donors)}), flush=True)


if __name__ == '__main__':
    main()
