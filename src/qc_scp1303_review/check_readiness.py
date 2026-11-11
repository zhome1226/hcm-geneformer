"""Read-only readiness and registry audit; never filters or rewrites counts."""
from datetime import datetime, timezone
from importlib import metadata, util
from pathlib import Path
import hashlib
import json
import sys

import numpy as np
import pandas as pd

root = Path.home() / 'xiaoe-internship-data'
work = root / '_processing/SCP1303_review_20261111'
out = root / 'processed/SCP1303_review_20261111'
ref = root / 'references/Litvinukova2020'
report = {'checked_utc': datetime.now(timezone.utc).isoformat(),
          'python': sys.version, 'executable': sys.executable}
report['versions'] = {}
for package in ['scanpy', 'anndata', 'numpy', 'scipy', 'pandas', 'h5py',
                'matplotlib', 'umap-learn', 'igraph', 'leidenalg']:
    try:
        report['versions'][package] = metadata.version(package)
    except metadata.PackageNotFoundError:
        report['versions'][package] = None
report['clustering_modules'] = {name: util.find_spec(name) is not None
                                for name in ['igraph', 'leidenalg', 'umap']}
report['singleR_smoke_test'] = (work / 'R_SMOKE_TEST.txt').read_text().strip()
summary = json.loads((out / 'REGISTRY_SUMMARY.json').read_text())
report['registry'] = {k: summary[k] for k in ['status', 'membership_checks',
    'context_nuclei', 'context_cell_types', 'context_donors', 'priority_light_nuclei']}
report['registry_file_checks'] = []
for entry in summary['files']:
    path = out / entry['file']
    with path.open('rb') as handle:
        digest = hashlib.file_digest(handle, 'sha256').hexdigest()
    report['registry_file_checks'].append({'file': entry['file'],
        'bytes_match': path.stat().st_size == entry['bytes'],
        'sha256_match': digest == entry['sha256']})
query = pd.read_csv(out / 'query_metadata.tsv', sep='\t', index_col='barcode', low_memory=False)
evidence = pd.read_csv(out / 'marker_program_evidence.tsv.gz', sep='\t', index_col='barcode')
report['registry_structure_checks'] = {
    'unique_barcodes': bool(query.index.is_unique),
    'evidence_same_order': bool(query.index.equals(evidence.index)),
    'nuclei_match': len(query) == summary['context_nuclei'],
    'pending_candidates': int(query['review_action'].eq('pending_review').sum()),
    'context_only': int(query['review_action'].eq('context_only').sum()),
    'sample_id_count': int(query['sample_id'].nunique()),
    'samples_with_multiple_cell_types': int((query.groupby('sample_id')['cell_type_original'].nunique() > 1).sum()),
}
groups = {
    'priority_26': evidence['priority_two_panels_absent_light'],
    'HCM_CM_III': query['disease'].eq('HCM') & query['cell_type_original'].eq('Cardiomyocyte_III'),
    'P1735_at_most_5000': query['review_p1735'] & query['n_genes_by_counts'].le(5000),
    'P1735_above_5000': query['review_p1735'] & query['n_genes_by_counts'].gt(5000),
}
report['initial_evidence'] = {}
for name, mask in groups.items():
    subset = query.loc[mask]
    ev = evidence.loc[mask]
    report['initial_evidence'][name] = {
        'nuclei': len(subset), 'donors': int(subset['donor_id'].nunique()),
        'diseases': subset['disease'].value_counts().to_dict(),
        'author_types': subset['cell_type_original'].value_counts().to_dict(),
        'detected_genes_median': float(subset['n_genes_by_counts'].median()),
        'author_doublet_score_median': float(subset['doublet_score'].median()),
        'cardiomyocyte_panel_detected_median': float(ev['cardiomyocyte_markers_detected'].median()),
        'any_expanded_cm_marker': int(ev['cardiomyocyte_markers_detected'].gt(0).sum()),
    }
catalog = json.loads((ref / 'collection.json').read_text())
dataset = next(item for item in catalog['datasets'] if item['dataset_id'] == 'd4e69e01-3ba2-4d6b-a15d-e7048f78f22e')
asset = next(item for item in dataset['assets'] if item['filetype'] == 'H5AD')
report['reference'] = {'doi': catalog.get('doi'),
    'collection_id': catalog['collection_id'], 'dataset_id': dataset['dataset_id'],
    'dataset_version_id': dataset['dataset_version_id'],
    'title': dataset.get('title'), 'expected_bytes': int(asset['filesize']),
    'source_url': asset['url']}
target = ref / 'heart_atlas_49312f8c.h5ad'
report['reference']['download_record'] = json.loads((ref / 'DOWNLOAD_VERIFIED.json').read_text()) if (ref / 'DOWNLOAD_VERIFIED.json').exists() else None
report['reference']['assembled_file_bytes'] = target.stat().st_size if target.exists() else 0
report['reference']['completed_range_bytes'] = sum(p.stat().st_size for p in (ref / 'range_parts').glob('*.part'))
report['reference']['incoming_range_bytes'] = sum(p.stat().st_size for p in (ref / 'range_parts').glob('*.incoming'))
if report['reference']['download_record']:
    import anndata as ad
    import h5py
    with h5py.File(target, 'r') as handle:
        obs = ad.io.read_elem(handle['obs'])
        var = ad.io.read_elem(handle['var'])
        structure = {'top_level_keys': list(handle.keys()),
                     'layers': list(handle['layers'].keys()) if 'layers' in handle else [],
                     'has_raw': 'raw' in handle}
    report['reference']['metadata'] = {'nuclei_and_cells': len(obs), 'genes': len(var),
        'obs_columns': list(obs.columns), 'var_columns': list(var.columns),
        'categories': {key: obs[key].astype(str).value_counts().to_dict()
                       for key in ['tissue', 'suspension_type', 'assay', 'disease', 'cell_type']
                       if key in obs.columns}}
    report['reference']['metadata']['matrix_structure'] = structure
    for name, tissues in {'left_ventricle_nuclei': ['heart left ventricle'],
        'ventricular_regions_nuclei': ['heart left ventricle', 'heart right ventricle',
                                       'interventricular septum', 'apex of heart']}.items():
        selected = obs['suspension_type'].eq('nucleus') & obs['tissue'].isin(tissues)
        subset = obs.loc[selected]
        report['reference']['metadata'][name] = {'nuclei': len(subset),
            'donors': int(subset['donor_id'].nunique()),
            'cell_types': subset['cell_type'].astype(str).value_counts().to_dict()}
    query_genes = pd.read_csv(out / 'query_genes.tsv', sep='\t')
    report['reference']['metadata']['query_ensembl_ids_in_reference_index'] = int(
        query_genes['gene_id'].isin(var.index).sum())
report['scope'] = 'Preparation only. No new identity labels, doublet calls, threshold changes, or H5AD exports.'
(out / 'READINESS_AUDIT.json').write_text(json.dumps(report, indent=2, ensure_ascii=False) + '\n')
print(json.dumps({'versions': report['versions'], 'registry': report['registry_structure_checks'],
    'reference': report['reference'], 'report': str(out / 'READINESS_AUDIT.json')}, indent=2))
