"""Build an adult-LV-nuclear annotation reference with held-out reference donors."""
from datetime import datetime, timezone
from pathlib import Path
import gzip
import hashlib
import json

import anndata as ad
import h5py
import numpy as np
import pandas as pd
from scipy.io import mmwrite

root = Path.home() / 'xiaoe-internship-data'
out = root / 'processed/SCP1303_review_20261111/identity'
out.mkdir(exist_ok=True)
refpath = root / 'references/Litvinukova2020/heart_atlas_49312f8c.h5ad'
seed = 20261111
mapping = {
 'regular ventricular cardiac myocyte': 'Cardiomyocyte',
 'regular atrial cardiac myocyte': 'Cardiomyocyte',
 'pericyte': 'Pericyte', 'fibroblast': 'Fibroblast',
 'capillary endothelial cell': 'Endothelial', 'endothelial cell of artery': 'Endothelial',
 'endothelial cell': 'Endothelial', 'vein endothelial cell': 'Endothelial',
 'endothelial cell of lymphatic vessel': 'Endothelial',
 'smooth muscle cell': 'Smooth_muscle', 'neural cell': 'Neural',
 'epicardial adipocyte': 'Adipocyte', 'mesothelial cell': 'Mesothelial',
 'macrophage': 'Immune', 'natural killer cell': 'Immune',
 'CD14-positive, CD16-positive monocyte': 'Immune',
 'CD4-positive, alpha-beta cytotoxic T cell': 'Immune',
 'activated CD8-positive, alpha-beta T cell': 'Immune',
 'CD8-positive, alpha-beta cytotoxic T cell': 'Immune',
 'CD14-positive monocyte': 'Immune', 'monocyte': 'Immune',
 'mast cell': 'Immune', 'mature NK T cell': 'Immune', 'B cell': 'Immune',
 'activated CD4-positive, alpha-beta T cell': 'Immune', 'dendritic cell': 'Immune',
}
with h5py.File(refpath) as f:
    obs = ad.io.read_elem(f['obs'])
    var = ad.io.read_elem(f['raw/var'])
    mask = obs['tissue'].eq('heart left ventricle') & obs['suspension_type'].eq('nucleus')
    mask &= obs['disease'].eq('normal') & obs['is_primary_data'].astype(bool)
    selected = obs.loc[mask].copy()
    selected['source_row'] = np.flatnonzero(mask)
    selected['broad_type'] = selected['cell_type'].astype(str).map(mapping)
    excluded_unknown = int(selected['broad_type'].isna().sum())
    selected = selected.loc[selected['broad_type'].notna()].copy()
    rng = np.random.default_rng(seed)
    donors = np.sort(selected['donor_id'].astype(str).unique())
    heldout = sorted(rng.choice(donors, size=4, replace=False))
    selected['partition'] = np.where(selected['donor_id'].astype(str).isin(heldout), 'calibration', 'reference')
    chosen = []
    for (partition, donor, label), group in selected.groupby(['partition', 'donor_id', 'broad_type'], observed=True, sort=True):
        cap = 100 if partition == 'calibration' else 60
        chosen.extend(rng.choice(group['source_row'].to_numpy(), min(cap, len(group)), replace=False))
    chosen = np.sort(chosen)
    metadata = selected.set_index('source_row').loc[chosen].copy()
    metadata.index.name = 'source_row'
    metadata['reference_barcode'] = obs.index[chosen]
    genes = pd.read_csv(out.parent / 'query_genes.tsv', sep='\t')
    common = genes.loc[genes['gene_id'].isin(var.index)].copy()
    positions = var.index.get_indexer(common['gene_id'])
    if not common['gene_id'].is_unique or (positions < 0).any():
        raise ValueError('Reference gene mapping must be unique and exact')
    counts = ad.io.sparse_dataset(f['raw/X'])[chosen, :][:, positions].tocsr()
    if (counts.data < 0).any() or not np.isfinite(counts.data).all() or not np.allclose(counts.data, np.round(counts.data)):
        raise ValueError('Reference raw layer is not nonnegative integer counts')
    common.to_csv(out / 'reference_genes.tsv', sep='\t', index=False)
    metadata.to_csv(out / 'reference_selection.tsv', sep='\t', index=False)
    for partition in ['reference', 'calibration']:
        keep = metadata['partition'].eq(partition).to_numpy()
        metadata.loc[keep].to_csv(out / f'{partition}_metadata.tsv', sep='\t', index=False)
        with gzip.open(out / f'{partition}_counts.mtx.gz', 'wb') as handle:
            mmwrite(handle, counts[keep].T, comment='Unnormalized raw adult left ventricular nuclear counts')

counts_table = metadata.groupby(['partition', 'broad_type'], observed=True).agg(
    nuclei=('reference_barcode', 'size'), donors=('donor_id', 'nunique')).reset_index()
counts_table.to_csv(out / 'reference_design.tsv', sep='\t', index=False)
design = {'created_utc': datetime.now(timezone.utc).isoformat(), 'seed': seed,
 'source': str(refpath), 'reference_donor_ids': sorted(set(donors)-set(heldout)),
 'heldout_reference_donor_ids': heldout, 'broad_label_mapping': mapping,
 'reference_regions': ['heart left ventricle'], 'suspension_type': 'nucleus',
 'reference_gene_count': len(common), 'excluded_unknown_nuclei': excluded_unknown,
 'selection': 'Fixed seed; max 60 reference or 100 calibration nuclei per donor and broad cell type',
 'pruning_plan': 'Report SingleR default nmads=3 pruning; additionally use the fifth percentile of correct held-out reference delta and delta.next by assigned type when at least 30 correct examples are available. These are empirical reference-tail flags, not posterior confidence.',
 'independence_note': 'Reference donors held out from reference training. Biological donor overlap between studies remains unestablished; this reference is not a disease evaluation set.',
 'hashes': {}}
for path in sorted(out.glob('*.tsv')):
    with path.open('rb') as handle:
        design['hashes'][path.name] = hashlib.file_digest(handle, 'sha256').hexdigest()
(out / 'REFERENCE_DESIGN.json').write_text(json.dumps(design, indent=2) + '\n')
print(counts_table.to_string(index=False), flush=True)
print('REFERENCE_INPUTS_READY', flush=True)
