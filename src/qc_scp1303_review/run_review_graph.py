"""Diagnostic graph and paired background evidence; no filtering or count export."""
from datetime import datetime, timezone
from pathlib import Path
import gc
import gzip
import json

import anndata as ad
import h5py
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc
from scipy.io import mmread

root = Path.home() / 'xiaoe-internship-data'
base = root / 'processed/SCP1303_review_20261111'
out = base / 'identity'
seed = 20261111
o = pd.read_csv(base / 'query_metadata.tsv', sep='\t', index_col='barcode', low_memory=False)
genes = pd.read_csv(base / 'query_genes.tsv', sep='\t')
with gzip.open(base / 'query_counts.mtx.gz', 'rb') as handle:
    counts = mmread(handle, spmatrix=True).T.tocsr()
assert counts.shape == (len(o), len(genes))
obj = ad.AnnData(counts.astype('float32'), obs=o.copy(),
    var=pd.DataFrame(index=genes['gene_id'].to_numpy()))
obj.var['gene_symbol'] = genes['gene_symbol'].to_numpy()
sc.pp.normalize_total(obj, target_sum=10000)
sc.pp.log1p(obj)
sc.pp.highly_variable_genes(obj, n_top_genes=3000, flavor='seurat', subset=True)
sc.pp.scale(obj, zero_center=False, max_value=10)
sc.tl.pca(obj, n_comps=50, random_state=seed, svd_solver='arpack')
sc.pp.neighbors(obj, n_neighbors=30, n_pcs=40, random_state=seed,
                method='gauss', transformer='sklearn')
for res in [0.5, 1.0]:
    sc.tl.leiden(obj, resolution=res, key_added=f'leiden_{res}',
                 flavor='igraph', directed=False, n_iterations=3, random_state=seed)
original = o['cell_type_original'].astype(str)
is_cm = original.str.startswith('Cardiomyocyte').to_numpy()
distance = obj.obsp['distances'].tocsr()
neighbors = np.array([is_cm[distance.indices[distance.indptr[i]:distance.indptr[i+1]]].mean()
                      for i in range(len(o))])
result = obj.obs[['leiden_0.5', 'leiden_1.0']].copy()
result['author_cm_neighbor_fraction'] = neighbors
result['PC1'] = obj.obsm['X_pca'][:, 0]
result['PC2'] = obj.obsm['X_pca'][:, 1]
result.to_csv(out / 'query_graph_evidence.tsv', sep='\t', index_label='barcode')
groups = np.where(is_cm, 'Author cardiomyocyte', 'Other author types')
fig, axes = plt.subplots(1, 2, figsize=(13, 5))
for label in np.unique(groups):
    k = groups == label
    axes[0].scatter(result.loc[k, 'PC1'], result.loc[k, 'PC2'], s=2, alpha=.35, label=label)
for label, k in [('CM_III', original.eq('Cardiomyocyte_III').to_numpy()),
    ('High-gene candidates', o['review_high_gene'].to_numpy()),
    ('Baseline warnings', o['review_baseline_warning'].to_numpy())]:
    axes[1].scatter(result.loc[k, 'PC1'], result.loc[k, 'PC2'], s=3, alpha=.5, label=label)
for ax in axes:
    ax.legend(markerscale=3, fontsize=8)
    ax.set_xlabel('PC1'); ax.set_ylabel('PC2')
fig.suptitle('Diagnostic review context; author labels are not independent ground truth')
fig.tight_layout(); fig.savefig(out / 'review_context_pca.png', dpi=160); plt.close(fig)
obj.uns['review_scope'] = 'Diagnostic normalized HVG graph; never use this object as raw training counts'
obj.write_h5ad(out / 'diagnostic_graph.h5ad', compression='gzip')
del obj
gc.collect()
panels = json.loads((base / 'REGISTRY_SUMMARY.json').read_text())['panels']
with h5py.File(root / '_processing/SCP1303_20261111/author_input.h5ad') as f:
    var = ad.io.read_elem(f['var'])
    symbols = list(dict.fromkeys(g for values in panels.values() for g in values if g in var.index))
    positions = var.index.get_indexer(symbols)
    # Explicit symbol lookup preserves the exact query gene order.
    qpos = pd.Index(genes['gene_symbol']).get_indexer(symbols)
    cb = counts[:, qpos].toarray()
    raw = ad.io.sparse_dataset(f['layers/cellranger_raw'])[o['original_row'].to_numpy(), :][:, positions].toarray()
paired = pd.DataFrame(index=o.index)
for label, panel in panels.items():
    ids = [symbols.index(g) for g in panel if g in symbols]
    paired[label + '_cellbender_fraction'] = cb[:, ids].sum(axis=1) / np.maximum(o['total_counts'].to_numpy(), 1)
    paired[label + '_cellranger_fraction'] = raw[:, ids].sum(axis=1) / np.maximum(o['total_counts_cellranger'].to_numpy(), 1)
    paired[label + '_raw_markers_detected'] = (raw[:, ids] > 0).sum(axis=1)
paired.to_csv(out / 'paired_background_evidence.tsv.gz', sep='\t', index_label='barcode', compression='gzip')
summary = {'completed_utc': datetime.now(timezone.utc).isoformat(), 'nuclei': len(o),
    'seed': seed, 'hvg': 3000, 'pcs': 50, 'neighbors': 30,
    'graph_method': 'Exact sklearn neighbors with Gaussian connectivity; PCA plot. Avoids installed UMAP/Numba layout incompatibility.',
    'note': 'Author-neighbor support is correlated evidence. Raw/CellBender differences alone do not prove contamination.'}
(out / 'GRAPH_SUCCESS.json').write_text(json.dumps(summary, indent=2) + '\n')
print('REVIEW_GRAPH_COMPLETE', flush=True)
