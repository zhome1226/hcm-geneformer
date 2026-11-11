"""Exercise installed clustering dependencies using synthetic counts only."""
from datetime import datetime, timezone
from importlib.metadata import version
from pathlib import Path
import json

import anndata as ad
import igraph
import leidenalg
import numpy as np
import scanpy as sc
from scipy import sparse

rng = np.random.default_rng(20261111)
obj = ad.AnnData(sparse.csr_matrix(rng.poisson(2, size=(100, 60)).astype('float32')))
sc.pp.normalize_total(obj, target_sum=10000)
sc.pp.log1p(obj)
sc.tl.pca(obj, n_comps=12, random_state=20261111)
sc.pp.neighbors(obj, n_neighbors=8, n_pcs=12, random_state=20261111)
sc.tl.leiden(obj, flavor='igraph', directed=False, n_iterations=2,
             random_state=20261111, key_added='test_leiden')
partition = leidenalg.find_partition(igraph.Graph.Ring(12),
    leidenalg.RBConfigurationVertexPartition, seed=20261111)
assert len(obj.obs['test_leiden']) == 100
assert obj.obsp['connectivities'].shape == (100, 100)
assert len(partition.membership) == 12
out = Path.home() / 'xiaoe-internship-data/processed/SCP1303_review_20261111'
result = {'status': 'passed', 'checked_utc': datetime.now(timezone.utc).isoformat(),
    'versions': {name: version(name) for name in ['scanpy', 'anndata', 'igraph', 'leidenalg']},
    'checks': ['synthetic_sparse_counts', 'normalization', 'PCA', 'neighbors',
               'scanpy_igraph_Leiden', 'leidenalg_partition'],
    'scope': 'Synthetic smoke test only; no biological data modified.'}
(out / 'PYTHON_SMOKE_TEST.json').write_text(json.dumps(result, indent=2) + '\n')
print(json.dumps(result, indent=2))
