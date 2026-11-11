"""Resolve concentrated uncertainty with reference-derived programs and donor/cluster profiles."""
from pathlib import Path
import gzip
import json

import numpy as np
import pandas as pd
from scipy import sparse
from scipy.io import mmread, mmwrite

root = Path.home() / 'xiaoe-internship-data'
base = root / 'processed/SCP1303_review_20261111'
out = base / 'identity/group_review'
out.mkdir(exist_ok=True)
reference_dir = out.parent
genes = pd.read_csv(reference_dir / 'reference_genes.tsv', sep='\t')
rmeta = pd.read_csv(reference_dir / 'reference_metadata.tsv', sep='\t')
cmeta = pd.read_csv(reference_dir / 'calibration_metadata.tsv', sep='\t')
qmeta = pd.read_csv(base / 'query_metadata.tsv', sep='\t', index_col='barcode', low_memory=False)
qgenes = pd.read_csv(base / 'query_genes.tsv', sep='\t')
graph = pd.read_csv(reference_dir / 'query_graph_evidence.tsv', sep='\t', index_col='barcode')
def read(path):
    with gzip.open(path, 'rb') as handle:
        return mmread(handle, spmatrix=True).tocsr()
ref = read(reference_dir / 'reference_counts.mtx.gz').T.tocsr()
norm = sparse.diags(10000/np.maximum(np.asarray(ref.sum(axis=1)).ravel(), 1)) @ ref
labels = rmeta['broad_type'].to_numpy()
means = {label: np.asarray(norm[labels==label].mean(axis=0)).ravel() for label in np.unique(labels)}
cm = means['Cardiomyocyte']
other_max = np.max(np.stack([value for label, value in means.items() if label!='Cardiomyocyte']), axis=0)
detected = np.asarray((ref[labels=='Cardiomyocyte']>0).mean(axis=0)).ravel()
fold = np.log2((cm+.1)/(other_max+.1))
eligible = (detected>=.2) & (cm>=.5) & (fold>=1)
positions = np.flatnonzero(eligible)
positions = positions[np.argsort(-fold[positions], kind='stable')[:100]]
if len(positions)<30:
    raise ValueError('No sufficiently broad reference-derived cardiomyocyte program')
signature = genes.iloc[positions].copy()
signature['mean_CM_per_10k'] = cm[positions]
signature['max_other_mean_per_10k'] = other_max[positions]
signature['log2_CM_over_max_other'] = fold[positions]
signature['fraction_CM_detected'] = detected[positions]
signature.to_csv(out / 'reference_CM_program.tsv', sep='\t', index=False)
query = read(base / 'query_counts.mtx.gz').T.tocsr()
sig_positions = pd.Index(qgenes['gene_id']).get_indexer(signature['gene_id'])
sig = query[:, sig_positions]
program = pd.DataFrame(index=qmeta.index)
program['reference_CM_program_genes_detected'] = np.asarray((sig>0).sum(axis=1)).ravel()
program['reference_CM_program_count_fraction'] = np.asarray(sig.sum(axis=1)).ravel()/np.maximum(np.asarray(query.sum(axis=1)).ravel(), 1)
program.to_csv(out / 'query_CM_program.tsv', sep='\t', index_label='barcode')
profiles, metadata = [], []
cluster = graph['leiden_1.0'].astype(str)
for label, indexes in cluster.groupby(cluster, sort=True).groups.items():
    rows = qmeta.index.get_indexer(indexes)
    profiles.append(sparse.csr_matrix(query[rows].sum(axis=0)))
    metadata.append({'profile_id': 'cluster_'+label, 'cluster': label, 'donor_id':'ALL', 'nuclei':len(rows)})
for (label, donor), part in qmeta.assign(cluster=cluster).groupby(['cluster','donor_id'], sort=True, observed=True):
    if len(part)<20:
        continue
    rows = qmeta.index.get_indexer(part.index)
    profiles.append(sparse.csr_matrix(query[rows].sum(axis=0)))
    metadata.append({'profile_id': 'cluster_'+label+'|donor_'+str(donor), 'cluster': label, 'donor_id':str(donor), 'nuclei':len(rows)})
with gzip.open(out / 'query_group_counts.mtx.gz', 'wb', compresslevel=1) as handle:
    mmwrite(handle, sparse.vstack(profiles).T)
pd.DataFrame(metadata).to_csv(out / 'query_group_metadata.tsv', sep='\t', index=False)
qgenes.to_csv(out / 'query_group_genes.tsv', sep='\t', index=False)
cal = read(reference_dir / 'calibration_counts.mtx.gz').T.tocsr()
profiles, metadata = [], []
for (donor, label), part in cmeta.groupby(['donor_id','broad_type'], sort=True, observed=True):
    if len(part)<20:
        continue
    rows = part.index.to_numpy()
    profiles.append(sparse.csr_matrix(cal[rows].sum(axis=0)))
    metadata.append({'profile_id':'reference_donor_'+str(donor)+'|'+label,
        'donor_id':str(donor), 'truth':label, 'nuclei':len(rows)})
with gzip.open(out / 'calibration_group_counts.mtx.gz', 'wb', compresslevel=1) as handle:
    mmwrite(handle, sparse.vstack(profiles).T)
pd.DataFrame(metadata).to_csv(out / 'calibration_group_metadata.tsv', sep='\t', index=False)
(out / 'METHOD.json').write_text(json.dumps({
    'program_derivation': 'Training-reference nuclei only; CM mean>=0.5 counts/10k, detection>=20%, log2 fold over highest other broad type>=1; up to 100 highest-fold genes.',
    'grouping': 'Diagnostic Leiden resolution 1.0 and donor-by-cluster aggregates with >=20 nuclei, unnormalized counts pooled; reference-donor holdout type aggregates provide an applicability check.',
    'purpose': 'Assess sparse individual nuclei in a coherent CM group; no disease labels or disease-classifier performance used.',
    'candidate_retention_rule': 'Individual best reference label CM, independent broad-program >=3 genes, global cluster CM and >=80% of donor-profile nuclei assigned CM, no applicable non-CM call or multi-evidence doublet exclusion. Group evidence is correlated, not an independent vote.'}, indent=2)+'\n')
print('GROUP_INPUTS_READY', flush=True)
