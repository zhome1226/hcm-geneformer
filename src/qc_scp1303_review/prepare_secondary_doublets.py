"""Export complete author-retained 10x libraries for triggered residual-doublet review."""
from datetime import datetime, timezone
from pathlib import Path
import gzip
import json

import anndata as ad
import h5py
import numpy as np
import pandas as pd
from scipy.io import mmwrite

root = Path.home() / 'xiaoe-internship-data'
review = root / 'processed/SCP1303_review_20261111/identity'
out = review / 'doublets'
out.mkdir(exist_ok=True)
evidence = pd.read_csv(review / 'joined_identity_evidence.tsv.gz', sep='\t', index_col='barcode', low_memory=False)
captures = sorted(evidence.loc[evidence['secondary_doublet_review_trigger'], 'capture_id'].unique())
with h5py.File(root / '_processing/SCP1303_20261111/author_input.h5ad') as f:
    obs = ad.io.read_elem(f['obs'])
    var = ad.io.read_elem(f['var'])
    per_library_donors = obs.groupby('biosample_id', observed=True)['donor_id'].nunique()
    if (per_library_donors > 1).any():
        raise ValueError('Library contains multiple donor identities; revisit capture assumptions')
    # Ample host memory permits one read instead of many slow backed random reads.
    matrix = ad.io.read_elem(f['layers/cellranger_raw'])
    records = []
    for number, capture in enumerate(captures, 1):
        mask = obs['biosample_id'].astype(str).eq(capture).to_numpy()
        rows = np.flatnonzero(mask)
        path = out / capture
        path.mkdir(exist_ok=True)
        complete = False
        if (path / 'counts.mtx.gz').exists() and (path / 'cells.tsv').exists():
            saved = pd.read_csv(path / 'cells.tsv', sep='\t')
            if list(saved['barcode']) == list(obs.index[rows]):
                try:
                    with gzip.open(path / 'counts.mtx.gz', 'rb') as handle:
                        while handle.read(8 * 1024**2):
                            pass
                    complete = True
                except (EOFError, OSError):
                    complete = False
        if not complete:
            counts = matrix[rows, :].tocsr()
            if not np.isfinite(counts.data).all() or (counts.data < 0).any():
                raise ValueError('Invalid CellRanger counts')
            with gzip.open(path / 'counts.mtx.gz', 'wb', compresslevel=1) as handle:
                mmwrite(handle, counts.T)
            pd.DataFrame({'barcode': obs.index[rows], 'donor_id': obs.iloc[rows]['donor_id'].astype(str).to_numpy(),
                'author_type': obs.iloc[rows]['cell_type_leiden0.6'].astype(str).to_numpy()}).to_csv(path / 'cells.tsv', sep='\t', index=False)
        records.append({'capture_id': capture, 'nuclei': len(rows),
            'author_types': int(obs.iloc[rows]['cell_type_leiden0.6'].nunique())})
        print(f'Prepared {number}/{len(captures)} {capture}: {len(rows)} nuclei (resumed={complete})', flush=True)
pd.DataFrame(records).to_csv(out / 'capture_inputs.tsv', sep='\t', index=False)
design = {'prepared_utc': datetime.now(timezone.utc).isoformat(), 'capture_count': len(records),
    'total_nuclei': sum(item['nuclei'] for item in records), 'input_layer': 'cellranger_raw',
    'capture_definition': 'Author README explicitly defines biosample_id as the sequencing-library identifier. Each biosample maps to one donor; use as the 10x library/capture unit.',
    'input_limitation': 'Complete cell-type composition of each AUTHOR-RETAINED library, not original unfiltered droplets. Cannot reconstruct the original doublet removal.',
    'rate_assumption': 'Residual dbr=0.01, dbr.sd=0.015; no experimentally known residual rate. Calls are supporting evidence and never a standalone exclusion.',
    'seed': 20261111}
(out / 'DOUBLET_DESIGN.json').write_text(json.dumps(design, indent=2) + '\n')
print('SECONDARY_DOUBLET_INPUTS_READY', flush=True)
