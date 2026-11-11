"""Check the installed environment and the writable-mask fix before resuming."""
from importlib.metadata import version
import json
import numpy as np
import pandas as pd

frame = pd.DataFrame({'target_gene': ['A', '', 'B'], 'guide_id': ['g1', 'g2', None], 'batch': ['b1'] * 3})
valid = frame.fillna('').ne('').all(axis=1).to_numpy(dtype=bool, copy=True)
assert valid.flags.writeable
valid &= ~frame.isna().any(axis=1).to_numpy()
assert np.array_equal(valid, [True, False, False])
print(json.dumps({'writable_mask_regression': 'passed', 'versions': {
    name: version(name) for name in ['scanpy', 'anndata', 'numpy', 'pandas', 'scipy', 'h5py']
}}, indent=2))
