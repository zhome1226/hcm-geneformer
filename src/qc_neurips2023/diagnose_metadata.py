"""Record why the first independent metadata comparison stopped."""
from pathlib import Path
import json
import anndata as ad
import h5py
import numpy as np
import pandas as pd
from verify_neurips2023 import native_columns_preserved

ROOT = Path.home() / 'xiaoe-internship-data'
WORK = ROOT / '_processing/all_data_cleaning_20261111/NeurIPS2023'
OUT = ROOT / 'processed/NeurIPS2023_qc_20261111'
with h5py.File(ROOT / 'benchmarks/competitions/neurips2023/sc_counts_reannotated_with_counts.h5ad', 'r') as handle:
    source = ad.io.read_elem(handle['var'])
saved = pd.read_csv(OUT / 'sc_reannotated_feature_QC.tsv.gz', sep='\t', index_col=0, keep_default_na=False, float_precision='round_trip')
tokens = pd.read_csv(OUT / 'sc_reannotated_feature_QC.tsv.gz', sep='\t', index_col=0, keep_default_na=False, dtype=str)
report = {'index_equal': source.index.equals(saved.index), 'columns': {}}
report['corrected_verifier_feature_comparison'] = native_columns_preserved(source, saved)
assert report['corrected_verifier_feature_comparison']
for name in source:
    column = source[name]
    info = {'source_dtype': str(column.dtype), 'saved_dtype': str(saved[name].dtype)}
    if pd.api.types.is_numeric_dtype(column) and not pd.api.types.is_bool_dtype(column):
        numeric = pd.to_numeric(saved[name], errors='coerce').to_numpy()
        original = column.to_numpy()
        mismatch = ~(np.equal(numeric, original) | (np.isnan(numeric) & np.isnan(original)))
        exact_dtype = np.asarray([float(value) if value else np.nan for value in tokens[name]], dtype=original.dtype)
        info.update(mismatches_after_float64_comparison=int(mismatch.sum()),
            roundtrip_to_source_dtype_exact=bool(np.array_equal(exact_dtype, original, equal_nan=True)))
        if mismatch.any():
            first = int(np.flatnonzero(mismatch)[0])
            info.update(example_feature=str(source.index[first]), source_value=float(original[first]),
                saved_value=float(numeric[first]), serialized_token=tokens[name].iloc[first])
    report['columns'][name] = info
(WORK / 'METADATA_DTYPE_DIAGNOSIS.json').write_text(json.dumps(report, indent=2) + '\n')
print(json.dumps(report, indent=2))
