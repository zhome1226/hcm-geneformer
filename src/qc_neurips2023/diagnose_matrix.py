"""Record exact differences in the baseline ATAC processed-X storage audit."""
import json
import h5py
import numpy as np
from verify_neurips2023 import SOURCE, WORK, read_table, numeric_audit

recorded = read_table('MATRIX_QC_SUMMARY.tsv').set_index(['unit', 'relative_slot']).loc[('baseline_atac', 'X')]
with h5py.File(SOURCE / 'multiome_counts.h5mu', 'r') as handle:
    actual, _ = numeric_audit(handle['mod/atac/X'], 'processed_nonnegative')
result = {'component': 'baseline_atac:X', 'fields': {}}
result['corrected_native_storage_comparison'] = bool(actual['has_noncanonical_storage'] == (recorded.noncanonical_blocks > 0))
assert result['corrected_native_storage_comparison']
for name, value in actual.items():
    previous = recorded.get(name, None)
    if isinstance(value, np.generic):
        value = value.item()
    if isinstance(previous, np.generic):
        previous = previous.item()
    result['fields'][name] = {'recorded': previous, 'independent': value,
        'equal': bool(value == previous) if previous is not None else None}
(WORK / 'MATRIX_SUMMARY_DIAGNOSIS.json').write_text(json.dumps(result, indent=2) + '\n')
print(json.dumps(result, indent=2))
