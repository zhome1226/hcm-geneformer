"""Meaningful boundary, leakage and limited-supply tests; no expression data."""
import unittest

import numpy as np
import pandas as pd

from qc_core import baseline_mask, balanced_donor_sample, make_donor_splits


class QCCoreTests(unittest.TestCase):
    def test_threshold_boundaries_and_missing(self):
        genes = [499, 500, 5000, 5001, 1000, 1000, np.nan]
        mt = [0, 4.999, 0, 0, 5, np.nan, 0]
        self.assertEqual(baseline_mask(genes, mt).tolist(), [False, True, True, False, False, False, False])

    def test_donors_never_cross_splits(self):
        o = pd.DataFrame([{'donor_id': f'{disease}_{d}', 'disease': disease}
                          for disease, count in [('NF', 16), ('HCM', 15)]
                          for d in range(count) for _ in range(4)])
        splits = make_donor_splits(o)
        o['split'] = o.donor_id.map(splits)
        self.assertEqual(int(o.groupby('donor_id').split.nunique().max()), 1)
        self.assertEqual(len(splits), 31)
        self.assertEqual(splits, make_donor_splits(o))
        self.assertEqual(set(o.groupby('disease').split.nunique()), {3})

    def test_donor_disease_conflict_raises(self):
        with self.assertRaises(ValueError):
            make_donor_splits(pd.DataFrame({'donor_id': ['A', 'A'], 'disease': ['NF', 'HCM']}))

    def test_sampling_preserves_small_donors_without_duplication(self):
        rows = []
        for label in [0, 1]:
            for donor, n in enumerate([2, 20, 100]):
                rows.extend({'label': label, 'donor_id': f'{label}_{donor}'} for _ in range(n))
        o = pd.DataFrame(rows)
        o.index = 'cell_' + o.index.astype(str)
        selected = balanced_donor_sample(o, total=100)
        self.assertEqual(len(selected), 100)
        self.assertEqual(len(set(selected)), 100)
        self.assertEqual(o.loc[selected, 'label'].value_counts().to_dict(), {0: 50, 1: 50})
        self.assertEqual(o.loc[selected, 'donor_id'].nunique(), 6)
        self.assertEqual(selected, balanced_donor_sample(o, total=100))

    def test_small_dataset_never_upsampled(self):
        o = pd.DataFrame({'label': [0, 0, 0, 1], 'donor_id': ['A', 'A', 'B', 'C']})
        self.assertEqual(len(balanced_donor_sample(o, total=20000)), 2)


if __name__ == '__main__':
    unittest.main()
