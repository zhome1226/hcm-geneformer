"""Deterministic filtering, donor splits and sampling for the HCM QC workflow."""
from __future__ import annotations

import numpy as np
import pandas as pd


def baseline_mask(genes, mt_pct, minimum=500, maximum=5000, mt_limit=5.0):
    genes, mt_pct = np.asarray(genes), np.asarray(mt_pct)
    return (np.isfinite(genes) & np.isfinite(mt_pct) &
            (genes >= minimum) & (genes <= maximum) & (mt_pct < mt_limit))


def make_donor_splits(obs: pd.DataFrame, seed=20261111):
    """Provisional 70/15/15 split, stratified on disease, exclusively by donor."""
    groups = obs.groupby('donor_id', observed=True)['disease'].agg(['first', 'nunique'])
    if (groups['nunique'] != 1).any() or obs['donor_id'].isna().any():
        raise ValueError('Missing donor or inconsistent donor disease')
    rng = np.random.default_rng(seed)
    split_map = {}
    names = ['train', 'validation', 'test']
    for disease in sorted(groups['first'].unique()):
        donors = np.sort(groups.index[groups['first'] == disease].astype(str).to_numpy())
        if len(donors) < 3:
            raise ValueError('Need at least 3 donors per disease for provisional split')
        rng.shuffle(donors)
        exact = len(donors) * np.array([.7, .15, .15])
        allocation = np.floor(exact).astype(int)
        for i in np.argsort(-(exact - allocation), kind='stable')[:len(donors)-sum(allocation)]:
            allocation[i] += 1
        for j in (1, 2):
            if allocation[j] == 0:
                allocation[0] -= 1
                allocation[j] += 1
        start = 0
        for name, n in zip(names, allocation):
            split_map.update({str(d): name for d in donors[start:start+n]})
            start += n
    return split_map


def balanced_donor_sample(obs: pd.DataFrame, total=20000, seed=20261111):
    """Balance diseases and spread each allocation across all eligible donors."""
    rng = np.random.default_rng(seed)
    labels = sorted(obs['label'].unique())
    if labels != [0, 1]:
        raise ValueError('Require both NF=0 and HCM=1')
    target_per_label = min(total // 2, int(obs['label'].value_counts().min()))
    selected = []
    for label in labels:
        part = obs.loc[obs['label'] == label]
        donor_cells = {str(d): g.index.to_numpy() for d, g in part.groupby('donor_id', observed=True)}
        donors = np.array(sorted(donor_cells))
        rng.shuffle(donors)
        allocation = {d: 0 for d in donors}
        remaining = target_per_label
        while remaining:
            available = [d for d in donors if allocation[d] < len(donor_cells[d])]
            if not available:
                raise ValueError('Insufficient distinct cells')
            quota = max(1, remaining // len(available))
            for d in available:
                take = min(quota, len(donor_cells[d]) - allocation[d], remaining)
                allocation[d] += take
                remaining -= take
                if not remaining:
                    break
        for donor in donors:
            selected.extend(rng.choice(donor_cells[donor], size=allocation[donor], replace=False).tolist())
    rng.shuffle(selected)
    if len(selected) != len(set(selected)):
        raise AssertionError('Sampling duplicated cells')
    return selected
