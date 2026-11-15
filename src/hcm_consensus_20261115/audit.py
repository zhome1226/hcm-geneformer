"""Independent reconciliation of full-gene consensus evidence; no model refitting.

Run after fit, coverage and integration complete: python code/audit.py WORK
The audit reads the preserved matrices and donor detail, independently recomputes
filters/BH/stability/eligibility/ranks, and preserves missing measurements as NA.
"""
import argparse
import hashlib
import json
import pickle
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import pandas as pd

OLD_LOCK_SHA = 'd2fc63ed14d76c8d291afd09c067c0339584f38550369c8691248e0126bfa09d'
PANEL = {'SCP1303': (114711, 31, 15), 'liu2023': (8547, 12, 10),
         'chin_GSE174691': (10150, 9, 9), 'chin_GSE161921': (3534, 4, 0),
         'GSE138262_whole_cell': (3595, 12, 5), 'GSE255296': (7516, 4, 4),
         'GSE292754': (3306, 2, 2)}
KEY = ['dataset', 'disease', 'gene_id']


def require(value, message):
    if not bool(value):
        raise ValueError(message)


def same(actual, expected, message):
    require(np.array_equal(np.asarray(actual), np.asarray(expected)), message)


def near(actual, expected, message, rtol=3e-9, atol=1e-12):
    a, b = (np.asarray(x, dtype=float) for x in [actual, expected])
    require(a.shape == b.shape and np.allclose(a, b, rtol=rtol, atol=atol, equal_nan=True), message)


def booleans(values):
    mapped = values.map({True: True, False: False, 'TRUE': True, 'FALSE': False,
                         'True': True, 'False': False})
    require(mapped.notna().all(), 'Unrecognized boolean in ' + str(values.name))
    return mapped.to_numpy(dtype=bool)


def sha(path):
    with Path(path).open('rb') as f:
        return hashlib.file_digest(f, 'sha256').hexdigest()


def bh(p):
    p = np.asarray(p, dtype=float)
    require(np.isfinite(p).all() and ((p >= 0) & (p <= 1)).all(), 'Invalid P values')
    order = np.argsort(p)
    ranked = np.minimum.accumulate((p[order] * len(p) / np.arange(1, len(p)+1))[::-1])[::-1]
    result = np.empty(len(p)); result[order] = np.minimum(1, ranked)
    return result


class Audit:
    def __init__(self, work):
        self.w = work
        self.phase = ''
        self.done = []
        self.metrics = {}

    def read(self, relative):
        return pd.read_csv(self.w/relative, sep='\t', keep_default_na=False, na_values=['NA', ''])

    def indexed(self, relative, index='gene_id'):
        t = self.read(relative)
        require(t[index].notna().all() and t[index].is_unique, str(relative)+' duplicate/missing IDs')
        return t.set_index(index)

    def json(self, relative):
        return json.loads((self.w/relative).read_text(encoding='utf-8'))

    def run(self, label, function):
        self.phase = label
        function()
        self.done.append(label)

    def provenance(self):
        for tag in ['prepare', 'fit', 'coverage', 'integrate']:
            require((self.w/f'logs/{tag}.exit').read_text().strip() == '0', tag+' did not finish successfully')
        preparation = self.json('input/PREPARATION.json')
        require(preparation['design_commit'] == 'aefa26c', 'Unexpected second-round design commit')
        require(preparation['prior_outcomes_seen'] is True and
                preparation['Liu_role'] == 'joint_exploratory_discovery_not_validation', 'Exploratory role label')
        require(preparation['design_sha256'] == sha(self.w/'code/DESIGN.md'), 'Design changed after snapshot')
        for r in preparation['files']:
            source, copied = Path(r['source']), self.w/r['path']
            require(source.is_file() and copied.is_file(), 'Snapshot/source missing: '+r['path'])
            require(copied.stat().st_size == r['bytes'] and sha(copied) == r['sha256'] == sha(source),
                    'Snapshot/source bytes changed: '+r['path'])
        previous = self.w/'input/previous/CANDIDATES.tsv'
        require(sha(previous) == OLD_LOCK_SHA, 'Original 200 candidates changed')
        self.old = self.indexed('input/previous/CANDIDATES.tsv')
        require(len(self.old) == 200, 'Original shortlist size')
        self.universe = self.indexed('input/GENE_UNIVERSE.tsv')
        self.ids = self.universe.index.tolist()
        same(self.ids, sorted(set(self.ids)), 'Full gene universe ordering')
        require(self.universe.index.str.fullmatch(r'ENSG[0-9]+').all(), 'Invalid Ensembl IDs')
        require(preparation['gene_universe'] == len(self.ids) and
                preparation['gene_universe_sha256'] == sha(self.w/'input/GENE_UNIVERSE.tsv'), 'Universe seal')
        self.metrics.update(genes=len(self.ids), original_candidates=200, original_candidate_sha256=OLD_LOCK_SHA,
                            design_sha256=preparation['design_sha256'], mode='exploratory_joint_discovery')

    def counts(self):
        self.ro = self.read('input/pb/DONOR_ROSTER.tsv')
        require(not self.ro.duplicated(['dataset', 'donor_id']).any() and
                len(self.ro) == 74 and self.ro.nuclei.sum() == 151359, 'Donor panel size')
        for dataset, (n, donors, hcm) in PANEL.items():
            d = self.ro[self.ro.dataset.eq(dataset)]
            require(len(d) == donors and d.disease.eq('HCM').sum() == hcm and d.nuclei.sum() == n,
                    dataset+' donor/disease/cell counts')
        same(self.ro.model_train_test_split.unique(), ['unassigned'], 'Model split changed')
        same(self.ro.cross_study_person_match.unique(), ['unresolved'], 'Cross-study independence overclaimed')
        self.x, self.md, self.kept = {}, {}, {}
        for mod in ['scp_raw', 'scp_cb', 'scp_hhatl', 'liu', 'chin', 'whole']:
            x = self.indexed('input/pb/'+mod+'.tsv.gz')
            m = self.indexed('input/pb/'+mod+'_donors.tsv', 'donor_id')
            same(x.columns, m.index, mod+' donor axis')
            same(self.read('input/pb/'+mod+'_genes.tsv').gene_id, x.index, mod+' gene axis')
            require(np.issubdtype(x.values.dtype, np.integer) and (x.values >= 0).all(), mod+' integer counts')
            lib = x.sum(axis=0)
            same(lib, m.pseudobulk_library_size, mod+' count conservation')
            minimum = m.disease.value_counts().min()
            if minimum > 10:
                minimum = 10 + .7 * (minimum-10)
            threshold = 10 / np.median(lib) * 1e6
            # edgeR 4.8 calc_cpm_raw multiplies count by 1e6 BEFORE dividing.
            # Algebraic rearrangement changes IEEE-754 rounding at exact cutoffs:
            # Chin_2879 count=10 sits one ULP below the filter cutoff in edgeR.
            # Preserve the installed implementation, without adding tolerance.
            cpm = x.mul(1e6).div(lib, axis=1)
            keep = ((cpm >= threshold).sum(axis=1) >= minimum-1e-14) & (x.sum(axis=1) >= 15-1e-14)
            self.x[mod], self.md[mod], self.kept[mod] = x, m, list(x.index[keep])
        union = set().union(*(set(self.x[m].index) for m in ['scp_raw', 'liu', 'chin', 'whole']))
        same(self.ids, sorted(union), 'Union includes every analyzed source gene')
        same(self.md['liu'].index[self.md['liu'].disease.eq('Control')],
             ['Liu2023_CL0412', 'Liu2023_CL2'], 'Two distinct Liu controls')
        for mod in ['liu', 'whole']:
            a = self.read('input/pb/'+mod+'_feature_audit.tsv.gz')
            native = booleans(a.measured_in_all_included_samples)
            mapped = a.gene_id.fillna('').str.fullmatch(r'ENSG[0-9]+').to_numpy()
            same(booleans(a.included_in_model), native & mapped, mod+' source measurement mask')
            require(set(a.loc[native & mapped, 'gene_id']) == set(self.x[mod].index), mod+' measured gene axis')
        self.metrics['expression_filtered_genes'] = {m:len(self.kept[m]) for m in ['scp_raw', 'liu', 'chin', 'whole']}
        self.metrics['CPM_arithmetic'] = '(counts * 1e6) / library_size; edgeR calc_cpm_raw order, no added tolerance'

    def fits(self):
        s = self.indexed('input/scp/discovery.tsv')
        same(s.index, self.kept['scp_raw'], 'SCP independent expression filter')
        near(s.FDR, bh(s.PValue), 'SCP genome-wide BH')
        folds = self.indexed('input/scp/scp_loo_logFC.tsv.gz')
        same(folds.index, s.index, 'SCP omitted-donor gene axis')
        same(folds.columns, self.x['scp_raw'].columns, 'SCP 31 omitted donor IDs')
        require(len(folds.columns) == 31 and np.isfinite(folds.values).all(), 'SCP finite 31 folds')
        main = s.logFC.to_numpy()[:,None]
        for col, value in [('loo_same_direction_fraction', (np.sign(folds.values) == np.sign(main)).mean(axis=1)),
                           ('loo_min_abs_logFC', np.abs(folds.values).min(axis=1)),
                           ('loo_min_logFC', folds.min(axis=1)), ('loo_max_logFC', folds.max(axis=1)),
                           ('loo_max_abs_delta', np.abs(folds.values-main).max(axis=1))]:
            near(s[col], value, 'SCP independently recalculated '+col)
        design = self.indexed('input/scp/scp_design.tsv', 'donor_id')
        m = self.md['scp_raw'].loc[design.index]
        near(design.values, np.column_stack([np.ones(len(m)), m.age-m.age.mean(),
             m.sex.str.lower().isin(['m', 'male']), m.disease.eq('HCM')]), 'SCP covariate design')
        require(np.linalg.matrix_rank(design.values) == 4, 'SCP design rank')
        self.s = s.reindex(self.ids)
        self.stats = {}
        for mod in ['liu', 'chin', 'whole']:
            d = self.indexed('stats/'+mod+'_all.tsv.gz')
            same(d.index, self.x[mod].index, mod+' full source output axis')
            same(booleans(d.source_feature_present), np.ones(len(d), bool), mod+' feature availability')
            keep = d.index.isin(self.kept[mod])
            same(booleans(d.passes_expression_filter), keep, mod+' independently recomputed filter')
            if mod == 'liu':
                omit = [c for c in d if c.startswith('logFC_omit_')]
                same(sorted(omit), ['logFC_omit_Liu2023_CL0412', 'logFC_omit_Liu2023_CL2'], 'Liu omissions')
                for c in ['logFC', 'logCPM', 'F', 'PValue', 'FDR', *omit]:
                    same(np.isfinite(d[c]), keep, 'Liu missingness '+c)
                near(d.loc[keep, 'FDR'], bh(d.loc[keep, 'PValue']), 'Liu full-genome BH')
                previous = self.indexed('input/previous/liu.tsv')
                for c in previous:
                    a, b = d.loc[previous.index, c], previous[c]
                    if pd.api.types.is_numeric_dtype(b):
                        near(a, b, 'Liu original 200 unchanged '+c)
                    else:
                        same(a.fillna('NA'), b.fillna('NA'), 'Liu original 200 unchanged '+c)
                ld = self.indexed('stats/liu_design.tsv', 'donor_id')
                lm = self.md['liu'].loc[ld.index]
                near(ld.values, np.column_stack([np.ones(len(lm)), lm.disease.eq('HCM')]), 'Liu design')
            else:
                require(not set(d)&{'FDR', 'PValue', 'logFC'}, mod+' must remain descriptive')
                same(np.isfinite(d.descriptive_log2_cpm_difference), keep, mod+' descriptive missingness')
                lib = self.indexed('stats/'+mod+'_libraries.tsv', 'donor_id').reindex(self.x[mod].columns)
                measured = self.x[mod].loc[self.kept[mod]]
                near(lib.analysis_library_size, measured.sum(axis=0), mod+' filtered library sizes')
                near(lib.effective_library_size, lib.analysis_library_size*lib.TMM_factor, mod+' normalized library sizes')
                signal = np.log2(measured.mul(1e6).div(lib.effective_library_size, axis=1)+.5)
                group = self.md[mod].disease
                h = signal.loc[:, group.eq('HCM')].median(axis=1)
                c = signal.loc[:, group.eq('Control')].median(axis=1)
                near(d.loc[self.kept[mod], 'median_log2CPM_HCM'], h, mod+' equal-donor HCM median')
                near(d.loc[self.kept[mod], 'median_log2CPM_Control'], c, mod+' equal-donor control median')
                near(d.loc[self.kept[mod], 'descriptive_log2_cpm_difference'], h-c, mod+' descriptive effect')
            self.stats[mod] = d.reindex(self.ids)
        self.l = self.stats['liu']
        self.metrics['SCP_LODO_folds'] = 31
        self.metrics['Liu_BH_recomputed_from_all_tested_PValues'] = len(self.kept['liu'])

    def coverage(self):
        receipt = self.json('coverage/COVERAGE_RECEIPT.json')
        require(receipt['status'] == 'passed' and receipt['genes'] == len(self.ids), 'Coverage receipt scope')
        for relative, digest in receipt['inputs'].items():
            require(sha(self.w/relative) == digest, 'Coverage input hash '+relative)
        for relative, rec in receipt['outputs'].items():
            f = self.w/'coverage'/relative
            require(f.stat().st_size == rec['bytes'] and sha(f) == rec['sha256'], 'Coverage output hash '+relative)
        require(sha(receipt['vocab']['path']) == receipt['vocab']['sha256'], 'Vocabulary hash')
        with Path(receipt['vocab']['path']).open('rb') as f:
            vocabulary = pickle.load(f)
        self.token_ids = np.array([vocabulary.get(g, np.nan) for g in self.ids], dtype=float)
        total = self.read('coverage/ALL_COVERAGE.tsv.gz')
        require(not total.duplicated(KEY).any() and len(total) == 10*len(self.ids), 'Full coverage aggregate grid')
        per_group, seen = {}, set()
        expected_donors = self.ro.set_index(['dataset', 'donor_id'])
        modes = {'SCP1303':'scp_raw', 'liu2023':'liu', 'chin_GSE174691':'chin',
                 'chin_GSE161921':'chin', 'GSE138262_whole_cell':'whole'}
        with pd.read_csv(self.w/'coverage/DONOR_COVERAGE.tsv.gz', sep='\t', keep_default_na=False,
                         na_values=['NA', ''], chunksize=len(self.ids)) as chunks:
            for d in chunks:
                require(len(d) == len(self.ids) and d.dataset.nunique() == d.donor_id.nunique() == d.disease.nunique() == 1,
                        'Each coverage block must be exactly one donor and all genes')
                ds, donor, disease = d.iloc[0][['dataset', 'donor_id', 'disease']]
                key = (ds, donor)
                require(key not in seen and key in expected_donors.index, 'Unexpected/duplicate donor coverage')
                seen.add(key)
                origin = expected_donors.loc[key]
                require(origin.disease == disease and d.donor_cells.eq(origin.nuclei).all(), 'Donor label/denominator')
                same(d.gene_id, self.ids, 'Donor coverage gene axis')
                near(d.token_id, self.token_ids, 'Donor vocabulary identity')
                same(booleans(d.in_v1_vocab), np.isfinite(self.token_ids), 'Vocabulary eligibility')
                n = d.donor_cells.to_numpy(); t = d.token_cells.to_numpy()
                h = d.expressing_cells.to_numpy(); nc = d.native_cells.to_numpy()
                require(((0 <= t)&(t <= nc)&(0 <= h)&(h <= nc)&(nc <= n)&(d.count_sum.to_numpy() >= h)).all(),
                        'Per-donor coverage arithmetic')
                if ds != 'SCP1303':
                    require((t <= h).all(), 'Token without measured expression')
                near(d.expr_fraction_native, np.divide(h, nc, out=np.full(len(nc), np.nan), where=nc>0), 'Native detection denominator')
                near(d.token_fraction_all_cells, t/n, 'Token denominator')
                if ds in modes:
                    x = self.x[modes[ds]]
                    aligned = d.set_index('gene_id').loc[x.index]
                    same(aligned.count_sum, x[donor], ds+' counts vs independently prepared donor matrix')
                row = pd.DataFrame({'gene_id': self.ids, 'donors': np.ones(len(d), int), 'cells': n,
                    'native_cells': nc, 'expressing_cells': h, 'count_sum': d.count_sum.to_numpy(),
                    'donors_native': (nc>0).astype(int), 'donors_expressing': (h>0).astype(int),
                    'token_cells': t, 'donors_token_positive': (t>0).astype(int),
                    'donors_token_ge5pct': (20*t>=n).astype(int), 'donors_token_ge10pct': (10*t>=n).astype(int)}).set_index('gene_id')
                group = (ds, disease)
                if group not in per_group:
                    per_group[group] = [row, [t/n]]
                else:
                    per_group[group][0] += row
                    per_group[group][1].append(t/n)
        require(seen == set(expected_donors.index), 'Coverage must include all 74 donor records')
        self.cov = {}
        for (ds, disease), (sums, fractions) in per_group.items():
            s = total[total.dataset.eq(ds)&total.disease.eq(disease)].set_index('gene_id')
            same(s.index, self.ids, 'Aggregate coverage gene axis')
            sums['expr_fraction'] = sums.expressing_cells/sums.native_cells.replace(0, np.nan)
            sums['token_fraction'] = sums.token_cells/sums.cells
            sums['median_donor_token_fraction'] = np.median(np.asarray(fractions), axis=0)
            for c in sums:
                near(s[c], sums[c], ds+' '+disease+' independently aggregated '+c)
            near(s.token_id, self.token_ids, 'Aggregate token identity')
            self.cov[(ds, disease)] = s
        require(len(self.cov) == 10, 'Disease-specific coverage groups')
        self.metrics['donor_gene_coverage_rows'] = len(seen)*len(self.ids)

    def final(self):
        f = self.indexed('final/ALL_GENE_EVIDENCE.tsv.gz')
        same(f.index, self.ids, 'Final union axis and order')
        for c in ['gene_symbol', 'gene_symbol_aliases', 'source_axes']:
            same(f[c].fillna(''), self.universe[c].fillna(''), 'Final annotation '+c)
        for mod, data in [('scp', self.s), *self.stats.items()]:
            for c in data:
                if c in ['source_feature_present', 'passes_expression_filter']:
                    continue
                col = mod+'_'+c
                require(col in f, 'Final missing original field '+col)
                if pd.api.types.is_numeric_dtype(data[c]):
                    near(f[col], data[c], 'Final original statistic '+col)
                else:
                    same(f[col].fillna('NA'), data[c].fillna('NA'), 'Final original metadata '+col)
            axis_mod = 'scp_raw' if mod == 'scp' else mod
            same(booleans(f[mod+'_source_feature_present']), f.index.isin(self.x[axis_mod].index), mod+' final source availability')
            same(booleans(f[mod+'_passes_expression_filter']), f.index.isin(self.kept[axis_mod]), mod+' final test availability')
        s, l = self.s, self.l
        scp_deg = (s.FDR < .05)&(s.logFC.abs() >= .5)
        liu_deg = (l.FDR < .05)&(l.logFC.abs() >= .5)
        same_dir = s.logFC*l.logFC > 0
        scp_stable = (s.logFC*s.cb_logFC > 0)&s.loo_same_direction_fraction.eq(1)&s.loo_min_abs_logFC.ge(.25)
        liu_stable = (l.logFC*l.logFC_omit_Liu2023_CL0412 > 0)&(l.logFC*l.logFC_omit_Liu2023_CL2 > 0)
        shared = scp_deg & liu_deg & same_dir
        a = shared & scp_stable & liu_stable
        b = (scp_deg ^ liu_deg) & same_dir & scp_stable & liu_stable
        sc = self.cov[('SCP1303', 'HCM')]
        lc = self.cov[('liu2023', 'HCM')]
        technical = np.isfinite(self.token_ids)&sc.donors_token_ge10pct.ge(8)&lc.donors_token_ge10pct.ge(5)
        flags = {'scp_deg':scp_deg, 'liu_deg':liu_deg, 'same_direction':same_dir,
                 'scp_stable':scp_stable, 'liu_stable':liu_stable, 'shared_deg':shared,
                 'token_eligible':technical, 'eligible_A':a&technical, 'eligible_B':b&technical}
        for col, expected in flags.items():
            same(booleans(f[col]), expected, 'Final independent decision '+col)
        tier = np.full(len(f), 'none', object)
        tier[np.asarray(a)] = 'A_shared_strong'; tier[np.asarray(b)] = 'B_single_study_directional'
        same(f.tier, tier, 'Final evidence tier')
        near(f.rank_score_worst_FDR, np.maximum(s.FDR, l.FDR), 'Ranking score is max individual FDR')
        near(f.rank_score_min_abs_logFC, np.minimum(s.logFC.abs(), l.logFC.abs()), 'Ranking score is smaller effect')
        near(f.token_id, self.token_ids, 'Final token identity')
        prefixes = {'SCP1303':'scp', 'liu2023':'liu', 'chin_GSE174691':'chin_hcm',
                    'chin_GSE161921':'chin_control', 'GSE138262_whole_cell':'whole',
                    'GSE255296':'gse255296', 'GSE292754':'gse292754'}
        for (ds, disease), cov in self.cov.items():
            prefix = 'cov_'+prefixes[ds]+'_'+disease+'_'
            for c in ['donors','cells','native_cells','expr_fraction','token_fraction',
                      'donors_native','donors_token_ge10pct','median_donor_token_fraction']:
                near(f[prefix+c], cov[c], 'Final per-study coverage '+prefix+c)
        near(f.scp_hcm_donors_token_ge10pct, sc.donors_token_ge10pct, 'SCP named coverage gate')
        near(f.liu_hcm_donors_token_ge10pct, lc.donors_token_ge10pct, 'Liu named coverage gate')
        same(booleans(f.in_v1_vocab), np.isfinite(self.token_ids), 'Final vocabulary eligibility')
        old_rank = pd.Series(self.old['rank'].to_dict()).reindex(self.ids)
        near(f.first_round_rank, old_rank, 'Old ranks vs novel candidates')
        def relation(x,y):
            if not np.isfinite(x) or not np.isfinite(y): return 'not_evaluable'
            if x == 0 or y == 0: return 'zero'
            return 'same' if x*y > 0 else 'opposite'
        same(f.comparison_status, [relation(x,y) for x,y in zip(s.logFC,l.logFC)], 'Main comparison direction labels')
        for mod in ['scp','liu','chin','whole']:
            axis_mod = 'scp_raw' if mod == 'scp' else mod
            source_set, kept_set = set(self.x[axis_mod].index), set(self.kept[axis_mod])
            availability = ['absent_from_analysis_axis' if g not in source_set else
                            'below_expression_filter' if g not in kept_set else 'evaluable' for g in self.ids]
            same(f[mod+'_availability'], availability, mod+' detailed missingness')
        support = np.where(scp_deg & liu_deg, 'both', np.where(scp_deg,'SCP1303',np.where(liu_deg,'Liu2023','neither')))
        same(f.main_support, support, 'Main supporting cohort')
        for mod in ['chin','whole']:
            same(f[mod+'_direction'], [relation(x,y) for x,y in zip(self.stats[mod].descriptive_log2_cpm_difference,s.logFC)], mod+' descriptive direction labels')
        for field, labels in [('descriptive_same_count',['same']),('descriptive_opposite_count',['opposite']),
                              ('descriptive_evaluable_count',['same','opposite','zero'])]:
            near(f[field], f[['chin_direction','whole_direction']].isin(labels).sum(axis=1), field)
        for ds, cutoff in [('gse255296',2),('gse292754',1)]:
            same(booleans(f[ds+'_half_donors_token_covered']), f['cov_'+ds+'_HCM_donors_token_ge10pct'].ge(cutoff), ds+' half-donor coverage')
        expected_lists = {}
        for name, mask in [('A_BIOLOGICAL',a),('B_BIOLOGICAL',b),('A_SIMULATABLE',a&technical),('B_SIMULATABLE',b&technical)]:
            order = sorted(np.flatnonzero(mask), key=lambda i: (max(s.FDR.iloc[i],l.FDR.iloc[i]),
                           -min(abs(s.logFC.iloc[i]),abs(l.logFC.iloc[i])), self.ids[i]))
            ids = [self.ids[i] for i in order]
            expected_lists[name] = ids
        expected_lists['PRIORITY_TOP200'] = expected_lists['A_SIMULATABLE'][:200]
        for name, ids in expected_lists.items():
            actual = self.indexed('final/'+name+'.tsv')
            same(actual.index, ids, 'Independently ordered '+name)
            same(actual.candidate_rank, np.arange(1,len(ids)+1), name+' sequential rank')
            for c in actual:
                if c == 'candidate_rank': continue
                require(c in f, name+' contains unknown field '+c)
                expected = f.loc[ids, c]
                if pd.api.types.is_numeric_dtype(expected):
                    near(actual[c], expected, name+' value '+c)
                else:
                    same(actual[c].fillna('NA'), expected.fillna('NA'), name+' value '+c)
        comparison = self.indexed('final/PREVIOUS_200_COMPARISON.tsv')
        same(comparison.index, self.old.index, 'Previous 200 order')
        same(comparison.old_rank, self.old['rank'], 'Previous ranks preserved')
        for c in f:
            expected = f.loc[comparison.index, c]
            if pd.api.types.is_numeric_dtype(expected):
                near(comparison[c], expected, 'Previous comparison value '+c)
            else:
                same(comparison[c].fillna('NA'), expected.fillna('NA'), 'Previous comparison value '+c)
        expected_status = []
        for g in comparison.index:
            if bool(a.loc[g]): expected_status.append('A_simulatable' if bool(technical.loc[g]) else 'A_low_coverage')
            elif bool(b.loc[g]): expected_status.append('B_simulatable' if bool(technical.loc[g]) else 'B_low_coverage')
            else: expected_status.append('not_selected')
        same(comparison.new_status, expected_status, 'Previous candidate disposition')
        for label, col in [('A_SIMULATABLE','A_simulatable_rank'),('B_SIMULATABLE','B_simulatable_rank')]:
            expected = pd.Series(dict(zip(expected_lists[label],range(1,len(expected_lists[label])+1)))).reindex(self.ids)
            near(f[col], expected, 'Full table '+col)
        same(booleans(f.priority_top200), f.index.isin(expected_lists['PRIORITY_TOP200']), 'Priority membership flag')
        categories = []
        for i in range(len(f)):
            if a.iloc[i]: category = 'A_shared_strong'
            elif b.iloc[i]: category = 'B_single_study_directional'
            elif f.comparison_status.iloc[i] != 'same': category = f.comparison_status.iloc[i]
            elif shared.iloc[i]: category = 'shared_direction_sensitive'
            elif scp_deg.iloc[i] != liu_deg.iloc[i]: category = 'single_study_direction_sensitive'
            else: category = 'same_neither_DEG'
            categories.append(category)
        same(f.evidence_category, categories, 'Complete evidence categories')
        summary = self.json('final/SUMMARY.json')
        require(summary['held_out_validation'] is False and summary['joint_FDR_claim'] is False and
                summary['selection_scope'] == 'full_gene_union_not_previous_200', 'Final interpretation boundaries')
        require(summary['gene_universe'] == len(f) and summary['shared_DEG'] == int(shared.sum()) and
                summary['shared_stability_failed'] == int((shared & ~a).sum()), 'Shared-gene summary')
        require(summary['table_rows'] == {k:len(v) for k,v in expected_lists.items()}, 'Candidate table counts')
        require(summary['prior_200_sha256'] == OLD_LOCK_SHA and
                summary['prior_200_status'] == pd.Series(expected_status).value_counts().to_dict(), 'Previous cohort summary')
        require(summary['comparison_counts'] == f.comparison_status.value_counts().to_dict() and
                summary['evidence_category_counts'] == pd.Series(categories).value_counts().to_dict(), 'All-gene evidence counts')
        require(summary['both_evaluable'] == int((np.isfinite(s.logFC)&np.isfinite(l.logFC)).sum()), 'Both-study estimability')
        for mod in ['scp','liu','chin','whole']:
            axis_mod = 'scp_raw' if mod == 'scp' else mod
            metrics = summary['studies'][mod]
            require(metrics['source_genes'] == len(self.x[axis_mod]) and
                    metrics['expression_filtered'] == len(self.kept[axis_mod]), mod+' summary gene counts')
            if mod in ['scp','liu']:
                require(metrics['DEG'] == int(flags[mod+'_deg'].sum()), mod+' summary DEG count')
        for name, ids in expected_lists.items():
            group = summary['candidate_groups'][name]
            t = f.loc[ids]
            expected = {'higher_HCM': int(t.scp_logFC.gt(0).sum()), 'lower_HCM':int(t.scp_logFC.lt(0).sum()),
                        'in_previous_200':int(t.first_round_rank.notna().sum()),
                        'new_vs_previous_200':int(t.first_round_rank.isna().sum()),
                        'main_support':t.main_support.value_counts().to_dict(),
                        'chin_directions':t.chin_direction.value_counts().to_dict(),
                        'whole_directions':t.whole_direction.value_counts().to_dict(),
                        'gse255296_half_donors_covered':int(t.gse255296_half_donors_token_covered.sum()),
                        'gse292754_half_donors_covered':int(t.gse292754_half_donors_token_covered.sum())}
            require(group == expected, name+' complete summary reconciliation')
        for path, digest in self.json('final/OUTPUT_HASHES.json').items():
            require(sha(self.w/'final'/path) == digest, 'Final output seal '+path)
        self.metrics['predicted_counts'] = {name:len(ids) for name,ids in expected_lists.items()}
        self.metrics['shared_DEG_before_sensitivity'] = int(shared.sum())
        self.metrics['priority_previous_overlap'] = len(set(expected_lists['PRIORITY_TOP200'])&set(self.old.index))
        self.metrics['final_evidence_sha256'] = sha(self.w/'final/ALL_GENE_EVIDENCE.tsv.gz')
        self.metrics['candidate_sha256'] = {name:sha(self.w/'final'/ (name+'.tsv')) for name in expected_lists}
        self.expected_lists = expected_lists
        self.f = f


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('work', type=Path)
    args = p.parse_args()
    a = Audit(args.work.resolve())
    out = a.w/'audit'
    out.mkdir(exist_ok=True)
    error = None
    try:
        for label, function in [('preserved source hashes and exploratory design',a.provenance),
                                ('donor counts, gene universe and expression filters',a.counts),
                                ('whole-genome BH, 31 omissions and original Liu values',a.fits),
                                ('seven-study coverage from all donor detail',a.coverage),
                                ('full evidence joins, candidate rules and complete rankings',a.final)]:
            a.run(label,function)
    except Exception as exc:
        error = a.phase+': '+type(exc).__name__+': '+str(exc)
    report = {'status':'failed' if error else 'passed', 'mode':'independent reconciliation, not a second model fit',
              'completed_checks':a.done, 'metrics':a.metrics, 'failure':error,
              'utc':datetime.now(timezone.utc).isoformat()}
    dest = out/'FINAL_AUDIT.json'
    temp = out/'FINAL_AUDIT.partial'
    temp.write_text(json.dumps(report,indent=2,allow_nan=False)+'\n',encoding='utf-8'); temp.replace(dest)
    print(json.dumps(report,allow_nan=False))
    if error: raise SystemExit(1)


if __name__ == '__main__':
    main()
