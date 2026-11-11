"""Independently reopen delivered objects and compare every count with original GEO matrices."""
from datetime import datetime, timezone
from pathlib import Path
import gzip
import hashlib
import json

import anndata as ad
import numpy as np
import pandas as pd
from scipy.io import mmread
from threadpoolctl import threadpool_limits

ROOT = Path.home() / 'xiaoe-internship-data'
WORK = ROOT / '_processing/all_data_cleaning_20261111/chin_pair'
samples = pd.read_csv(WORK / 'samples.tsv', sep='\t')
official = json.loads((ROOT / 'processed/all_data_cleaning_20261111/GEO_SAMPLES.json').read_text())


def check(checks, name, result):
    checks[name] = bool(result)


def exact(left, right):
    return left.shape == right.shape and (left != right).nnz == 0


summary = []
for dataset, group in samples.groupby('dataset', sort=True):
    output = ROOT / f'processed/{dataset}_qc_20261111'
    checks = {}
    exports = json.loads((output / 'EXPORT_SUCCESS.json').read_text())
    policy = json.loads((output / 'qc_policy.json').read_text())
    names = ['all_cells_qc.h5ad', 'cm_qc_full.h5ad', 'cm_qc_light.h5ad',
             'cm_review_pending.h5ad', 'suspected_doublets.h5ad']
    objects = {name: ad.read_h5ad(output / name) for name in names}
    allcell = objects['all_cells_qc.h5ad']
    dbl = objects['suspected_doublets.h5ad']
    core = objects['cm_qc_full.h5ad']
    pending = objects['cm_review_pending.h5ad']
    light = objects['cm_qc_light.h5ad']
    identity = pd.read_csv(output / 'cell_identity_audit.tsv.gz', sep='\t', index_col='barcode', low_memory=False)
    check(checks, 'policy_strict_limits', policy['min_genes'] == 500 and policy['max_genes'] == 5000 and policy['mt_percent_exclusive_max'] == 5)
    check(checks, 'core_all_identity_retained', core.obs.review_action.eq('retain_CM').all())
    check(checks, 'core_best_and_pruned_reference_CM', core.obs.reference_label.eq('Cardiomyocyte').all() and core.obs.reference_pruned.eq('Cardiomyocyte').all())
    check(checks, 'core_multigene_and_program_support', core.obs.Cardiomyocyte_markers.ge(2).all() and core.obs.cm_program_dominates.all())
    check(checks, 'core_reference_tail_or_context_support', (core.obs.reference_tail_pass | (core.obs.Cardiomyocyte_markers.ge(3) & core.obs.cm_anchor_neighbor_fraction.ge(.5))).all())
    check(checks, 'pending_explicit_and_disjoint_core', pending.obs.review_action.eq('pending_CM').all() and not pending.obs_names.isin(core.obs_names).any())
    check(checks, 'doublet_exclusions_have_two_evidence_sources', dbl.obs.scDblFinder_class.eq('doublet').all() and dbl.obs.coherent_mixed_panel_flag.all())
    check(checks, 'allcell_no_excluded_doublets', not allcell.obs.exclude_suspected_doublet.any())
    check(checks, 'core_matches_audit_membership', set(core.obs_names) == set(identity.index[identity.review_action.eq('retain_CM')]))
    check(checks, 'pending_matches_audit_membership', set(pending.obs_names) == set(identity.index[identity.review_action.eq('pending_CM')]))
    check(checks, 'all_strict_audit_members_exported', set(identity.index) == set(allcell.obs_names) | set(dbl.obs_names))
    check(checks, 'doublets_disjoint_allcell', not dbl.obs_names.isin(allcell.obs_names).any())
    diagnostic = ad.read_h5ad(output / 'diagnostic_normalized.h5ad')
    check(checks, 'diagnostic_explicitly_analysis_only', diagnostic.uns.get('analysis_only', False) and 'Diagnostic log1p CP10k' in diagnostic.uns.get('count_source', ''))
    check(checks, 'diagnostic_full_strict_axes', set(diagnostic.obs_names) == set(identity.index) and diagnostic.var_names.equals(allcell.var_names))
    original_diagnostic = ad.concat([allcell, dbl], join='inner', merge='same')[diagnostic.obs_names].X.astype(float).tocsr()
    size = np.asarray(original_diagnostic.sum(axis=1)).ravel()
    expected_diagnostic = original_diagnostic.multiply((10000 / np.maximum(size, 1))[:, None]).tocsr()
    expected_diagnostic.data = np.log1p(expected_diagnostic.data)
    diagnostic_difference = expected_diagnostic - diagnostic.X
    check(checks, 'diagnostic_values_match_declared_log_CP10k', np.max(np.abs(diagnostic_difference.data), initial=0) < 1e-5)
    check(checks, 'light_unique_no_padding', light.obs_names.is_unique and light.n_obs == min(20000, core.n_obs))
    check(checks, 'light_subset_core', light.obs_names.isin(core.obs_names).all())
    for name in ['cm_qc_full.h5ad', 'cm_qc_light.h5ad', 'cm_review_pending.h5ad']:
        obj = objects[name]
        check(checks, name + '_matches_allcell_X', exact(obj.X, allcell[obj.obs_names].X))
    for name, obj in objects.items():
        check(checks, name + '_unique_axes', obj.obs_names.is_unique and obj.var_names.is_unique)
        check(checks, name + '_strict_bounds', obj.obs.n_genes_by_counts.between(500,5000).all() and obj.obs.pct_counts_mt.lt(5).all())
        check(checks, name + '_called_cells_only', obj.obs.emptydrops_called.all() and obj.obs.emptydrops_FDR.le(.001).all())
        check(checks, name + '_finite_integer_nonnegative', np.isfinite(obj.X.data).all() and (obj.X.data >= 0).all() and (obj.X.data == np.round(obj.X.data)).all())
        check(checks, name + '_X_equals_counts', exact(obj.X, obj.layers['counts']))
        check(checks, name + '_no_training_split', 'split' not in obj.obs)
    strict_expected = set()
    for row in group.itertuples():
        sample = row.sample_id
        stage = json.loads((WORK / sample / 'STAGED.json').read_text())
        for feature, record in stage['sources'].items():
            source = Path(record['path'])
            with source.open('rb') as handle:
                sha = hashlib.file_digest(handle, 'sha256').hexdigest()
            check(checks, sample + '_' + feature + '_source_unchanged', source.stat().st_size == record['bytes'] and sha == record['sha256'])
        prefix = row.source_prefix
        var = pd.read_csv(prefix + '_features.tsv.gz', sep='\t', header=None, dtype=str, keep_default_na=False)
        barcodes = pd.read_csv(prefix + '_barcodes.tsv.gz', sep='\t', header=None, dtype=str)[0]
        with gzip.open(prefix + '_matrix.mtx.gz', 'rb') as handle, threadpool_limits(limits=2):
            raw = mmread(handle).T.tocsr()
        raw.sum_duplicates(); raw.eliminate_zeros()
        counts = np.asarray(raw.sum(axis=1)).ravel()
        genes = np.diff(raw.indptr)
        mt_counts = np.asarray(raw[:, var[1].str.startswith('MT-').to_numpy()].sum(axis=1)).ravel()
        mt = 100 * mt_counts / np.maximum(counts, 1)
        ed = pd.read_csv(WORK / sample / 'emptydrops.tsv.gz', sep='\t', index_col='barcode',
                         usecols=['barcode', 'Total', 'FDR', 'emptydrops_called'])
        prefixed = pd.Index(sample + ':' + barcodes)
        check(checks, sample + '_emptydrops_original_axes_and_totals', ed.index.equals(prefixed) and np.array_equal(ed.Total, counts))
        called = ed.FDR.le(.001).to_numpy()
        check(checks, sample + '_cell_calling_matches_FDR', np.array_equal(called, ed.emptydrops_called))
        strict = called & (genes >= 500) & (genes <= 5000) & (mt < 5)
        strict_expected.update(prefixed[strict])
        model = pd.read_csv(WORK / sample / 'doublets_v2.tsv', sep='\t')
        check(checks, sample + '_doublet_model_complete_ge500_capture', set(model.barcode) == set(prefixed[called & (genes >= 500)]) and model.barcode.is_unique)
        check(checks, sample + '_doublet_scores_finite', np.isfinite(model.scDblFinder_score).all() and model.scDblFinder_class.isin(['singlet', 'doublet']).all())
        audit = pd.read_csv(output / f'{sample}_droplet_audit.tsv.gz', sep='\t', index_col='barcode')
        check(checks, sample + '_complete_droplet_audit', audit.index.equals(prefixed) and np.array_equal(audit.total_counts, counts) and np.array_equal(audit.n_genes_by_counts, genes) and np.allclose(audit.pct_counts_mt, mt))
        original_record = [r for r in official if r['accession'] == sample]
        check(checks, sample + '_unique_official_record', len(original_record) == 1)
        original_donor = original_record[0]['title'][0].split(':')[0].strip()
        check(checks, sample + '_donor_matches_official', str(row.donor_original) == original_donor)
        expected_disease = 'Control' if row.series == 'GSE161921' else 'HCM'
        expected_label = 0 if expected_disease == 'Control' else 1
        for name, obj in objects.items():
            part = obj[obj.obs.sample_id.eq(sample)]
            positions = barcodes.to_numpy()
            positions = pd.Index(positions).get_indexer(part.obs.barcode_original.astype(str))
            stem = sample + '_' + name
            check(checks, stem + '_all_genes_original_order', np.array_equal(part.var_names.astype(str), var[0]) and np.array_equal(part.var.gene_symbol.astype(str), var[1]))
            check(checks, stem + '_all_barcodes_original', (positions >= 0).all() and np.array_equal(part.obs_names.astype(str), sample + ':' + part.obs.barcode_original.astype(str)))
            check(checks, stem + '_all_counts_original', exact(part.X, raw[positions]))
            check(checks, stem + '_recomputed_QC', np.array_equal(part.obs.total_counts, counts[positions]) and np.array_equal(part.obs.n_genes_by_counts, genes[positions]) and np.allclose(part.obs.pct_counts_mt, mt[positions]))
            check(checks, stem + '_source_metadata', part.obs.donor_id.eq('Chin_' + original_donor).all() and part.obs.disease.eq(expected_disease).all() and part.obs.label.eq(expected_label).all() and part.obs.capture_id.eq(sample).all())
        print('Verified original count values: ' + sample, flush=True)
    check(checks, 'every_original_strict_candidate_accounted_for', strict_expected == set(allcell.obs_names) | set(dbl.obs_names))
    check(checks, 'study_not_merged', all(obj.obs.dataset.eq(dataset).all() for obj in objects.values()))
    check(checks, 'expected_original_sample_count', len(group) == (4 if dataset == 'chin_GSE161921' else 9))
    hashes = {}
    for name, record in exports['hashes'].items():
        path = output / name
        md5 = hashlib.md5(); sha = hashlib.sha256()
        with path.open('rb') as handle:
            while chunk := handle.read(32 * 1024**2):
                md5.update(chunk); sha.update(chunk)
        hashes[name] = {'bytes': path.stat().st_size, 'md5': md5.hexdigest(), 'sha256': sha.hexdigest()}
        check(checks, name + '_size_and_export_hash', hashes[name]['bytes'] == record['bytes'] and hashes[name]['sha256'] == record['sha256'])
    status = 'passed' if all(checks.values()) else 'failed'
    result = {'status': status, 'dataset': dataset, 'verified_utc': datetime.now(timezone.utc).isoformat(),
              'checks': checks, 'hashes': hashes, 'all_cells_qc': allcell.n_obs, 'core_CM': core.n_obs,
              'light_CM': light.n_obs, 'pending_CM': pending.n_obs, 'suspected_doublets': dbl.n_obs,
              'input_donors': len(group), 'core_donors': int(core.obs.donor_id.nunique()), 'genes': core.n_vars,
              'verifier_sha256': hashlib.sha256(Path(__file__).read_bytes()).hexdigest()}
    (output / 'INDEPENDENT_VERIFICATION.json').write_text(json.dumps(result, indent=2) + '\n')
    if status != 'passed':
        raise AssertionError([k for k, v in checks.items() if not v])
    for algorithm in ['md5', 'sha256']:
        (output / ('MD5SUMS' if algorithm == 'md5' else 'SHA256SUMS')).write_text(''.join(f"{r[algorithm]}  {n}\n" for n, r in hashes.items()))
    coverage = pd.read_csv(output / 'sample_qc_summary.tsv', sep='\t')
    donor_warning = '; '.join(coverage.loc[coverage.core_CM.eq(0), 'donor_id']) or '无'
    report = f'''# {dataset} 清洗与验收完成

验收时间：{result['verified_utc']}；{len(checks)}项检查全部通过。服务器：`{output}/`。

| 对象 | 核数 |
| --- | ---: |
| 全细胞QC（扣除联合证据疑似双细胞） | {allcell.n_obs} |
| 心肌核心 | {core.n_obs} |
| 心肌轻量 | {light.n_obs} |
| 心肌身份待定，单独保留 | {pending.n_obs} |
| 多证据疑似双细胞，单独保留 | {dbl.n_obs} |

全部{core.n_vars}基因及未归一化计数保留；核心覆盖{result['core_donors']}/{len(group)}位原供体。核心无细胞的供体：{donor_warning}。不复制细胞凑20k。

原始完整液滴矩阵先经EmptyDrops FDR≤0.001细胞识别，再按500≤检测基因≤5000、MT<5%执行项目严格QC。原始矩阵及完整逐液滴审计可追溯。每个GEO样本单独建scDblFinder模型，包含各细胞类型；采用10x默认加载率假设，而非SCP1303过滤后的残余率。工具调用加互斥多基因程序共同支持时才排除疑似双细胞；不能保证识别所有同型双细胞。

心肌身份由成人左心室单核参考SingleR、多标记程序及必要的邻域支持联合判断，参考分数不当作概率。参考训练和留出供体不混用；参考对HCM状态的适配仍有限，冲突者保存为待定。各项操作阈值见qc_policy.json；不把项目阈值称为HCM行业统一标准。

sample_qc_summary.tsv报告每供体细胞识别、基因/MT筛选、心肌保留、待定及双细胞影响；该严格版本可能损失低复杂度或高MT核。原始液滴数据未做去背景，稀疏背景仍可能存在；交付计数未插补、归一化、整合或HVG裁切。diagnostic_normalized.h5ad明确为分析副本，不能当原始计数输入。

独立验收重新读取每个原始GEO矩阵，对全部交付对象逐值核对原始计数、全基因及其顺序、条码、QC与疾病/供体元数据；所有严格合格细胞均在全细胞或疑似双细胞对象中有记录。重新计算文件MD5/SHA256，并确认输入文件摘要不变。

疾病标签来自GEO系列设计：GSE161921为未使用供体心脏，Control=0；GSE174691为9位HCM切肌标本，HCM=1。标签是来源供体诊断，不是逐核病理确证；不作为与其他研究生物供体完全独立的证明。未指定训练/测试划分，不进行模型训练、tokenization或疾病差异分析。
'''
    (output / 'REPORT.md').write_text(report)
    summary.append({k: v for k, v in result.items() if k not in ['checks', 'hashes']})
    print(f'VERIFIED {dataset}: {len(checks)} checks, core={core.n_obs}', flush=True)
(WORK / 'VERIFICATION_SUMMARY.json').write_text(json.dumps(summary, indent=2) + '\n')
