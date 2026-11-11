"""Independently verify every bulk value, source mapping and explicit source limitation."""
from collections import Counter
from datetime import datetime, timezone
from pathlib import Path
import csv
import gzip
import hashlib
import json
import re
import traceback

import anndata as ad
import numpy as np
import pandas as pd

ROOT = Path.home() / 'xiaoe-internship-data'
SOURCE = ROOT / 'validation/bulk_GSE249925'
WORK = ROOT / '_processing/all_data_cleaning_20261111/GSE249925'
OUT = ROOT / 'processed/bulk_GSE249925_qc_20261111'
CHECKS = []


def check(name, passed):
    CHECKS.append({'name': name, 'passed': bool(passed)})
    if not passed:
        raise AssertionError(name)


def digest(path, algorithm='sha256'):
    with path.open('rb') as handle:
        return hashlib.file_digest(handle, algorithm).hexdigest()


def run():
    baseline = json.loads((WORK / 'ORIGINAL_DOWNLOAD_BASELINE.json').read_text())
    check('three_original_payloads_accounted_for', len(baseline) == 3)
    original_checksums = {}
    for row in baseline:
        path = Path(row['server_path'])
        check('original_bytes_' + path.name, path.stat().st_size == int(row['bytes']))
        check('original_MD5_' + path.name, digest(path, 'md5') == row['md5'])
        measured_sha256 = digest(path)
        historical_sha256 = row.get('sha256', '')
        if historical_sha256:
            check('original_historical_SHA256_' + path.name, measured_sha256 == historical_sha256)
        original_checksums[path.name] = {'bytes': path.stat().st_size, 'MD5': row['md5'],
            'SHA256': measured_sha256, 'historical_SHA256_available': bool(historical_sha256),
            'historical_size_and_MD5_matched': True}
    # Use the original CSV stream, not the exported AnnData or preparation tables.
    with gzip.open(SOURCE / 'GSE249925_20231208-For-GEO-Raw-Counts.tsv.gz', 'rt') as handle:
        reader = csv.reader(handle, delimiter='\t')
        header = next(reader)
        rows = list(reader)
    check('native_header_and_complete_row_widths', header[:2] == ['Gene_Index', 'Length'] and len(header) == 122 and all(len(row) == len(header) for row in rows))
    check('no_duplicate_native_column_headers', len(set(header)) == len(header))
    native_ids = pd.Index([row[0] for row in rows], name='source_feature_id')
    lengths = [row[1] for row in rows]
    x = np.array([list(map(int, row[2:])) for row in rows], dtype=np.int64)
    samples = pd.Index(header[2:], name='sample_id')
    check('source_matrix_is_real_nonnegative_integer_release', x.shape == (64253, 120) and np.isfinite(x).all() and (x >= 0).all())
    check('native_feature_ids_are_unique', native_ids.is_unique)
    full = ad.read_h5ad(OUT / 'all_published_expression.h5ad')
    core = ad.read_h5ad(OUT / 'bulk_qc_expression.h5ad')
    check('full_axes_exactly_preserved', full.var_names.equals(native_ids) and full.obs_names.equals(samples))
    check('every_full_published_value_matches_original', np.array_equal(full.X.toarray(), x.T) and np.array_equal(full.layers['published_values'].toarray(), x.T))
    check('source_length_and_composite_ID_strings_preserved', full.var.source_Gene_Index.tolist() == native_ids.tolist() and full.var.source_Length.astype(str).tolist() == lengths)
    check('original_feature_row_order_preserved', np.array_equal(full.var.source_row, np.arange(len(x))))
    pieces = [identifier.split('@') for identifier in native_ids]
    check('every_native_composite_feature_has_five_parts', all(len(parts) == 5 for parts in pieces))
    for position, name in enumerate(['chromosome', 'gene_id', 'gene_symbol', 'start_original', 'end_original']):
        check('native_component_' + name, full.var[name].astype(str).tolist() == [parts[position] for parts in pieces])
    valid_rows = []
    for parts, length in zip(pieces, lengths):
        try:
            valid_rows.append(bool(re.fullmatch(r'ENSG\d{11}(?:\.\d+)?', parts[1])) and bool(parts[0]) and bool(parts[2])
                              and int(parts[3]) >= 1 and int(parts[4]) >= int(parts[3]) and int(length) > 0)
        except (ValueError, TypeError):
            valid_rows.append(False)
    valid = np.asarray(valid_rows, dtype=bool)
    totals = x.sum(axis=0)
    detected = np.count_nonzero(x, axis=0)
    zero = np.count_nonzero(x, axis=1) == 0
    feature_mask = valid & ~zero
    sample_mask = (totals > 0) & (detected > 0)
    check('feature_identity_and_zero_QC_recomputed', np.array_equal(full.var.identity_valid, valid) and np.array_equal(full.var.all_zero, zero) and np.array_equal(full.var.qc_usable, feature_mask))
    check('sample_information_QC_recomputed', np.array_equal(full.obs.qc_usable, sample_mask))
    check('every_feature_QC_total_and_detection_matches', np.array_equal(full.var.published_expression_total, x.sum(axis=1)) and np.array_equal(full.var.n_samples_detected, np.count_nonzero(x, axis=1)))
    check('every_sample_QC_total_and_detection_matches', np.array_equal(full.obs.published_expression_total, totals) and np.array_equal(full.obs.n_features_detected, detected))
    mito = np.array([parts[0] in ['chrM', 'chrMT', 'MT', 'M'] for parts in pieces])
    check('bulk_MT_is_a_reported_metric_only', np.array_equal(full.var.mt_by_source_chromosome, mito) and np.allclose(full.obs.pct_published_expression_mt, 100 * x[mito].sum(axis=0) / np.maximum(totals, 1), rtol=0, atol=1e-10) and full.obs.qc_usable.all())
    expected_reasons = np.where(~valid, 'invalid_feature_identity', np.where(zero, 'all_samples_zero', 'retained'))
    check('every_feature_exclusion_reason_recomputed', np.array_equal(full.var.qc_reason, expected_reasons))
    gene_ids = pd.Series([parts[1] for parts in pieces])
    symbols = pd.Series([parts[2] for parts in pieces])
    check('duplicate_native_gene_IDs_and_symbols_flagged_without_aggregation', np.array_equal(full.var.gene_id_nonunique, gene_ids.duplicated(keep=False)) and np.array_equal(full.var.gene_symbol_nonunique, symbols.duplicated(keep=False)))
    check('QC_axes_exact_retained_subset', core.var_names.equals(native_ids[feature_mask]) and core.obs_names.equals(samples[sample_mask]))
    expected = x[feature_mask][:, sample_mask].T
    check('every_QC_value_matches_original_without_normalization', np.array_equal(core.X.toarray(), expected) and np.array_equal(core.layers['published_values'].toarray(), expected))
    check('QC_metadata_is_exact_subset', core.var.astype(str).equals(full.var.loc[core.var_names].astype(str)) and core.obs.astype(str).equals(full.obs.loc[core.obs_names].astype(str)))
    qc_tsv = pd.read_csv(OUT / 'bulk_qc_matrix.tsv.gz', sep='\t', dtype=str, keep_default_na=False)
    check('QC_TSV_native_order_and_annotations', qc_tsv.columns.tolist() == header and qc_tsv.Gene_Index.tolist() == native_ids[feature_mask].tolist() and qc_tsv.Length.tolist() == np.array(lengths)[feature_mask].tolist())
    check('every_QC_TSV_value_matches_original', np.array_equal(qc_tsv.iloc[:, 2:].to_numpy(dtype=np.int64), expected.T))
    check('all_zero_removals_really_have_no_information', zero.sum() == 15624 and x[zero].sum() == 0)
    feature_audit = pd.read_csv(OUT / 'feature_audit.tsv.gz', sep='\t', index_col=0, keep_default_na=False)
    feature_exclusions = pd.read_csv(OUT / 'feature_exclusions.tsv.gz', sep='\t', index_col=0, keep_default_na=False)
    sample_exclusions = pd.read_csv(OUT / 'sample_exclusions.tsv', sep='\t', index_col=0, keep_default_na=False)
    check('complete_gene_audit_and_exclusion_indices', feature_audit.index.equals(native_ids) and feature_exclusions.index.equals(native_ids[~feature_mask]) and np.array_equal(feature_audit.qc_usable, feature_mask))
    check('sample_exclusions_are_real_empty_result', len(sample_exclusions) == 0 and sample_mask.all())
    # Parse source SOFT independently by complete SAMPLE blocks.
    source_records = {}
    for block in (SOURCE / 'sample_metadata.soft').read_text().split('^SAMPLE = ')[1:]:
        accession, *lines = block.splitlines()
        fields = {}
        for line in lines:
            if line.startswith('!Sample_') and ' = ' in line:
                field, value = line.split(' = ', 1)
                fields.setdefault(field, []).append(value)
        source_records[accession] = fields
    reference_geo = json.loads((WORK / 'GEO_SAMPLES_GSE249925.json').read_text())
    reference_geo = {row['accession']: row for row in reference_geo}
    check('all120_original_and_previously_retrieved_GEO_records_accounted_for', len(source_records) == len(reference_geo) == 120 and set(full.obs.GSM) == set(source_records) == set(reference_geo))
    source_processing = set()
    for sample, row in full.obs.iterrows():
        record = source_records[row.GSM]
        characteristics = dict(value.split(': ', 1) for value in record['!Sample_characteristics_ch1'])
        relation = record['!Sample_relation']
        biosample = next(value.rsplit('/', 1)[1] for value in relation if value.startswith('BioSample: '))
        sra = next(value.split('term=', 1)[1] for value in relation if value.startswith('SRA: '))
        check(row.GSM + '_exact_source_title_and_registry_mapping', record['!Sample_title'] == [sample] and row.source_title == sample and row.BioSample == biosample and row.SRA_experiment == sra)
        check(row.GSM + '_disease_and_tissue_from_GEO', row.disease == characteristics['group'] and row.label == {'HCM': 1, 'Control': 0}[characteristics['group']] and row.tissue_original == characteristics['tissue'])
        check(row.GSM + '_second_source_GEO_crosscheck', reference_geo[row.GSM]['title'] == [sample] and set(reference_geo[row.GSM]['characteristics_ch1']) == set(record['!Sample_characteristics_ch1']))
        source_processing.update(record['!Sample_data_processing'])
    check('only_true_source_groups_labeled', full.obs.disease.value_counts().to_dict() == {'HCM': 97, 'Control': 23} and full.obs.label.isin([0, 1]).all())
    check('donor_age_sex_batch_not_invented', full.obs.donor_id.astype(str).eq('').all() and full.obs.donor_mapping_status.eq('not_provided_per_sample').all() and not full.obs[['age_available', 'sex_available', 'batch_available']].any().any())
    check('all_source_GEO_specimens_are_human_bulk_RNA', full.obs.modality.eq('bulk_RNA').all() and full.obs.organism.eq('Homo sapiens').all() and not full.obs.training_eligible_single_cell.any())
    check('numeric_sample_titles_and_R_suffixes_preserved', samples.tolist() == header[2:] and all(isinstance(value, str) for value in full.obs_names))
    check('registry_duplicates_recomputed', np.array_equal(full.obs.biosample_nonunique, full.obs.BioSample.duplicated(keep=False)) and np.array_equal(full.obs.sra_experiment_nonunique, full.obs.SRA_experiment.duplicated(keep=False)))
    expression_hashes = [hashlib.sha256(np.ascontiguousarray(x[:, i]).tobytes()).hexdigest() for i in range(len(samples))]
    frequency = Counter(expression_hashes)
    check('exact_expression_duplicates_recomputed_not_used_as_donor_map', full.obs.expression_vector_sha256.tolist() == expression_hashes and np.array_equal(full.obs.exact_expression_duplicate, [frequency[value] > 1 for value in expression_hashes]))
    for values, column in [(totals, 'log_total_5MAD_flag_only'), (detected, 'log_detected_5MAD_flag_only')]:
        transformed = np.log1p(values)
        center = np.median(transformed)
        mad = np.median(np.abs(transformed - center))
        flagged = np.abs(transformed - center) > 5 * 1.4826 * mad if mad > 0 else np.zeros(len(samples), dtype=bool)
        check('review_only_' + column, np.array_equal(full.obs[column], flagged))
    policy = json.loads((OUT / 'qc_policy.json').read_text())
    check('source_conflict_is_real_and_preserved', any('Timmed means of M-values' in text for text in source_processing) and any('Raw gene-count matrix' in text for text in source_processing) and set(policy['published_processing']) == source_processing)
    check('no_false_raw_count_certification', policy['certified_unnormalized_counts'] is False and not full.obs.certified_unnormalized_counts.any() and 'counts' not in full.layers and 'counts' not in core.layers and 'TMM' in full.uns['X_semantics'])
    check('bulk_specific_policy_and_no_silent_transformation', set(policy['not_applied']) == {'single-cell500to5000', 'MT_lt5pct', 'EmptyDrops', 'doublet_calling', 'CM_extraction', 'normalization', 'imputation', 'batch_correction', 'differential_expression', 'tokenization', 'training'})
    check('training_split_not_assigned_without_subject_mapping', full.obs.split_status.eq('unassigned_subject_and_cross_study_overlap_review_required').all())
    sample_audit = pd.read_csv(OUT / 'sample_metadata_QC.tsv', sep='\t', index_col=0, keep_default_na=False)
    check('sample_audit_native_string_order', sample_audit.index.astype(str).tolist() == samples.tolist())
    check('sample_audit_disease_and_donor_limits', np.array_equal(sample_audit.disease, full.obs.disease) and np.array_equal(sample_audit.label, full.obs.label) and sample_audit.donor_id.eq('').all() and not sample_audit.certified_unnormalized_counts.any())
    grouped = pd.read_csv(OUT / 'group_QC_summary.tsv', sep='\t', index_col=0)
    for disease in ['HCM', 'Control']:
        part = full.obs.loc[full.obs.disease.eq(disease)]
        row = grouped.loc[disease]
        check('group_summary_' + disease, row.samples == len(part) and row.usable == part.qc_usable.sum() and row.total_expression_min == part.published_expression_total.min() and row.total_expression_max == part.published_expression_total.max() and row.detected_min == part.n_features_detected.min() and row.detected_max == part.n_features_detected.max())
    export = json.loads((OUT / 'EXPORT_SUCCESS.json').read_text())
    check('export_shape_summary_matches_real_outputs', export['source_features'] == full.n_vars and export['QC_features'] == core.n_vars and export['source_samples'] == full.n_obs and export['QC_samples'] == core.n_obs)
    manifest = []
    for name, row in export['hashes'].items():
        path = OUT / name
        check(name + '_bytes_MD5_SHA256', path.stat().st_size == row['bytes'] and digest(path, 'md5') == row['MD5'] and digest(path) == row['SHA256'])
        manifest.append({'dataset': 'bulk_GSE249925', 'file': name, 'server_path': str(path), **row,
                         'use': 'QC_of_published_expression_raw_vs_TMM_unresolved'})
    check('all_nine_required_deliverables_present', set(export['hashes']) == {'all_published_expression.h5ad', 'bulk_qc_expression.h5ad', 'bulk_qc_matrix.tsv.gz', 'sample_metadata_QC.tsv', 'feature_audit.tsv.gz', 'feature_exclusions.tsv.gz', 'sample_exclusions.tsv', 'group_QC_summary.tsv', 'qc_policy.json'})
    code_hashes = {}
    for line in (WORK / 'CODE_EXPECTED_SHA256SUMS.txt').read_text().splitlines():
        expected_hash, name = line.split(maxsplit=1)
        check('executed_code_SHA256_' + name, digest(WORK / name) == expected_hash)
        code_hashes[name] = expected_hash
    check('prepare_process_exited_successfully', (WORK / 'prepare_exit_code.txt').read_text().strip() == '0')
    completed = datetime.now(timezone.utc).isoformat()
    summary = {'dataset': 'bulk_GSE249925', 'status': 'complete_verified', 'checked_utc': completed,
               'checks_passed': len(CHECKS), 'source_features': full.n_vars, 'QC_features': core.n_vars,
               'source_samples': full.n_obs, 'QC_samples': core.n_obs, 'HCM_samples': 97, 'Control_samples': 23,
               'zero_features_removed_from_QC_only': int(zero.sum()), 'invalid_feature_identity': int((~valid).sum()),
               'duplicate_Ensembl_rows': int(gene_ids.duplicated(keep=False).sum()),
               'duplicate_symbol_rows': int(symbols.duplicated(keep=False).sum()),
               'MAD_total_flags': int(full.obs.log_total_5MAD_flag_only.sum()), 'MAD_detected_flags': int(full.obs.log_detected_5MAD_flag_only.sum()),
               'certified_unnormalized_counts': False, 'single_cell_training_eligible': False,
               'original_input_checksums': original_checksums,
               'all_checks': CHECKS, 'verified_code_sha256': code_hashes}
    (OUT / 'INDEPENDENT_VERIFICATION.json').write_text(json.dumps(summary, indent=2) + '\n')
    report = f'''# GSE249925：bulk 表达矩阵清洗与逐值验收

{completed}；{len(CHECKS)} 项检查通过。服务器产物：`{OUT}/`。

## 完成内容

原发布 64253 行特征×120 个样本，97 HCM、23 Control，全数通过样本对应与非零信息检查。可用矩阵保留 {core.n_vars} 行，移除全队列恒为零的 15624 行；结构无效特征 {summary['invalid_feature_identity']} 行。原件和完整发布矩阵全部保留，所有数值与原始 gzip 内表格逐值一致，未归一化、插补、重建或汇总任何表达值。清洗仅表示该发布矩阵、源标签与 QC 产物已验收；不是原始计数语义认证或训练许可。

两个真实 AnnData：all_published_expression.h5ad 保存所有源特征与样本；bulk_qc_expression.h5ad 保存通过特征 QC 的矩阵；bulk_qc_matrix.tsv.gz 保存相同可用矩阵及原始 Gene_Index/Length。X 与 published_values 层完全相同；没有伪称 raw counts 的 counts 层。sample_metadata_QC.tsv、feature_audit.tsv.gz 和排除索引可追溯每个样本及原特征行。sample_exclusions.tsv 的空结果表示真实没有被排除样本，不是占位数据。

## 来源与必须保留的语义限制

GEO GSE249925（https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE249925）；Garmany et al., Journal of Molecular and Cellular Cardiology 2024，PMID39255898，DOI10.1016/j.yjmcc.2024.09.005。原研究报告梗阻性 HCM 切肌病例和供体心脏对照。GEO 的 hg38、Ensembl v78、STAR/HT-Seq/edgeR 处理说明已存档。

GEO 文件名及补充文件描述称其为原始基因计数，但同一处理字段又明确提到所提供读数经过 TMM 标准化。全部数值虽为非负整数，这一点不足以排除已处理或取整来源。已检查 GEO、论文元数据与可访问的出版商摘要/方法片段，尚未取得能解决“该具体文件是否经过 TMM”的权威说明；出版商全文返回403、预印本全文在本次查询中不可访问。EuropePMC在服务器超时，本地成功取得期刊论文元数据；该期刊记录没有可直接读取的PMC全文。保留这些证据与检索限制，**certified_unnormalized_counts=False**。文件可用于受此限制的后续表达验证设计，不能冒充已确认的未归一化 HT-Seq 原计数直接投入计数模型。

## 基因与样本 QC

原始 chromosome@Ensembl@symbol@start@end、原 Length、顺序和所有行均保留；Length 是源注释，不按起止跨度重新推定。重复 Ensembl 行 {summary['duplicate_Ensembl_rows']}，重复 symbol 行 {summary['duplicate_symbol_rows']}，只标记、保留原始特征，不按符号或 Ensembl 擅自求和。恒零行从可用矩阵移出，完整矩阵及 feature_exclusions.tsv.gz 仍含这些原始零行。未按低表达或疾病差异筛基因。

所有标题与 GSM/BioSample/SRA 精确对应，标签直接取 GEO group：HCM=1、Control=0。数字字符串和 R 等前后缀原样保留，不据此猜性别、供体或技术重复。每个 GSM 的源关系均核对；逐样本 donor、age、sex、batch 未提供，donor_id 留空并显式标为未知，不将120个GSM自动认证为120位可独立划分患者，也不将论文性别总数猜配到样本。同一病例跨研究或多标本重叠在训练前仍需审查。

发布表达总量 {int(totals.min())}–{int(totals.max())}；检测特征数 {int(detected.min())}–{int(detected.max())}。log总量/检出特征的5MAD标记分别 {summary['MAD_total_flags']}/{summary['MAD_detected_flags']} 个，仅用于复核、不作删除规则。bulk MT仅按来源染色体计算描述性比例，未套用MT<5%。未执行液滴识别、双细胞检测、心肌提取、单细胞500–5000阈值、差异表达或模型训练。该矩阵是组织混合表达，不能把每个样本当作一枚心肌细胞。

## 验收与使用边界

三个原件大小和此前记录的本地 MD5 均一致；旧下载清单没有 SHA256，本轮新增测得值写入验收记录，不能说与不存在的旧 SHA256 比对通过。该 MD5 是下载基线复查，不宣称获得出版方独立 MD5。首轮验收误把空SHA256当作旧校验值导致退出1，失败日志、源码和原清单已归档；修正后只重跑独立验收，不重跑已经成功的清洗。独立脚本从 gzip 原表及 SOFT 重新解析，核对所有表达值、基因轴、样本轴、120份GEO分组和注册关系、全部元数据/QC选择、文件摘要及执行源码。FINAL_COMPLETION_GATE.json 只在这些检查全部通过后生成；还须核实 verify_exit_code.txt=0 才能记完成。

bulk清洗完成不消除 raw/TMM 语义疑点，也不证明测序质量、性别/批次无混杂或跨队列供体独立；FASTQ/比对级QC本次未执行。训练划分未分配。所有矩阵、逐特征审计保留服务器，Git仅收代码和文本汇总。
'''
    (OUT / 'REPORT.md').write_text(report)
    with (OUT / 'MANIFEST.tsv').open('w', newline='') as handle:
        writer = csv.DictWriter(handle, fieldnames=list(manifest[0]), delimiter='\t')
        writer.writeheader()
        writer.writerows(manifest)
    (OUT / 'MD5SUMS').write_text(''.join(row['MD5'] + '  ' + row['file'] + '\n' for row in manifest))
    gate = {'dataset': 'bulk_GSE249925', 'status': 'complete_verified', 'verified_utc': completed,
            'checks_passed': len(CHECKS), 'verification_sha256': digest(OUT / 'INDEPENDENT_VERIFICATION.json'),
            'report_sha256': digest(OUT / 'REPORT.md'), 'output_directory': str(OUT), 'code_sha256': code_hashes,
            'certified_unnormalized_counts': False, 'single_cell_training_eligible': False,
            'limitation': 'Published integer expression validated; original GEO raw-count label conflicts with TMM statement; no per-sample donor/sex/batch map.'}
    for directory in [WORK, OUT]:
        (directory / 'FINAL_COMPLETION_GATE.json').write_text(json.dumps(gate, indent=2) + '\n')
    print(json.dumps({key: value for key, value in summary.items() if key not in ['all_checks', 'verified_code_sha256']}, indent=2), flush=True)


if __name__ == '__main__':
    try:
        run()
    except BaseException as error:
        (WORK / 'VERIFICATION_FAILURE.json').write_text(json.dumps({'status': 'failed', 'error': str(error),
            'traceback': traceback.format_exc(), 'checks': CHECKS}, indent=2) + '\n')
        raise
