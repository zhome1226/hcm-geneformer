"""Independently compare the prepared object with original compressed HDF5 bytes."""
from datetime import datetime, timezone
from pathlib import Path
import gzip
import hashlib
import io
import json

import anndata as ad
import h5py
import numpy as np
import pandas as pd

ROOT = Path.home() / 'xiaoe-internship-data'
SOURCE = ROOT / 'validation/cox7b_GSE227734'
OUT = ROOT / 'processed/cox7b_GSE227734_qc_20261111_pending_annotations'
checks = {}


def check(name, condition):
    checks[name] = bool(condition)


with gzip.open(SOURCE / 'GSE227734_hcm_raw.h5ad.gz', 'rb') as stream:
    native_bytes = stream.read()
with h5py.File(io.BytesIO(native_bytes)) as handle:
    raw = np.asarray(handle['X'])
    obs = ad.io.read_elem(handle['obs'])
    var = ad.io.read_elem(handle['var'])
    check('official_HDF5_has_no_hidden_gene_annotation', set(handle.keys()) == {'X', 'obs', 'var'} and len(var.columns) == 0)
obj = ad.read_h5ad(OUT / 'author_counts_annotation_pending.h5ad')
preparation = json.loads((OUT / 'PREPARATION.json').read_text())
policy = json.loads((OUT / 'qc_policy.json').read_text())
check('all_published_rows_and_features_retained', obj.shape == raw.shape == (338, 9760))
check('all_counts_exactly_equal_original_compressed_source', np.array_equal(obj.X.toarray(), raw))
with h5py.File(OUT / 'author_counts_annotation_pending.h5ad') as handle:
    stored_layers = set(handle['layers'])
check('counts_layer_exact_and_only_stored_layer', stored_layers == {'counts'} and (obj.X != obj.layers['counts']).nnz == 0)
check('None_key_is_only_API_alias_of_X', None not in obj.layers or (obj.X != obj.layers[None]).nnz == 0)
check('count_domain_is_nonnegative_finite_integer', np.isfinite(raw).all() and (raw >= 0).all() and (raw == np.rint(raw)).all())
check('original_cell_IDs_unique_and_preserved', obj.obs_names.is_unique and np.array_equal(obj.obs.source_cell_id.astype(str), obs.index.astype(str)) and np.array_equal(obj.obs_names, 'GSE227734:' + obs.index.astype(str)))
check('source_row_order_preserved', np.array_equal(obj.obs.source_row, np.arange(338)))
check('anonymous_axis_preserved_without_guessing_gene_names', np.array_equal(obj.var_names, var.index.astype(str)) and np.array_equal(obj.var_names, np.arange(9760).astype(str)) and not any(key in obj.var for key in ['gene_symbol', 'gene_id', 'mt']))
check('anonymous_axis_positions_traceable', np.array_equal(obj.var.source_feature_position, np.arange(9760)))
check('source_obs_columns_declared', json.loads(obj.uns['author_obs_columns_json']) == list(obs.columns))
for key in obs:
    original, delivered = obs[key], obj.obs[key]
    equal = np.array_equal(original.to_numpy(), delivered.to_numpy())
    check('original_obs_' + key + '_unchanged', equal)
totals = raw.sum(axis=1, dtype=np.float64)
detected = (raw > 0).sum(axis=1)
numeric = (detected >= 500) & (detected <= 5000)
check('recomputed_counts_match_original_and_author', np.array_equal(obj.obs.total_counts, totals) and np.array_equal(totals, obs.nCount_RNA))
check('recomputed_detected_features_match_original_and_author', np.array_equal(obj.obs.n_anonymous_features_by_counts, detected) and np.array_equal(detected, obs.nFeature_RNA))
check('gene_bound_flags_recomputed', np.array_equal(obj.obs.pass_gene_count_bounds, numeric) and numeric.sum() == 240)
check('MT_missing_not_fabricated_zero', obj.obs.pct_counts_mt.isna().all() and obj.obs.mt_qc_status.eq('not_computable_without_gene_mapping').all())
check('all_gene_annotations_explicitly_missing', obj.obs.gene_annotation_status.eq('missing_in_official_published_matrix').all() and obj.var.feature_identity_status.eq('anonymous_unmapped_column_not_gene_ID').all())
check('no_false_identity_or_doublet_decisions', obj.obs.cell_identity_status.eq('author_CM_only_not_independently_verified').all() and obj.obs.doublet_status.eq('not_assessed_no_gene_mapping_or_full_capture_population').all())
check('all_training_eligibility_false', not obj.obs.eligible_for_training.any())
check('no_false_full_QC_or_training_split', not obj.uns['dataset_cleaning_complete'] and obj.obs.cleaning_status.eq('blocked_missing_gene_annotation').all() and 'split' not in obj.obs and 'pass_strict_numeric' not in obj.obs)
check('policy_and_preparation_not_marked_clean', policy['dataset_cleaning_complete'] is False and preparation['cleaning_complete'] is False and policy['dataset_status'] == 'blocked_missing_gene_annotation')
check('no_placeholder_clean_objects', not any((OUT / name).exists() for name in ['cm_qc_full.h5ad', 'cm_qc_light.h5ad', 'all_cells_qc.h5ad']))
check('published_morphology_not_recoded_as_disease', obj.obs.disease.eq('HCM').all() and obj.obs.label.eq(1).all() and np.array_equal(obj.obs.area_class.astype(str), obs.area_class.astype(str)))

records = []
for line in (SOURCE / 'sample_metadata.soft').read_text().splitlines():
    if line.startswith('^SAMPLE = '):
        record = {'accession': line.split(' = ', 1)[1]}; records.append(record)
    elif line.startswith('!Sample_') and ' = ' in line:
        key, value = line.removeprefix('!Sample_').split(' = ', 1)
        record.setdefault(key, []).append(value)
summary = pd.read_csv(OUT / 'sample_preparation_summary.tsv', sep='\t')
check('all_five_samples_in_summary', len(summary) == 5 and summary.sample_id.is_unique and set(summary.sample_id) == {r['accession'] for r in records})
for record in records:
    alias = record['description'][0]
    original_mask = obs['sample'].astype(str).eq(alias).to_numpy()
    part = obj.obs.loc[original_mask]
    characteristics = dict(v.split(': ', 1) for v in record['characteristics_ch1'])
    key = record['accession']
    check(key + '_sample_and_patient_trace_to_GEO', part.sample_id.eq(key).all() and part.donor_id.eq('GSE227734:' + alias).all() and part.source_title.eq(record['title'][0]).all())
    check(key + '_age_sex_and_tissue_source', part.age_years.eq(int(characteristics['age'])).all() and part.tissue.eq(characteristics['tissue']).all() and np.array_equal(part.sex.astype(str), obs.loc[original_mask, 'gender'].astype(str)))
    row = summary.loc[summary.sample_id.eq(key)].iloc[0]
    check(key + '_summary_matches_source_and_QC', row.published_cells == original_mask.sum() and row.gene_bound_candidates_only == numeric[original_mask].sum() and row.above_5000_features == (detected[original_mask] > 5000).sum() and row.fully_QC_passed_cells == 0 and row.MT_assessable_cells == 0 and row.training_eligible_cells == 0)
check('summary_whole_dataset_accounting', summary.published_cells.sum() == 338 and summary.gene_bound_candidates_only.sum() == 240 and summary.above_5000_features.sum() == 98)
check('repeated_barcode_sequences_not_deduplicated', obj.obs.Barcode_sequence.duplicated().sum() == obs.Barcode_sequence.duplicated().sum() == 328)
check('whole_cell_and_feature_audits_present', len(pd.read_csv(OUT / 'cell_preparation_audit.tsv.gz', sep='\t')) == 338 and len(pd.read_csv(OUT / 'feature_annotation_audit.tsv.gz', sep='\t')) == 9760)
for name, record in preparation['source_files'].items():
    path = SOURCE / name
    with path.open('rb') as stream:
        sha = hashlib.file_digest(stream, 'sha256').hexdigest()
    check(name + '_original_unchanged', path.stat().st_size == record['bytes'] and sha == record['sha256'])
check('decompressed_source_hash_exact', hashlib.sha256(native_bytes).hexdigest() == preparation['decompressed_sha256'])
hashes = {}
for path in sorted(OUT.iterdir()):
    if path.is_file() and path.name not in ['PREPARATION_VERIFICATION.json', 'MD5SUMS', 'SHA256SUMS', 'REPORT.md']:
        data = path.read_bytes()
        hashes[path.name] = {'bytes': len(data), 'md5': hashlib.md5(data).hexdigest(), 'sha256': hashlib.sha256(data).hexdigest()}
check('prepared_object_hash_matches_export', hashes['author_counts_annotation_pending.h5ad']['sha256'] == preparation['counts_file_sha256'] and hashes['author_counts_annotation_pending.h5ad']['bytes'] == preparation['counts_file_bytes'])
result = {'audit_status': 'passed' if all(checks.values()) else 'failed',
          'dataset_cleaning_status': 'blocked_missing_gene_annotation', 'dataset_cleaning_complete': False,
          'verified_utc': datetime.now(timezone.utc).isoformat(), 'checks': checks, 'hashes': hashes,
          'published_cells': 338, 'anonymous_features': 9760, 'source_patients': 5,
          'gene_bound_candidates_only': int(numeric.sum()), 'full_QC_passed_cells': 0,
          'training_eligible_cells': 0, 'verifier_sha256': hashlib.sha256(Path(__file__).read_bytes()).hexdigest()}
(OUT / 'PREPARATION_VERIFICATION.json').write_text(json.dumps(result, indent=2) + '\n')
if not all(checks.values()):
    raise AssertionError([k for k, passed in checks.items() if not passed])
for algorithm, filename in [('md5', 'MD5SUMS'), ('sha256', 'SHA256SUMS')]:
    (OUT / filename).write_text(''.join(f'{record[algorithm]}  {name}\n' for name, record in hashes.items()))
report = f'''# GSE227734：准备审计通过，清洗尚未完成

{result['verified_utc']}；{len(checks)} 项准备审计检查通过。服务器：`{OUT}/`。

338 个已发表成人原代心肌细胞，来自 5 位 HCM 切肌患者；所有 9760 列原始计数与作者元数据保留。逐值重读原始 gzip 内的 HDF5，计数、细胞顺序、匿名列顺序及作者计数/检出数全部一致，原件摘要未变化。这里的通过仅指准备及来源审计，不能解释为已完成分子质控。

## 必须解决的缺失

GEO 官方 HDF5 只有 X、obs、var；var 无字段，列名仅 0–9759，没有 gene symbol、Ensembl ID、raw 或其他基因映射。不能按字母、参考基因组顺序或其他队列的基因顺序猜配。因而线粒体百分比为 NaN（未测得），不是 0%；心肌标记、参考身份及合适的混合谱系复核也无法执行。240 个细胞符合匿名检出列数 500–5000，98 个高于 5000；这不等于 240 个细胞通过 MT<5% 和身份质控。全部338个细胞 training eligibility=false。

技术来自 STRT-seq 分离成人原代心肌细胞，不能套用 10x 液滴 EmptyDrops 或按装载量估计双细胞。8 碱基条码反复使用，338 个细胞只有 10 种条码；以原始患者前缀细胞 ID 保持唯一，不按条码序列去重。作者 Large/Small 是细胞形态分类，保留原值及面积、长宽、年龄和性别；不转成 Control/HCM 或扰动标签。可访问的元数据未确认形态数值单位，未臆定单位。

GEO 系列设计与 PMID37979444 摘要说明这338个原代心肌来自HCM患者，五个GSM与作者样本别名按原记录绑定；本发布对象没有健康对照。论文包含小鼠AAV9实验，不把人细胞标为COX7B干预。作者心肌标签仅作为来源注释，尚未独立验证；跨研究供体重叠和训练划分未确认。

产物 author_counts_annotation_pending.h5ad 是真实计数与元数据的待注释整理版，不是可用于基因模型的清洗对象；未生成虚假的核心或轻量对象。cell_preparation_audit.tsv.gz、feature_annotation_audit.tsv.gz、sample_preparation_summary.tsv、qc_policy.json 和 PREPARATION_VERIFICATION.json 保存证据。

下一步必须取得与原9760列顺序对应的权威基因表或计数可对应的修正版矩阵，核对来源、顺序与计数后补做MT及身份复核。若只能从原始测序重建，将属于另一个需明确文库/细胞映射、参考版本和计算流程的步骤，不能把新计数冒充原作者矩阵。
'''
(OUT / 'REPORT.md').write_text(report)
print(f'PREPARATION_AUDIT_PASSED: {len(checks)} checks; DATASET_CLEANING_STILL_BLOCKED', flush=True)
