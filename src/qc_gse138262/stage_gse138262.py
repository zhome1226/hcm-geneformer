"""Stage full corrected count tables, retaining native axes and author provenance."""
from datetime import datetime, timezone
from pathlib import Path
import hashlib
import json
import re

import anndata as ad
import h5py
import numpy as np
import pandas as pd
from scipy import sparse
from scipy.io import mmwrite
from threadpoolctl import threadpool_limits

ROOT = Path.home() / 'xiaoe-internship-data'
WORK = ROOT / '_processing/all_data_cleaning_20261111/gse138262'
SOURCE = ROOT / 'validation/hcm_GSE138262/author_corrected_zenodo18206554'
REF = ROOT / 'processed/SCP1303_review_20261111/identity'


def sha256(path):
    with Path(path).open('rb') as handle:
        return hashlib.file_digest(handle, 'sha256').hexdigest()


def load_table(path):
    # Author write.table uses spaces despite .tsv suffix; header has no row-name label.
    table = pd.read_csv(path, sep=r'\s+', index_col=0, keep_default_na=False)
    table.index = table.index.astype(str)
    table.columns = table.columns.astype(str)
    assert table.index.is_unique and table.columns.is_unique
    values = table.to_numpy()
    assert np.isfinite(values).all() and (values >= 0).all()
    assert np.equal(values, np.rint(values)).all(), path
    # Missing author symbols are genuine empty suffixes; stable Ensembl IDs remain usable.
    assert all(re.fullmatch(r'ENSG\d+_.*', x) for x in table.index), path
    return table


def read_h5(path):
    with h5py.File(path) as f:
        m = f['matrix']
        matrix = sparse.csc_matrix((m['data'][:], m['indices'][:], m['indptr'][:]), shape=m['shape'][:]).T.tocsr()
        decode = lambda d: pd.Index(d.asstr()[:])
        return matrix, decode(m['barcodes']), decode(m['features/id'])


def sample_definitions():
    notes = (WORK / 'AUTHOR_DATA_NOTES.md').read_text()
    main = (WORK / 'AUTHOR_MAIN.R').read_text()
    hcm_block = main.split('HCM_SCS_ids =', 1)[1].split('dataset_list_paths_HCM', 1)[0]
    hu_block = main.split('WANG_ids =', 1)[1].split('dataset_list_paths_WANG', 1)[0]
    hcm_ids = dict((value, key) for key, value in re.findall(r'(\w+)\s*=\s*\x27([^\x27]+)\x27', hcm_block))
    hu_ids = dict((value, key) for key, value in re.findall(r'([\w.]+)\s*=\s*\x27([^\x27]+)\x27', hu_block)
                  if value.startswith('p.N'))
    assert len(hcm_ids) == 16 and len(hu_ids) == 16
    geo = {r['accession']: r for r in json.loads((WORK / 'GEO_SAMPLES_GSE138262.json').read_text(encoding='utf-8-sig'))}
    hcm_rows = re.findall(r'\|\s*(GSM\d+)\s*\|\s*Patient\s+(\d+)\s*\|\s*Plate\s+(\d+)\s*\|\s*([^|\s]+\.tsv)\s*\|', notes)
    assert len(hcm_rows) == 16
    rows = []
    for gsm, donor, plate, filename in hcm_rows:
        prefix = filename.split('_pT.nonRibo_')[0]
        path = list((SOURCE / 'counttables-ROOIJ').rglob(filename))
        assert len(path) == 1 and prefix in hcm_ids
        assert f'patient id: {donor}' in geo[gsm]['characteristics_ch1']
        assert 'patient diagnosis: hypertrophic cardiomyopathy' in geo[gsm]['characteristics_ch1']
        rows.append({'sample_id': gsm, 'author_prefix': hcm_ids[prefix], 'donor_id': 'Wehrens_Patient' + donor,
                     'donor_original': donor, 'plate': plate, 'source_path': str(path[0]), 'disease': 'HCM', 'label': 1,
                     'study_origin': 'Wehrens2022', 'source_accession': 'GSE138262', 'technology': 'CEL-seq2 FACS-sorted whole cells',
                     'tissue': 'left_ventricular_myectomy', 'clinical_stage': 'HCM_septal_myectomy',
                     'donor_overlap_status': 'unresolved_across_HCM_studies',
                     'disease_source_url': 'https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=' + gsm})
    hu_rows = re.findall(r'\|\s*(GSM\d+)\s*\|\s*(p\.N\d+\.plate\.[^|\s]+\.tsv)\s*\|', notes)
    assert len(hu_rows) == 16
    for gsm, filename in hu_rows:
        prefix = filename.split('_nc.nonRibo_')[0]
        author_prefix = hu_ids[prefix]
        donor = re.fullmatch(r'p\.(N\d+)\.plate\.(\d+)\.part\.(\d+)_cat', prefix)
        assert donor
        path = list((SOURCE / 'counttables-HU').rglob(filename))
        assert len(path) == 1
        rows.append({'sample_id': 'Wang_' + author_prefix, 'author_prefix': author_prefix,
                     'donor_id': 'Wang_' + donor[1], 'donor_original': donor[1],
                     'plate': donor[2] + '_part' + donor[3], 'source_path': str(path[0]), 'disease': 'Control', 'label': 0,
                     'study_origin': 'Wang2020_reprocessed_by_Wehrens', 'source_accession': gsm,
                     'technology': 'iCell8 whole cells; selected LV clusters', 'tissue': 'left_ventricle',
                     'clinical_stage': 'healthy_donor', 'donor_overlap_status': 'bind_Wang_donor_across_GSE121893_GSE109816',
                     'disease_source_url': 'https://doi.org/10.5281/zenodo.18206554'})
    assert len(rows) == 32 and len({r['sample_id'] for r in rows}) == 32
    assert len({r['donor_id'] for r in rows}) == 12
    return rows


def main():
    assert (WORK / 'AUTHOR_EXPORT_SUCCESS.txt').exists()
    assert json.loads((WORK / 'DOWNLOAD_VERIFICATION.json').read_text())['record_id'] == 18206554
    rows = sample_definitions()
    metadata = pd.read_csv(WORK / 'author_ALL/metadata.tsv.gz', sep='\t', index_col='barcode', keep_default_na=False)
    assert metadata.index.is_unique
    tables, audits, author_comparisons = [], [], []
    author_assays = []
    for kind in ['ROOIJonly', 'ALL']:
        for path in sorted((WORK / ('author_' + kind)).glob('*_counts.h5')):
            matrix, cells, genes = read_h5(path)
            assert cells.is_unique and genes.is_unique
            author_assays.append((kind, path.stem, matrix, cells, genes))
    for row in rows:
        directory = WORK / row['sample_id']
        directory.mkdir(exist_ok=True)
        source = Path(row['source_path'])
        table = load_table(source)
        x = sparse.csr_matrix(table.to_numpy(dtype=np.int64).T)
        gene_ids = table.index.str.split('_', n=1).str[0]
        symbols = table.index.str.split('_', n=1).str[1]
        assert gene_ids.is_unique
        author_cells = pd.Index([row['author_prefix'] + '_' + str(c) for c in table.columns])
        obs = pd.DataFrame(index=pd.Index(row['sample_id'] + ':' + table.columns, name='cell_id'))
        obs['source_cell_id'] = table.columns.to_numpy()
        obs['author_cell_id'] = author_cells.to_numpy()
        for key in ['sample_id','author_prefix','donor_id','donor_original','plate','disease','label',
                    'study_origin','source_accession','technology','tissue','clinical_stage','donor_overlap_status','disease_source_url']:
            obs[key] = row[key]
        author_meta = metadata.reindex(author_cells).reset_index(drop=True).set_axis(obs.index)
        obs['author_analysis_included'] = author_cells.isin(metadata.index)
        obs = obs.join(author_meta.add_prefix('author_'))
        for col in obs.select_dtypes(include=['object']).columns:
            obs[col] = obs[col].fillna('').astype(str)
        var = pd.DataFrame({'gene_id': gene_ids.to_numpy(), 'gene_symbol': symbols.to_numpy()}, index=table.index.copy())
        var.index.name = 'source_feature_id'
        var['mt'] = var.gene_symbol.str.startswith('MT-')
        totals = np.asarray(x.sum(axis=1)).ravel()
        obs['total_counts'] = totals
        obs['n_genes_by_counts'] = np.diff(x.indptr)
        obs['pct_counts_mt'] = np.divide(np.asarray(x[:,var.mt.to_numpy()].sum(axis=1)).ravel() * 100., totals,
                                         out=np.full(len(obs), np.nan), where=totals>0)
        obs['pass_strict_numeric'] = obs.n_genes_by_counts.between(500,5000) & obs.pct_counts_mt.lt(5)
        obs['numeric_exclusion_reason'] = np.select([totals==0, obs.n_genes_by_counts.lt(500),
                                                     obs.n_genes_by_counts.gt(5000), obs.pct_counts_mt.ge(5)],
                                                    ['zero_counts','genes_lt500','genes_gt5000','MT_ge5pct'],default='pass')
        obs['doublet_assessment'] = 'not_modeled_plate_sorted_no_10x_loading_assumption'
        obs['split_status'] = 'unassigned_requires_donor_overlap_and_design_review'
        a = ad.AnnData(x, obs=obs, var=var)
        a.uns['count_source'] = str(source)
        a.uns['source_sha256'] = sha256(source)
        a.uns['axis_policy'] = 'Native plate feature/cell axes retained; no normalization or gene filtering'
        a.write_h5ad(directory / 'native_counts.h5ad', compression='gzip')
        obs.to_csv(directory / 'cell_numeric_audit.tsv.gz', sep='\t')
        mapped_features = pd.Index([g.replace('_', ':') for g in table.index])
        for kind, assay, matrix, cells, genes in author_assays:
            if kind == 'ROOIJonly' and row['disease'] == 'Control':
                continue
            selected = np.flatnonzero(author_cells.isin(cells))
            shared = np.flatnonzero(mapped_features.isin(genes))
            if not len(selected) or not len(shared):
                continue
            expected = x[selected][:,shared]
            observed = matrix[cells.get_indexer(author_cells[selected])][:,genes.get_indexer(mapped_features[shared])]
            difference = (expected - observed).tocsr()
            difference.eliminate_zeros()
            assert difference.nnz == 0, (row['sample_id'],kind,assay,difference.nnz)
            author_comparisons.append({'sample_id':row['sample_id'],'object':kind,'assay':assay,
                                       'shared_cells':len(selected),'shared_genes':len(shared),'differing_values':difference.nnz})
        row.update({'native_cells':a.n_obs,'native_genes':a.n_vars,'strict_cells':int(obs.pass_strict_numeric.sum()),
                    'author_analysis_cells':int(obs.author_analysis_included.sum()),'mt_gene_count':int(var.mt.sum()),
                    'source_bytes':source.stat().st_size,'source_sha256':sha256(source)})
        tables.append(a)
        audits.append(obs)
        print(f'{row["sample_id"]}: {a.n_obs} cells, {a.n_vars} genes, strict={row["strict_cells"]}, MT genes={row["mt_gene_count"]}', flush=True)
    all_counts = ad.concat(tables, join='outer', fill_value=0, merge='first')
    all_counts.X = all_counts.X.tocsr()
    assert all_counts.obs_names.is_unique and all_counts.var_names.is_unique
    native_vars = pd.concat([a.var for a in tables])
    for column in ['gene_id','gene_symbol','mt']:
        assert native_vars.groupby(level=0)[column].nunique().le(1).all()
    all_counts.var = native_vars.loc[~native_vars.index.duplicated()].reindex(all_counts.var_names).copy()
    feature_mask = np.column_stack([all_counts.var_names.isin(a.var_names) for a in tables])
    all_counts.varm['native_feature_presence_by_sample'] = feature_mask
    all_counts.uns['feature_presence_sample_order'] = [r['sample_id'] for r in rows]
    all_counts.uns['missing_feature_policy'] = 'Union stores zeros; per-sample native axis/presence mask retained. Source tables were generated from observed UMI records; gene omissions are not evidence of independent detection.'
    all_counts.uns['count_source'] = 'Zenodo version DOI 10.5281/zenodo.18206554 corrected full UMI count tables; includes MT before author RNA-assay removal'
    all_counts.write_h5ad(WORK / 'corrected_native_union_counts.h5ad', compression='gzip')
    all_counts.obs.to_csv(WORK / 'cell_numeric_audit.tsv.gz', sep='\t')
    samples = pd.DataFrame(rows)
    samples.to_csv(WORK / 'samples.tsv', sep='\t', index=False)
    pd.DataFrame(author_comparisons).to_csv(WORK / 'AUTHOR_COUNT_COMPARISONS.tsv', sep='\t', index=False)
    metadata.loc[~metadata.index.isin(all_counts.obs.author_cell_id)].to_csv(WORK / 'author_external_or_unmatched_cells.tsv.gz', sep='\t')
    reference = pd.read_csv(REF / 'reference_genes.tsv', sep='\t')
    assert all_counts.var.gene_id.is_unique
    union_ids = pd.Index(all_counts.var.gene_id)
    common = reference.gene_id[reference.gene_id.isin(union_ids)].tolist()
    assert len(common)>5000
    positions = union_ids.get_indexer(common)
    strict = all_counts.obs.pass_strict_numeric.to_numpy()
    with threadpool_limits(limits=2):
        mmwrite(WORK / 'strict_reference_query.mtx', all_counts.X[strict][:,positions].T)
    pd.Series(common).to_csv(WORK / 'reference_common_genes.txt', index=False, header=False)
    pd.Series(all_counts.obs_names[strict]).to_csv(WORK / 'strict_query_cells.txt', index=False, header=False)
    summary = {'checked_utc': datetime.now(timezone.utc).isoformat(), 'cells': all_counts.n_obs,
               'genes': all_counts.n_vars, 'samples':len(rows),'donors':12,'HCM_donors':5,'Control_donors':7,
               'strict_cells':int(strict.sum()),'common_reference_genes':len(common),
               'author_assay_comparisons':len(author_comparisons),
               'source_hashes':{r['sample_id']:r['source_sha256'] for r in rows},
               'scope':'Primary HCM plus separately identified Wang reprocessed control companion; author ALL Litvinukova cells inventoried separately, not appended or treated as independent disease validation'}
    (WORK / 'STAGE_SUCCESS.json').write_text(json.dumps(summary,indent=2)+'\n')
    print(json.dumps(summary),flush=True)


if __name__ == '__main__':
    main()
