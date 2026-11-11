"""Inspect published Liu counts and donor metadata before any additional filtering."""
from pathlib import Path
from datetime import datetime, timezone
import hashlib
import json
import subprocess
import urllib.request
import xml.etree.ElementTree as ET

import anndata as ad
import numpy as np
import pandas as pd
from scipy.io import mmread
from scipy.sparse import csc_matrix
from threadpoolctl import threadpool_limits
import h5py

PANELS = {
    'Cardiomyocyte': ['TNNT2','TTN','TNNI3','ACTC1','MYH7','MYBPC3','MYL2','TNNC1','RYR2','PLN','ACTN2','MYH6'],
    'Fibroblast': ['COL1A1','COL1A2','DCN','LUM','PDGFRA'],
    'Endothelial': ['PECAM1','VWF','KDR','EMCN','CLDN5'],
    'Immune': ['PTPRC','LST1','TYROBP','CD3D','CD79A'],
    'Mural': ['RGS5','PDGFRB','CSPG4','MYH11','ACTA2','TAGLN'],
}
PAIRS = [('Cardiomyocyte','Fibroblast'),('Cardiomyocyte','Endothelial'),('Cardiomyocyte','Immune'),
         ('Cardiomyocyte','Mural'),('Fibroblast','Endothelial'),('Fibroblast','Immune'),
         ('Endothelial','Immune'),('Immune','Mural')]

ROOT = Path.home() / 'xiaoe-internship-data'
WORK = ROOT / '_processing/all_data_cleaning_20261111/liu2023'


def digest(path):
    with path.open('rb') as stream:
        return hashlib.file_digest(stream, 'sha256').hexdigest()


def read_mtx(directory):
    with threadpool_limits(limits=2):
        x = mmread(directory / 'matrix.mtx.gz').T.tocsr()
    x.sum_duplicates(); x.eliminate_zeros(); x.sort_indices()
    return x


def read_author():
    with h5py.File(WORK / 'liu_RNA_counts.h5') as f:
        g = f['matrix']
        x = csc_matrix((g['data'][:], g['indices'][:], g['indptr'][:]), shape=g['shape'][:]).T.tocsr()
        genes = g['features/id'].asstr()[:]
        barcodes = g['barcodes'].asstr()[:]
    x.eliminate_zeros(); x.sort_indices()
    return ad.AnnData(x, obs=pd.DataFrame(index=barcodes), var=pd.DataFrame(index=genes))


def main():
    processes = subprocess.run(['pgrep', '-af', '[p]reflight_liu|[r]un_liu|[f]inalize_liu'], capture_output=True, text=True).stdout
    print('Preflight process context: ' + processes, flush=True)
    xml_path = WORK / 'PMC9842679.xml'
    if not xml_path.exists():
        request = urllib.request.Request('https://www.ebi.ac.uk/europepmc/webservices/rest/PMC9842679/fullTextXML', headers={'User-Agent': 'HCM-validation-data-QC/1.0'})
        with urllib.request.urlopen(request, timeout=40) as stream:
            data = stream.read()
        ET.fromstring(data)
        xml_path.write_bytes(data)
    tree = ET.parse(xml_path)
    sections = []
    for sec in tree.iter('sec'):
        title = sec.find('title')
        if title is not None and ('quality control' in ''.join(title.itertext()).lower() or 'preprocessing' in ''.join(title.itertext()).lower()):
            sections.append(' '.join(sec.itertext()))
    (WORK / 'AUTHOR_METHODS.txt').write_text('\n\n'.join(sections) + '\n')
    meta = pd.read_csv(WORK / 'author_metadata.tsv.gz', sep='\t', index_col='barcode', keep_default_na=False)
    author = read_author()
    assert author.obs_names.equals(meta.index)
    assert author.obs_names.is_unique and author.var_names.is_unique
    assert len(meta) == 55122
    expected_donors = {'CL0412': 'Control', 'CL2': 'Control', **{s: 'HCM' for s in ['HCM300','HCM456','HCM483','HCM516','HCM534','HCM621','HCM752','HCM763','HCM770','HCM823']}}
    rows, inspection = [], []
    accounted = []
    originals = {}
    for path in sorted((ROOT / 'validation/liu2023').glob('*')):
        if path.is_file():
            md5 = hashlib.md5()
            sha = hashlib.sha256()
            with path.open('rb') as stream:
                while chunk := stream.read(32 * 1024**2):
                    md5.update(chunk); sha.update(chunk)
            originals[path.name] = {'bytes': path.stat().st_size, 'md5': md5.hexdigest(), 'sha256': sha.hexdigest()}
    for name, size, md5 in [('HCM_snRNAseq_seuratObj.Rdata',2376598237,'22153f776c3e667d87feead247272799'),('snRNA-seq.tar.gz',450043610,'3d850cd3017a9c8623febdc4269e71fa')]:
        assert originals[name]['bytes'] == size and originals[name]['md5'] == md5
    for sample, frame in meta.groupby('orig.ident', sort=True):
        raw_dir = WORK / 'snRNA-seq' / sample
        features = pd.read_csv(raw_dir / 'features.tsv.gz', sep='\t', header=None, dtype=str, keep_default_na=False)
        barcodes = pd.read_csv(raw_dir / 'barcodes.tsv.gz', header=None, dtype=str)[0]
        x = read_mtx(raw_dir)
        assert x.shape == (len(barcodes), len(features))
        assert barcodes.is_unique and features[0].is_unique
        assert np.array_equal(features[0], features[1]) and np.array_equal(features[0], author.var_names)
        assert set(barcodes) == set(frame.index)
        assert (x != author[barcodes].X).nnz == 0
        assert np.isfinite(x.data).all() and (x.data >= 0).all() and (x.data == np.round(x.data)).all()
        frame = frame.loc[barcodes]
        totals = np.asarray(x.sum(axis=1)).ravel()
        genes = np.diff(x.indptr)
        mt = 100 * np.asarray(x[:, features[1].str.startswith('MT-').to_numpy()].sum(axis=1)).ravel() / np.maximum(totals,1)
        assert np.array_equal(totals, frame.nCount_RNA) and np.array_equal(genes, frame.nFeature_RNA)
        assert np.allclose(mt, 100 * pd.to_numeric(frame['percent.mito']))
        assert frame.Subject.nunique() == 1 and frame.Group.nunique() == 1 and frame.group.nunique() == 1
        donor = frame.Subject.iloc[0]
        disease = expected_donors[donor]
        assert frame.Group.eq('HEALTHY' if disease == 'Control' else 'HCM').all()
        assert frame.group.eq('CTRL' if disease == 'Control' else 'CASE').all()
        sample_dir = WORK / sample
        sample_dir.mkdir(exist_ok=True)
        obs = frame.rename(columns={name: 'author_' + name.replace('.', '_') for name in frame.columns}).copy()
        obs['barcode_original'] = obs.index
        obs.index = sample + ':' + obs.index
        obs['dataset'] = 'liu2023'
        obs['study'] = 'Liu_CellDiscovery_2023'
        obs['sample_id'] = sample
        obs['capture_id'] = sample
        obs['donor_original'] = donor
        obs['donor_id'] = 'Liu2023_' + donor
        obs['disease'] = disease
        obs['label'] = 0 if disease == 'Control' else 1
        obs['tissue'] = 'interventricular septum'
        obs['assay'] = '10x 3prime v3.1 snRNA-seq'
        obs['source_pre_filtered'] = True
        obs['total_counts'] = totals
        obs['n_genes_by_counts'] = genes
        obs['pct_counts_mt'] = mt
        obs['pass_strict_numeric'] = (genes >= 500) & (genes <= 5000) & (mt < 5)
        marker_n, programs = {}, {}
        for label, panel in PANELS.items():
            values = x[:, features[1].isin(panel).to_numpy()].toarray()
            marker_n[label] = (values > 0).sum(axis=1)
            programs[label] = np.log1p(values * (10000 / np.maximum(totals,1))[:,None]).mean(axis=1)
        mixed = np.zeros(len(obs),dtype=bool)
        for left,right in PAIRS:
            ratio = np.minimum(programs[left],programs[right]) / np.maximum(np.maximum(programs[left],programs[right]),1e-12)
            mixed |= (marker_n[left]>=3) & (marker_n[right]>=3) & (ratio>=.25)
        obs['residual_mixed_screen'] = mixed
        n_mixed = int((mixed & obs.pass_strict_numeric.to_numpy()).sum())
        strict_n = int(obs.pass_strict_numeric.sum())
        trigger = n_mixed>=10 and n_mixed/max(strict_n,1)>=.01
        var = pd.DataFrame({'gene_symbol': features[1].to_numpy(), 'mt': features[1].str.startswith('MT-').to_numpy()},index=features[0].to_numpy())
        a = ad.AnnData(x, obs=obs, var=var)
        a.write_h5ad(sample_dir / 'author_counts.h5ad', compression='gzip')
        obs.to_csv(sample_dir / 'input_qc.tsv.gz', sep='\t', index_label='barcode')
        sources = {name: {'path': str(raw_dir / name), 'bytes': (raw_dir / name).stat().st_size,'sha256': digest(raw_dir / name)} for name in ['matrix.mtx.gz','barcodes.tsv.gz','features.tsv.gz']}
        (sample_dir / 'STAGED.json').write_text(json.dumps({'sample_id':sample,'sources':sources,'nuclei':a.n_obs,'genes':a.n_vars},indent=2)+'\n')
        rows.append({'dataset':'liu2023','sample_id':sample,'donor_id':'Liu2023_'+donor,'donor_original':donor,'disease':disease,'raw_directory':str(raw_dir),'residual_mixed_strict':n_mixed,'residual_mixed_strict_fraction':n_mixed/max(strict_n,1),'trigger_residual_doublet_screen':trigger})
        inspection.append({'sample_id':sample,'donor':donor,'disease':disease,'input_nuclei':a.n_obs,'strict_nuclei':int(obs.pass_strict_numeric.sum()),'below_500':int((genes<500).sum()),'above_5000':int((genes>5000).sum()),'MT_ge5':int((mt>=5).sum()),'author_CM':int(frame.cellType.eq('CM').sum()),'matrix_equal_author_RNA_counts':True,'author_MT_unit':'fraction; multiply by 100 to obtain percent'})
        accounted.extend(barcodes)
        print(json.dumps(inspection[-1]), flush=True)
    assert len(accounted) == len(set(accounted)) == len(meta) and set(accounted) == set(meta.index)
    samples = pd.DataFrame(rows)
    assert samples.sample_id.nunique() == 13 and samples.donor_id.nunique() == 12
    assert set(samples.donor_original) == set(expected_donors)
    assert samples[samples.sample_id.isin(['CL0412A','CL0412B'])].donor_id.eq('Liu2023_CL0412').all()
    samples.to_csv(WORK / 'samples.tsv', sep='\t', index=False)
    pd.DataFrame(inspection).to_csv(WORK / 'PREFLIGHT_SAMPLE_SUMMARY.tsv',sep='\t',index=False)
    refdir=ROOT / 'processed/SCP1303_review_20261111/identity'
    refgenes=pd.read_csv(refdir/'reference_genes.tsv',sep='\t')
    mapping=refgenes[~refgenes.gene_symbol.duplicated(keep=False) & refgenes.gene_symbol.isin(author.var_names)].copy()
    assert len(mapping)>25000 and mapping.gene_symbol.is_unique
    mapping.to_csv(WORK/'reference_gene_mapping.tsv',sep='\t',index=False)
    result={'status':'passed_preflight_only','checked_utc':datetime.now(timezone.utc).isoformat(),'nuclei':len(meta),'genes':author.n_vars,'libraries':13,'biological_donors':12,'Control_donors':2,'HCM_donors':10,'CL0412A_B':'Author Subject assigns both samples to CL0412. Treat as two samples of one biological donor; exact technical versus tissue replicate mechanism not established. Never split donor across partitions.','counts_provenance':'All 13 released matrices equal corresponding author RNA count layer exactly. They contain exactly the 55122 published post-QC nuclei, not full raw droplets. No new EmptyDrops or claim to recover removed nuclei/genes.','mitochondrial_units':'author percent.mito is a fraction; recalculate 100*MT_UMI/total_UMI','originals':originals,'reference_common_unique_symbols':len(mapping),'source_methods_xml':str(xml_path),'source_methods_xml_sha256':digest(xml_path),'source_methods_url':'https://www.ebi.ac.uk/europepmc/webservices/rest/PMC9842679/fullTextXML','samples':inspection}
    (WORK/'PREFLIGHT.json').write_text(json.dumps(result,indent=2)+'\n')
    print('PREFLIGHT_PASSED',flush=True)


if __name__ == '__main__':
    main()
