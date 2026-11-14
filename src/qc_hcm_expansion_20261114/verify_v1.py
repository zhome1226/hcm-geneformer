"""Independent validation from count matrices; does not call Geneformer rank functions."""
import argparse, hashlib, json, pickle, re
from datetime import datetime, timezone
from pathlib import Path
from collections import Counter
import anndata as ad
import h5py
import numpy as np
import pandas as pd
import scipy.sparse as sp
from anndata._io.specs import read_elem
from datasets import load_from_disk

W=Path('/home/zhome/xiaoe-internship-data/_processing/geneformer_v1_tokenization_20261112')
O=Path('/home/zhome/hcm-geneformer-data/tokenized/geneformer_v1_20261114')
D=W/'official_source/geneformer/gene_dictionaries_30m'

def digest(p):
    h=hashlib.sha256()
    with open(p,'rb') as f:
        for b in iter(lambda:f.read(8*1024*1024),b''):h.update(b)
    return h.hexdigest()

def main():
    parser=argparse.ArgumentParser();parser.add_argument('--mode',choices=['smoke','full'],required=True)
    mode=parser.parse_args().mode
    token=pickle.loads((D/'token_dictionary_gc30M.pkl').read_bytes())
    median=pickle.loads((D/'gene_median_dictionary_gc30M.pkl').read_bytes())
    mapping=pickle.loads((D/'ensembl_mapping_dict_gc30M.pkl').read_bytes())
    config=json.loads((Path('/home/zhome/hcm-geneformer-analysis/hcm_expansion_20261114/nuclear_qc')/'INPUT_CONFIG.json').read_text())
    prepared={r['dataset']:r for r in json.loads((O/'prepared_inputs/PREPARATION_VERIFICATION.json').read_text())['datasets']}
    assert all(k==v for k,v in mapping.items()),'Revisit aggregation: this fixed V1 mapping was expected to be identity'
    assert len(set(token.values()))==len(token)
    results=[];all_ids=set();manifest=[]
    for ds,record in config.items():
        inp=O/('smoke_inputs' if mode=='smoke' else 'prepared_inputs')/(ds+'.h5ad')
        a=ad.read_h5ad(inp)
        if mode=='full':
            assert inp.stat().st_size==prepared[ds]['prepared_bytes'] and digest(inp)==prepared[ds]['prepared_sha256'],(ds,'prepared input changed')
            assert Path(record['path']).stat().st_size==record['bytes'] and digest(record['path'])==record['sha256']
            with h5py.File(record['path'],'r') as f:
                original=read_elem(f['X']).tocsr();obs=read_elem(f['obs']);var=read_elem(f['var'])
            # CSR column order is representation, not biological values: preparation
            # canonicalizes storage while retaining every cell x gene coordinate.
            assert a.X.shape==original.shape and (a.X!=original).nnz==0
            assert np.array_equal(a.obs_names,obs.index) and np.array_equal(a.var_names,var.index)
            for c in ['donor_id','sample_id','label','barcode_original']:
                if c=='barcode_original' and c not in obs:
                    assert np.array_equal(a.obs[c].astype(str),obs.author_cell_id.astype(str)),(ds,c)
                else:assert np.array_equal(a.obs[c].astype(str),obs[c].astype(str)),(ds,c)
            assert np.array_equal(a.obs.source_disease.astype(str),obs.disease.astype(str))
            if 'split' in obs:assert np.array_equal(a.obs.source_split.astype(str),obs.split.astype(str))
            if 'ensembl_id' in var:assert np.array_equal(a.var.ensembl_id.astype(str),var.ensembl_id.astype(str))
            elif 'gene_id' in var:assert np.array_equal(a.var.ensembl_id.astype(str),var.gene_id.astype(str))
            else:assert np.array_equal(a.var.ensembl_id.astype(str),var.index.astype(str))
            with h5py.File(record['path'],'r') as f:
                for key in a.varm:assert np.array_equal(a.varm[key],read_elem(f['varm'][key]))
                original_counts=read_elem(f['layers']['counts']).tocsr()
                assert (a.layers['counts']!=original_counts).nnz==0 and (a.X!=a.layers['counts']).nnz==0
            del original
        assert np.array_equal(a.obs.n_counts,np.asarray(a.X.sum(axis=1)).ravel())
        assert (a.obs.split.astype(str)=='unassigned').all()
        assert np.array_equal(a.obs.disease.astype(str),np.where(a.obs.label==0,'Control','HCM'))
        ids=a.obs.cell_id.astype(str).tolist();assert len(set(ids))==len(ids)
        assert not all_ids.intersection(ids);all_ids.update(ids)
        assert ids==[ds+'::'+str(v) for v in a.obs_names]
        # Aggregate matching gene counts independently with a sparse projection.
        mapped=[mapping.get(str(g).upper()) for g in a.var.ensembl_id]
        accepted=sorted(set(g for g in mapped if g in median and g in token))
        cols={g:i for i,g in enumerate(accepted)}
        pos=[i for i,g in enumerate(mapped) if g in cols]
        projection=sp.csr_matrix((np.ones(len(pos),dtype=a.X.dtype),
                    (pos,[cols[mapped[i]] for i in pos])),shape=(a.n_vars,len(cols)))
        counts=(a.X@projection).tocsr();counts.sum_duplicates();counts.eliminate_zeros();counts.sort_indices()
        genes=np.array([int(token[g]) for g in accepted]);med=np.array([median[g] for g in accepted])
        destination=O/('smoke_datasets' if mode=='smoke' else 'datasets')/(ds+'.dataset')
        hf=load_from_disk(str(destination));assert len(hf)==a.n_obs
        required=['cell_id','source_cell_id','barcode_original','donor_id','sample_id','study','dataset','disease','source_disease','label','cell_type','split','split_status','n_counts']
        assert set(required+['input_ids','length'])<=set(hf.column_names)
        metadata=[c for c in hf.column_names if c not in ['input_ids','length']]
        expressed=np.diff(counts.indptr);assert (expressed>0).all()
        lengths=[];offset=0
        for batch in hf.iter(batch_size=256):
            n=len(batch['input_ids']);end=offset+n
            for c in metadata:
                assert np.array_equal(np.asarray(batch[c]).astype(str),a.obs[c].iloc[offset:end].to_numpy().astype(str)),(ds,c,offset)
            # Arithmetic follows the documented total-count and per-gene median normalization.
            scores=(counts[offset:end]/a.obs.n_counts.iloc[offset:end].to_numpy()[:,None]*10000)/med
            scores=sp.csr_matrix(scores);scores.sort_indices()
            for j,sequence in enumerate(batch['input_ids']):
                v=np.asarray(sequence,dtype=np.int64);i=offset+j;expected=min(int(expressed[i]),2048)
                assert len(v)==batch['length'][j]==expected and len(set(v))==len(v),(ds,i,'length/duplicates')
                assert (v>=2).all() and (v<25426).all(),(ds,i,'vocabulary/special tokens')
                sl=slice(scores.indptr[j],scores.indptr[j+1]);s=scores.data[sl];tokens=genes[scores.indices[sl]]
                token_scores=dict(zip(tokens.tolist(),s.tolist()))
                assert set(v)<=set(token_scores),(ds,i,'unexpected genes')
                ranked=np.array([token_scores[int(t)] for t in v])
                assert (ranked[:-1]>=ranked[1:]).all(),(ds,i,'rank order')
                if len(s)>len(v):
                    chosen=set(v.tolist());omitted=[value for t,value in token_scores.items() if t not in chosen]
                    assert max(omitted)<=ranked[-1],(ds,i,'truncation cutoff')
                # Only exactly equal computed scores may appear in either tie order.
                lengths.append(len(v))
            offset=end
        assert offset==a.n_obs
        coverage=np.asarray(counts.sum(axis=1)).ravel()/a.obs.n_counts.to_numpy()
        bydonor=pd.DataFrame({'donor_id':a.obs.donor_id.astype(str).to_numpy(),'coverage':coverage,'expressed':expressed}).groupby('donor_id').agg(cells=('coverage','size'),median_count_fraction=('coverage','median'),min_expressed=('expressed','min'))
        bydonor.to_csv(O/(mode+'_'+ds+'_DONOR_COVERAGE.tsv'),sep='\t')
        result={'dataset':ds,'cells':a.n_obs,'donors':a.obs.donor_id.nunique(),
                'disease_counts':{str(k):int(v) for k,v in a.obs.disease.value_counts().items()},
                'length_min':min(lengths),'length_median':float(np.median(lengths)),'length_max':max(lengths),
                'expressed_vocab_genes_min':int(expressed.min()),'expressed_vocab_genes_median':float(np.median(expressed)),
                'truncated_cells':int((expressed>2048).sum()),'truncated_fraction':float((expressed>2048).mean()),
                'source_features':a.n_vars,'matched_features':len(pos),'unique_vocab_genes':len(cols),
                'duplicate_matched_features_collapsed':len(pos)-len(cols),
                'count_fraction_min':float(coverage.min()),'count_fraction_median':float(np.median(coverage)),
                'all_cells_all_tokens_verified':True,'metadata_exact':True,'source_unchanged':mode=='full',
                'path':str(destination)}
        results.append(result)
        for p in sorted(destination.rglob('*')):
            if p.is_file():manifest.append({'dataset':ds,'path':str(p.relative_to(O)),'bytes':p.stat().st_size,'sha256':digest(p)})
        (O/(mode.upper()+'_VERIFICATION_PROGRESS.json')).write_text(json.dumps(results,indent=2)+'\n')
        print(json.dumps(result),flush=True)
    report={'status':'passed','mode':mode,'utc':datetime.now(timezone.utc).isoformat(),'total_cells':sum(r['cells'] for r in results),'datasets':results,
            'checks':'Every cell: source counts, metadata, vocabulary, uniqueness, lengths, descending median-scaled ranks and exact cutoff ties; no Geneformer rank function used.'}
    pd.DataFrame(manifest).to_csv(O/(mode.upper()+'_MANIFEST.tsv'),sep='\t',index=False)
    (O/(mode.upper()+'_VERIFICATION.json')).write_text(json.dumps(report,indent=2)+'\n')
    print('INDEPENDENT_VERIFICATION_PASSED '+mode,flush=True)

if __name__=='__main__':main()
