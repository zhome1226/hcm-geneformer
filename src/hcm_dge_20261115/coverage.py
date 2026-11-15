"""Count detected expression and actual V1 sequence presence per source donor.

V1 environment: python coverage.py --work WORK --phase discovery|support
Discovery reads only SCP1303 counts and actual V1 sequences. Support requires
WORK/lock/CANDIDATES.tsv (gene_id) and restricts output to the frozen IDs.
SCP expression uses Cell Ranger raw; all other studies use source counts, with
the Ensembl mappings used for encoding and native missing-feature masks preserved.
"""
import argparse
import hashlib
import json
import pickle
import re
from pathlib import Path

import h5py
import numpy as np
import pandas as pd
from scipy import sparse
from anndata._io.specs import read_elem
from datasets import load_from_disk

OLD = Path('/home/zhome/xiaoe-internship-data')
NEW = Path('/home/zhome/hcm-geneformer-data')
OLD_TOKEN = OLD / 'tokenized/geneformer_v1_20261112'
NEW_TOKEN = NEW / 'tokenized/geneformer_v1_20261114'
ROOTS = {**dict.fromkeys(('SCP1303', 'liu2023', 'chin_GSE174691', 'chin_GSE161921'), OLD_TOKEN),
         **dict.fromkeys(('GSE255296', 'GSE292754', 'GSE138262_whole_cell'), NEW_TOKEN)}
VOCAB = OLD / '_processing/geneformer_v1_tokenization_20261112/official_source/geneformer/gene_dictionaries_30m/token_dictionary_gc30M.pkl'
VOCAB_SHA = 'ab9dc40973fa5224d77b793e2fd114cacf3d08423ed9c4c49caf0ba9c7f218f1'
MASKS = [('feature_available_by_capture', 'feature_availability_captures', 'capture_id'),
         ('native_feature_presence_by_sample', 'feature_presence_sample_order', 'sample_id')]


def digest(p):
    with Path(p).open('rb') as handle:
        return hashlib.file_digest(handle, 'sha256').hexdigest()


def valid_gene(g):
    return bool(re.fullmatch(r'ENSG[0-9]+', str(g)))


def table_out(df, p):
    df.to_csv(p, sep='\t', index=False, na_rep='NA', float_format='%.10g')


def fr(num, den):
    num, den = np.broadcast_arrays(np.asarray(num), np.asarray(den))
    return np.divide(num, den, out=np.full(num.shape, np.nan), where=den > 0)


def count(ds, candidates, vocab, keep_groups):
    root = ROOTS[ds]
    records = json.loads((root/'prepared_inputs/PREPARATION_VERIFICATION.json').read_text())
    rec = next(r for r in records['datasets'] if r['dataset'] == ds)
    source = Path(rec['source'])
    prepared = Path(rec['prepared_path'])
    assert source.stat().st_size == rec['source_bytes'] and digest(source) == rec['source_sha256']
    assert prepared.stat().st_size == rec['prepared_bytes'] and digest(prepared) == rec['prepared_sha256']
    dpath = root/'datasets'/(ds+'.dataset')
    man = pd.read_csv(root/'FULL_MANIFEST.tsv', sep='\t')
    derived = man[man.dataset == ds]
    assert {p.relative_to(root).as_posix() for p in dpath.iterdir() if p.is_file()} == set(derived.path)
    for r in derived.itertuples():
        p = root/r.path
        assert p.stat().st_size == r.bytes and digest(p) == r.sha256
    d = load_from_disk(str(dpath))
    pa = d.data.table
    with h5py.File(prepared) as fp, h5py.File(source) as f:
        features = read_elem(fp['var'])
        origin = read_elem(f['var'])
        assert np.array_equal(features.index, origin.index)
        ens = features.ensembl_id.astype(str).to_numpy()
        obs = read_elem(f['obs'])
        assert np.array_equal(pa['source_cell_id'].to_pylist(), obs.index.astype(str))
        assert np.array_equal(pa['donor_id'].to_pylist(), obs.donor_id.astype(str))
        assert np.array_equal(pa['sample_id'].to_pylist(), obs.sample_id.astype(str))
        assert np.array_equal(pa['label'].to_pylist(), obs.label.astype(int))
        assert np.array_equal(pa['dataset'].to_pylist(), np.repeat(ds, len(obs)))
        assert set(obs.label.astype(int)) <= {0, 1}
        assert (obs.groupby('donor_id',observed=True).label.nunique()==1).all()
        disease = np.where(obs.label.to_numpy(dtype=int)==1, 'HCM', 'Control')
        donors = obs.donor_id.astype(str).to_numpy()
        used = np.isin(disease, keep_groups)
        donorids = sorted(set(donors[used]))
        dlookup = {v:i for i,v in enumerate(donorids)}
        dindex = np.array([dlookup[x] if use else -1 for x,use in zip(donors, used)])
        ncell = np.bincount(dindex[used], minlength=len(donorids))
        labels = np.array([disease[donors == donor][0] for donor in donorids])
        genes = np.array(candidates if candidates is not None else list(dict.fromkeys(g for g in ens if valid_gene(g))))
        n = len(genes)
        positions = {g:i for i,g in enumerate(genes)}
        token_ids = np.array([vocab.get(g, -1) for g in genes])
        present_in_vocab = token_ids >= 0
        feature_pos = np.array([positions.get(g, -1) for g in ens])
        indices = np.flatnonzero(feature_pos >= 0)
        projection = sparse.csr_matrix((np.ones(len(indices),dtype=np.int64),(indices,feature_pos[indices])),shape=(len(ens),n))
        included = np.bincount(feature_pos[indices],minlength=n)>0
        nd = len(donorids)
        has = np.zeros((nd,n), dtype=np.int64)
        counts = np.zeros_like(has)
        tok = np.zeros_like(has)
        avail = np.zeros_like(has)
        masks = [m for m in MASKS if m[0] in f['varm']]
        assert len(masks)<=1
        if masks:
            mask, order, col = masks[0]
            avail_mask = np.asarray(read_elem(f['varm'][mask]), dtype=bool)
            units = [str(s) for s in read_elem(f['uns'][order])]
            assert avail_mask.shape==(len(ens),len(units))
            assert set(obs[col].astype(str))<=set(units)
            for j,unit in enumerate(units):
                native = np.unique(feature_pos[avail_mask[:,j] & (feature_pos>=0)])
                freq = np.bincount(dindex[(obs[col].astype(str).to_numpy()==unit)&used],minlength=nd)
                avail[:,native] += freq[:,None]
        else:
            avail[:,included] = ncell[:,None]

        layer = 'cellranger_raw' if ds=='SCP1303' else 'counts'
        enc = f['layers'][layer]
        assert enc.attrs['encoding-type']=='csr_matrix'
        ptr = enc['indptr'][:]
        assert len(ptr)==len(obs)+1
        for start in range(0,len(obs),512):
            end = min(len(obs),start+512)
            dx = dindex[start:end]
            if not (dx>=0).any():
                continue
            lo,hi = ptr[start],ptr[end]
            mat = sparse.csr_matrix((enc['data'][lo:hi],enc['indices'][lo:hi],ptr[start:end+1]-lo),shape=(end-start,len(ens)))
            mat.sum_duplicates();mat.eliminate_zeros()
            assert np.all(mat.data>=0) and np.all(mat.data==np.floor(mat.data))
            expr = (mat@projection).tocsr()
            for donoridx in np.unique(dx[dx>=0]):
                block = expr[dx==donoridx]
                counts[donoridx] += np.asarray(block.sum(axis=0)).ravel().astype(np.int64)
                has[donoridx] += np.bincount(block.indices,minlength=n)
        toklookup=np.full(max(vocab.values())+1,-1,dtype=np.int32)
        toklookup[token_ids[present_in_vocab]]=np.flatnonzero(present_in_vocab)
        step=0
        for b in pa.select(['input_ids','length']).to_batches(max_chunksize=512):
            seqs=b.column(0)
            lengths=b.column(1).to_numpy()
            ids=seqs.flatten().to_numpy()
            assert np.array_equal(np.diff(seqs.offsets.to_numpy()),lengths)
            assert (lengths<=2048).all() and (lengths>0).all() and (ids>=2).all() and (ids<len(toklookup)).all()
            dx=np.repeat(dindex[step:step+len(b)],lengths)
            gx=toklookup[ids]
            active=(dx>=0)&(gx>=0)
            if active.any():
                sums=np.bincount(dx[active]*n+gx[active],minlength=nd*n).reshape(nd,n)
                tok+=sums
            step+=len(b)
        assert step==len(obs)
        assert np.all(has<=avail) and np.all(tok<=avail) and np.all(counts>=has)
        if ds!='SCP1303':assert np.all(tok<=has)
        frames=[]
        for i,donor in enumerate(donorids):
            frames.append(pd.DataFrame({'dataset':ds,'disease':labels[i],'donor_id':donor,'gene_id':genes,'donor_cells':ncell[i],'native_cells':avail[i],
                                        'expressing_cells':has[i],'count_sum':counts[i],'expr_fraction_native':fr(has[i],avail[i]),
                                        'token_cells':tok[i],'token_fraction_all_cells':fr(tok[i],ncell[i]),'in_v1_vocab':present_in_vocab,
                                        'token_id':pd.array(np.where(present_in_vocab,token_ids,np.nan),dtype='Int64')}))
        detail=pd.concat(frames,ignore_index=True)
        summaries=[]
        for label in keep_groups:
            ix=labels==label;num=int(ix.sum());cells=int(ncell[ix].sum())
            if not num:continue
            a=avail[ix].sum(0);x=has[ix].sum(0);t=tok[ix].sum(0)
            donor_frac=fr(tok[ix],ncell[ix,None])
            summaries.append(pd.DataFrame({'dataset':ds,'disease':label,'gene_id':genes,'in_v1_vocab':present_in_vocab,
                 'token_id':pd.array(np.where(present_in_vocab,token_ids,np.nan),dtype='Int64'),
                 'donors':num,'cells':cells,'native_cells':a,'expressing_cells':x,'count_sum':counts[ix].sum(0),
                 'expr_fraction':fr(x,a),'donors_native':(avail[ix]>0).sum(0),'donors_expressing':(has[ix]>0).sum(0),
                 'token_cells':t,'token_fraction':fr(t,cells),'donors_token_positive':(tok[ix]>0).sum(0),
                 'donors_token_ge5pct':(donor_frac>=.05).sum(0),'donors_token_ge10pct':(donor_frac>=.10).sum(0),
                 'median_donor_token_fraction':np.median(donor_frac,axis=0),'count_layer':layer}))
        return pd.concat(summaries,ignore_index=True),detail


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--work',type=Path,required=True)
    parser.add_argument('--phase',choices=['discovery','support'],required=True)
    args=parser.parse_args()
    out=args.work/'coverage';out.mkdir(exist_ok=True,parents=True)
    assert digest(VOCAB)==VOCAB_SHA
    with VOCAB.open('rb') as f:vocab=pickle.load(f)
    assert len(vocab)==25426 and sorted(vocab.values())==list(range(25426))
    if args.phase=='discovery':
        aggregate,detail=count('SCP1303',None,vocab,['HCM'])
        aggregate.rename(columns={c:'hcm_'+c for c in aggregate if c not in ['gene_id','token_id','in_v1_vocab','dataset','disease','count_layer']},inplace=True)
        table_out(aggregate,out/'discovery.tsv')
        table_out(detail,out/'discovery_donors.tsv.gz')
        print(f'Discovery coverage: {len(aggregate)} genes, {detail.donor_id.nunique()} HCM donors.',flush=True)
    else:
        lock=args.work/'lock/CANDIDATES.tsv'
        signature=digest(lock)
        genes=pd.read_csv(lock,sep='\t',dtype=str,keep_default_na=False).gene_id.tolist()
        assert len(set(genes))==len(genes) and all(valid_gene(g) for g in genes)
        if not genes:raise ValueError('No locked genes: skip support, do not substitute candidates.')
        agg=[];det=[]
        for ds in ROOTS:
            a,d=count(ds,genes,vocab,['HCM','Control']);agg.append(a);det.append(d)
            print(f'{ds}: {len(genes)} locked genes, {d.donor_id.nunique()} donors.',flush=True)
        assert digest(lock)==signature
        ao=pd.concat(agg,ignore_index=True);det=pd.concat(det,ignore_index=True)
        table_out(ao,out/'support.tsv')
        table_out(det,out/'support_donors.tsv.gz')
        (out/'support.lock.sha256').write_text(signature+'  CANDIDATES.tsv\n')


if __name__=='__main__':
    main()
