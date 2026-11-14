"""Run the unmodified official TranscriptomeTokenizer, V1, separately per study."""
import argparse, hashlib, importlib.metadata, inspect, json, os, sys
from datetime import datetime, timezone
from pathlib import Path
import anndata as ad
import numpy as np
from geneformer import TranscriptomeTokenizer
ad.settings.allow_write_nullable_strings=True

W=Path('/home/zhome/xiaoe-internship-data/_processing/geneformer_v1_tokenization_20261112')
O=Path('/home/zhome/hcm-geneformer-data/tokenized/geneformer_v1_20261114')
ATTRS=['cell_id','source_cell_id','barcode_original','donor_id','sample_id','study','dataset',
       'disease','source_disease','label','cell_type','split','split_status','source_split',
       'source_split_status','n_counts','age_stratum','age_years','sex','technology','tissue','clinical_stage','training_use','reference_domain','donor_independence','count_processing_history','study_origin','original_source_cell_id','source_accession','donor_original','capture_id']

def main():
    mode=argparse.ArgumentParser();mode.add_argument('--mode',choices=['smoke','full'],required=True)
    mode=mode.parse_args().mode
    prepared=json.loads((O/'prepared_inputs/PREPARATION_VERIFICATION.json').read_text())
    assert prepared['status']=='passed' and prepared['total_cells']>0
    source=Path(inspect.getsourcefile(TranscriptomeTokenizer))
    expected=json.loads((W/'LINUX_SOURCE_ASSETS.json').read_text())
    expected=next(i['sha256'] for i in expected if i['path']=='geneformer/tokenizer.py')
    assert hashlib.sha256(source.read_bytes()).hexdigest()==expected
    environment={'utc':datetime.now(timezone.utc).isoformat(),'executable':sys.executable,
                 'tokenizer_source':str(source),'tokenizer_sha256':expected,
                 'official_revision':'1f7fbae4e469a5f4f1af8c111a529cfe1b3829f5',
                 'versions':{p:importlib.metadata.version(p) for p in ['geneformer','numpy','scipy','pandas','anndata','scanpy','datasets','torch','transformers']}}
    (O/'ENCODING_ENVIRONMENT.json').write_text(json.dumps(environment,indent=2)+'\n')
    out=O/('smoke_datasets' if mode=='smoke' else 'datasets');out.mkdir(exist_ok=True)
    for r in prepared['datasets']:
        ds=r['dataset'];inp=Path(r['prepared_path'])
        a=ad.read_h5ad(inp,backed='r');attrs={k:k for k in ATTRS if k in a.obs}
        assert {'cell_id','donor_id','sample_id','disease','label'}<=set(attrs)
        if mode=='smoke':
            sample=O/'smoke_inputs';sample.mkdir(exist_ok=True)
            inds=np.linspace(0,a.n_obs-1,min(32,a.n_obs),dtype=int)
            copy=a[inds].to_memory();copy.write_h5ad(sample/(ds+'.h5ad'))
            inp=sample/(ds+'.h5ad')
        a.file.close()
        destination=out/(ds+'.dataset')
        assert not destination.exists(),str(destination)+' exists; inspect before rerunning'
        tk=TranscriptomeTokenizer(custom_attr_name_dict=attrs,model_version='V1',nproc=1,chunk_size=256)
        assert tk.model_input_size==2048 and tk.special_token is False
        assert len(tk.gene_token_dict)==25426
        print(json.dumps({'mode':mode,'dataset':ds,'input':str(inp),'output':str(destination)}),flush=True)
        # input_identifier excludes all other studies in prepared_inputs.
        tk.tokenize_data(data_directory=inp.parent,output_directory=out,output_prefix=ds,
                         file_format='h5ad',input_identifier=ds)
    print('OFFICIAL_TOKENIZATION_FINISHED '+mode,flush=True)

if __name__=='__main__':main()
