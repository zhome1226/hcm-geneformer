"""Cross-check official supplementary DOCX against original HDF5 metadata.

This creates source evidence only; it never changes counts, resolves an age
disagreement by guessing, or declares the anonymous gene axis restored.
"""
from datetime import datetime, timezone
from pathlib import Path
from zipfile import ZipFile
from xml.etree import ElementTree as ET
import hashlib
import json
import anndata as ad
import h5py
import numpy as np
import pandas as pd

ROOT = Path.home() / "xiaoe-internship-data"
REF = ROOT / "references/GSE227734_publisher_20261111"
WORK = ROOT / "_processing/all_data_cleaning_20261111/GSE227734"
OUT = ROOT / "processed/cox7b_GSE227734_qc_20261111_pending_annotations/source_evidence"
OUT.mkdir(exist_ok=True)
source = REF / "mmc1.docx"
with source.open("rb") as stream:
    source_sha = hashlib.file_digest(stream, "sha256").hexdigest()
assert source_sha == "6b6ef2ca9c0737b0f300fa87244bf9ba9acd9bcd048839e8951b046bc2c5de54"
ns = {"w":"http://schemas.openxmlformats.org/wordprocessingml/2006/main"}
with ZipFile(source) as archive:
    doc = ET.fromstring(archive.read("word/document.xml"))
    tables = [[["".join(t.text or "" for t in cell.findall(".//w:t", ns)) for cell in row.findall("w:tc", ns)] for row in tab.findall("w:tr", ns)] for tab in doc.findall(".//w:tbl", ns)]
    embedded = [n for n in archive.namelist() if "/embeddings/" in n]
assert tables == json.loads((REF / "SUPPLEMENT_TABLES.json").read_text())
barcodes = {row[1] for row in tables[0][1:] if row[0].startswith("Barcode-")}
assert len(barcodes) == 10 and all(len(b)==8 and set(b)<={"A","C","G","T"} for b in barcodes)
clinical = tables[2]
assert clinical[0] == ["", "P1", "P2", "P3", "P4", "P5"]
clinical_rows = {r[0]:r[1:] for r in clinical[1:]}
age_key = next(k for k in clinical_rows if k.startswith("Age at surgery"))
library = tables[3]
assert library[0] == ["Sample ID","P1","P2","P3","P4","P5","Total"]
pre = list(map(int, library[1][1:]))
post = list(map(int, library[2][1:]))
with h5py.File(WORK / "GSE227734_hcm_raw.h5ad") as f:
    obs=ad.io.read_elem(f["obs"])
    var=ad.io.read_elem(f["var"])
    groups=list(f.keys())
assert set(obs.Barcode_sequence.astype(str)) == barcodes
assert len(obs)==post[-1]==sum(post[:-1])==338
assert pre[-1]==sum(pre[:-1])==390
mapping = {"P1":"CZQ","P2":"ZZW","P3":"GCZ","P4":"SCL","P5":"LGH"}
rows=[]
for i,(title,alias) in enumerate(mapping.items()):
    part=obs.loc[obs["sample"].astype(str).eq(alias)]
    ages=sorted(map(int,part.age.unique()))
    sex=sorted(part.gender.astype(str).unique())
    assert len(part)==post[i] and len(ages)==1
    assert sex == [clinical_rows["Sex"][i]] == ["male"]
    supplemental_age=int(clinical_rows[age_key][i])
    rows.append({"source_title":title,"source_sample_alias":alias,"pre_filter_cells_supplement":pre[i],"post_filter_cells_supplement":post[i],"published_cells":len(part),"age_years_original_GEO":ages[0],"age_years_supplement":supplemental_age,"age_conflict":ages[0]!=supplemental_age,"sex_original_and_supplement":"male","age_resolution":"unresolved_keep_both_sources" if ages[0]!=supplemental_age else "sources_agree"})
summary=pd.DataFrame(rows)
assert summary.age_conflict.sum()==1
assert summary.loc[summary.age_conflict,"source_title"].tolist()==["P5"]
checks={"supplement_sha256_matches_local_source":True,"tables_reparsed_independently_match_local":True,"all_ten_published_barcodes_match_supplement":True,"five_donor_postfilter_counts_match_original":True,"published_input_already_author_filtered_390_to_338":True,"all_five_sexes_agree":True,"P5_age_conflict_preserved_63_GEO_64_supplement":bool(rows[-1]["age_years_original_GEO"]==63 and rows[-1]["age_years_supplement"]==64),"no_hidden_gene_mapping_in_original_HDF5":groups==["X","obs","var"] and var.shape==(9760,0) and np.array_equal(var.index.astype(str),np.arange(9760).astype(str)),"supplement_has_six_small_tables_and_no_embedded_files":len(tables)==6 and max(map(len,tables))==21 and not embedded}
summary.to_csv(OUT/"SUPPLEMENT_DONOR_CROSSCHECK.tsv",sep="\t",index=False)
result={"checked_utc":datetime.now(timezone.utc).isoformat(),"source_url":"https://ars.els-cdn.com/content/image/1-s2.0-S0022282823001876-mmc1.docx","source_path":str(source),"source_bytes":source.stat().st_size,"source_sha256":source_sha,"checks":checks,"source_audit_passed":all(checks.values()),"dataset_cleaning_complete":False,"dataset_status":"blocked_missing_gene_annotation","training_eligible_cells":0,"tables_inspected":len(tables),"table_dimensions":[[len(t),max(map(len,t))] for t in tables],"embedded_files":embedded,"ordered_9760_feature_mapping_found":False,"raw_reconstruction_constraints":"Table S1 provides ten 8-base cell barcodes and primers, but not an authoritative raw-library/author-cell key, complete UMI/read-layout specification or exact reference annotation release; do not guess these from run ordering.","age_conflicts":summary.loc[summary.age_conflict].to_dict(orient="records"),"original_data_modified":False,"audit_source_sha256":hashlib.sha256(Path(__file__).read_bytes()).hexdigest()}
(OUT/"SUPPLEMENT_SOURCE_AUDIT.json").write_text(json.dumps(result,indent=2)+"\n")
assert result["source_audit_passed"],checks
print(json.dumps({"source_audit":"passed","checks":len(checks),"published_cells":338,"source_age_conflicts":1,"gene_mapping_restored":False,"dataset_cleaning_complete":False}))
