"""Measure native expression and actual V1 sequence coverage for the complete gene universe.

Run in the existing V1 Python environment:
    python code/all_coverage.py --work WORK

The source count implementation is the byte-identical, previously checked
coverage_base.py. This wrapper changes scope and output streaming, not counting.
Every gene is retained at every study/group/donor; unmeasured expression stays NA.
"""
import argparse
import gc
import gzip
import hashlib
import importlib.util
import io
import json
import pickle
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import pandas as pd


BASE_SHA256 = 'edf30d660d793a1f89318d8fff90a3607cee233fae0112fc97997e6443246d29'
EXPECTED = {
    'SCP1303': (114711, 31, {'HCM': 15, 'Control': 16}),
    'liu2023': (8547, 12, {'HCM': 10, 'Control': 2}),
    'chin_GSE174691': (10150, 9, {'HCM': 9}),
    'chin_GSE161921': (3534, 4, {'Control': 4}),
    'GSE255296': (7516, 4, {'HCM': 4}),
    'GSE292754': (3306, 2, {'HCM': 2}),
    'GSE138262_whole_cell': (3595, 12, {'HCM': 5, 'Control': 7}),
}


def sha(path):
    with Path(path).open('rb') as handle:
        return hashlib.file_digest(handle, 'sha256').hexdigest()


def read(path):
    return pd.read_csv(path, sep='\t', keep_default_na=False, na_values=['NA'])


def compressed_writer(path):
    """Single gzip member per table, reproducible metadata and exclusive creation."""
    raw = path.open('xb')
    compressed = gzip.GzipFile(filename='', mode='wb', fileobj=raw, compresslevel=6, mtime=0)
    return raw, io.TextIOWrapper(compressed, encoding='utf-8', newline='')


def source_receipt(base, ds):
    root = base.ROOTS[ds]
    prep = root/'prepared_inputs/PREPARATION_VERIFICATION.json'
    preparation = json.loads(prep.read_text())
    record = next(r for r in preparation['datasets'] if r['dataset'] == ds)
    man_path = root/'FULL_MANIFEST.tsv'
    entries = read(man_path)
    entries = entries[entries.dataset.eq(ds)]
    assert len(entries) > 0
    files = []
    for prefix in ['source', 'prepared']:
        name = prefix if prefix == 'source' else 'prepared_path'
        path = Path(record[name])
        files.append({'path': str(path), 'bytes': int(record[prefix+'_bytes']),
                      'sha256': record[prefix+'_sha256']})
    for row in entries.itertuples():
        files.append({'path': str(root/row.path), 'bytes': int(row.bytes), 'sha256': row.sha256})
    return {'dataset': ds, 'files': files,
            'preparation_manifest': str(prep), 'preparation_manifest_sha256': sha(prep),
            'dataset_manifest': str(man_path), 'dataset_manifest_sha256': sha(man_path)}


def verify_study(aggregate, detail, ds, genes, roster):
    n_cells, n_donors, groups = EXPECTED[ds]
    assert set(aggregate.dataset) == set(detail.dataset) == {ds}
    assert len(aggregate) == len(genes)*len(groups)
    assert len(detail) == len(genes)*n_donors
    assert not aggregate.duplicated(['disease', 'gene_id']).any()
    assert not detail.duplicated(['donor_id', 'gene_id']).any()
    actual = detail.groupby('donor_id', sort=True).agg(
        records=('gene_id', 'size'), cells=('donor_cells', 'first'),
        min_cells=('donor_cells', 'min'), max_cells=('donor_cells', 'max'),
        disease=('disease', 'first'), labels=('disease', 'nunique'))
    expected = roster[roster.dataset.eq(ds)].set_index('donor_id').sort_index()
    assert actual.index.equals(expected.index) and len(expected) == n_donors
    assert actual.records.eq(len(genes)).all() and actual.labels.eq(1).all()
    assert actual.cells.equals(actual.min_cells) and actual.cells.equals(actual.max_cells)
    assert np.array_equal(actual.cells, expected.nuclei)
    assert np.array_equal(actual.disease, expected.disease)
    assert actual.cells.sum() == n_cells and actual.disease.value_counts().to_dict() == groups
    for _, block in aggregate.groupby('disease', sort=False):
        assert block.gene_id.tolist() == genes
    for _, block in detail.groupby('donor_id', sort=False):
        assert block.gene_id.tolist() == genes
    assert aggregate.expr_fraction.notna().equals(aggregate.native_cells.gt(0))
    assert detail.expr_fraction_native.notna().equals(detail.native_cells.gt(0))
    return {'dataset': ds, 'cells': n_cells, 'donors': n_donors,
            'disease_donors': groups, 'genes': len(genes),
            'aggregate_rows': len(aggregate), 'donor_gene_rows': len(detail)}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--work', type=Path, required=True)
    w = parser.parse_args().work.resolve()
    universe_path = w/'input/GENE_UNIVERSE.tsv'
    roster_path = w/'input/pb/DONOR_ROSTER.tsv'
    base_path = w/'code/coverage_base.py'
    assert sha(base_path) == BASE_SHA256, 'The frozen counting implementation changed'
    spec = importlib.util.spec_from_file_location('coverage_base', base_path)
    base = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(base)
    assert dict.fromkeys(base.ROOTS) == dict.fromkeys(EXPECTED)
    universe = read(universe_path)
    genes = universe.gene_id.tolist()
    assert genes and genes == sorted(set(genes)) and all(base.valid_gene(g) for g in genes)
    roster = read(roster_path)
    assert len(roster) == 74 and roster.nuclei.sum() == 151359
    assert not roster.duplicated(['dataset', 'donor_id']).any()
    inputs = {str(p.relative_to(w)): sha(p) for p in [universe_path, roster_path, base_path]}
    if (w/'code/DESIGN.md').exists():
        inputs['code/DESIGN.md'] = sha(w/'code/DESIGN.md')
    assert sha(base.VOCAB) == base.VOCAB_SHA
    with base.VOCAB.open('rb') as handle:
        vocab = pickle.load(handle)
    assert len(vocab) == 25426 and sorted(vocab.values()) == list(range(25426))
    out = w/'coverage'
    out.mkdir(parents=True, exist_ok=True)
    names = ['ALL_COVERAGE.tsv.gz', 'DONOR_COVERAGE.tsv.gz', 'COVERAGE_RECEIPT.json']
    assert not any((out/n).exists() for n in names), 'Do not replace completed coverage'
    sources = [source_receipt(base, ds) for ds in EXPECTED]
    a_partial, d_partial = [out/(n+'.partial') for n in names[:2]]
    a_raw, a_writer = compressed_writer(a_partial)
    d_raw, d_writer = compressed_writer(d_partial)
    records = []
    try:
        for i, ds in enumerate(EXPECTED):
            aggregate, detail = base.count(ds, genes, vocab, ['HCM', 'Control'])
            records.append(verify_study(aggregate, detail, ds, genes, roster))
            for frame, writer in [(aggregate, a_writer), (detail, d_writer)]:
                frame.to_csv(writer, sep='\t', index=False, header=i == 0,
                             na_rep='NA', float_format='%.10g')
                writer.flush()
            print(json.dumps(records[-1]), flush=True)
            del aggregate, detail
            gc.collect()
    finally:
        a_writer.close()
        d_writer.close()
        a_raw.close()
        d_raw.close()
    # Verify inputs again after counting; a manifest is not a substitute for file hashes.
    for relative, digest in inputs.items():
        assert sha(w/relative) == digest, relative+' changed while counting'
    for study in sources:
        for item in study['files']:
            p = Path(item['path'])
            assert p.stat().st_size == item['bytes'] and sha(p) == item['sha256'], str(p)
        for key in ['preparation_manifest', 'dataset_manifest']:
            assert sha(study[key]) == study[key+'_sha256']
    assert sum(r['aggregate_rows'] for r in records) == 10*len(genes)
    assert sum(r['donor_gene_rows'] for r in records) == 74*len(genes)
    for temp, name in zip([a_partial, d_partial], names[:2]):
        temp.rename(out/name)
    receipt = {
        'status': 'passed', 'completed_utc': datetime.now(timezone.utc).isoformat(),
        'scope': 'all genes in frozen union; no previous candidate filtering',
        'genes': len(genes), 'datasets': len(EXPECTED), 'cells': 151359,
        'donor_study_records': 74, 'aggregate_rows': 10*len(genes),
        'donor_gene_rows': 74*len(genes), 'inputs': inputs,
        'vocab': {'path': str(base.VOCAB), 'sha256': base.VOCAB_SHA},
        'studies': records, 'source_files_verified_before_and_after': sources,
        'outputs': {n: {'bytes': (out/n).stat().st_size, 'sha256': sha(out/n)} for n in names[:2]},
    }
    (out/names[2]).write_text(json.dumps(receipt, indent=2, allow_nan=False)+'\n', encoding='utf-8')
    print(json.dumps({k: receipt[k] for k in ['status', 'genes', 'datasets', 'cells', 'donor_gene_rows']}), flush=True)


if __name__ == '__main__':
    main()
