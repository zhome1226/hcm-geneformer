"""Export analysis evidence and receipts; original cell and pseudobulk matrices stay on Linux."""
import argparse
import hashlib
import json
import tarfile
from datetime import datetime, timezone
from pathlib import Path

p = argparse.ArgumentParser()
p.add_argument('work', type=Path)
w = p.parse_args().work.resolve()
assert json.loads((w/'audit/FINAL_AUDIT.json').read_text())['status'] == 'passed'
for tag in ['fit_d', 'lock', 'pb_s', 'coverage_s', 'fit_s', 'integrate', 'audit']:
    assert (w/f'logs/{tag}.exit').read_text().strip() == '0', tag
skip = {f'pb/{n}.tsv.gz' for n in ['scp_raw', 'scp_cb', 'scp_hhatl', 'liu', 'chin', 'whole']}
paths = []
for d in ['code', 'pb', 'stats', 'lock', 'coverage', 'final', 'audit', 'logs']:
    for path in sorted((w/d).glob('*')):
        name = path.relative_to(w).as_posix()
        if path.is_file() and not path.is_symlink() and name not in skip:
            if d == 'logs' and path.stem.startswith('export'):
                continue
            assert path.suffix in {'.tsv', '.gz', '.json', '.sha256', '.txt', '.md', '.py', '.R', '.log', '.exit'}, name
            paths.append(path)
records = [{'path': str(path.relative_to(w)), 'bytes': path.stat().st_size,
            'sha256': hashlib.sha256(path.read_bytes()).hexdigest()} for path in paths]
out = w/'export'
out.mkdir(exist_ok=True)
manifest = out/'EVIDENCE_MANIFEST.json'
assert not manifest.exists(), 'Do not silently replace a completed export'
manifest.write_text(json.dumps({'created_utc': datetime.now(timezone.utc).isoformat(),
    'server_workspace': str(w), 'excludes': sorted(skip), 'files': records},
    indent=2)+'\n', encoding='utf8')
archive = out/'evidence.tar.gz'
with tarfile.open(archive, 'x:gz') as tf:
    for path in paths:
        tf.add(path, arcname=str(path.relative_to(w)), recursive=False)
    tf.add(manifest, arcname='EVIDENCE_MANIFEST.json', recursive=False)
receipt = {'file': archive.name, 'bytes': archive.stat().st_size,
    'sha256': hashlib.sha256(archive.read_bytes()).hexdigest(), 'files': len(records),
    'manifest_sha256': hashlib.sha256(manifest.read_bytes()).hexdigest()}
(out/'EXPORT_RECEIPT.json').write_text(json.dumps(receipt, indent=2)+'\n', encoding='utf8')
print(json.dumps(receipt))
