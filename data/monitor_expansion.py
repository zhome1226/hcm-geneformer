"""Write factual progress for this one-off download batch until all jobs finish."""
import csv
import hashlib
import json
from pathlib import Path
import time
from datetime import datetime, timezone


def main():
    root=Path.home()/'xiaoe-internship-data'
    work=root/'_download_expansion'
    names=['vcc','geo138262','geo2021_203274','causalbench','neurips2023','acdc']
    while True:
        rows=[]
        for name in names:
            path=work/(name+'_assets.json')
            if not path.exists():
                continue
            for a in json.loads(path.read_text(encoding='utf-8-sig')):
                target=root/a['path']
                current=target.stat().st_size if target.exists() else 0
                record=work/(hashlib.sha256(a['path'].encode()).hexdigest()[:20]+'.verified.json')
                state='verified' if record.exists() else 'downloading' if current else 'not_started'
                rows.append(dict(dataset=a['dataset'],relative_path=a['path'],expected_bytes=a['bytes'],
                                 current_bytes=current,status=state,source_url=a['url']))
        tmp=work/'STATUS.tsv.tmp'
        with tmp.open('w',newline='') as out:
            writer=csv.DictWriter(out,fieldnames=['dataset','relative_path','expected_bytes','current_bytes','status','source_url'],delimiter='\t')
            writer.writeheader();writer.writerows(rows)
        tmp.replace(work/'STATUS.tsv')
        summaries=[]
        for name in names:
            p=work/(name+'_assets.summary.json')
            if p.exists():summaries.append(json.loads(p.read_text()))
        state=dict(updated_utc=datetime.now(timezone.utc).isoformat(),files_planned=len(rows),
                   files_verified=sum(r['status']=='verified' for r in rows),
                   planned_bytes=sum(r['expected_bytes'] for r in rows),
                   current_bytes=sum(min(r['current_bytes'],r['expected_bytes']) for r in rows),
                   jobs_finished=len(summaries),errors=[e for s in summaries for e in s['errors']])
        state['status']=('complete' if len(summaries)==len(names) and not state['errors'] else
                         'partial_failure' if len(summaries)==len(names) else 'downloading')
        (work/'STATUS.json').write_text(json.dumps(state,indent=2)+'\n')
        print(json.dumps(state),flush=True)
        if state['status']!='downloading':
            break
        time.sleep(60)


if __name__=='__main__':
    main()
