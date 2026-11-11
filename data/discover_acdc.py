"""List the official public Girder ACDC collection; download no image contents."""
import concurrent.futures
import json
from pathlib import Path
import time
import urllib.parse
import urllib.request

API='https://humanheart-project.creatis.insa-lyon.fr/database/api/v1/'
PROXY='http://127.0.0.1:39843'
WORK=Path.home()/'xiaoe-internship-data/_download_expansion'
ROOT='637218e573e9f0047faa00fc'


def request(endpoint):
    for attempt in range(5):
        try:
            opener=urllib.request.build_opener(urllib.request.ProxyHandler({'http':PROXY,'https':PROXY}))
            with opener.open(API+endpoint,timeout=30) as r:
                return json.load(r)
        except Exception:
            if attempt==4:
                raise
            time.sleep(2**attempt)


def main():
    snapshots=[]
    todo=[('benchmarks/competitions/acdc2017',ROOT)]
    splits=request('folder?'+urllib.parse.urlencode(dict(parentType='folder',parentId=ROOT,limit=300)))
    snapshots.append(dict(endpoint='root folders',data=splits))
    for split in splits:
        prefix='benchmarks/competitions/acdc2017/'+split['name']
        todo.append((prefix,split['_id']))
        patients=request('folder?'+urllib.parse.urlencode(dict(parentType='folder',parentId=split['_id'],limit=300)))
        snapshots.append(dict(endpoint=prefix+' folders',data=patients))
        todo.extend((prefix+'/'+p['name'],p['_id']) for p in patients)

    def list_items(item):
        prefix,folder_id=item
        data=request('item?'+urllib.parse.urlencode(dict(folderId=folder_id,limit=300)))
        assets=[]
        for obj in data:
            name=obj['name']
            if name.endswith('.ipynb') or name.startswith('STATE-OF-THE-ART'):
                continue
            if '/patient' in prefix and name=='MANDATORY_CITATION.md':
                continue
            assets.append(dict(dataset='ACDC2017',path=prefix+'/'+name,bytes=obj['size'],
                               url=API+'item/'+obj['_id']+'/download',proxy=True))
        return dict(endpoint=prefix,data=data),assets

    assets=[]
    with concurrent.futures.ThreadPoolExecutor(8) as pool:
        for snapshot,rows in pool.map(list_items,todo):
            snapshots.append(snapshot)
            assets.extend(rows)
    WORK.mkdir(exist_ok=True)
    (WORK/'acdc_discovery.json').write_text(json.dumps(snapshots,indent=2))
    (WORK/'acdc_assets.json').write_text(json.dumps(assets,indent=2))
    print(json.dumps(dict(folders=len(todo),files=len(assets),bytes=sum(a['bytes'] for a in assets))))


if __name__=='__main__':
    main()
