"""Freeze the best exact-six-cluster lambda on strong-mix seed2025."""
import argparse
import json
from pathlib import Path
import pandas as pd
from common import *

def main():
    p=argparse.ArgumentParser();p.add_argument('--root',type=Path,default=BASE);root=p.parse_args().root
    rows=[];missing=[]
    for task in json.loads((root/'sweep_tasks.json').read_text()):
        path=Path(task['directory'])/'metrics.json'
        if not path.exists():missing.append(task['task']);continue
        m=json.loads(path.read_text());assert not m['harmony_applied']
        rows.append(dict(modality=task['modality'],**{'lambda':task['lambda']},domain_ari=m['domain_ari'],slice_ari=m['slice_ari'],achieved_clusters=m['achieved_clusters'],target_clusters=m['target_clusters'],path=str(path)))
    pd.DataFrame(rows).to_csv(root/'summary/sweep_metrics.csv',index=False)
    if missing:raise RuntimeError(f'Incomplete sweep tasks: {missing}; no parameters selected')
    choices={}
    for modality in MODALITIES:
        eligible=[r for r in rows if r['modality']==modality and r['achieved_clusters']==6]
        if not eligible:raise RuntimeError(f'No eligible six-domain run for {modality}; inspect the sweep')
        best=sorted(eligible,key=lambda r:(-r['domain_ari'],r['lambda']))[0]
        choices[modality]=dict(best,k_geom=K_GEOM[modality])
    result=dict(selection_seed=2025,selection_mix='strong',criterion='Highest domain ARI among exact-six-cluster runs; exact ties choose smaller lambda.',modalities=choices)
    path=root/'selected_parameters.json'
    if path.exists():assert json.loads(path.read_text())==result,'Refusing to change frozen selection'
    else:save(path,result)
    print(json.dumps(result,indent=2),flush=True)

if __name__=='__main__':main()
