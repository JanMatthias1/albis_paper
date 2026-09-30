import json
from pathlib import Path
import numpy as np
import pandas as pd
OUT=Path('/dcs04/hicks/data/Jan/sim_project/sim_paper/data/figure_5/clustering/cell_type/final_600k_noharmony_20260928')
def aggregate():
 import matplotlib
 matplotlib.use('Agg')
 import matplotlib.pyplot as plt
 rows=[];status=[]
 for r in json.loads((OUT/'runs.json').read_text()):
  p=Path(r['directory'])/'metrics.json';d=json.loads(p.read_text()) if p.exists() else None
  status.append(dict(r,status='missing' if d is None else ('complete' if d['achieved_clusters']==8 else 'cluster_count_mismatch')))
  if d:rows.append(d)
 pd.DataFrame(status).to_csv(OUT/'status.csv',index=False)
 if rows:pd.DataFrame(rows).drop(columns=['trials']).to_csv(OUT/'metrics.csv',index=False)
 if len(rows)!=15 or any(d['achieved_clusters']!=8 for d in rows):raise RuntimeError('Incomplete or count-mismatched trial; inspect status.csv')
 frame=pd.DataFrame(rows);conditions=[('albis',556),('sccube',556),('sccube',2000),('spider',556),('spider',2000)]
 fig,ax=plt.subplots(figsize=(9,5))
 ax.boxplot([frame[(frame.method==m)&(frame.genes==g)].cell_type_ari for m,g in conditions],showfliers=False)
 for j,seed in enumerate([2025,101,202]):
  vals=[frame[(frame.method==m)&(frame.genes==g)&(frame.seed==seed)].cell_type_ari.item() for m,g in conditions]
  ax.scatter(np.arange(1,6)+(j-1)*.08,vals,label=f'Seed {seed}',zorder=3)
 ax.set(xticks=range(1,6),xticklabels=['ALBIS\n556','scCube\n556','scCube\n2000','SPIDER\n556','SPIDER\n2000'],ylabel='Cell-type ARI',ylim=(-.05,1.05),title='No Harmony · Gaussian KNN (K=256)');ax.legend(frameon=False)
 fig.text(.5,.015,'Gene counts shown; 600k cells per method/seed; ALBIS batch σ=0.7; no Harmony',ha='center');fig.tight_layout(rect=(0,.04,1,1))
 for ext in ['png','pdf','svg']:fig.savefig(OUT/f'cell_type_ari_boxplot.{ext}',dpi=220)

if __name__ == "__main__":
 aggregate()
