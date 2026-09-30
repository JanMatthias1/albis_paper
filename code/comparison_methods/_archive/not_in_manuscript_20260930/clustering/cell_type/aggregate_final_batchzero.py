import json
from pathlib import Path
import numpy as np
import pandas as pd
OUT=Path('/dcs04/hicks/data/Jan/sim_project/sim_paper/data/figure_5/clustering/cell_type/final_600k_batchzero_noharmony_20260928')
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
 # Reuse the Figure3 purple/teal/orange palette, mapped here to simulator.
 colors={'albis':'#8E63C7','sccube':'#43B7A5','spider':'#F2A65A'}
 plt.rcParams.update({'font.family':'DejaVu Sans','font.size':12,
                      'axes.titlesize':17,'axes.labelsize':16,
                      'xtick.labelsize':12,'ytick.labelsize':12,
                      'pdf.fonttype':42,'svg.fonttype':'none'})
 fig,ax=plt.subplots(figsize=(10,5.6))
 values=[frame[(frame.method==m)&(frame.genes==g)].cell_type_ari.to_numpy() for m,g in conditions]
 means=np.array([v.mean() for v in values])
 sds=np.array([v.std(ddof=1) for v in values])
 ax.bar(np.arange(1,6),means,width=.5,color=[colors[m] for m,g in conditions],edgecolor='white',linewidth=1.2,zorder=2)
 ax.errorbar(np.arange(1,6),means,yerr=sds,fmt='none',ecolor='#000000',elinewidth=1.1,capsize=3,zorder=4)
 for x,mean,sd,v in zip(range(1,6),means,sds,values):
  ax.text(x,max(mean+sd,v.max())+.02,f'{mean:.2f}',ha='center',va='bottom',fontsize=10,fontweight='bold')
 for j,seed in enumerate([2025,101,202]):
  vals=[frame[(frame.method==m)&(frame.genes==g)&(frame.seed==seed)].cell_type_ari.item() for m,g in conditions]
  ax.scatter(np.arange(1,6)+(j-1)*.08,vals,color='#3d3d3d',s=20,
             edgecolors='white',linewidths=.5,zorder=5)
 ax.set(xticks=range(1,6),xticklabels=['ALBIS\n556 genes','scCube\n556 genes','scCube\n2,000 genes','SPIDER\n556 genes','SPIDER\n2,000 genes'],
        ylabel='Cell-type ARI',ylim=(-.05,1.08),yticks=np.arange(0,1.01,.2),
        title='Cell-type recovery · No Harmony')
 ax.title.set_fontweight('bold')
 for side in ['top','right','left']:ax.spines[side].set_visible(False)
 ax.spines['bottom'].set_color('#c3c2b7')
 ax.tick_params(colors='#3d3d3d')
 ax.yaxis.grid(True,color='#DDDDDD',linewidth=.9)
 ax.set_axisbelow(True)
 ax.axhline(0,color='#c3c2b7',linewidth=1)
 fig.text(.5,.025,'600k cells per method/seed · ALBIS batch σ = 0 · Gaussian K256\nBars: mean ± sample SD · Dots: simulation seeds 2025, 101, 202',ha='center',fontsize=11)
 fig.tight_layout(rect=(0,.09,1,1))
 for stem in ['cell_type_ari_barplot', 'cell_type_ari_boxplot']:
  for ext in ['png','pdf','svg']:fig.savefig(OUT/f'{stem}.{ext}',dpi=300,bbox_inches='tight',facecolor='white')
 plt.close(fig)
 (OUT/'plot_style.json').write_text(json.dumps({'method_colors':colors,'reference':'Figure3 condensed purple/teal/orange palette; colors identify methods in this figure','bars':'Mean ARI; error bars are sample SD across three seeds', 'legacy_filename':'cell_type_ari_boxplot is an alias of the barplot for existing links','points':'Three simulation seeds; dark gray with white edges'},indent=2)+'\n')

if __name__ == "__main__":
 aggregate()
