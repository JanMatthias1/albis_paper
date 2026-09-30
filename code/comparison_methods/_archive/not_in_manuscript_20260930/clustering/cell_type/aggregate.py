"""Collect run status; draw the final boxplot only when every run succeeds."""
import json,os
from pathlib import Path
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
ROOT=Path(__file__).resolve().parents[5]
OUT=ROOT/'sim_paper/data/figure_5/clustering/cell_type/experiment_600k'
def main():
    analysis=os.environ.get('CELLTYPE_ANALYSIS','clustering')
    destination=OUT if analysis=='clustering' else OUT/analysis
    destination.mkdir(exist_ok=True)
    rows=[];status=[]
    for run in json.loads((OUT/'runs.json').read_text()):
        run_analysis='clustering_no_harmony' if analysis=='clustering' and run['method']!='albis' else analysis
        path=Path(run['directory'])/run_analysis/'metrics.json'
        metric=json.loads(path.read_text()) if path.exists() else None
        state='incomplete' if metric is None else ('complete' if metric['achieved_clusters']==metric['target_clusters'] else 'cluster_count_mismatch')
        status.append(dict(run,status=state))
        if metric is not None:rows.append(metric)
    pd.DataFrame(status).to_csv(destination/'status.csv',index=False)
    if rows:pd.DataFrame(rows).drop(columns=['trials','shared_code_sha256']).to_csv(destination/'metrics.csv',index=False)
    if len(rows)!=15:raise SystemExit(f'{len(rows)}/15 complete; see status.csv and per-run logs. Final plot not produced.')
    if any(r['achieved_clusters']!=r['target_clusters'] for r in rows):
        raise SystemExit('Cluster-count mismatch: metrics retained, final comparison plot withheld; inspect status.csv.')
    frame=pd.DataFrame(rows);conditions=[('albis',556),('sccube',556),('sccube',2000),('spider',556),('spider',2000)]
    fig,ax=plt.subplots(figsize=(8,4.5),layout='constrained')
    vals=[frame[(frame.method==m)&(frame.genes==g)].cell_type_ari.to_numpy() for m,g in conditions]
    ax.boxplot(vals,showfliers=False,widths=.5)
    for j,seed in enumerate([2025,101,202]):
        y=[frame[(frame.method==m)&(frame.genes==g)&(frame.simulation_seed==seed)].cell_type_ari.item() for m,g in conditions]
        ax.scatter([i+1+(j-1)*.09 for i in range(5)],y,label=f'Seed {seed}',zorder=3,s=35)
    ax.set_xticks(range(1,6),['ALBIS\n556 genes','scCube\n556 genes','scCube\n2,000 genes','SPIDER\n556 genes','SPIDER\n2,000 genes'])
    ax.set_ylabel('Cell-type ARI');ax.set_ylim(min(-.05,min(map(min,vals))-.03),1.03)
    ax.set_title('600,000 cells · 10 slices · PCA–Harmony–Leiden');ax.legend(frameon=False,fontsize=8)
    ax.spines[['top','right']].set_visible(False)
    for suffix in ['png','pdf','svg']:fig.savefig(destination/f'cell_type_ari_boxplot.{suffix}',dpi=300)
    print('Completed all 15 runs and comparison plot')
if __name__=='__main__':main()
