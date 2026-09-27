"""Summarize completed native-count RCTD runs; one point per simulation seed."""
import json
from pathlib import Path
import pandas as pd,numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
ROOT=Path(__file__).resolve().parents[5]
BASE=ROOT/'sim_paper/data/comparison_methods/spatial_deconvolution'
OUT=BASE/'spider600k_summary'
def main():
 OUT.mkdir(exist_ok=True);rows=[]
 for genes in [556,2000]:
  for seed in [2025,101,202]:
   run=BASE/f'spider600k_genes{genes}_seed{seed}'
   d=json.loads((run/'results/metrics_summary.json').read_text());audit=json.loads((run/'inputs/input_audit.json').read_text())
   assert audit['underlying_cells']==600000 and audit['genes']==genes
   assert d['n_spots_scored']==d['n_spots_total']
   rows.append(dict(genes=genes,seed=seed,pearson_r=d['overall_pearson_r'],rmse=d['overall_rmse'],mean_per_spot_r=d['mean_per_spot_pearson_r'],mean_per_spot_rmse=d['mean_per_spot_rmse'],spots_scored=d['n_spots_scored'],input_audit=str(run/'inputs/input_audit.json')))
 df=pd.DataFrame(rows);df.to_csv(OUT/'metrics_by_seed.csv',index=False)
 df.groupby('genes')[['pearson_r','rmse']].agg(['mean','std','min','max']).to_csv(OUT/'summary.csv')
 fig,axes=plt.subplots(1,2,figsize=(9,4.5))
 for ax,key,title in zip(axes,['pearson_r','rmse'],['Overall Pearson correlation','Overall RMSE']):
  ax.boxplot([df[df.genes==g][key] for g in [556,2000]],widths=.4,showfliers=False)
  for j,seed in enumerate([2025,101,202]):
   ax.scatter(np.array([1,2])+(j-1)*.07,[df[(df.genes==g)&(df.seed==seed)][key].item() for g in [556,2000]],label=f'Seed {seed}',zorder=3)
  ax.set(xticks=[1,2],xticklabels=['556 genes','2,000 genes'],ylabel=title)
  ax.spines[['top','right']].set_visible(False)
 axes[0].set_ylim(0,1.03);axes[1].set_ylim(bottom=0);axes[0].legend(frameon=False,fontsize=8)
 fig.suptitle('SPIDER · 600,000 cells · Native-count RCTD')
 fig.text(.5,.01,'Source-pool reference; not an independent-reference benchmark',ha='center',fontsize=9)
 fig.tight_layout(rect=(0,.04,1,.95))
 for ext in ['png','pdf','svg']:fig.savefig(OUT/f'rctd_seed_overview.{ext}',dpi=220)
 (OUT/'README.md').write_text('# SPIDER 600k RCTD overview\n\nSix completed runs: two native gene counts, three simulation seeds each. All supplied spots scored. Native integer expression and native circular aggregation; no rounding/rescaling. The original 10k Splatter pool is the reference, so reference/query are not independent. Truth is contributing-cell fractions; ALBIS molecule-fraction truth must not be silently pooled with these scores. scCube legacy rounded-log pilots are excluded.\n')
 print(df.to_string(index=False));print(OUT/'rctd_seed_overview.png')
if __name__=='__main__':main()
