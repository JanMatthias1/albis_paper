"""Descriptive RCTD error panels; native simulator-specific truth, no refitting."""
import json
from pathlib import Path
import numpy as np,pandas as pd,anndata as ad
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
ROOT=Path(__file__).resolve().parents[5]
PAPER=ROOT/'sim_paper'
BASE=PAPER/'data/comparison_methods/spatial_deconvolution'
OUT=BASE/'albis_spider_error_comparison'
TYPES=[f'type{i}' for i in range(1,9)]
COLORS=['#6D87DA','#DE8B5A','#79BB70','#A881CF','#D86596','#8D9C63','#D3B747','#63B7BE']
def main():
 OUT.mkdir(exist_ok=True)
 runs=[('ALBIS · 556 genes · molecule-fraction truth',PAPER/'data/figure_3/spatial_deconvolution/RCTD/strong_mix'),('SPIDER · 556 genes · contributing-cell-fraction truth',BASE/'spider600k_genes556_seed2025/results'),('SPIDER · 2,000 genes · contributing-cell-fraction truth',BASE/'spider600k_genes2000_seed2025/results')]
 fig,axes=plt.subplots(3,1,figsize=(12,11),sharex=True,sharey=True)
 summaries=[];rows=[]
 for ax,(label,run) in zip(axes,runs):
  summary=json.loads((run/'metrics_summary.json').read_text());source=Path(summary['spot_h5ad'])
  if not source.exists():source=Path(str(source).replace('/figure_3/cellbin_batch_sigma_slide/','/figure_3/strong_domain_mix/batch_sigma_slide/'))
  a=ad.read_h5ad(source,backed='r');pred=pd.read_csv(run/'estimated_fractions_wide.csv')
  idx=pred.spot_id.str.removeprefix('spot_').astype(int).to_numpy()-1
  assert pred.spot_id.is_unique and (idx>=0).all() and (idx<a.n_obs).all()
  true=np.asarray(a.obsm['cell_type_frac_true'])[idx];est=pred[TYPES].to_numpy();a.file.close()
  assert true.shape==est.shape and np.allclose(true.sum(1),1,atol=1e-5) and np.allclose(est.sum(1),1,atol=1e-5)
  error=np.abs(est-true);rmse=float(np.sqrt(np.mean(error**2)))
  assert abs(rmse-summary['overall_rmse'])<.0001
  values=[error[:,i] for i in range(8)]
  vp=ax.violinplot(values,positions=np.arange(8),widths=.8,showextrema=False)
  for i,body in enumerate(vp['bodies']):body.set_facecolor(COLORS[i]);body.set_edgecolor(COLORS[i]);body.set_alpha(.55)
  ax.boxplot(values,positions=np.arange(8),widths=.13,showfliers=False,patch_artist=True,manage_ticks=False,boxprops={'facecolor':'white','edgecolor':'#444444'},medianprops={'color':'#333333'},whiskerprops={'color':'#444444'},capprops={'color':'#444444'})
  for j,t in enumerate(TYPES):
   r=float(np.corrcoef(true[:,j],est[:,j])[0,1]);e=float(np.sqrt(np.mean(error[:,j]**2)))
   ax.text(j,1.10,f'r = {r:.2f}',ha='center',fontsize=9);ax.text(j,1.025,f'RMSE = {e:.3f}',ha='center',fontsize=8)
   rows.append({'panel':label,'cell_type':t,'n_spots':len(idx),'pearson_r':r,'rmse':e,'median_absolute_error':float(np.median(error[:,j]))})
  ax.set_title(f'{label}   (n = {len(idx):,} scored spots)',loc='left',fontsize=12,pad=10)
  ax.set_ylim(0,1.18);ax.set_yticks([0,.25,.5,.75,1]);ax.set_ylabel('Spot absolute error');ax.grid(axis='y',alpha=.2);ax.set_axisbelow(True)
  ax.spines[['top','right']].set_visible(False)
  summaries.append({'panel':label,'run':str(run),'truth_file':str(source.resolve()),'recorded_truth_file':summary['spot_h5ad'],'n_scored':len(idx),'overall_rmse':rmse})
 axes[-1].set_xticks(range(8),[f'Type {i}' for i in range(1,9)])
 fig.suptitle('RCTD deconvolution errors · Seed 2025',fontsize=16,y=.995)
 fig.text(.5,.013,'All scored spots across slices · Full error distributions, no percentile clipping\nCell-type indices are simulator-specific; fraction targets and reference settings differ.',ha='center',fontsize=9)
 fig.tight_layout(rect=(0,.055,1,.975))
 for ext in ['png','pdf','svg']:fig.savefig(OUT/f'rctd_error_comparison_seed2025.{ext}',dpi=230,bbox_inches='tight')
 pd.DataFrame(rows).to_csv(OUT/'per_celltype_metrics.csv',index=False)
 (OUT/'provenance.json').write_text(json.dumps({'panels':summaries,'error':'absolute(estimated fraction - native truth fraction)','violin_clipping':False,'shared_y_limits':[0,1.18],'seed':2025,'limitations':['ALBIS molecule fractions versus SPIDER contributing-cell fractions.','Cell-type labels are not matched biological identities across simulators.','Different tissue sizes, expression models, batch settings and reference designs.','This is a descriptive comparison, not a controlled simulator ranking.']},indent=2)+'\n')
 print(OUT/'rctd_error_comparison_seed2025.png')
if __name__=='__main__':main()
