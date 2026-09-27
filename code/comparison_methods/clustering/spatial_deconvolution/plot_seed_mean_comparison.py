"""One-axis native-truth error comparison: mean and SD of three run-level MAEs."""
import json
from pathlib import Path
import numpy as np,pandas as pd,anndata as ad
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
ROOT=Path(__file__).resolve().parents[5]
PAPER=ROOT/'sim_paper';BASE=PAPER/'data/comparison_methods/spatial_deconvolution'
OUT=BASE/'albis_spider_error_comparison';TYPES=[f'type{i}' for i in range(1,9)];SEEDS=[2025,101,202]
CONDITIONS=[('ALBIS · 556 genes','#3572B0'),('SPIDER · 556 genes','#D98029'),('SPIDER · 2,000 genes','#7253A3')]
def main():
 OUT.mkdir(exist_ok=True);rows=[];provenance=[]
 for ci,(condition,color) in enumerate(CONDITIONS):
  for seed in SEEDS:
   if ci==0:run=PAPER/'data/figure_3/spatial_deconvolution/RCTD'/('strong_mix' if seed==2025 else f'strong_mix_seed{seed}')
   else:run=BASE/f'spider600k_genes{556 if ci==1 else 2000}_seed{seed}/results'
   meta=json.loads((run/'metrics_summary.json').read_text());src=Path(meta['spot_h5ad'])
   if not src.exists():src=Path(str(src).replace('/figure_3/cellbin_batch_sigma_slide/','/figure_3/strong_domain_mix/batch_sigma_slide/'))
   pred=pd.read_csv(run/'estimated_fractions_wide.csv');idx=pred.spot_id.str.removeprefix('spot_').astype(int).to_numpy()-1
   a=ad.read_h5ad(src,backed='r');truth=np.asarray(a.obsm['cell_type_frac_true'])[idx];a.file.close();est=pred[TYPES].to_numpy()
   assert pred.spot_id.is_unique and np.allclose(truth.sum(1),1,atol=1e-5) and np.allclose(est.sum(1),1,atol=1e-5)
   err=np.abs(est-truth);assert abs(np.sqrt(np.mean(err**2))-meta['overall_rmse'])<.0001
   for j,t in enumerate(TYPES):rows.append({'condition':condition,'seed':seed,'cell_type':t,'n_spots':len(idx),'mean_absolute_error':float(err[:,j].mean()),'rmse':float(np.sqrt(np.mean(err[:,j]**2)))})
   provenance.append({'condition':condition,'seed':seed,'run':str(run),'truth':str(src.resolve()),'reference':meta['cell_h5ad'],'replication':'reference seed; fixed query' if ci==0 else 'simulation seed; reference and query vary'})
 df=pd.DataFrame(rows);assert len(df)==72
 stats=df.groupby(['condition','cell_type'],sort=False).agg(mean_absolute_error=('mean_absolute_error','mean'),sd_absolute_error=('mean_absolute_error','std'),n_runs=('seed','size')).reset_index();assert (stats.n_runs==3).all()
 df.to_csv(OUT/'absolute_error_by_type_and_seed.csv',index=False);stats.to_csv(OUT/'absolute_error_seed_mean_sd.csv',index=False)
 fig,ax=plt.subplots(figsize=(11,5.5));positions=np.arange(8)
 for ci,(condition,color) in enumerate(CONDITIONS):
  sub=stats[stats.condition==condition].set_index('cell_type').loc[TYPES];x=positions+(ci-1)*.23
  for si,seed in enumerate(SEEDS):
   vals=df[(df.condition==condition)&(df.seed==seed)].set_index('cell_type').loc[TYPES,'mean_absolute_error']
   ax.scatter(x+(si-1)*.035,vals,s=16,color=color,alpha=.35,zorder=2)
  ax.errorbar(x,sub.mean_absolute_error,yerr=sub.sd_absolute_error,fmt='o',color=color,markersize=5,capsize=3,linewidth=1.5,label=condition,zorder=3)
 ax.set(xticks=positions,xticklabels=[f'Type {i}' for i in range(1,9)],ylabel='Mean spot absolute error',title='RCTD deconvolution · Mean ± SD across three runs',ylim=(0,None))
 ax.legend(frameon=False,ncol=3,loc='upper center',bbox_to_anchor=(.5,1.0));ax.margins(y=.2);ax.spines[['top','right']].set_visible(False);ax.grid(axis='y',alpha=.2);ax.set_axisbelow(True)
 fig.tight_layout()
 for ext in ['png','pdf','svg']:fig.savefig(OUT/f'rctd_error_seed_mean_comparison.{ext}',dpi=230,bbox_inches='tight')
 (OUT/'seed_mean_provenance.json').write_text(json.dumps({'metric':'Mean absolute estimated-minus-true fraction across scored spots for each type, then unweighted mean and sample SD across three runs. Not pooled spots or averaged predictions.','sources':provenance,'descriptive_only':True},indent=2)+'\n')
 print(OUT/'rctd_error_seed_mean_comparison.png')
if __name__=='__main__':main()
