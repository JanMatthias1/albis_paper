"""Final SPIDER RCTD figures from six completed native-count runs; no refitting."""
import json,hashlib
from pathlib import Path
import numpy as np,pandas as pd,anndata as ad
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
ROOT=Path(__file__).resolve().parents[5]
BASE=ROOT/'sim_paper/data/figure_5/spatial_deconvolution'
OUT=BASE/'final_spider600k'
SEEDS=[2025,101,202];GENES=[556,2000];TYPES=[f'type{i}' for i in range(1,9)]
COLORS=['#365CC0','#D66B24','#30934A']
def save(fig,name):
 for ext in ['png','pdf','svg']:fig.savefig(OUT/f'{name}.{ext}',dpi=250,bbox_inches='tight')
 plt.close(fig)
def main():
 OUT.mkdir(exist_ok=True);metrics=[];type_rows=[];spatial={};provenance=[]
 plt.rcParams.update({'font.size':10,'axes.spines.top':False,'axes.spines.right':False,'pdf.fonttype':42,'svg.fonttype':'none'})
 for genes in GENES:
  for seed in SEEDS:
   run=BASE/f'spider600k_genes{genes}_seed{seed}';results=run/'results'
   summary=json.loads((results/'metrics_summary.json').read_text());audit=json.loads((run/'inputs/input_audit.json').read_text())
   assert audit['underlying_cells']==600000 and audit['genes']==genes
   pred=pd.read_csv(results/'estimated_fractions_wide.csv');assert pred.spot_id.is_unique
   a=ad.read_h5ad(run/'inputs/spot.h5ad',backed='r')
   idx=pred.spot_id.str.removeprefix('spot_').astype(int).to_numpy()-1
   assert (idx>=0).all() and (idx<a.n_obs).all()
   truth=np.asarray(a.obsm['cell_type_frac_true'])[idx];estimate=pred[TYPES].to_numpy()
   assert np.allclose(truth.sum(1),1) and np.allclose(estimate.sum(1),1)
   rmse=float(np.sqrt(np.mean((truth-estimate)**2)));corr=float(np.corrcoef(truth.ravel(),estimate.ravel())[0,1])
   assert abs(rmse-summary['overall_rmse'])<.0001 and abs(corr-summary['overall_pearson_r'])<.0001
   metrics.append(dict(genes=genes,seed=seed,pearson_r=corr,rmse=rmse,spots_scored=len(idx),spots_total=a.n_obs))
   for j,t in enumerate(TYPES):type_rows.append(dict(genes=genes,seed=seed,cell_type=t,rmse=float(np.sqrt(np.mean((truth[:,j]-estimate[:,j])**2)))))
   keep=a.obs.iloc[idx].slice_id.astype(int).to_numpy()==4
   spatial[(genes,seed)]=(np.asarray(a.obsm['spatial'])[idx][keep,:2],truth[keep],estimate[keep])
   a.file.close()
   provenance.append({'run':str(run),'audit':audit,'metrics_sha256':hashlib.sha256((results/'metrics_summary.json').read_bytes()).hexdigest()})
 df=pd.DataFrame(metrics);ct=pd.DataFrame(type_rows);assert len(df)==6 and len(ct)==48
 df.to_csv(OUT/'metrics_by_seed.csv',index=False);ct.to_csv(OUT/'cell_type_rmse_by_seed.csv',index=False)
 df.groupby('genes')[['pearson_r','rmse']].agg(['mean','std','min','max']).to_csv(OUT/'metrics_summary.csv')
 fig,axes=plt.subplots(1,2,figsize=(8,4))
 for ax,key,title in zip(axes,['pearson_r','rmse'],['Pearson correlation','RMSE']):
  ax.boxplot([df[df.genes==g][key] for g in GENES],widths=.35,showfliers=False,medianprops={'color':'black'})
  for j,seed in enumerate(SEEDS):ax.scatter(np.array([1,2])+(j-1)*.065,[df[(df.genes==g)&(df.seed==seed)][key].item() for g in GENES],color=COLORS[j],s=35,zorder=3,label=f'Seed {seed}')
  ax.set(xticks=[1,2],xticklabels=['556 genes','2,000 genes'],ylabel=title)
 axes[0].set_ylim(0,1.02);axes[1].set_ylim(0,.22);axes[0].legend(frameon=False,fontsize=8)
 fig.suptitle('SPIDER deconvolution · 600,000 cells · 10 slices');fig.tight_layout();save(fig,'deconvolution_overview')
 fig,axes=plt.subplots(1,2,figsize=(11,4),sharey=True)
 for ax,genes in zip(axes,GENES):
  for j,seed in enumerate(SEEDS):
   vals=ct[(ct.genes==genes)&(ct.seed==seed)].set_index('cell_type').loc[TYPES,'rmse']
   ax.scatter(np.arange(8)+(j-1)*.12,vals,color=COLORS[j],s=25,label=f'Seed {seed}')
  ax.set(xticks=range(8),xticklabels=TYPES,title=f'{genes:,} genes',ylim=(0,None));ax.tick_params(axis='x',rotation=45)
 axes[0].set_ylabel('Cell-type fraction RMSE');axes[1].legend(frameon=False,fontsize=8);fig.tight_layout();save(fig,'cell_type_errors')
 for seed in SEEDS:
  fig,axes=plt.subplots(4,8,figsize=(17,9),layout='constrained')
  for gi,genes in enumerate(GENES):
   xy,true,est=spatial[(genes,seed)]
   for j,t in enumerate(TYPES):
    for off,values in enumerate([true,est]):
     ax=axes[gi*2+off,j];im=ax.scatter(xy[:,0],xy[:,1],c=values[:,j],vmin=0,vmax=1,cmap='viridis',s=7,marker='o',linewidths=0,rasterized=True)
     ax.set_aspect('equal');ax.set_xlim(0,4100);ax.set_ylim(0,4100);ax.set_xticks([]);ax.set_yticks([])
     if gi==0 and off==0:ax.set_title(t)
     if j==0:ax.set_ylabel(f'{genes:,} genes\n'+('True' if off==0 else 'RCTD'))
  fig.colorbar(im,ax=axes.ravel().tolist(),shrink=.6,label='Cell-type fraction')
  fig.suptitle(f'SPIDER · Slice 5 · Simulation seed {seed}')
  save(fig,f'spatial_true_vs_estimated_seed{seed}')
 (OUT/'provenance.json').write_text(json.dumps({'runs':provenance,'plot_code_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),'spatial_slice_id':4,'spatial_color_limits':[0,1],'spatial_marker':'Circular display glyph; native capture is circular','replication':'Three simulation seeds; no spot-level pseudoreplication in summary boxes'},indent=2)+'\n')
 (OUT/'README.md').write_text('''# Final SPIDER deconvolution figures

- deconvolution_overview: overall Pearson correlation and RMSE; each point is one simulation seed, boxes summarize three seeds per gene panel.
- cell_type_errors: per-cell-type fraction RMSE for every seed.
- spatial_true_vs_estimated_seed*: all scored spots from Slice 5, all eight types, both gene panels; common 0–1 color scale. Circular glyphs represent spots; marker size is for display, not a physical scale.

All figures are available as PNG/PDF/SVG. Tables contain recomputed full-precision scores, validated against RCTD summaries. Six native-count runs use 600,000 cells, ten slices, 556/2000 genes, seeds 2025/101/202, and native SPIDER circular aggregation (55 µm diameter, 100 µm spacing). Empty/zero-expression spots were excluded before RCTD; scored/total supplied counts are reported. RCTD full mode, fixed R seed 0, default internal filtering. No count rounding or rescaling.

Reference is the original 10k-cell Splatter expression pool sampled by SPIDER; it is not independent of the queries. Truth is contributing-cell fractions. Scores do not establish equivalence to ALBIS molecule-fraction recovery. The two gene conditions were independently generated rather than obtained by subsetting one matrix. No scCube results are included; previous scCube runs are archived separately.
''')
 print('Validated six runs; wrote overview, per-type errors and three spatial figures to',OUT,flush=True)
if __name__=='__main__':main()
