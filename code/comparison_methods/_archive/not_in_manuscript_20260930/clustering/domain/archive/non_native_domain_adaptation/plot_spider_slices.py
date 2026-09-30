"""Show all native domain-labelled cells in slices 1–10; no clustering or subsampling."""
import argparse,json
from pathlib import Path
import anndata as ad
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

def main():
 p=argparse.ArgumentParser();p.add_argument('--input',type=Path,required=True);a=p.parse_args()
 x=ad.read_h5ad(a.input,backed='r');obs=x.obs.copy();xy=np.asarray(x.obsm['spatial']).copy();x.file.close()
 domains=sorted(obs.domain_true.astype(str).unique());colors=['#3366AA','#E68A2E','#32905C','#9854A6'];palette=dict(zip(domains,colors))
 slices=obs.slice_id.astype(int).to_numpy();labels=obs.domain_true.astype(str).to_numpy()
 assert set(slices)==set(range(10)) and len(obs)==600000 and len(domains)==4
 lo=xy.min(0);hi=xy.max(0);pad=(hi-lo)*.015
 fig,axes=plt.subplots(2,5,figsize=(18,8),sharex=True,sharey=True)
 counts={}
 for sid,ax in enumerate(axes.flat):
  idx=np.flatnonzero(slices==sid);idx=np.random.default_rng(sid).permutation(idx)
  ax.scatter(xy[idx,0],xy[idx,1],c=[palette[v] for v in labels[idx]],s=.65,linewidths=0,alpha=1,rasterized=True)
  ax.set_title(f'Slice {sid+1}',fontsize=13);ax.set_aspect('equal');ax.set_xlim(lo[0]-pad[0],hi[0]+pad[0]);ax.set_ylim(lo[1]-pad[1],hi[1]+pad[1]);ax.set_xticks([]);ax.set_yticks([])
  counts[str(sid+1)]=len(idx)
 fig.suptitle('SPIDER gyrus · Generated domains · 600,000 cells',fontsize=18)
 fig.legend(handles=[Line2D([],[],marker='o',linestyle='',color=palette[d],label=d,markersize=7) for d in domains],loc='lower center',ncol=4,frameon=False)
 fig.subplots_adjust(left=.015,right=.99,bottom=.09,top=.87,wspace=.08,hspace=.28)
 prefix=a.input.parent/'domain_slices_1_to_10'
 for ext in ['png','pdf']:fig.savefig(prefix.with_suffix('.'+ext),dpi=250,bbox_inches='tight')
 (prefix.with_suffix('.json')).write_text(json.dumps({'source':str(a.input.resolve()),'truth_column':'domain_true','slice_numbering':'displayed 1–10 = stored slice_id 0–9','all_cells_plotted':True,'counts_by_slice':counts,'colors':palette},indent=2)+'\n')
 print(prefix.with_suffix('.png'))
if __name__=='__main__':main()
