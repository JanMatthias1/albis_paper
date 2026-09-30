"""Figure supp 6: existing geometry-only occupancy sensitivity on scCube slice 5."""
import json,hashlib
import argparse
from pathlib import Path
import numpy as np,pandas as pd,anndata as ad
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
ROOT=Path(__file__).resolve().parents[5]
# Retain the input tissue recorded in the promoted figure's original provenance.
BASE=ROOT/'sim_paper/data/comparison_methods/_archive/figure_5A_superseded_20260929/overview_figure2_600k_35881874'
OUT=ROOT/'sim_paper/data/comparison_methods/figure_supp_6'
OUT.mkdir(exist_ok=True)
parser=argparse.ArgumentParser(description=__doc__)
parser.add_argument('--redraw-only',action='store_true',help='Read saved spot tables; update figure exports only.')
args=parser.parse_args()
cfg=json.loads((BASE/'settings.json').read_text());scale=cfg['extent_um']/cfg['sccube_grid_size']
a=ad.read_h5ad(BASE/'sccube/cell.h5ad',backed='r');mask=a.obs.slice_id.to_numpy().astype(int)==4
xy=np.asarray(a.obsm['spatial'])[mask]/scale
labels=a.obs.cell_type_true.astype(str).to_numpy()[mask];ids=a.obs_names.to_numpy()[mask];a.file.close()
types=[f'type{i}' for i in range(1,9)];codes=pd.Categorical(labels,categories=types).codes
assert (codes>=0).all()
# The archived tissue settings predate this palette field; retain the colors
# recorded in the promoted figure's SVG legend.
colors=cfg.get('celltype_colors',dict(zip(types,[
 '#365CC0','#D66B24','#30934A','#8750B5','#CA3976','#687632','#C4A414','#16999D'])))
N=len(ids)
def aggregate(target):
 # Exact native scCube grid assignment (select_min_grid: first grid value > x).
 gap=np.floor(np.sqrt(N/target));grid_range=int(np.ceil(xy.max()))
 grid=np.array([i/gap for i in range(0,int(grid_range*(gap+1)),grid_range)])
 ix=np.searchsorted(grid,xy[:,0],side='right');iy=np.searchsorted(grid,xy[:,1],side='right')
 assert ix.max()<len(grid) and iy.max()<len(grid)
 ux=np.unique(ix);uy=np.unique(iy)
 xrank=np.searchsorted(ux,ix);yrank=np.searchsorted(uy,iy)
 native_id=xrank*len(uy)+yrank+1
 unique,inverse=np.unique(native_id,return_inverse=True)
 first=np.unique(native_id,return_index=True)[1]
 centers=np.column_stack([grid[ix[first]],grid[iy[first]]])
 centers[:,1]+=(xrank[first]%2)*grid_range/(gap*2)
 counts=np.zeros((len(unique),8),dtype=int);np.add.at(counts,(inverse,codes),1)
 occupancy=counts.sum(1);assert occupancy.sum()==N
 return centers*scale,counts,occupancy,inverse,grid_range/gap*scale,unique
if not args.redraw_only:
 # Validate vectorized geometry against the original native n_cell=10 output.
 centers,counts,occupancy,_,_,_=aggregate(10)
 b=ad.read_h5ad(BASE/'sccube/spot.h5ad',backed='r');sel=b.obs.slice_id.to_numpy().astype(int)==4
 saved=np.asarray(b.obsm['spatial'])[sel];saved_labels=b.obs.cell_type_true.astype(str).to_numpy()[sel];b.file.close()
 assert len(saved)==len(centers)
 sort=lambda x:np.lexsort((x[:,1],x[:,0]))
 i,j=sort(centers),sort(saved)
 assert np.allclose(centers[i],saved[j],atol=1e-3)
 assert np.array_equal(np.array(types)[counts.argmax(1)][i],saved_labels[j])
else:
 summary=pd.read_csv(OUT/'aggregation_summary.csv').set_index('target_cells_per_spot')
plt.rcParams.update({'font.family':'DejaVu Sans','font.size':11,'pdf.fonttype':42,'svg.fonttype':'none'})
fig,axes=plt.subplots(1,6,figsize=(19,4.7))
axes[0].scatter(*(xy*scale).T,c=[colors[t] for t in labels],s=.15,linewidths=0,rasterized=True)
axes[0].set_title(f'Cell input\n{N:,} cells',fontweight='bold',fontsize=12)
rows=[]
panel_centers=[((xy*scale).min(axis=0)+(xy*scale).max(axis=0))/2]
for ax,target in zip(axes[1:],[3,100,1000,5000,10000]):
 if args.redraw_only:
  table=pd.read_csv(OUT/f'spots_target{target}.csv')
  centers=table[['x_um','y_um']].to_numpy();counts=table[types].to_numpy()
  occupancy=table.n_cells.to_numpy();names=table.spot.to_numpy()
  step=summary.loc[target,'grid_step_um']
  assert occupancy.sum()==N
 else:
  centers,counts,occupancy,inverse,step,names=aggregate(target)
 panel_centers.append((centers.min(axis=0)+centers.max(axis=0))/2)
 purity=counts.max(1)/occupancy;dominant=np.array(types)[counts.argmax(1)]
 # Display glyph diameter scales with grid spacing, not a physical capture radius.
 size=max(1.2,(step/4100*165*.68)**2)
 ax.scatter(*centers.T,c=[colors[t] for t in dominant],s=size,linewidths=0,rasterized=True)
 ax.set_title(f'Target: {target:,} cells/spot\n{len(centers):,} occupied spots',fontsize=12,fontweight='bold')
 ax.text(.5,-.08,f'Mean occupancy: {occupancy.mean():.1f}\nMean dominant fraction: {purity.mean():.2f}',transform=ax.transAxes,ha='center',va='top',fontsize=10)
 row=dict(target_cells_per_spot=target,n_spots=len(centers),n_cells=N,mean_occupancy=occupancy.mean(),median_occupancy=float(np.median(occupancy)),min_occupancy=int(occupancy.min()),max_occupancy=int(occupancy.max()),mean_dominant_fraction=purity.mean(),grid_step_um=step)
 rows.append(row)
 table=pd.DataFrame(dict(spot=names,x_um=centers[:,0],y_um=centers[:,1],n_cells=occupancy,dominant_type=dominant,dominant_fraction=purity))
 for k,t in enumerate(types):table[t]=counts[:,k]
 if not args.redraw_only:
  table.to_csv(OUT/f'spots_target{target}.csv',index=False)
  pd.DataFrame(dict(cell=ids,spot=names[inverse])).to_csv(OUT/f'membership_target{target}.csv.gz',index=False)
for ax,center in zip(axes,panel_centers):
 ax.set(xlim=(center[0]-3450,center[0]+3450),ylim=(-400,6500),aspect='equal');ax.set_xticks([]);ax.set_yticks([])
 for spine in ax.spines.values():spine.set_visible(False)
 ax.plot([100,1100],[100,100],color='black',lw=2)
axes[0].text(600,190,'1 mm',ha='center',fontsize=9)
fig.suptitle('scCube spot aggregation · Slice 5',fontsize=17,fontweight='bold',y=1.05)
fig.legend(handles=[Patch(color=colors[t],label=f'Type {i+1}') for i,t in enumerate(types)],loc='lower center',bbox_to_anchor=(.5,.035),ncol=8,frameon=False,fontsize=10)
fig.subplots_adjust(left=.01,right=.99,top=.84,bottom=.28,wspace=.12)
for ext in ['png','pdf','svg']:fig.savefig(OUT/f'slice5_spot_aggregation.{ext}',dpi=300,bbox_inches='tight',facecolor='white')
plt.close(fig)
if args.redraw_only:
 print(OUT)
 raise SystemExit(0)
pd.DataFrame(rows).to_csv(OUT/'aggregation_summary.csv',index=False)
(OUT/'provenance.json').write_text(json.dumps(dict(input=str((BASE/'sccube/cell.h5ad').resolve()),slice_id=4,slice_number=5,n_cells=N,targets=[3,100,1000,5000,10000],platform='Visium',scope='Geometry and cell-type membership only; no expression regenerated or transformed',implementation='Vectorized exact native grid assignment and alternate-column Visium display shift',validation='n_cell10 centers and dominant labels match saved Figure5A native scCube spots',current_figure5a=dict(bin_target=3,spot_target=10),symbol_size='Display glyph diameter proportional to grid spacing, not physical capture footprint',native_source_sha256=hashlib.sha256((ROOT/'comparison_methods/env/scCube_src/scCube/utils.py').read_bytes()).hexdigest()),indent=2)+'\n')
print(pd.DataFrame(rows).to_string(index=False));print(OUT)
