"""Figure 5D: actual Figure5A simulator coordinates under an explicitly custom hex overlay.

ALBIS and scCube come from figure_5A_600k. SPIDER comes from a separate run with its
own seed (distinct_seeds/, see submit_distinct_seeds.sh): with Figure 5A's shared seed,
scCube and SPIDER draw identical uniform cell positions, which this panel would expose.
"""
from pathlib import Path
import json,hashlib
import numpy as np,anndata as ad
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Circle,Patch
from matplotlib.colors import to_rgb
from matplotlib.collections import PatchCollection
import argparse
import albis as ab
import albis.simulation_sphere
ROOT=Path('/dcs04/hicks/data/Jan/sim_project');PAPER=ROOT/'sim_paper'
BASE=PAPER/'data/figure_5/figure_5A_600k'
ALBIS_DIR=PAPER/'data/figure_5/figure_5D_actual_coordinates'  # shared ALBIS molecule cache
PANEL=ALBIS_DIR/'distinct_seeds'  # manuscript panel
parser=argparse.ArgumentParser();parser.add_argument('--spider-base',type=Path,default=PANEL);parser.add_argument('--out',type=Path,default=PANEL)
args=parser.parse_args();SPIDER_BASE=args.spider_base
OUT=args.out;OUT.mkdir(exist_ok=True)
SOURCE=BASE/'albis_generation/data/spot/generation_manifest.json'  # Figure 5A's own ALBIS generation record
manifest=json.loads(SOURCE.read_text());cfg=manifest['simulator_parameters'];palette=json.loads((BASE/'settings.json').read_text())['celltype_colors']
angles=np.arange(6)*np.pi/3;centers=np.vstack([[0,0],100*np.column_stack([np.cos(angles),np.sin(angles)])]);radius=27.5
source_code=Path(albis.simulation_sphere.__file__);digest=hashlib.sha256(source_code.read_bytes()).hexdigest()
expected=manifest['source_sha256'].get('albis/simulation_sphere.py')
assert digest==expected,f'Installed albis {ab.__version__} simulation_sphere.py does not match the one that generated {BASE}; regenerate Figure 5A with this albis'
cache=ALBIS_DIR/'albis_molecules_roi.npz'  # ALBIS is identical in every variant
record=ALBIS_DIR/'albis_regeneration.json'
if cache.exists():
 cached=json.loads(record.read_text()).get('source_sha256') if record.exists() else None
 assert cached==digest,f'{cache} was made by a different albis simulation_sphere.py; move it and {record.name} to an archive folder to regenerate'
if not cache.exists():
 print('Regenerating native ALBIS molecular realization',flush=True)
 base=ab.simulate_3d_molecule_sphere_base(**cfg)
 saved=ad.read_h5ad(BASE/'albis/cell.h5ad',backed='r');truth=base['adata_cell_obs']
 assert saved.n_obs==truth.n_obs
 assert np.allclose(saved.obsm['spatial_3d'],truth.obsm['spatial'])
 assert np.array_equal(saved.obs.cell_type_true.astype(str),truth.obs.cell_type_true.astype(str));saved.file.close()
 mol=base['molecules'];xyz=mol['full_xyz']
 # Slice5: fixed native Z slab [-410,0); central field selected before inspecting types.
 zlo=-cfg['sphere_R_um']+4*(2*cfg['sphere_R_um']/cfg['n_slices']);zhi=zlo+2*cfg['sphere_R_um']/cfg['n_slices']
 # Cache +-200 um (native_spots.py centres its ALBIS field on an ALBIS spot); this panel uses +-145.
 keep=(xyz[:,2]>=zlo)&(xyz[:,2]<zhi)&(np.abs(xyz[:,0])<=200)&(np.abs(xyz[:,1])<=200)
 np.savez_compressed(cache,xyz=xyz[keep],gene=mol['full_gene'][keep],source_type=mol['full_src_celltype'][keep])
 record.write_text(json.dumps(dict(source_sha256=digest,seed=cfg['seed'],config=cfg,validated='All600k original cell coordinates and type labels match saved Figure5A',z_bounds=[zlo,zhi],xy_half_width_um=200,molecules='Native pre-batch molecule instances; not post-resampling molecule coordinates'),indent=2)+'\n')
 del base,mol,xyz,truth
points={};types={};ids={}
m=np.load(cache);roi=(np.abs(m['xyz'][:,0])<=145)&(np.abs(m['xyz'][:,1])<=145)
points['ALBIS']=m['xyz'][roi,:2];types['ALBIS']=np.array([f'type{i+1}' for i in m['source_type'][roi]]);ids['ALBIS']=np.arange(len(points['ALBIS']))
# Project the saved spherical cells from the same slice onto the XY field.
# Include circles intersecting the field even when their centers are outside it.
a=ad.read_h5ad(BASE/'albis/cell.h5ad',backed='r')
cell_xyz=np.asarray(a.obsm['spatial_3d']);cell_radii=a.obs.cell_radius.to_numpy()
cell_keep=(a.obs.slice_id.to_numpy().astype(int)==4)&(np.abs(cell_xyz[:,0])<=145+cell_radii)&(np.abs(cell_xyz[:,1])<=145+cell_radii)
cell_xy=cell_xyz[cell_keep,:2];cell_radii=cell_radii[cell_keep]
cell_types=a.obs.cell_type_true.astype(str).to_numpy()[cell_keep]
cell_colors=[tuple(.55+.45*np.asarray(to_rgb(palette[t]))) for t in cell_types]
a.file.close()
for method,name in [('sccube','scCube'),('spider','SPIDER')]:
 a=ad.read_h5ad((SPIDER_BASE if method=='spider' else BASE)/method/'cell.h5ad',backed='r');xyz=np.asarray(a.obsm['spatial_3d']);xy=xyz[:,:2]-2050
 keep=(a.obs.slice_id.to_numpy().astype(int)==4)&(np.abs(xy[:,0])<=145)&(np.abs(xy[:,1])<=145)
 points[name]=xy[keep];types[name]=a.obs.cell_type_true.astype(str).to_numpy()[keep];ids[name]=a.obs_names.to_numpy()[keep];a.file.close()
 np.savez_compressed(OUT/f'{method}_cells_roi.npz',xy=points[name],types=types[name],cell_ids=ids[name])
plt.rcParams.update({'font.family':'DejaVu Sans','font.size':11,'pdf.fonttype':42,'svg.fonttype':'none'})
fig,axes=plt.subplots(1,3,figsize=(13,5.8));reports=[]
for ax,name in zip(axes,['ALBIS','scCube','SPIDER']):
 xy=points[name];labels=types[name];dist=((xy[:,None,:]-centers[None,:,:])**2).sum(2);assignment=dist.argmin(1);inside=dist.min(1)<=radius**2
 if name=='ALBIS':
  outlines=PatchCollection([Circle(center,float(r)) for center,r in zip(cell_xy,cell_radii)],facecolors=cell_colors,edgecolors=cell_colors,linewidths=.65,alpha=.45,zorder=.5)
  ax.add_collection(outlines)
  ax.scatter(*cell_xy.T,s=9,c=cell_colors,alpha=.65,marker='o',linewidths=0,zorder=2)
 ax.scatter(*xy[~inside].T,s=3 if name=='ALBIS' else 10,c='#BDBDBD',alpha=.35,marker='x' if name=='ALBIS' else 'o',linewidths=.35 if name=='ALBIS' else 0,rasterized=True)
 ax.scatter(*xy[inside].T,s=5 if name=='ALBIS' else 20,c=[palette[t] for t in labels[inside]],marker='x' if name=='ALBIS' else 'o',linewidths=.45 if name=='ALBIS' else 0,rasterized=True)
 counts=np.bincount(assignment[inside],minlength=7)
 for i,center in enumerate(centers):
  ax.add_patch(Circle(center,radius,fill=False,edgecolor='#333333' if i==0 else '#888888',lw=1.6 if i==0 else 1))
  ax.text(center[0],center[1]-radius-4,str(counts[i]),ha='center',va='top',fontsize=9)
 ax.text(.5,1.10,name,transform=ax.transAxes,ha='center',va='bottom',fontweight='bold',fontsize=15)
 ax.set_title('Native mRNA instances' if name=='ALBIS' else 'Native cell positions',fontweight='bold',fontsize=11,pad=12)
 ax.set(xlim=(-145,145),ylim=(-145,145),aspect='equal');ax.set_axis_off()
 ax.plot([-135,-85],[-137,-137],color='black',lw=2);ax.text(-110,-133,'50 µm',ha='center',fontsize=9)
 reports.append(dict(method=name,roi_points=len(xy),inside_overlay=int(inside.sum()),counts_by_overlay_spot=counts.tolist()))
fig.legend(handles=[Patch(color=palette[f'type{i}'],label=f'Type {i}') for i in range(1,9)],loc='lower center',bbox_to_anchor=(.5,.10),ncol=8,frameon=False,fontsize=10)
fig.subplots_adjust(left=.025,right=.975,top=.80,bottom=.20,wspace=.15)
for ext in ['png','pdf','svg']:fig.savefig(OUT/f'figure5d_actual_coordinates.{ext}',dpi=300,bbox_inches='tight',facecolor='white')
(OUT/'provenance.json').write_text(json.dumps(dict(source=str(BASE),spider_source=str(SPIDER_BASE/'spider'),slice_id=4,roi_rule='Central XY field of each tissue, selected independently of types; same physical scale',geometry='Custom seven-spot hexagonal overlay, not native captures',albis_cell_overlay=dict(source='albis/cell.h5ad',radius_field='obs.cell_radius',coordinates='obsm.spatial_3d',n_cells=len(cell_xy),selection='Centers in slice_id 4; projected circles intersecting the XY field',interpretation='Full-radius XY projections of spherical cells, not thin-plane cross-sections',style='Pastel translucent cell-type-colored circles and centers'),reports=reports),indent=2)+'\n');print(reports,flush=True)
