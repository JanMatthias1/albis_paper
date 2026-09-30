"""3D view of the cached Figure 5D ALBIS molecular field; no regeneration."""
from pathlib import Path
import json
import numpy as np
import anndata as ad
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import to_rgb
from matplotlib.patches import Patch
from mpl_toolkits.mplot3d.art3d import Poly3DCollection

ROOT=Path(__file__).resolve().parents[4]
OUT=ROOT/'sim_paper/data/figure_5/figure_5D_actual_coordinates'
BASE=ROOT/'sim_paper/data/figure_5/figure_5A_600k'
palette=json.loads((BASE/'settings.json').read_text())['celltype_colors']
m=np.load(OUT/'albis_molecules_roi.npz');xyz=m['xyz']
labels=np.array([f'type{i+1}' for i in m['source_type']])
zlo,zhi=json.loads((OUT/'albis_regeneration.json').read_text())['z_bounds']
a=ad.read_h5ad(BASE/'albis/cell.h5ad',backed='r')
c=np.asarray(a.obsm['spatial_3d']);r=a.obs.cell_radius.to_numpy()
keep=(a.obs.slice_id.to_numpy().astype(int)==4)&(np.abs(c[:,0])+r<=145)&(np.abs(c[:,1])+r<=145)
ct=a.obs.cell_type_true.astype(str).to_numpy()[keep];c=c[keep];r=r[keep];a.file.close()
angles=np.arange(6)*np.pi/3
centers=np.vstack([[0,0],100*np.column_stack([np.cos(angles),np.sin(angles)])])
inside=(((xyz[:,None,:2]-centers[None,:,:])**2).sum(2).min(1)<=27.5**2)
plt.rcParams.update({'font.family':'DejaVu Sans','pdf.fonttype':42,'svg.fonttype':'none'})
fig=plt.figure(figsize=(9,10));ax=fig.add_subplot(projection='3d',computed_zorder=False)
# Translucent spherical surfaces using saved radii and centers.
u=np.linspace(0,2*np.pi,13);v=np.linspace(0,np.pi,9)
unit=np.stack([np.outer(np.cos(u),np.sin(v)),np.outer(np.sin(u),np.sin(v)),np.outer(np.ones_like(u),np.cos(v))],axis=-1)
faces=[];colors=[]
for center,radius,label in zip(c,r,ct):
 surface=center+radius*unit
 pastel=.55+.45*np.asarray(to_rgb(palette[label]))
 for i in range(len(u)-1):
  for j in range(len(v)-1):
   face=surface[[i,i+1,i+1,i],[j,j,j+1,j+1]]
   if np.all((face[:,2]>=zlo)&(face[:,2]<=zhi)):
    faces.append(face);colors.append((*pastel,.45))
ax.add_collection3d(Poly3DCollection(faces,facecolors=colors,edgecolors='none',zorder=1,rasterized=True))
ax.scatter(*xyz[~inside].T,s=1,marker='x',c='#BBBBBB',alpha=.12,linewidths=.25,depthshade=False,rasterized=True,zorder=2)
ax.scatter(*xyz[inside].T,s=4,marker='x',c=[palette[t] for t in labels[inside]],alpha=.85,linewidths=.45,depthshade=False,rasterized=True,zorder=3)
ax.scatter(*c.T,s=7,c=[.55+.45*np.asarray(to_rgb(palette[t])) for t in ct],depthshade=False,alpha=.8,zorder=4)
t=np.linspace(0,2*np.pi,121)
for i,(x,y) in enumerate(centers):
 ax.plot(x+27.5*np.cos(t),y+27.5*np.sin(t),np.full_like(t,zlo),color='#333333' if i==0 else '#888888',lw=1.5,zorder=5)
ax.set(xlim=(-145,145),ylim=(-145,145),zlim=(zlo,zhi),xlabel='X (µm)',ylabel='Y (µm)',zlabel='Z (µm)')
ax.set_box_aspect((290,290,zhi-zlo));ax.view_init(elev=24,azim=-55)
ax.set_title('ALBIS · mRNA and cells in 3D',fontsize=18,fontweight='bold',pad=18)
for axis in (ax.xaxis,ax.yaxis,ax.zaxis):axis.pane.fill=False
ax.grid(False)
fig.legend(handles=[Patch(color=palette[f'type{i}'],label=f'Type {i}') for i in range(1,9)],loc='lower center',ncol=4,frameon=False)
fig.subplots_adjust(bottom=.10,top=.93,left=0,right=.95)
for ext in ['png','pdf','svg']:fig.savefig(OUT/f'figure5d_albis_3d.{ext}',dpi=250,bbox_inches='tight')
(OUT/'figure5d_albis_3d_provenance.json').write_text(json.dumps(dict(molecules=len(xyz),cells=len(c),z_bounds=[zlo,zhi],cell_selection='Slice 5 centers with full XY circle inside field; sphere mesh clipped at slab Z boundaries',geometry='Saved XYZ and cell radii; illustrative capture circles on lower Z plane; color membership uses XY distance through entire slab',rendering='Translucent spheres and mRNA crosses; matplotlib transparency is approximate'),indent=2)+'\n')
print('Saved',OUT/'figure5d_albis_3d.png')
