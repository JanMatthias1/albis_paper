"""Explicit conceptual schematic; synthetic positions, not measured simulator output."""
from pathlib import Path
import json
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Circle,Patch
ROOT=Path('/dcs04/hicks/data/Jan/sim_project/sim_paper')
OUT=ROOT/'data/comparison_methods/figure_5D_hexagonal_capture'
OUT.mkdir(exist_ok=True)
colors=json.loads((ROOT/'data/comparison_methods/figure_5A_600k/settings.json').read_text())['celltype_colors']
rng=np.random.default_rng(20260928)
angles=np.arange(6)*np.pi/3
centers=np.vstack([[0,0],np.column_stack([np.cos(angles),np.sin(angles)])*2.55])
# Identical illustrative cell centers/types across columns isolate representation.
positions=[];types=[];owners=[]
for spot,center in enumerate(centers):
 for k in range(9):
  angle=rng.uniform(0,2*np.pi);radius=.83*np.sqrt(rng.random())
  positions.append(center+radius*np.array([np.cos(angle),np.sin(angle)]))
  types.append(int(rng.choice(8,p=np.roll(np.array([.32,.24,.12,.08,.07,.06,.06,.05]),spot))))
  owners.append(spot)
positions=np.array(positions);types=np.array(types);owners=np.array(owners)
# Illustrative molecule clouds around source cells, retained inside capture spots.
mrna=[];mt=[];mo=[]
for pos,t,spot in zip(positions,types,owners):
 points=pos+rng.normal(0,.115,(70,2))
 points=points[np.linalg.norm(points-centers[spot],axis=1)<=1]
 mrna.extend(points);mt.extend([t]*len(points));mo.extend([spot]*len(points))
mrna=np.array(mrna);mt=np.array(mt)
plt.rcParams.update({'font.family':'DejaVu Sans','font.size':12,'pdf.fonttype':42,'svg.fonttype':'none'})
fig,axes=plt.subplots(1,3,figsize=(13,5.8))
for ax,name in zip(axes,['ALBIS','scCube','SPIDER']):
 for i,center in enumerate(centers):
  ax.add_patch(Circle(center,1,facecolor='#FAFAFA',edgecolor='#333333' if i==0 else '#999999',lw=1.8 if i==0 else 1,zorder=0))
 if name=='ALBIS':
  ax.scatter(*mrna.T,c=[colors[f'type{t+1}'] for t in mt],s=2.4,linewidths=0,zorder=2,rasterized=True)
  subtitle='Individual mRNA instances'
 else:
  # Glyphs denote contributing cells, not reconstructed physical cell boundaries.
  for pos,t in zip(positions,types):
   ax.add_patch(Circle(pos,.12,facecolor=colors[f'type{t+1}'],edgecolor='white',lw=.6,zorder=2))
  subtitle='Contributing cells'
 ax.set_title(name+'\n'+subtitle,fontsize=15,fontweight='bold',pad=15)
 ax.set(xlim=(-3.8,3.8),ylim=(-3.5,3.5),aspect='equal');ax.set_axis_off()
fig.suptitle('Hexagonal capture example · Center spot and six neighbors',fontsize=17,fontweight='bold',y=.98)
fig.legend(handles=[Patch(facecolor=colors[f'type{i}'],label=f'Type {i}') for i in range(1,9)],loc='lower center',bbox_to_anchor=(.5,.105),ncol=8,frameon=False,fontsize=10)
fig.text(.5,.064,'Colors indicate source cell type for mRNA and cell identity for contributing cells.',ha='center',fontsize=10)
fig.text(.5,.023,'Illustrative schematic: synthetic positions and common hexagonal capture geometry; not native simulator outputs or a quantitative comparison.',ha='center',fontsize=9,color='#444444')
fig.subplots_adjust(left=.025,right=.975,top=.79,bottom=.19,wspace=.15)
for ext in ['png','pdf','svg']:fig.savefig(OUT/f'figure5d_hexagonal_capture.{ext}',dpi=300,bbox_inches='tight',facecolor='white')
plt.close(fig)
np.savez_compressed(OUT/'schematic_positions.npz',spot_centers=centers,cell_positions=positions,cell_types=types,cell_spots=owners,mrna_positions=mrna,mrna_source_types=mt,mrna_spots=np.array(mo))
(OUT/'provenance.json').write_text(json.dumps(dict(kind='Conceptual schematic, not simulator-derived data',seed=20260928,spots=7,geometry='Common illustrative hexagonal arrangement; unit capture radius',cell_positions='Same synthetic cell positions and types for all methods',mrna='Synthetic Gaussian clouds for illustration only; not regenerated ALBIS molecules',interpretation='Illustrates molecule-level versus cell-profile aggregation; does not claim native hexagonal/circular capture support in scCube or reproduce Figure5A',scCube='Circular boundaries illustrate the common capture example, not scCube native grid memberships'),indent=2)+'\n')
print(OUT/'figure5d_hexagonal_capture.png')
