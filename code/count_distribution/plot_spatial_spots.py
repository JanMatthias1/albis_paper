"""Plot saved ALBIS spots at physical capture radius; slice IDs are zero-based."""
import argparse
from pathlib import Path
import sys
import json
import anndata as ad
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.collections import PatchCollection
from matplotlib.patches import Circle, Patch
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from manuscript_style import CELLTYPE_COLORS

p = argparse.ArgumentParser(description=__doc__)
p.add_argument('--input', type=Path, required=True)
p.add_argument('--output-dir', type=Path, required=True)
p.add_argument('--slice-id', type=int, default=4)  # 5th slice from the bottom
a = p.parse_args()
x = ad.read_h5ad(a.input, backed='r')
mask = x.obs.slice_id.to_numpy() == a.slice_id
obs = x.obs.loc[mask].copy()
xy = np.asarray(x.obsm['spatial'])[mask, :2]
x.file.close()
if not len(obs):
    p.error('Slice ID has no observations')
keep = (~obs.is_empty.to_numpy(dtype=bool)) & (obs.cell_type_true.astype(str).to_numpy() != 'unassigned')
labels = obs.cell_type_true.astype(str).to_numpy()[keep]
colors = {f'type{i}': CELLTYPE_COLORS[f'Cell Type {i}'] for i in range(1, 9)}
plt.rcParams.update({'font.family': 'DejaVu Sans', 'font.size': 12, 'pdf.fonttype': 42, 'svg.fonttype': 'none'})
fig, ax = plt.subplots(figsize=(7, 7))
patches = [Circle(pos, radius) for pos, radius in zip(xy[keep], obs.spot_radius_um.to_numpy()[keep])]
ax.add_collection(PatchCollection(patches, facecolors=[colors[t] for t in labels], edgecolors='none', alpha=.75))
# Retain the complete capture window, including empty locations.
radius = float(obs.spot_radius_um.max())
low, high = xy.min(axis=0)-radius, xy.max(axis=0)+radius
pad = .03 * max(high-low)
ax.set(xlim=(low[0]-pad,high[0]+pad), ylim=(low[1]-pad,high[1]+pad), aspect='equal', xlabel='X (µm)', ylabel='Y (µm)', title=f'ALBIS spots — slice_id={a.slice_id}')
ax.spines[['top', 'right']].set_visible(False)
fig.legend(handles=[Patch(facecolor=c, alpha=.75, label=f'Cell Type {i}') for i,c in enumerate(colors.values(),1)], loc='lower center', ncol=4, frameon=False, fontsize=10)
fig.subplots_adjust(left=.14, right=.97, top=.91, bottom=.17)
a.output_dir.mkdir(parents=True, exist_ok=True)
stem = a.output_dir / f'spatial_spots_slice{a.slice_id}'
for ext in ('png','pdf','svg'):
    fig.savefig(stem.with_suffix('.'+ext), dpi=500)
plt.close(fig)
summary = dict(input=str(a.input.resolve()), slice_id=a.slice_id, total_spots=len(obs), displayed_spots=int(keep.sum()), hidden_spots=int((~keep).sum()), spot_radius_um=radius, spot_spacing_um=sorted(obs.spot_spacing_um.unique().tolist()), coordinates='spatial', color='cell_type_true')
stem.with_suffix('.json').write_text(json.dumps(summary, indent=2)+'\n')
print(json.dumps(summary, indent=2))
