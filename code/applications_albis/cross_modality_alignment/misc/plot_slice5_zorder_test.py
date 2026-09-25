"""Diagnostic comparison of overlay drawing order using manuscript colors."""
from pathlib import Path
import numpy as np
import matplotlib.pyplot as plt
from plot_cross_tech_stair import load, TECH_COLORS, BASE


def main():
    adata, _ = load(5)
    coords = np.asarray(adata.obsm['spatial'])[:, :2]
    tech = adata.obs['technology'].astype(str).to_numpy()
    fig, axes = plt.subplots(1, 3, figsize=(13, 4.6))
    titles = ['spot drawn LAST', 'spot drawn FIRST (bottom layer)',
              'spot drawn on top, zorder=100, size=20']
    for ax, order, title in zip(axes,
            [('cell','bin16um','spot'), ('spot','cell','bin16um'), ('cell','bin16um','spot')], titles):
        for key in order:
            xy = coords[tech == key]
            emphasized = ax is axes[2] and key == 'spot'
            ax.scatter(xy[:,0], xy[:,1], c=TECH_COLORS[key],
                       s=20 if emphasized else (6 if key=='spot' else 1.5),
                       alpha=1 if key=='spot' else .85,
                       linewidths=.3 if emphasized else 0,
                       edgecolors='black' if emphasized else 'none',
                       zorder=100 if emphasized else 2)
        ax.set_aspect('equal'); ax.set_axis_off(); ax.set_title(title,fontsize=11)
    fig.tight_layout()
    out=Path(BASE)/'plots/slice5_zorder_test.png'
    fig.savefig(out,dpi=200);plt.close(fig)
    print('wrote',out)


if __name__=='__main__':
    main()
