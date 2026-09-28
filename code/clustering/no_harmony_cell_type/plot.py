"""Separate cell-type recovery figure; existing tuned-batch figures stay intact."""
import json
import sys
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
from manuscript_style import apply_style, MODALITY_LOOKUP


def plot(root):
    df = pd.read_csv(root/'summary/metrics_by_seed.csv')
    assert len(df) == 36 and (df.status == 'complete').all()
    apply_style()
    fig, axes = plt.subplots(1, 2, figsize=(12, 5.5), sharey=True)
    for ax, mix in zip(axes, ['strong', 'weak']):
        for x, mod in enumerate(['cell', 'bin16um', 'spot']):
            for j, method in enumerate(['banksy', 'expression']):
                rows = df[(df.mix == mix) & (df.modality == mod) & (df.pipeline == method)]
                assert set(rows.seed) == {2025, 101, 202}
                vals = rows.cell_type_ari.to_numpy()
                mean, sd = vals.mean(), vals.std(ddof=1)
                xp, color = x+(j-.5)*.36, MODALITY_LOOKUP[mod]
                ax.bar(xp, mean, width=.34, color=color if j == 0 else 'white', edgecolor=color,
                       hatch=None if j == 0 else '///', zorder=2)
                ax.errorbar(xp, mean, yerr=sd, fmt='none', color='#333333', capsize=3, zorder=3)
                ax.scatter(np.full(3,xp), vals, s=16, color='#555555', edgecolors='white', linewidths=.5, zorder=4)
                ax.text(xp, max(mean+sd, vals.max(), 0)+.025, f'{mean:.2f}', ha='center', fontsize=11, fontweight='bold')
        ax.set(title=f'{mix.title()} domain mix', xticks=[0,1,2], xticklabels=['Cell','Bin (16 µm)','Spot'], ylim=(-.05,1.12))
        ax.set_axisbelow(True)
        ax.grid(axis='y', color='#DDDDDD', linewidth=.7)
        ax.spines[['top','right']].set_visible(False)
    axes[0].set_ylabel('Cell-type Adjusted Rand Index')
    fig.suptitle('Cell-type recovery · Batch = 0 · No Harmony', fontweight='bold')
    fig.legend(handles=[Patch(facecolor='#888888', label='BANKSY → PCA → Leiden'),
                        Patch(facecolor='white',edgecolor='#888888',hatch='///',label='Genes → PCA → Leiden')],
               loc='lower center', ncol=2, frameon=False)
    fig.tight_layout(rect=(0,.10,1,.95))
    for ext in ['png','pdf','svg']:
        fig.savefig(root/f'summary/cell_type_recovery_no_harmony.{ext}', dpi=300, bbox_inches='tight')
    plt.close(fig)
    (root/'summary/caption.txt').write_text('Cell-type recovery at batch_sigma=0 without Harmony, using BANKSY or expression-only PCA and K15/Leiden. Bars: three-seed mean; error bars: sample SD; dots: seeds 2025,101,202. BANKSY lambda selected by cell-type ARI on strong-mix seed2025, also included in this summary. Bin and spot labels denote dominant cell type, not mixture recovery.\n')
    (root/'summary/figure_provenance.json').write_text(json.dumps(dict(source=str(root/'summary/metrics_by_seed.csv'),
        selected_parameters=str(root/'selected_parameters.json'),code=str(Path(__file__).resolve())),indent=2)+'\n')
