"""Wide barplots: batch-zero domain recovery and tuned-batch cell-type recovery.

Run with sim_paper/env/albis-tutorial/bin/python from the project root.
Existing figures remain unchanged; outputs use banksy_vs_pca_recovery_condensed.
"""
import csv
import json
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
import plot_banksy_vs_pca_recovery as source


def main():
    data, lambdas = source.load()
    out = source.FIG3_DIR / 'ari_recovery_summary'
    stem = out / 'banksy_vs_pca_recovery_condensed'
    styles = [
        ('banksy', '#0072B2', 'o', r'BANKSY, $\sigma_{\mathrm{batch}}=0$'),
        ('genes_bs0', '#D55E00', 's', r'Expression, $\sigma_{\mathrm{batch}}=0$'),
        ('genes_bs', '#555555', '^', r'Expression, modality-specific $\sigma_{\mathrm{batch}}$'),
    ]
    groups = [(mix, mod) for mix in ['strong', 'weak'] for mod in source.MODALITIES]
    positions = np.array([0., 1., 2., 3.5, 4.5, 5.5])
    rows, summary = [], []
    source.apply_style()
    plt.rcParams.update({'axes.titlesize': source.TITLE_SIZE, 'axes.labelsize': source.LABEL_SIZE,
                         'xtick.labelsize': source.TICK_SIZE, 'ytick.labelsize': source.TICK_SIZE,
                         'pdf.fonttype': 42, 'svg.fonttype': 'none'})
    fig, axes = plt.subplots(1, 2, figsize=(16, 5.6), sharey=True)
    for ax, target, title in zip(axes, source.TARGETS, ['Spatial domain recovery\nσ = 0', 'Cell-type recovery\nTuned σ']):
        active = styles[:1] if target == 'domain_true' else styles[2:]
        for j, (pipeline, color, marker, label) in enumerate(active):
            for x, (mix, mod) in zip(positions, groups):
                values = data.get((mix, pipeline, mod, target), {})
                if set(values) != set(source.SEEDS):
                    raise RuntimeError(f'Expected three completed seeds: {mix}/{pipeline}/{mod}/{target}: {values}')
                vals = np.array([values[s] for s in source.SEEDS])
                mean, sd = float(vals.mean()), float(vals.std(ddof=1))
                color = source.MODALITY_LOOKUP[mod]
                xp = x + (j-(len(active)-1)/2)*.38
                ax.bar(xp, mean, width=.38 if len(active)>1 else .5,
                       color=color if pipeline == 'banksy' else 'white',
                       edgecolor='white' if pipeline == 'banksy' else color,
                       hatch=None if pipeline == 'banksy' else '///',
                       linewidth=1.2, zorder=2)
                ax.errorbar(xp, mean, yerr=sd, fmt='none', ecolor=source.INK,
                            elinewidth=1.1, capsize=3, zorder=4)
                ax.scatter(np.full(len(vals), xp), vals, color=source.INK_MUTED,
                           marker='o', s=12, edgecolors='white', linewidths=.5, zorder=5)
                ax.text(xp, max(mean+sd, float(vals.max()), 0)+.02, f'{round(mean,2)+0.:.2f}',
                        ha='center', va='bottom', fontsize=source.ANNOT_SIZE,
                        fontweight='bold', color=source.INK)
                lam = lambdas[(mod, target)] if pipeline == 'banksy' else ''
                sigma = float(source.BATCH_SIGMA[mod]) if pipeline == 'genes_bs' else 0.
                common = dict(mix=mix, modality=mod, target=target, pipeline=pipeline,
                              batch_sigma=sigma, banksy_lambda=lam)
                summary.append(dict(common, mean_ari=mean, sample_sd=sd, n_seeds=3))
                for seed, value in zip(source.SEEDS, vals):
                    rows.append(dict(common, seed=seed, ari=float(value)))
        ax.set(title=title, xticks=positions,
               xticklabels=[source.MODALITY_DISPLAY[m] + (f'\nσ={source.BATCH_SIGMA[m]}' if target == 'cell_type_true' else '') for _, m in groups],
               xlim=(-.55, 6.05), ylim=(-.05, 1.08), yticks=np.arange(0, 1.01, .2))
        ax.axvline(2.75, color='#CCCCCC', lw=.8)
        source.style_axis(ax)
        for x, text in [(1., 'Strong domain mix'), (4.5, 'Weak domain mix')]:
            ax.text(x, -.25, text, transform=ax.get_xaxis_transform(),
                    ha='center', va='top', fontsize=source.TICK_SIZE, fontweight='bold')
    axes[0].set_ylabel('Adjusted Rand Index')
    handles = [Patch(facecolor='#808080', edgecolor='white', label='BANKSY → PCA → Harmony'),
               Patch(facecolor='white', edgecolor='#808080', hatch='///', label='Genes → PCA → Harmony')]
    fig.legend(handles=handles, loc='lower center', bbox_to_anchor=(.5, .005),
               ncol=2, frameon=False, fontsize=source.TICK_SIZE, columnspacing=2.0, handletextpad=.6)
    fig.subplots_adjust(left=.055, right=.992, top=.83, bottom=.30, wspace=.08)
    for ext in ['png', 'pdf', 'svg']:
        fig.savefig(f'{stem}.{ext}', dpi=300, facecolor='white', bbox_inches='tight')
    plt.close(fig)
    for suffix, records in [('by_seed', rows), ('summary', summary)]:
        with open(f'{stem}_{suffix}.csv', 'w') as f:
            writer = csv.DictWriter(f, fieldnames=list(records[0]))
            writer.writeheader()
            writer.writerows(records)
    assert len(rows) == 36 and len(summary) == 12
    metadata = dict(source_plot_code=str(source.__file__), n_groups=12, n_seed_scores=36,
                    retained='Domain: BANKSY only at batch0. Cell type: expression-only at tuned batch. Both mixes and all modalities.',
                    aggregation='Bars: mean; error bars: sample SD; dots: simulation seeds 2025, 101, 202.',
                    downstream='Every series uses PCA, Harmony by slice, and Leiden with resolution selected for known category count.',
                    banksy='Domain lambda: cell1.0/bin0.5/spot0.3. All BANKSY inputs have batch_sigma0; no tuned-batch BANKSY series is shown.',
                    batch='Nonzero expression series: cell1.5/bin0.7/spot0.3.',
                    spatial='k_geom cell60/bin100/spot8.',
                    styling='Original Figure 3 modality palette and typography. Solid: BANKSY; hatched: expression-only. Numeric labels are means.',
                    aggregate_labels='For bins and spots, cell-type recovery scores the stored dominant cell-type label, not mixture recovery.')
    with open(f'{stem}_provenance.json', 'w') as f:
        json.dump(metadata, f, indent=2)
    caption = ('Spatial-domain recovery uses BANKSY embeddings at '
               'batch_sigma=0; cell-type recovery uses expression-only embeddings with tuned '
               'batch_sigma=1.5, 0.7, and 0.3 for cells, 16um bins, and spots, respectively. '
               'Both panels show strong and weak domain mixing. All series use PCA, Harmony '
               'by slice, and Leiden with resolution selected to match the known category count. '
               'Bars show mean ARI across simulation seeds 2025, 101, 202; '
               'error bars show sample SD and dots show all three seeds. Numeric labels show means. '
               'BANKSY domain lambda is 1.0, 0.5, 0.3 for cells, bins, spots; k_geom is 60, 100, 8. '
               'Parameters were selected on strong-mix data and reused for weak-mix data. '
               'Bin/spot cell-type recovery scores dominant labels, not mixture fractions.\n')
    with open(f'{stem}_caption.txt', 'w') as f:
        f.write(caption)
    print(f'Validated {len(rows)} seed scores; saved {stem}.png')


if __name__ == '__main__':
    main()
