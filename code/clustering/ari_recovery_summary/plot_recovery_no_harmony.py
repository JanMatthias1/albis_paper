"""No-Harmony domain recovery with unchanged tuned-batch cell-type recovery.

Run with albis_paper/env/albis-tutorial/bin/python from the project root.
Writes a separate banksy_vs_pca_recovery_condensed_no_harmony figure.
"""
import csv
from pathlib import Path
import json
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
import plot_banksy_vs_pca_recovery as source


def main():
    data, lambdas = source.load()
    experiment = source.FIG3_DIR / 'no_harmony_domain'
    selected = json.loads((experiment/'selected_parameters.json').read_text())
    records = list(csv.DictReader((experiment/'summary/metrics_by_seed.csv').open()))
    assert len(records) == 36 and all(r['status'] in ('complete', 'not_converged') for r in records)
    # runs whose Leiden search missed 6 domains within 15 trials: left out of the bar (labelled n=2)
    excluded = [r for r in records if r['status'] == 'not_converged']
    for m in source.MODALITIES:
        lambdas[(m,'domain_true')] = str(selected['modalities'][m]['lambda'])
        for mix in ['strong','weak']:
            for pipeline, key in [('banksy','banksy'),('expression','genes_bs0')]:
                matches = [r for r in records if r['mix']==mix and r['modality']==m and r['pipeline']==pipeline]
                assert {int(r['seed']) for r in matches} == set(source.SEEDS)
                data[(mix,key,m,'domain_true')] = {int(r['seed']):float(r['domain_ari']) for r in matches if r['status'] == 'complete'}
    out = source.FIG3_DIR / 'ari_recovery_summary'
    stem = out / 'banksy_vs_pca_recovery_condensed_no_harmony'
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
    for ax, target, title in zip(axes, source.TARGETS, ['Spatial domain recovery\nNo Harmony · σ = 0', 'Cell-type recovery\nHarmony · Tuned σ']):
        active = styles[:2] if target == 'domain_true' else styles[2:]
        for j, (pipeline, color, marker, label) in enumerate(active):
            for x, (mix, mod) in zip(positions, groups):
                values = data.get((mix, pipeline, mod, target), {})
                dropped = {int(r['seed']) for r in excluded if target == 'domain_true' and (r['mix'], r['modality']) == (mix, mod)
                           and pipeline == ('banksy' if r['pipeline'] == 'banksy' else 'genes_bs0')}
                if set(values) != set(source.SEEDS) - dropped or len(values) < 2:
                    raise RuntimeError(f'Expected completed seeds: {mix}/{pipeline}/{mod}/{target}: {values}')
                vals = np.array([values[s] for s in sorted(values)])
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
                ax.text(xp, max(mean+sd, float(vals.max()), 0)+.02, f'{round(mean,2)+0.:.2f}' + (f'\nn={len(vals)}' if dropped else ''),
                        ha='center', va='bottom', fontsize=source.ANNOT_SIZE,
                        fontweight='bold', color=source.INK)
                lam = lambdas[(mod, target)] if pipeline == 'banksy' else ''
                sigma = float(source.BATCH_SIGMA[mod]) if pipeline == 'genes_bs' else 0.
                common = dict(mix=mix, modality=mod, target=target, pipeline=pipeline,
                              batch_sigma=sigma, banksy_lambda=lam)
                summary.append(dict(common, mean_ari=mean, sample_sd=sd, n_seeds=len(vals)))
                for seed, value in sorted(values.items()):
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
    handles = [Patch(facecolor='#808080', edgecolor='white', label='BANKSY-augmented expression'),
               Patch(facecolor='white', edgecolor='#808080', hatch='///', label='Expression only')]
    fig.legend(handles=handles, loc='lower center', bbox_to_anchor=(.5, .005),
               ncol=2, frameon=False, fontsize=source.TICK_SIZE, columnspacing=2.0, handletextpad=.6)
    fig.subplots_adjust(left=.055, right=.992, top=.83, bottom=.30, wspace=.08)
    for ext in ['png', 'pdf', 'svg']:
        fig.savefig(f'{stem}.{ext}', dpi=500, facecolor='white', bbox_inches='tight')
    plt.close(fig)
    for suffix, records in [('by_seed', rows), ('summary', summary)]:
        with open(f'{stem}_{suffix}.csv', 'w') as f:
            writer = csv.DictWriter(f, fieldnames=list(records[0]))
            writer.writeheader()
            writer.writerows(records)
    assert len(rows) == 54 - len(excluded) and len(summary) == 18
    not_converged = [f"{r['mix']} mix / {source.MODALITY_DISPLAY[r['modality']]} / seed {r['seed']} / "
                     f"{'BANKSY' if r['pipeline'] == 'banksy' else 'expression-only'} ({r['achieved_clusters']} clusters)" for r in excluded]
    metadata = dict(source_plot_code=str(Path(__file__).resolve()), selected_parameters=str(experiment/'selected_parameters.json'), n_groups=18, n_seed_scores=len(rows),
                    not_converged_excluded=not_converged,
                    retained='Domain: BANKSY and expression-only at batch0. Cell type: expression-only at tuned batch. Both mixes and all modalities.',
                    aggregation='Bars: mean; error bars: sample SD; dots: simulation seeds 2025, 101, 202.',
                    downstream='Domain: PCA/Leiden without Harmony. Cell type: existing PCA/Harmony/Leiden. Resolution selected for known category count.',
                    banksy={m:selected['modalities'][m]['lambda'] for m in source.MODALITIES}, selection_note='Chosen on strong seed2025, which remains in the displayed three-seed summary; not held out.',
                    batch='Nonzero expression series: cell1.5/bin0.7/spot0.3.',
                    spatial='k_geom cell60/bin100/spot8.',
                    styling='Original Figure 3 modality palette and typography. Solid: BANKSY; hatched: expression-only. Numeric labels are means.',
                    aggregate_labels='For bins and spots, cell-type recovery scores the stored dominant cell-type label, not mixture recovery.')
    with open(f'{stem}_provenance.json', 'w') as f:
        json.dump(metadata, f, indent=2)
    caption = ('Spatial-domain recovery compares BANKSY and expression-only embeddings at '
               'batch_sigma=0; cell-type recovery uses expression-only embeddings with tuned '
               'batch_sigma=1.5, 0.7, and 0.3 for cells, 16um bins, and spots, respectively. '
               'Both panels show strong and weak domain mixing. Domain recovery uses PCA and Leiden without Harmony; '
               'cell-type recovery retains the existing PCA, Harmony-by-slice, and Leiden results. '
               'Leiden resolution targets the known category count. '
               'Bars show mean ARI across simulation seeds 2025, 101, 202; '
               'error bars show sample SD and dots show all three seeds. Numeric labels show means. '
               f"BANKSY domain lambda (cell/bin/spot): {[selected['modalities'][m]['lambda'] for m in source.MODALITIES]}; k_geom is 60, 100, 8. "
               'Parameters were selected on strong-mix seed2025 and frozen across final runs; seed2025 is not held out. '
               'Bin/spot cell-type recovery scores dominant labels, not mixture fractions.'
               + (' Bars labelled n=2 average the remaining seeds: Leiden did not reach exactly 6 domains within 15 '
                  'resolution trials for ' + '; '.join(not_converged) + ', which is excluded.' if not_converged else '')
               + '\n')
    with open(f'{stem}_caption.txt', 'w') as f:
        f.write(caption)
    print(f'Validated {len(rows)} seed scores; saved {stem}.png')


if __name__ == '__main__':
    main()
