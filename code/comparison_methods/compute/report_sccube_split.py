"""Validate exact native equivalence and plot the in-memory scCube pilot."""
import argparse
import json
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import pandas as pd

parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument('--root', type=Path, required=True)
parser.add_argument('--source', type=Path, required=True)
args = parser.parse_args()
root = args.root
combined = json.loads((root / 'combined/measurements.json').read_text())[0]
split = json.loads((root / 'split/measurements.json').read_text())
setup = json.loads((root / 'split/setup.json').read_text())
assert [r['n_cells'] for r in split] == [10000, 100000]
assert all(r['status'] == 'ok' for r in split)
assert split[0]['output_hashes'] == combined['output_hashes'], 'Native output equivalence failed'
validation = dict(status='passed', n_cells=10000, epochs=200,
                  exact_expression_and_metadata_hash_match=True,
                  output_hashes=combined['output_hashes'])
(root / 'validation.json').write_text(json.dumps(validation, indent=2))
rows = pd.DataFrame(split).drop(columns='output_hashes')
for name in ['input_seconds', 'reference_generation_seconds', 'vae_training_seconds']:
    rows[name] = setup[name]
rows['simulation_with_preparation_seconds'] = rows.simulation_seconds + rows.preparation_seconds
rows['first_use_seconds'] = rows[['simulation_with_preparation_seconds', 'input_seconds',
                                  'reference_generation_seconds', 'vae_training_seconds']].sum(axis=1)
rows.to_csv(root / 'measurements.csv', index=False)
baseline = pd.read_csv(args.source / 'measurements.csv')
baseline = baseline.loc[baseline.status.eq('ok') & baseline.method.isin(['albis', 'spider'])].copy()
baseline['first_use_seconds'] = baseline[['generation_seconds', 'input_seconds',
                                         'reference_generation_seconds']].sum(axis=1, min_count=3)
out = root / 'figures'
out.mkdir(exist_ok=True)
plt.rcParams.update({'font.size': 9, 'axes.titlesize': 10, 'axes.labelsize': 9,
                     'legend.fontsize': 8, 'axes.spines.top': False,
                     'axes.spines.right': False, 'savefig.dpi': 300,
                     'svg.fonttype': 'none', 'pdf.fonttype': 42})
colors = {'albis': '#4C8FD5', 'spider': '#8E63C7', 'sccube': '#8A949E'}
fig, axes = plt.subplots(2, 2, figsize=(10, 7.5))
for ax, total in zip(axes[0], [False, True]):
    for method, label in [('albis', 'ALBIS'), ('spider', 'SPIDER')]:
        d = baseline.loc[baseline.method.eq(method)].sort_values('n_cells')
        ax.plot(d.n_cells, d['first_use_seconds' if total else 'generation_seconds'],
                'o-', color=colors[method], label=label)
    ax.plot(rows.n_cells, rows['first_use_seconds' if total else 'simulation_with_preparation_seconds'],
            'o-', color=colors['sccube'], label='scCube (in-memory VAE)')
    ax.set(xscale='log', yscale='log', xlabel='Requested tissue cells', ylabel='Seconds (log scale)',
           title='First-use cost: setup + simulation' if total else 'Simulation with prepared inputs / VAE')
    ax.set_xticks(rows.n_cells, ['10k', '100k'])
    ax.legend(frameon=False)
    ax.grid(axis='y', color='#DDDDDD', lw=0.5)
ax = axes[1, 0]
stage_names = ['Splatter reference', 'Input / preprocessing', 'VAE training']
stage_values = [setup['reference_generation_seconds'], setup['input_seconds'], setup['vae_training_seconds']]
bars = ax.barh(stage_names, stage_values, color=['#B9C0C7', '#7393B3', colors['sccube']])
ax.bar_label(bars, labels=[f'{v:.2f} s' for v in stage_values], padding=4, fontsize=8)
ax.invert_yaxis()
ax.set(xlim=(0, max(stage_values)*1.3), xlabel='Seconds', title='scCube one-time setup · seed 2025')
ax = axes[1, 1]
bottom = pd.Series(0., index=rows.index)
for column, label, color in [('preparation_seconds', 'Native input preparation', '#B9C0C7'),
                             ('expression_seconds', 'Expression generation', '#7393B3'),
                             ('spatial_seconds', 'Spatial placement', '#8E63C7')]:
    ax.bar(['10k', '100k'], rows[column], bottom=bottom, label=label, color=color, width=0.55)
    bottom += rows[column]
for i, value in enumerate(bottom):
    ax.text(i, value, f'{value:.2f} s', ha='center', va='bottom', fontsize=8)
ax.set(ylim=(0, bottom.max()*1.45), xlabel='Requested tissue cells', ylabel='Seconds', title='scCube simulation stages')
ax.legend(frameon=False)
fig.suptitle('Preliminary native compute · VAE training measured separately', y=0.99)
fig.text(0.04, 0.015,
         '556 genes; one CPU thread; seed 2025. scCube trains once for 200 epochs; no model files saved.\n'
         'First-use totals charge setup once per size (accounting view, not independent retraining). Imports/startup excluded.\n'
         'Native representations differ. ALBIS/SPIDER points are from the original pilot. Lines are descriptive; no error bars.',
         fontsize=8)
fig.tight_layout(rect=(0, 0.11, 1, 0.96))
for ext in ['png', 'pdf', 'svg']:
    fig.savefig(out / f'native_split_runtime.{ext}', bbox_inches='tight')
plt.close(fig)

fig, axes = plt.subplots(1, 2, figsize=(9, 4.5))
for method, label in [('albis', 'ALBIS'), ('spider', 'SPIDER')]:
    d = baseline.loc[baseline.method.eq(method)].sort_values('n_cells')
    axes[0].plot(d.n_cells, d.process_peak_rss_gib, 'o-', color=colors[method], label=label)
axes[0].set(title='Independent method process peaks', ylabel='Peak RSS (GiB)')
axes[0].legend(frameon=False)
axes[1].plot(rows.n_cells, rows.cumulative_process_peak_rss_gib, 'o-', color=colors['sccube'])
axes[1].set(title='scCube cumulative process high-water mark', ylabel='Cumulative peak RSS (GiB)')
for ax in axes:
    ax.set_xscale('log')
    ax.set_xticks(rows.n_cells, ['10k', '100k'])
    ax.set(xlabel='Requested tissue cells', ylim=(0, None))
    ax.grid(axis='y', color='#DDDDDD', lw=0.5)
fig.text(0.04, 0.02,
         'scCube uses one process: training → 10k → 100k, retaining the VAE. Peaks include prior work and validation.\n'
         'These are not independent per-size scCube memory measurements. Separate Splatter process excluded.', fontsize=8)
fig.tight_layout(rect=(0, 0.14, 1, 1))
for ext in ['png', 'pdf', 'svg']:
    fig.savefig(out / f'native_split_memory.{ext}', bbox_inches='tight')
plt.close(fig)
(root / 'RESULTS.md').write_text(
    '# scCube split-timing pilot\n\n'
    'Exact output-hash comparison against the combined native API passed at 10k cells, '
    '200 epochs, seed 2025. Expression, metadata, axes and dtypes were checked.\n\n'
    f"VAE training: {setup['vae_training_seconds']:.3f} seconds. "
    f"Input/preprocessing: {setup['input_seconds']:.3f} seconds. "
    f"Existing shared Splatter reference generation: {setup['reference_generation_seconds']:.3f} seconds.\n\n"
    + '\n'.join(f"- {r['n_cells']:,} cells: expression {r['expression_seconds']:.3f} s; "
                f"spatial placement {r['spatial_seconds']:.3f} s; "
                f"native preparation {r['preparation_seconds']:.3f} s." for r in split)
    + '\n\nThe VAE was trained once and retained in memory; no model files were saved. '
    'The 10k generation follows the training RNG state for equivalence validation; '
    '100k generation resets RNGs to seed 2025. First-use totals charge setup once '
    'per size but do not represent independent fresh training runs. '
    'Only generation-stage timings are compared in the simulation panel. '
    'Memory is cumulative process high-water RSS, including earlier stages and '
    'validation allocations, not isolated generation or independent per-size peaks. '
    'ALBIS/SPIDER measurements are reused from the direct-native pilot only. '
    'Imports/startup, aggregation and exports are excluded. No uncertainty estimate '
    'is available from one seed. Larger sizes have not yet been run.\n')
print((root / 'RESULTS.md').read_text())
