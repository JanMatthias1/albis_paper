"""Combine direct-native pilots and the single-seed extension; plot measurements only."""
import argparse
import json
from pathlib import Path
import subprocess
import sys

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import pandas as pd

parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument('--root', type=Path, required=True)
parser.add_argument('--extra-plots', action='store_true',
                    help='Also regenerate the detailed stage, setup, memory and molecule-count figures.')
args = parser.parse_args()
root = args.root
protocol = json.loads((root / 'protocol.json').read_text())
source = Path(protocol['source_pilot'])
pilot = Path(protocol['validated_split'])
assert json.loads((pilot / 'validation.json').read_text())['status'] == 'passed'
base = pd.read_csv(source / 'measurements.csv')
base = base.loc[base.status.eq('ok') & base.method.isin(['albis', 'spider'])].copy()
records = []
for row in base.to_dict('records'):
    records.append(dict(method=row['method'], n_cells=row['n_cells'], seed=row['seed'], status='ok',
                        simulation_seconds=row['generation_seconds'], input_seconds=row['input_seconds'],
                        reference_seconds=row['reference_generation_seconds'], training_seconds=0.,
                        process_peak_rss_gib=row['process_peak_rss_gib'], n_molecules=row['n_molecules'],
                        batch='original pilot', memory_scope='independent method process', source=row['source']))
ref_seconds = json.loads((source / 'references/seed2025/measurement.json').read_text())['generation_seconds']
for n in protocol['sizes']:
    for method in ['albis', 'spider']:
        path = root / 'raw' / f'{method}_n{n}_seed2025' / 'measurement.json'
        record = dict(method=method, n_cells=n, seed=2025, batch='extension', source=str(path),
                      memory_scope='independent method process', status='missing_measurement')
        if path.exists():
            d = json.loads(path.read_text())
            assert d['n_genes'] == 556 and d['n_cells'] == n
            record.update(status=d['status'], simulation_seconds=d['generation_seconds'],
                          input_seconds=0. if method == 'albis' else d['input_seconds'],
                          reference_seconds=ref_seconds if method == 'spider' else 0.,
                          training_seconds=0., process_peak_rss_gib=d['process_peak_rss_bytes'] / 2**30,
                          n_molecules=d.get('n_molecules'))
        records.append(record)
setups = []
for folder, batch, sizes in [(pilot, '10k–100k pilot', [10000, 100000]),
                              (root, '200k–600k extension', protocol['sizes'])]:
    path = folder / 'split/measurements.json'
    measured = {d['n_cells']: d for d in json.loads(path.read_text())} if path.exists() else {}
    setup_path = folder / 'split/setup.json'
    setup = json.loads(setup_path.read_text()) if setup_path.exists() else None
    if setup:
        setups.append(dict(batch=batch, **setup))
    for n in sizes:
        record = dict(method='sccube', n_cells=n, seed=2025, batch=batch, source=str(path),
                      memory_scope='cumulative training and generation process', status='missing_measurement')
        if n in measured and setup:
            d = measured[n]
            assert d['n_genes'] == 556
            record.update(status=d['status'], simulation_seconds=d['simulation_seconds'] + d['preparation_seconds'],
                          expression_seconds=d['expression_seconds'], spatial_seconds=d['spatial_seconds'],
                          preparation_seconds=d['preparation_seconds'], input_seconds=setup['input_seconds'],
                          reference_seconds=setup['reference_generation_seconds'],
                          training_seconds=setup['vae_training_seconds'],
                          process_peak_rss_gib=d['cumulative_process_peak_rss_gib'])
        records.append(record)
data = pd.DataFrame(records).sort_values(['method', 'n_cells'])
data['first_use_seconds'] = data[['simulation_seconds', 'input_seconds', 'reference_seconds',
                                 'training_seconds']].sum(axis=1, min_count=4)
data.to_csv(root / 'measurements.csv', index=False)
pd.DataFrame(setups).to_csv(root / 'sccube_setup.csv', index=False)
ok = data.loc[data.status.eq('ok')]
out = root / 'figures'
out.mkdir(exist_ok=True)
methods = {'albis': ('ALBIS', '#4C8FD5'), 'spider': ('SPIDER', '#8E63C7'), 'sccube': ('scCube', '#8A949E')}
# Optional figures use the same retained detailed measurements.
if args.extra_plots:
    plt.rcParams.update({'font.size': 9, 'axes.titlesize': 10, 'axes.labelsize': 9,
                         'legend.fontsize': 8, 'axes.spines.top': False,
                         'axes.spines.right': False, 'savefig.dpi': 300,
                         'svg.fonttype': 'none', 'pdf.fonttype': 42})
    sizes = [10000, 100000] + protocol['sizes']
    ticklabels = [f'{n//1000}k' for n in sizes]


    def save(fig, name):
        for ext in ['png', 'pdf', 'svg']:
            fig.savefig(out / f'{name}.{ext}', bbox_inches='tight')
        plt.close(fig)


    fig, axes = plt.subplots(2, 2, figsize=(11, 8))
    for ax, column, title in [(axes[0, 0], 'simulation_seconds', 'Simulation with prepared inputs / VAE'),
                               (axes[0, 1], 'first_use_seconds', 'First-use cost: setup + simulation')]:
        for method, (label, color) in methods.items():
            d = ok.loc[ok.method.eq(method)]
            ax.plot(d.n_cells, d[column], 'o-', label=label, color=color)
        ax.set(xscale='log', yscale='log', xlabel='Requested tissue cells', ylabel='Seconds (log scale)', title=title)
        ax.set_xticks(sizes, ticklabels)
        ax.tick_params(axis='x', labelrotation=35)
        ax.legend(frameon=False)
        ax.grid(axis='y', color='#DDDDDD', lw=0.5)
    ax = axes[1, 0]
    for i, setup in enumerate(setups):
        values = [setup['reference_generation_seconds'], setup['input_seconds'], setup['vae_training_seconds']]
        ax.barh(i, values[0], color='#B9C0C7', label='Splatter' if i == 0 else None)
        ax.barh(i, values[1], left=values[0], color='#7393B3', label='Preprocessing' if i == 0 else None)
        ax.barh(i, values[2], left=sum(values[:2]), color='#8A949E', label='VAE training' if i == 0 else None)
        ax.text(sum(values), i, f"  VAE: {values[2]/60:.2f} min", va='center', fontsize=8)
    ax.set_yticks(range(len(setups)), [s['batch'] for s in setups])
    if setups:
        ax.set_xlim(0, max(sum(s[k] for k in ['reference_generation_seconds', 'input_seconds', 'vae_training_seconds']) for s in setups)*1.38)
    ax.set(title='scCube setup: one training run per batch', xlabel='Seconds')
    ax.legend(frameon=False, loc='upper left', bbox_to_anchor=(0, -0.2), ncol=3)
    ax = axes[1, 1]
    sc = ok.loc[ok.method.eq('sccube')]
    bottom = pd.Series(0., index=sc.index)
    for column, label, color in [('preparation_seconds', 'Native input preparation', '#B9C0C7'),
                                 ('expression_seconds', 'Expression generation', '#7393B3'),
                                 ('spatial_seconds', 'Spatial placement', '#8E63C7')]:
        ax.bar([f'{n//1000}k' for n in sc.n_cells], sc[column], bottom=bottom, label=label, color=color)
        bottom += sc[column]
    for i, value in enumerate(bottom):
        ax.text(i, value, f'{value:.1f}', ha='center', va='bottom', fontsize=8)
    ax.set(title='scCube simulation stages', xlabel='Requested tissue cells', ylabel='Seconds',
           ylim=(0, bottom.max()*1.4 if len(bottom) else 1))
    ax.legend(frameon=False)
    fig.suptitle('Native compute scaling · 556 genes · one CPU thread · seed 2025', y=0.99)
    fig.text(0.03, 0.015,
             'VAE training is separate; no model files saved. First-use totals charge the batch setup once per size; imports/startup excluded.\n'
             'Two scCube training batches ran on different nodes; setup-time differences are not tissue-size scaling. No seed replicates.\n'
             'Native outputs differ (ALBIS molecules, scCube cell expression, SPIDER native view). No aggregation/export; lines are descriptive.', fontsize=8)
    fig.tight_layout(rect=(0, 0.13, 1, 0.96))
    save(fig, 'native_extension_runtime')

    fig, axes = plt.subplots(1, 3, figsize=(12, 4.6))
    for method in ['albis', 'spider']:
        label, color = methods[method]
        d = ok.loc[ok.method.eq(method)]
        axes[0].plot(d.n_cells, d.process_peak_rss_gib, 'o-', color=color, label=label)
    axes[0].set(title='Independent method process peaks', ylabel='Peak RSS (GiB)')
    axes[0].legend(frameon=False)
    for batch, d in sc.groupby('batch', sort=False):
        axes[1].plot(d.n_cells, d.process_peak_rss_gib, 'o-' if 'pilot' in batch else 's--',
                     color=methods['sccube'][1], label=batch)
    axes[1].set(title='scCube cumulative process peaks', ylabel='Cumulative peak RSS (GiB)')
    axes[1].legend(frameon=False)
    d = ok.loc[ok.method.eq('albis')]
    axes[2].plot(d.n_cells, d.n_molecules / 1e6, 'o-', color=methods['albis'][1])
    axes[2].set(title='ALBIS molecular output', ylabel='Molecules (millions)')
    for ax in axes:
        ax.set_xscale('log')
        ax.set_xticks(sizes, ticklabels)
        ax.tick_params(axis='x', labelrotation=45)
        ax.set(xlabel='Requested tissue cells', ylim=(0, None))
        ax.grid(axis='y', color='#DDDDDD', lw=0.5)
    fig.text(0.03, 0.02,
             'scCube peaks include training, prior generation and validation in each batch; they are not independent per-size memory measurements.\n'
             'Separate Splatter process excluded. One seed only; failed or missing measurements are not plotted as successes.', fontsize=8)
    fig.tight_layout(rect=(0, 0.17, 1, 1))
    save(fig, 'native_extension_memory')

lines = ['# Native compute extension: 200k, 400k, 600k', '',
         f"Successful measurements: {len(ok)} / {len(data)} (including previous 10k/100k pilots).", '',
         '| Method | Cells | Simulation (s) | First-use accounting (s) | Process peak RSS (GiB) |',
         '|---|---:|---:|---:|---:|']
for r in ok.to_dict('records'):
    lines.append(f"| {methods[r['method']][0]} | {r['n_cells']:,} | {r['simulation_seconds']:.3f} | {r['first_use_seconds']:.3f} | {r['process_peak_rss_gib']:.3f} |")
lines += ['', 'scCube setup:', '']
for s in setups:
    lines.append(f"- {s['batch']}: VAE {s['vae_training_seconds']:.3f} s; input/preprocessing {s['input_seconds']:.3f} s; Splatter {s['reference_generation_seconds']:.3f} s.")
lines += ['', 'Seed 2025 only. scCube trains once per batch, with the VAE held in memory and no model files saved. '
          'The new batch starts at 200k; its first generation follows training RNG state and later sizes reset generation RNGs to 2025. '
          'The prior 10k native-equivalence validation passed; method calls are unchanged. '
          'First-use accounting charges one setup cost to each size, not an independent training run. '
          'scCube memory is cumulative within each batch, including earlier stages and validation. '
          'ALBIS/SPIDER memory is from independent processes. Native outputs differ; SPIDER retains its expression view. '
          'Reference generation is reused and its measured cost reported separately. '
          'Runtime excludes imports/startup, aggregation and export. Training timings from different nodes are not a tissue-size trend.', '']
missing = data.loc[~data.status.eq('ok')]
if not missing.empty:
    lines += ['Missing or failed measurements (inspect scheduler/logs):', missing[['method', 'n_cells', 'status']].to_string(index=False)]
(root / 'RESULTS.md').write_text('\n'.join(lines))
subprocess.run([sys.executable, str(Path(__file__).with_name('plot_compute_three_panel.py')),
                '--root', str(root)], check=True)
print('\n'.join(lines))
