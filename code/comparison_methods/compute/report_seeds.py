"""Collect a multi-seed native compute run into tables, then draw the figure.

Reads, under --root (a run made by submit_seeds.sh):
  raw/<method>_n<N>_seed<S>/measurement.json   ALBIS and SPIDER, one per size x seed
  references/seed<S>/measurement.json           Splatter reference (SPIDER/scCube input)
  sccube/seed<S>_<batch>/{measurements,setup}.json and
  sccube/seed<S>_<batch>_rss_samples.csv        scCube, one process (own VAE) per size, 100 ms RSS samples
Writes measurements.csv (one row per method x size x seed), summary.csv (mean, SD,
min, max over seeds) and RESULTS.md, then runs plot_compute.py.

Time: simulation = the method's calls from 3D tissue to cell-level data in Z sections
(10 slices); setup = Splatter reference +
input loading + scCube VAE training (measured in each size's own run). Memory: ALBIS
and SPIDER run one size per process, so their process peak RSS is their generation peak;
scCube reports its process peak (includes VAE training) and the peak of the 100 ms samples taken during each size's generation.
"""
import argparse
import json
from pathlib import Path
import subprocess
import sys

import pandas as pd

GIB = 2 ** 30
METRICS = ['simulation_seconds', 'setup_seconds', 'total_seconds',
           'generation_peak_rss_gib', 'process_peak_rss_gib']

parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
parser.add_argument('--root', type=Path, required=True)
parser.add_argument('--no-plot', action='store_true')
args = parser.parse_args()
root = args.root.resolve()
protocol = json.loads((root / 'protocol.json').read_text())
sizes, seeds = protocol['sizes'], protocol['seeds']


def load(path):
    return json.loads(path.read_text()) if path.exists() else None


rows = []
for seed in seeds:
    reference = load(root / 'references' / f'seed{seed}' / 'measurement.json')
    reference_seconds = reference['generation_seconds'] if reference and reference['status'] == 'ok' else None
    for n in sizes:
        for method in ['albis', 'spider']:
            path = root / 'raw' / f'{method}_n{n}_seed{seed}' / 'measurement.json'
            row = dict(method=method, n_cells=n, seed=seed, status='missing', source=str(path.relative_to(root)))
            d = load(path)
            if d:
                assert (d['n_cells'], d['seed'], d['n_genes']) == (n, seed, protocol['n_genes']), path
                setup = 0.0 if method == 'albis' else d['input_seconds'] + (reference_seconds or float('nan'))
                row.update(status=d['status'], simulation_seconds=d['generation_seconds'], setup_seconds=setup,
                           process_peak_rss_gib=d['process_peak_rss_bytes'] / GIB,
                           generation_peak_rss_gib=d['process_peak_rss_bytes'] / GIB,
                           generation_memory='process peak (one size per process)')
            rows.append(row)
    for batch, batch_sizes in protocol['sccube_batches'].items():
        folder = root / 'sccube' / f'seed{seed}_{batch}'
        measured = {d['n_cells']: d for d in (load(folder / 'measurements.json') or [])}
        setup = load(folder / 'setup.json')
        samples_path = root / 'sccube' / f'seed{seed}_{batch}_rss_samples.csv'
        samples = pd.read_csv(samples_path) if samples_path.exists() else pd.DataFrame(columns=['phase', 'rss_bytes'])
        for n in batch_sizes:
            row = dict(method='sccube', n_cells=n, seed=seed, status='missing', batch=batch,
                       source=str((folder / 'measurements.json').relative_to(root)))
            d = measured.get(n)
            if d and setup:
                assert (d['seed'], d['n_genes']) == (seed, protocol['n_genes']), folder
                phase = samples.loc[samples.phase.eq(f'simulation_n{n}'), 'rss_bytes']
                row.update(status=d['status'], simulation_seconds=d['simulation_seconds'] + d['preparation_seconds'],
                           setup_seconds=setup['input_seconds'] + setup['vae_training_seconds']
                           + setup['reference_generation_seconds'],
                           vae_training_seconds=setup['vae_training_seconds'],
                           process_peak_rss_gib=d['cumulative_process_peak_rss_gib'],
                           generation_peak_rss_gib=phase.max() / GIB if len(phase) else float('nan'),
                           generation_rss_samples=len(phase), generation_memory='100 ms samples')
            rows.append(row)

data = pd.DataFrame(rows).sort_values(['method', 'n_cells', 'seed']).reset_index(drop=True)
data['total_seconds'] = data['simulation_seconds'] + data['setup_seconds']
data.to_csv(root / 'measurements.csv', index=False)

ok = data.loc[data.status.eq('ok')]
summary = ok.groupby(['method', 'n_cells'])[METRICS].agg(['mean', 'std', 'min', 'max'])
summary.columns = [f'{metric}_{stat}' for metric, stat in summary.columns]
summary.insert(0, 'n_seeds', ok.groupby(['method', 'n_cells']).size())
summary = summary.reset_index()
summary.to_csv(root / 'summary.csv', index=False)

lines = [f'# Native compute benchmark: {len(seeds)} seeds ({", ".join(map(str, seeds))})', '',
         f'Successful measurements: {len(ok)} / {len(data)}. {protocol["n_genes"]} genes, '
         f'{protocol["threads"]} CPU thread. Values are mean ± SD over seeds.', '',
         '| Method | Cells | Seeds | Simulation (s) | Including setup (s) | Generation peak RAM (GiB) | Process peak RAM (GiB) |',
         '|---|---:|---:|---:|---:|---:|---:|']
for r in summary.to_dict('records'):
    cell = lambda m, f: f"{r[m + '_mean']:{f}} ± {r[m + '_std']:{f}}"
    lines.append(f"| {r['method']} | {r['n_cells']:,} | {r['n_seeds']} | {cell('simulation_seconds', '.1f')} | "
                 f"{cell('total_seconds', '.1f')} | {cell('generation_peak_rss_gib', '.2f')} | {cell('process_peak_rss_gib', '.2f')} |")
lines += ['', __doc__.split('\n\n', 2)[2].strip(), '']
missing = data.loc[~data.status.eq('ok'), ['method', 'n_cells', 'seed', 'status']]
if len(missing):
    lines += ['Missing or failed measurements:', '', missing.to_string(index=False), '']
(root / 'RESULTS.md').write_text('\n'.join(lines))
print('\n'.join(lines))
if not args.no_plot:
    subprocess.run([sys.executable, str(Path(__file__).with_name('plot_compute.py')), '--root', str(root)], check=True)
