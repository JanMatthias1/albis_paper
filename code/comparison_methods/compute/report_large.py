"""Preserve prior detail, add the staged large runs and sampled scCube generation RSS."""
import argparse
import json
from pathlib import Path
import subprocess
import sys

import pandas as pd

parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument('--root', type=Path, required=True)
args = parser.parse_args()
root = args.root.resolve()
protocol = json.loads((root / 'protocol.json').read_text())
previous = Path(protocol['previous_run'])
source = Path(protocol['source_pilot'])
base = pd.read_csv(previous / 'measurements.csv')
# Resolve source provenance after archival without modifying historical tables.
relocation = previous.parent / 'relocation.json'
if relocation.exists():
    for old, new in json.loads(relocation.read_text())['paths'].items():
        base['source'] = base['source'].str.replace(old, new, regex=False)
        relative_old = old[old.index('sim_paper/'):] if 'sim_paper/' in old else old
        relative_new = new[new.index('sim_paper/'):] if 'sim_paper/' in new else new
        base['source'] = base['source'].str.replace(relative_old+'/', relative_new+'/', regex=False)
setup_table = pd.read_csv(previous / 'sccube_setup.csv')
ref_seconds = json.loads((source / 'references/seed2025/measurement.json').read_text())['generation_seconds']
receipt = pd.read_csv(root / 'submissions.tsv', sep='\t', dtype=str)
job_ids = receipt.loc[receipt.method.ne('report'), 'job_id'].tolist()
accounting = {}
try:
    result = subprocess.run(['sacct', '-j', ','.join(job_ids), '--format=JobID,State,ExitCode,Elapsed,MaxRSS,NodeList',
                             '-P'], check=True, capture_output=True, text=True, timeout=20)
    (root / 'slurm_accounting.psv').write_text(result.stdout)
    import csv
    import io
    accounting = {r['JobID']: r for r in csv.DictReader(io.StringIO(result.stdout), delimiter='|') if '.' not in r['JobID']}
except (OSError, subprocess.SubprocessError) as exc:
    (root / 'accounting_unavailable.txt').write_text(str(exc) + '\n')
sampling = pd.read_csv(root / 'rss_samples.csv') if (root / 'rss_samples.csv').exists() else pd.DataFrame()
rows = []
for n in protocol['sizes']:
    for method in ['albis', 'spider']:
        p = root / 'raw' / f'{method}_n{n}_seed2025' / 'measurement.json'
        submitted = receipt.loc[receipt.method.eq(method) & receipt.n_cells.eq(str(n))]
        row = dict(method=method, n_cells=n, seed=2025, batch='1M–5M staged run',
                   source=str(p), memory_scope='independent method process',
                   status='submitted_without_measurement' if len(submitted) else 'not_submitted',
                   job_id=submitted.job_id.iloc[-1] if len(submitted) else None)
        if p.exists():
            d = json.loads(p.read_text())
            assert d['n_cells'] == n and d['n_genes'] == 556
            row.update(status=d['status'], simulation_seconds=d['generation_seconds'],
                       input_seconds=0. if method == 'albis' else d['input_seconds'],
                       reference_seconds=ref_seconds if method == 'spider' else 0., training_seconds=0.,
                       process_peak_rss_gib=d['process_peak_rss_bytes']/2**30,
                       n_molecules=d.get('n_molecules'), expression_is_view=d.get('expression_is_view'))
        rows.append(row)
p = root / 'split/measurements.json'
sc = {d['n_cells']: d for d in json.loads(p.read_text())} if p.exists() else {}
p = root / 'split/setup.json'
setup = json.loads(p.read_text()) if p.exists() else None
if setup:
    setup_table = pd.concat([setup_table, pd.DataFrame([dict(batch='1M–5M staged run', **setup)])], ignore_index=True)
for n in protocol['sizes']:
    row = dict(method='sccube', n_cells=n, seed=2025, batch='1M–5M staged run',
               source=str(root / 'split/measurements.json'),
               memory_scope='cumulative training and generation process; additional sampled generation RSS',
               status='submitted_without_measurement',
               job_id=receipt.loc[receipt.method.eq('sccube'), 'job_id'].iloc[-1])
    if n in sc and setup:
        d = sc[n]
        row.update(status=d['status'], simulation_seconds=d['simulation_seconds'] + d['preparation_seconds'],
                   expression_seconds=d['expression_seconds'], spatial_seconds=d['spatial_seconds'],
                   preparation_seconds=d['preparation_seconds'], input_seconds=setup['input_seconds'],
                   reference_seconds=setup['reference_generation_seconds'], training_seconds=setup['vae_training_seconds'],
                   process_peak_rss_gib=d['cumulative_process_peak_rss_gib'])
        if not sampling.empty:
            samples = sampling.loc[sampling.phase.eq(f'simulation_n{n}')]
            if len(samples):
                row.update(generation_sampled_peak_rss_gib=samples.rss_bytes.max()/2**30,
                           generation_sampled_start_rss_gib=samples.rss_bytes.iloc[0]/2**30,
                           generation_rss_sample_count=len(samples), generation_rss_interval_seconds=0.1)
    rows.append(row)
new = pd.DataFrame(rows)
new['first_use_seconds'] = new[['simulation_seconds', 'input_seconds', 'reference_seconds',
                              'training_seconds']].sum(axis=1, min_count=4)
data = pd.concat([base, new], ignore_index=True).sort_values(['method', 'n_cells'])
for index, row in data.iterrows():
    record = accounting.get(str(row.get('job_id')))
    if record:
        data.loc[index, 'scheduler_state'] = record['State']
        data.loc[index, 'scheduler_exit_code'] = record['ExitCode']
        data.loc[index, 'scheduler_elapsed'] = record['Elapsed']
        data.loc[index, 'scheduler_node'] = record['NodeList']
        if row['status'] == 'submitted_without_measurement' and record['State'] not in ['PENDING', 'RUNNING', 'COMPLETING']:
            data.loc[index, 'status'] = record['State'].lower() + '_without_measurement'
data.to_csv(root / 'measurements.csv', index=False)
setup_table.to_csv(root / 'sccube_setup.csv', index=False)
ok = data.loc[data.status.eq('ok')]
lines = ['# Native single-seed scaling through 5M cells', '',
         f'Successful measurements: {len(ok)} / {len(data)} (including previous runs).', '',
         '| Method | Cells | Simulation (s) | Setup-inclusive (s) | Cumulative/process RSS (GiB) | Sampled generation RSS (GiB) |',
         '|---|---:|---:|---:|---:|---:|']
for r in ok.to_dict('records'):
    sampled = r.get('generation_sampled_peak_rss_gib')
    sampled_text = f'{sampled:.3f}' if pd.notna(sampled) else '—'
    lines.append(f"| {r['method']} | {r['n_cells']:,} | {r['simulation_seconds']:.3f} | {r['first_use_seconds']:.3f} | {r['process_peak_rss_gib']:.3f} | {sampled_text} |")
lines += ['', 'Setup by batch:', '']
for s in setup_table.to_dict('records'):
    lines.append(f"- {s['batch']}: VAE {s['vae_training_seconds']:.3f} s; preprocessing {s['input_seconds']:.3f} s; Splatter {s['reference_generation_seconds']:.3f} s.")
lines += ['', 'Seed 2025 only, 556 genes, one computational CPU thread. scCube retains one VAE in memory per batch; no model files. '
          'Native representations differ; SPIDER retains its reference-backed expression view. ALBIS generates explicit molecules. '
          'Timing excludes imports/startup, validation, aggregation and export. Setup totals charge each size one batch setup cost. '
          'The first generation in a scCube batch follows training RNG state; subsequent sizes reset RNGs to 2025.', '',
          'Memory: ALBIS/SPIDER use independent whole-process high-water RSS. Prior scCube batches report cumulative high-water RSS. '
          'The 1M–5M batch additionally samples worker RSS externally every 100 ms during generation, excluding training/validation intervals '
          'but including the resident VAE, reference, and allocations retained from earlier work. This is a sampled process peak, not '
          'isolated incremental memory or a process-tree total; short peaks may be missed. Raw samples and cumulative peaks are retained. '
          'No subtraction of high-water marks is used.', '']
missing = data.loc[~data.status.eq('ok')]
if len(missing):
    lines += ['Incomplete measurements:', missing[['method', 'n_cells', 'status', 'job_id']].to_string(index=False), '']
(root / 'RESULTS.md').write_text('\n'.join(lines))
subprocess.run([sys.executable, str(Path(__file__).with_name('plot_compute_four_panel_preliminary.py')),
                '--root', str(root)], check=True)
print('\n'.join(lines))
