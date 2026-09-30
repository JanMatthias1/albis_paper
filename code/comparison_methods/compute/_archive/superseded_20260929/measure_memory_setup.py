"""External process RSS measurement around unchanged direct native API workers."""
import argparse
import csv
import json
import os
from pathlib import Path
import subprocess
import time

parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument('--root', type=Path, required=True)
parser.add_argument('--task', type=int, required=True)
args = parser.parse_args()
root = args.root.resolve()
protocol = json.loads((root / 'protocol.json').read_text())
method = protocol['methods'][args.task % 3]
n = protocol['sizes'][args.task // 3]
work = root / 'raw' / f'{method}_n{n}_seed2025'
work.mkdir(exist_ok=False)
project = Path('/dcs04/hicks/data/Jan/sim_project')
envs = project / 'comparison_methods/env'
phase = work / 'phase.txt'
settings = root / 'settings' / f'n{n}_seed2025.json'
reference = work / 'references/seed2025'
env = dict(os.environ, SETTINGS=str(settings), OUTPUT=str(work/'native'),
           REFERENCE=str(reference), PHASE_FILE=str(phase),
           MPLCONFIGDIR=str(work/'cache/mpl'), XDG_CACHE_HOME=str(work/'cache'),
           NUMBA_CACHE_DIR=str(work/'cache/numba'))
commands = []
if method != 'albis':
    commands.append(('reference', [str(envs/'splatter/bin/Rscript'), str(root/'code/reference_native.R')]))
if method == 'sccube':
    command = [str(envs/'sccube/bin/python'), '-u', str(root/'code/sccube_split.py'),
               '--root', str(work), '--source', str(work), '--mode', 'split', '--sizes', str(n),
               '--settings-dir', str(root/'settings'), '--skip-output-hashes', '--phase-file', str(phase)]
else:
    environment = 'our_method' if method == 'albis' else 'spider'
    command = [str(envs/environment/'bin/python'), '-u', str(root/'code'/f'{method}_native.py')]
commands.append((method, command))
peaks = {}
counts = {}
page = os.sysconf('SC_PAGE_SIZE')
with (work/'rss_samples.csv').open('w', buffering=1) as output:
    writer = csv.writer(output)
    writer.writerow(['monotonic_seconds','process','phase','rss_bytes'])
    for stage, command in commands:
        phase.write_text('imports')
        with (work/f'{stage}.log').open('w') as log:
            child = subprocess.Popen(command, env=env, stdout=log, stderr=subprocess.STDOUT)
            while child.poll() is None:
                try:
                    current = phase.read_text().strip()
                    rss = int(Path(f'/proc/{child.pid}/statm').read_text().split()[1])*page
                except (FileNotFoundError, ProcessLookupError):
                    continue
                writer.writerow([time.monotonic(), stage, current, rss])
                peaks[current] = max(peaks.get(current, 0), rss)
                counts[current] = counts.get(current, 0)+1
                time.sleep(protocol['interval_seconds'])
            if child.returncode:
                raise SystemExit(f'{stage} failed ({child.returncode}); see {work}/{stage}.log')
generation_phase = f'simulation_n{n}' if method == 'sccube' else 'generation'
assert counts.get(generation_phase, 0), 'No generation samples captured'
reference_seconds = json.loads((reference/'measurement.json').read_text())['generation_seconds'] if method != 'albis' else 0.
if method == 'sccube':
    result = json.loads((work/'split/measurements.json').read_text())[0]
    setup = json.loads((work/'split/setup.json').read_text())
    generation_seconds = result['simulation_seconds']+result['preparation_seconds']
    input_seconds = setup['input_seconds']
    training_seconds = setup['vae_training_seconds']
else:
    result = json.loads((work/'native/measurement.json').read_text())
    generation_seconds = result['generation_seconds']
    input_seconds = result.get('input_seconds', 0.)
    training_seconds = 0.
relevant = ['reference','setup','training',generation_phase]
record = dict(status='ok',method=method,n_cells=n,seed=2025,
    simulation_seconds=generation_seconds, input_seconds=input_seconds,
    training_seconds=training_seconds, reference_seconds=reference_seconds,
    total_seconds=generation_seconds+input_seconds+training_seconds+reference_seconds,
    generation_peak_rss_gib=peaks[generation_phase]/2**30,
    including_setup_peak_rss_gib=max(peaks.get(p,0) for p in relevant)/2**30,
    phase_peak_rss_bytes=peaks, phase_sample_counts=counts,
    sample_interval_seconds=protocol['interval_seconds'], job_id=os.environ.get('SLURM_JOB_ID'),
    source=str(work), native_measurement=result)
(work/'measurement.json').write_text(json.dumps(record,indent=2)+'\n')
print(json.dumps(record),flush=True)
