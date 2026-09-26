#!/usr/bin/env python
"""Submit shared Figure 3 inputs, with exact generation-array dependencies.

Default is dry-run. --submit writes an incremental job receipt.
"""
import argparse
import csv
import json
import subprocess
from datetime import datetime
from pathlib import Path
from stagate_inputs import DATASETS, PROJECT


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--seed-array', required=True)
    ap.add_argument('--control-array', required=True)
    ap.add_argument('--submit', action='store_true')
    args = ap.parse_args()
    here = Path(__file__).resolve().parent
    tasks = PROJECT / 'sim_paper/code/clustering/strong_domain_mix/batch_sigma_slide/batch_sigma_slide_seed_tasks.tsv'
    with tasks.open() as f:
        rows = list(csv.DictReader(f, delimiter='\t'))
    plan = []
    for name, spec in DATASETS.items():
        path = Path(spec['h5ad'])
        seed = spec['simulation_seed']
        if seed == 2025:
            if not path.is_file():
                raise FileNotFoundError(path)
            dependency = 'none'
        else:
            mod = 'bin' if spec['technology'] == 'bin16um' else spec['technology']
            matches = [i for i, r in enumerate(rows) if r['modality'] == mod
                       and float(r['batch_sigma']) == spec['batch_sigma'] and r['seed'] == str(seed)]
            if len(matches) != 1:
                raise ValueError((name, matches))
            i = matches[0]
            dependency = f'{args.seed_array if i < 68 else args.control_array}_{i}'
        plan.append(dict(dataset=name, input=str(path), dependency=dependency))
    receipt = here / 'logs_stagate_3D' / f'strongmix_grid_{datetime.now():%Y%m%d_%H%M%S}.json'
    receipt.parent.mkdir(exist_ok=True)
    submitted = []
    for run in plan:
        if args.submit:
            result = subprocess.run(['bash', str(here / 'run_stagate_3D.sh'), run['dependency'], run['dataset']],
                                    check=True, text=True, capture_output=True)
            run['job_id'] = result.stdout.strip().split()[-1]
            submitted.append(run)
            receipt.write_text(json.dumps(submitted, indent=2) + '\n')
        print(json.dumps(run), flush=True)
    if args.submit:
        print(f'Receipt: {receipt}')


if __name__ == '__main__':
    main()
