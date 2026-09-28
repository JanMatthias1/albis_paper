"""Reconcile interrupted compute stages with Slurm; retain successful measurements."""
import argparse
import json
from pathlib import Path
import sys
sys.path.insert(0, str(Path(__file__).resolve().parents[4] / 'comparison_methods/code'))
from common.contract import derive_seed, load_contract


def finalize(root):
    plan = json.loads((root/'stress_plan.json').read_text())
    accounting = {}
    for line in (root/'slurm_accounting.psv').read_text().splitlines():
        fields = line.split('|')
        job = fields[0]
        if '.' in job or '_' not in job:
            continue
        index = job.rsplit('_', 1)[1]
        if index.isdigit():
            accounting[int(index)] = dict(state=fields[1], exit_code=fields[2], elapsed_seconds=fields[3])
    for task in plan['tasks']:
        directory = root/'raw'/task['directory']
        directory.mkdir(exist_ok=True)
        path = directory/'scaling_summary.json'
        names = (['splatter_expression'] if task['method'] != 'albis' else []) + [task['step']]
        summary = json.loads(path.read_text()) if path.exists() else dict(
            n_cells=task['n_cells'], n_genes=556, technology=task['technology'], scenario='random_null',
            geometry_version=plan.get('geometry_version', 'legacy_bin20_square100'),
            seed=derive_seed(load_contract().base_seed, 'scaling', task['n_cells'], task['technology']),
            steps=[dict(name=n, manifest_status='pending', returncode=None, time_v={}, out_dir=str(directory/n)) for n in names])
        resource = accounting.get(task['index'], dict(state='ACCOUNTING_MISSING'))
        state = resource['state'].split()[0].rstrip('+')
        category = {'OUT_OF_MEMORY':'out_of_memory', 'TIMEOUT':'timeout', 'CANCELLED':'cancelled',
                    'NODE_FAIL':'infrastructure_failure', 'PREEMPTED':'infrastructure_failure'}.get(state)
        summary['slurm'] = resource
        bad = [s for s in summary['steps'] if s['manifest_status'] != 'ok']
        for step in bad:
            previous = step.get('failure_reason')
            status = category or ('accounting_missing' if state == 'ACCOUNTING_MISSING' else
                                 'algorithm_error' if previous and step['manifest_status'] == 'failed' else 'incomplete')
            if step['manifest_status'] == 'skipped' and not category:
                status = 'prerequisite_failed'
            step['manifest_status'] = status
            step['failure_reason'] = f'Slurm {state}; ' + (previous or 'Stage did not complete successfully')
        if bad:
            # Explain shared-pool failures in the method row exported by the plotter.
            target = summary['steps'][-1]
            target['failure_reason'] = '; '.join(f"{s['name']}: {s['failure_reason']}" for s in bad)
        tmp = path.with_suffix('.json.tmp')
        tmp.write_text(json.dumps(summary, indent=2)+'\n')
        tmp.replace(path)

if __name__ == '__main__':
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--root', type=Path, required=True)
    finalize(p.parse_args().root)
