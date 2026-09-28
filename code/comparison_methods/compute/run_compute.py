#!/usr/bin/env python3
"""Run one scaling point with existing isolated workflow environments.

Defaults use 556 genes, random_null, CPU only, and volume scaled with cell
count; updated capture is 16um bins and 55um circular spots for ALBIS/SPIDER.
scCube retains its native occupancy-based aggregation.
Outputs are compatible with plot_compute.py and the original scaling summaries.
"""
from __future__ import annotations
import argparse
import json
import subprocess
import sys
import time
from pathlib import Path

PROJECT = Path(__file__).resolve().parents[4]
BENCH = PROJECT / 'comparison_methods'
sys.path.insert(0, str(BENCH / 'code'))
from common.contract import derive_seed, load_contract
from scaling.run_scaling_point import scaled_tissue_dimensions_um, parse_time_v

METHODS = {'albis': ('our_method', 'our_method', 'our_method.py'),
           'spider': ('splatter_spider', 'spider', 'spider_adapter.py'),
           'sccube': ('splatter_sccube', 'sccube', 'sccube_adapter.py')}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--n-cells', type=int, required=True)
    parser.add_argument('--technology', choices=['cell', 'bin', 'spot'], required=True)
    parser.add_argument('--method', nargs='+', choices=list(METHODS), default=list(METHODS))
    parser.add_argument('--seed', type=int, help='Default: deterministic original sweep seed')
    parser.add_argument('--out-dir', type=Path)
    parser.add_argument('--dry-run', action='store_true', help='Print plan without creating outputs')
    args = parser.parse_args()
    if args.n_cells <= 0:
        parser.error('--n-cells must be positive')
    out = (args.out_dir or PROJECT / 'sim_paper/data/comparison_methods/compute/raw' /
           f'n{args.n_cells}_{args.technology}').resolve()
    bc = load_contract()
    seed = args.seed if args.seed is not None else derive_seed(bc.base_seed, 'scaling', args.n_cells, args.technology)
    contract = dict(n_cells=args.n_cells, n_genes=556, n_cell_types=bc.n_cell_types,
                    cell_type_proportions=list(bc.cell_type_proportions),
                    spatial_dimensions=bc.spatial_dimensions,
                    tissue_dimensions_um=scaled_tissue_dimensions_um(args.n_cells, list(bc.tissue_dimensions_um)),
                    coordinate_units=bc.coordinate_units, technology=args.technology,
                    technology_params=bc.technologies[args.technology], sequencing_depth=bc.sequencing_depth,
                    scenario='random_null', organization_level=None)
    contract['technology_params'] = dict(contract['technology_params'])
    if args.technology == 'bin':
        contract['technology_params']['bin_size_um'] = 16.0
    if args.technology == 'spot':
        contract['technology_params'].update(spot_radius_um=27.5, spot_spacing_um=100.0)
        contract['spider_capture_implementation'] = 'pinned_circle_v1'
    contract['compute_geometry_version'] = 'bin16_circle55_spacing100'
    steps = []
    if any(m != 'albis' for m in args.method):
        steps.append(('splatter_expression', 'splatter', 'splatter_expression.R'))
    steps.extend(METHODS[m] for m in dict.fromkeys(args.method))
    commands = []
    for name, env, script in steps:
        executable = BENCH / 'env' / env / 'bin' / ('Rscript' if env == 'splatter' else 'python')
        cmd = [str(executable), str(BENCH / 'code/workflows' / script),
               '--contract', str(out / 'contract.json'), '--seed', str(seed), '--out-dir', str(out / name)]
        if name in ('splatter_spider', 'splatter_sccube'):
            cmd += ['--splatter-dir', str(out / 'splatter_expression')]
        commands.append((name, cmd))
    if args.dry_run:
        print(json.dumps({'out_dir': str(out), 'contract': contract, 'seed': seed, 'commands': commands}, indent=2))
        return
    if out.exists() and any(out.iterdir()):
        parser.error(f'Output directory is not empty: {out}. Choose a fresh --out-dir.')
    for _, cmd in commands:
        if not Path(cmd[0]).is_file() or not Path(cmd[1]).is_file():
            parser.error(f'Missing environment or adapter: {cmd[:2]}')
    out.mkdir(parents=True, exist_ok=True)
    (out / 'contract.json').write_text(json.dumps(contract, indent=2) + '\n')
    summary = dict(n_cells=args.n_cells, n_genes=556, technology=args.technology,
                   geometry_version=contract['compute_geometry_version'],
                   technology_params=contract['technology_params'],
                   scenario='random_null', seed=seed, steps=[])
    # Persist planned stages before starting: Slurm may kill the whole job.
    summary['steps'] = [dict(name=name, returncode=None, manifest_status='pending',
                             failure_reason=None, time_v={}, out_dir=str(out / name))
                        for name, _ in commands]
    def save_summary():
        temporary = out / 'scaling_summary.json.tmp'
        temporary.write_text(json.dumps(summary, indent=2) + '\n')
        temporary.replace(out / 'scaling_summary.json')
    save_summary()
    pool_ok = False
    for index, (name, cmd) in enumerate(commands):
        step_dir = out / name
        step_dir.mkdir(exist_ok=True)
        if name in ('splatter_spider', 'splatter_sccube') and not pool_ok:
            result = dict(name=name, returncode=None, manifest_status='skipped',
                          failure_reason='Splatter prerequisite failed', time_v={}, out_dir=str(step_dir))
        else:
            summary['steps'][index]['manifest_status'] = 'running'
            save_summary()
            print(f'Running {name}: {args.n_cells} cells, {args.technology}', flush=True)
            started = time.monotonic()
            with (step_dir / 'stdout.log').open('w') as stdout, (step_dir / 'stderr.log').open('w') as stderr:
                proc = subprocess.run(['/usr/bin/time', '-v', '-o', str(step_dir / 'time.txt')] + cmd,
                                      cwd=BENCH, stdout=stdout, stderr=stderr)
            timing_path = step_dir / 'time.txt'
            timing = parse_time_v(timing_path.read_text()) if timing_path.exists() else {}
            timing.setdefault('wall_clock_seconds', time.monotonic() - started)
            manifest_path = step_dir / 'manifest.json'
            manifest = json.loads(manifest_path.read_text()) if manifest_path.exists() else {}
            ok = proc.returncode == 0 and manifest.get('status') == 'ok'
            result = dict(name=name, returncode=proc.returncode, manifest_status='ok' if ok else 'failed',
                          failure_reason=manifest.get('failure_reason') or (None if ok else 'See stderr.log / missing successful manifest'),
                          adapter_runtime_seconds=manifest.get('runtime_seconds'),
                          adapter_peak_memory_bytes=manifest.get('peak_memory_bytes'),
                          time_v=timing, out_dir=str(step_dir))
            if name == 'splatter_expression':
                pool_ok = ok
        summary['steps'][index] = result
        save_summary()
        print(f"  {name}: {result['manifest_status']}", flush=True)
    if any(s['manifest_status'] != 'ok' for s in summary['steps']):
        sys.exit(1)


if __name__ == '__main__':
    main()
