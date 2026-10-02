"""Approved direct upstream calls; train once, reuse in memory, save no models. Cells are
placed in 3D and sectioned along Z natively (generate_pattern_random, is_split=True)."""
import argparse
import gc
import hashlib
import json
import os
from pathlib import Path
import random
import resource
import time

import numpy as np
import pandas as pd
from scipy.io import mmread
import torch
import scCube.sccube as native_module
import scCube.utils as native_utils
from scCube.sccube import scCube
from scCube.utils import train_vae, generate_vae

parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument('--root', type=Path, required=True)
parser.add_argument('--source', type=Path, required=True)
parser.add_argument('--mode', choices=['combined', 'split'], required=True)
parser.add_argument('--sizes', type=int, nargs='+', default=[10000, 100000])
parser.add_argument('--settings-dir', type=Path)
parser.add_argument('--skip-output-hashes', action='store_true',
                    help='For scaling after equivalence validation; retain shape/finite checks.')
parser.add_argument('--phase-file', type=Path,
                    help='Optional phase labels for an external read-only RSS sampler.')
parser.add_argument('--memory-budget-gib', type=float,
                    help='Check measured high-water RSS before advancing to each larger size.')
parser.add_argument('--seed', type=int, default=2025, help='Simulation seed of the settings files to read')
parser.add_argument('--out', type=Path, help='Output folder (default: ROOT/MODE)')
args = parser.parse_args()
out = args.out or args.root / args.mode
out.mkdir(parents=True, exist_ok=False)
settings_dir = args.settings_dir or args.source / 'settings'
settings = [json.loads((settings_dir / f'n{n}_seed{args.seed}.json').read_text())
            for n in args.sizes]
seed = settings[0]['seed']
random.seed(seed); np.random.seed(seed); torch.manual_seed(seed)
torch.set_num_threads(1); torch.set_num_interop_threads(1)
reference_dir = args.source / 'references' / f'seed{seed}'
provenance = dict(mode=args.mode, seed=seed, settings=settings,
                  job_id=os.environ.get('SLURM_JOB_ID'), host=os.uname().nodename,
                  torch_version=torch.__version__, model_files_saved=False,
                  output_hashes_enabled=not args.skip_output_hashes,
                  source_hashes={str(p): hashlib.sha256(p.read_bytes()).hexdigest()
                                 for p in [Path(__file__), Path(native_module.__file__),
                                           Path(native_utils.__file__), reference_dir / 'counts.mtx']})
(out / 'protocol.json').write_text(json.dumps(provenance, indent=2))
if args.phase_file:
    args.phase_file.write_text('setup')
started = time.perf_counter()
counts = mmread(reference_dir / 'counts.mtx')
genes = (reference_dir / 'genes.tsv').read_text().splitlines()
metadata = pd.read_csv(reference_dir / 'cells.tsv', sep='\t')
assert counts.shape == (556, 10000)
model = scCube()
reference = model.pre_process(pd.DataFrame(counts.toarray(), index=genes, columns=metadata.Cell),
                              metadata, is_normalized=False)
input_seconds = time.perf_counter() - started
reference_seconds = json.loads((reference_dir / 'measurement.json').read_text())['generation_seconds']
records = []
for i, s in enumerate(settings[:1] if args.mode == 'combined' else settings):
    if i and args.memory_budget_gib:
        projected = records[-1]['cumulative_process_peak_rss_gib'] * s['n_cells'] / settings[i-1]['n_cells'] * 1.75
        gate = dict(previous_cells=settings[i-1]['n_cells'], next_cells=s['n_cells'],
                    projected_with_headroom_gib=projected, budget_gib=args.memory_budget_gib,
                    proceed=projected < args.memory_budget_gib * 0.8)
        (out / f"gate_n{s['n_cells']}.json").write_text(json.dumps(gate, indent=2))
        if not gate['proceed']:
            raise SystemExit('Stopping before next size: memory projection exceeds conservative budget')
    if args.mode == 'combined':
        started = time.perf_counter()
        generated_metadata, expression = model.train_vae_and_generate_cell(
            reference, celltype_key='Cell_type', cell_key='Cell',
            target_num=s['target_cells_per_type'].copy(), epoch_num=s['sccube_epochs'],
            used_device='cpu', save_model=False)
        combined_seconds = time.perf_counter() - started
    else:
        started = time.perf_counter()
        inputs = model._scCube__parpare_generate(
            reference, celltype_key='Cell_type', cell_key='Cell',
            target_num=s['target_cells_per_type'].copy())
        preparation_seconds = time.perf_counter() - started
        if i == 0:
            if args.phase_file:
                args.phase_file.write_text('training')
            started = time.perf_counter()
            vae = train_vae(inputs, batch_size=512, epoch_num=s['sccube_epochs'],
                            lr=0.0001, hidden_size=128, used_device='cpu')
            training_seconds = time.perf_counter() - started
            setup = dict(input_seconds=input_seconds, vae_training_seconds=training_seconds,
                         reference_generation_seconds=reference_seconds,
                         cumulative_peak_rss_after_training_gib=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 2**20)
            (out / 'setup.json').write_text(json.dumps(setup, indent=2))
        # For subsequent sizes, reset generation RNGs; the first follows training
        # exactly, allowing an exact comparison against the combined public API.
        if i > 0:
            random.seed(seed); np.random.seed(seed); torch.manual_seed(seed)
        if args.phase_file:
            args.phase_file.write_text(f"simulation_n{s['n_cells']}")
        started = time.perf_counter()
        generated_metadata, expression = generate_vae(vae, inputs, used_device='cpu')
        expression_seconds = time.perf_counter() - started
    started = time.perf_counter()
    expression, generated_metadata = model.generate_pattern_random(
        expression, generated_metadata, spatial_dim=3, spatial_size=8,
        delta=2.0, lamda=0.75, is_split=True, split_coord='point_z', slice_num=s['albis']['n_slices'],
        set_seed=True, seed=seed)
    spatial_seconds = time.perf_counter() - started
    # Take high-water mark before validation/hashing allocations.
    peak = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 2**20
    if args.phase_file:
        args.phase_file.write_text(f"validation_n{s['n_cells']}")
    assert expression.shape == (556, s['n_cells'])
    assert len(generated_metadata) == s['n_cells']
    assert generated_metadata['slice'].nunique() == s['albis']['n_slices']
    # Read-only chunked checks avoid a whole-output boolean temporary at 5M.
    for start in range(0, s['n_cells'], 10000):
        assert np.isfinite(expression.iloc[:, start:start+10000].to_numpy()).all()
    assert np.isfinite(generated_metadata[['point_x', 'point_y', 'point_z']].to_numpy()).all()
    # Read-only digests include axes and dtypes as well as all returned values.
    digests = {}
    for name, frame in ([] if args.skip_output_hashes else
                        [('expression', expression), ('metadata', generated_metadata)]):
        digest = hashlib.sha256(pd.util.hash_pandas_object(frame, index=True).values.tobytes())
        digest.update(repr((frame.columns.tolist(), frame.dtypes.astype(str).tolist())).encode())
        digests[name] = digest.hexdigest()
    record = dict(status='ok', n_cells=s['n_cells'], seed=seed, n_genes=556, n_slices=s['albis']['n_slices'],
                  spatial_seconds=spatial_seconds, cumulative_process_peak_rss_gib=peak,
                  output_hashes=digests)
    if args.mode == 'combined':
        record['combined_expression_seconds'] = combined_seconds
    else:
        record.update(preparation_seconds=preparation_seconds, expression_seconds=expression_seconds,
                      simulation_seconds=expression_seconds + spatial_seconds)
    records.append(record)
    (out / 'measurements.json').write_text(json.dumps(records, indent=2))
    print(json.dumps(record), flush=True)
    del expression, generated_metadata
    if args.mode == 'split':
        del inputs
    gc.collect()
if args.phase_file:
    args.phase_file.write_text('complete')
