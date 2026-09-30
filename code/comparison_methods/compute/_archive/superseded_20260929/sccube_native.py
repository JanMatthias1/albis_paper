"""Direct scCube public APIs. Keep expression and coordinates exactly as returned."""
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
from scCube.sccube import scCube

settings_path = Path(os.environ['SETTINGS'])
settings = json.loads(settings_path.read_text())
out = Path(os.environ['OUTPUT'])
out.mkdir(parents=True, exist_ok=False)
(out/'settings.json').write_text(settings_path.read_text())
reference_dir = Path(os.environ['REFERENCE'])
random.seed(settings['seed']); np.random.seed(settings['seed']); torch.manual_seed(settings['seed'])
torch.set_num_threads(1); torch.set_num_interop_threads(1)
started = time.perf_counter()
counts = mmread(reference_dir/'counts.mtx')
genes = (reference_dir/'genes.tsv').read_text().splitlines()
metadata = pd.read_csv(reference_dir/'cells.tsv', sep='\t')
assert counts.shape == (556, settings['reference_cells'])
model = scCube()
reference = model.pre_process(
    pd.DataFrame(counts.toarray(), index=genes, columns=metadata.Cell),
    metadata, is_normalized=False)
input_seconds = time.perf_counter() - started

started = time.perf_counter()
generated_metadata, expression = model.train_vae_and_generate_cell(
    reference, celltype_key='Cell_type', cell_key='Cell',
    target_num=settings['target_cells_per_type'].copy(), epoch_num=settings['sccube_epochs'],
    used_device='cpu', save_model=False)
expression, generated_metadata = model.generate_pattern_random(
    expression, generated_metadata, spatial_dim=3, spatial_size=8,
    delta=2.0, lamda=0.75, is_split=False, set_seed=True, seed=settings['seed'])
generation_seconds = time.perf_counter() - started

assert expression.shape == (556, settings['n_cells'])
assert len(generated_metadata) == settings['n_cells']
assert {'point_x','point_y','point_z'}.issubset(generated_metadata.columns)
report = dict(status='ok', method='sccube', n_cells=settings['n_cells'], seed=settings['seed'],
    n_genes=expression.shape[0], input_seconds=input_seconds, generation_seconds=generation_seconds,
    process_peak_rss_bytes=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss*1024,
    output_policy='native DataFrames unchanged; native coordinate units; no simulated-data export',
    endpoint='train_vae_and_generate_cell plus generate_pattern_random returns; training included',
    script_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest())
(out/'measurement.json').write_text(json.dumps(report, indent=2)+'\n')
print(json.dumps(report), flush=True)
