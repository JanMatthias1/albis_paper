"""Direct SPIDER APIs; preserve its returned labels, XYZ and expression object."""
import hashlib
import importlib.metadata
import json
import os
from pathlib import Path
import random
import resource
import time
import anndata as ad
import numpy as np
import pandas as pd
from scipy.io import mmread
from spider import simulate_10X_3d, get_sim_cell_level_expr

assert importlib.metadata.version('st-spider') == '1.2.0'
phase_file = Path(os.environ['PHASE_FILE']) if os.environ.get('PHASE_FILE') else None
if phase_file:
    phase_file.write_text('setup')
settings_path = Path(os.environ['SETTINGS'])
settings = json.loads(settings_path.read_text())
out = Path(os.environ['OUTPUT'])
out.mkdir(parents=True, exist_ok=False)
(out/'settings.json').write_text(settings_path.read_text())
reference_dir = Path(os.environ['REFERENCE'])
random.seed(settings['seed']); np.random.seed(settings['seed'])
started = time.perf_counter()
counts = mmread(reference_dir/'counts.mtx').T.tocsr()
genes = (reference_dir/'genes.tsv').read_text().splitlines()
metadata = pd.read_csv(reference_dir/'cells.tsv', sep='\t')
assert counts.shape == (settings['reference_cells'],556)
reference = ad.AnnData(counts, obs=metadata.set_index('Cell'), var=pd.DataFrame(index=genes))
input_seconds = time.perf_counter() - started
requested_counts = np.array(list(settings['target_cells_per_type'].values()), dtype=int)

if phase_file:
    phase_file.write_text('generation')
started = time.perf_counter()
labels, xyz = simulate_10X_3d(
    cell_num=settings['n_cells'], Num_celltype=8,
    Num_ct_sample=requested_counts, prior=np.array(settings['proportions']),
    target_trans=np.array(settings['spider_transition']),
    image_width=settings['extent_um'], image_height=settings['extent_um'],
    image_depth=settings['extent_um'], smallsample_max_iter=settings['spider_iterations'])
expression = get_sim_cell_level_expr(
    celltype_assignment=labels, adata=reference, Num_celltype=8,
    Num_ct_sample=requested_counts, match_list=list(settings['target_cells_per_type']), ct_key='Cell_type')
generation_seconds = time.perf_counter() - started
if phase_file:
    phase_file.write_text('validation')

assert xyz.shape == (settings['n_cells'],3)
assert expression.shape == (settings['n_cells'],556)
assert len(labels) == settings['n_cells']
report = dict(status='ok', method='spider', n_cells=settings['n_cells'], seed=settings['seed'],
    n_genes=expression.n_vars, input_seconds=input_seconds, generation_seconds=generation_seconds,
    process_peak_rss_bytes=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss*1024,
    output_policy='native labels, coordinates and expression object unchanged; no simulated-data export',
    endpoint='simulate_10X_3d and get_sim_cell_level_expr returns',
    expression_is_view=bool(expression.is_view),
    script_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest())
(out/'measurement.json').write_text(json.dumps(report, indent=2)+'\n')
print(json.dumps(report), flush=True)
