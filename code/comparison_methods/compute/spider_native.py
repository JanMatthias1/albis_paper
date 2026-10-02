"""Direct SPIDER APIs: 3D cells and cell-level expression, then sections along Z with
spider.slice_anndata_by_z (equal z-bins over the cube) applied to SPIDER's returned expression
object, as in Figure 5A (overview/generate.py). Expression values are unchanged. If a step
raises (at 5M cells, sectioning the expression exceeds scipy's 32-bit sparse indices), the
run is recorded as failed with its error message instead of a timing."""
import hashlib
import importlib.metadata
import json
import os
from pathlib import Path
import random
import resource
import sys
import time
import traceback
import anndata as ad
import numpy as np
import pandas as pd
from scipy.io import mmread
from spider import simulate_10X_3d, get_sim_cell_level_expr, slice_anndata_by_z

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
n_slices = settings['albis']['n_slices']


def peak_bytes():
    return resource.getrusage(resource.RUSAGE_SELF).ru_maxrss * 1024


def write_report(**fields):
    report = dict(method='spider', n_cells=settings['n_cells'], seed=settings['seed'], n_genes=556,
                  n_slices=n_slices, input_seconds=input_seconds, process_peak_rss_bytes=peak_bytes(),
                  script_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(), **fields)
    (out/'measurement.json').write_text(json.dumps(report, indent=2)+'\n')
    print(json.dumps(report), flush=True)


stage = 'simulate_10X_3d'
started = time.perf_counter()
try:
    labels, xyz = simulate_10X_3d(
        cell_num=settings['n_cells'], Num_celltype=8,
        Num_ct_sample=requested_counts, prior=np.array(settings['proportions']),
        target_trans=np.array(settings['spider_transition']),
        image_width=settings['extent_um'], image_height=settings['extent_um'],
        image_depth=settings['extent_um'], smallsample_max_iter=settings['spider_iterations'])
    stage = 'get_sim_cell_level_expr'
    expression = get_sim_cell_level_expr(
        celltype_assignment=labels, adata=reference, Num_celltype=8,
        Num_ct_sample=requested_counts, match_list=list(settings['target_cells_per_type']), ct_key='Cell_type')
    returned_as_view = bool(expression.is_view)
    # SPIDER's slicer reads coordinates from obsm (z_key, and 'spatial' for x/y), so they are
    # attached to SPIDER's returned object; AnnData then materializes the view it returned.
    stage = 'slice_anndata_by_z'
    expression.obsm['spatial_3d'] = np.asarray(xyz, dtype=float)
    expression.obsm['spatial'] = expression.obsm['spatial_3d'][:, :2].copy()
    sections = slice_anndata_by_z(expression, z_key='spatial_3d',
                                  z_bins=np.linspace(0, settings['extent_um'], n_slices + 1))
    generation_seconds = time.perf_counter() - started
except Exception as error:  # record the failure; nothing is retried or worked around
    write_report(status='failed', failed_stage=stage, seconds_until_failure=time.perf_counter() - started,
                 error=f'{type(error).__name__}: {error}', traceback=traceback.format_exc())
    sys.exit(1)
if phase_file:
    phase_file.write_text('validation')

assert xyz.shape == (settings['n_cells'], 3)
assert expression.shape == (settings['n_cells'], 556)
assert len(labels) == settings['n_cells']
assert len(sections) == n_slices and sum(s.n_obs for s in sections) == settings['n_cells']
assert all(s.n_vars == 556 for s in sections)
write_report(status='ok', generation_seconds=generation_seconds,
             output_policy='native labels, coordinates and expression unchanged; sections from slice_anndata_by_z; '
                           'no simulated-data export',
             endpoint='simulate_10X_3d, get_sim_cell_level_expr and slice_anndata_by_z (with expression) returns',
             expression_returned_as_view=returned_as_view)
