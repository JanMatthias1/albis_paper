"""Direct ALBIS base API call. No method wrappers or output transformations."""
import json
import os
from pathlib import Path
import resource
import time
import hashlib
import albis as ab

REQUIRED_ALBIS_VERSION = '0.1.2'
if ab.__version__ != REQUIRED_ALBIS_VERSION:
    raise RuntimeError(f'albis {ab.__version__} found at {ab.__file__}; this script requires albis '
                       f'{REQUIRED_ALBIS_VERSION} (pip install albis=={REQUIRED_ALBIS_VERSION}).')

phase_file = Path(os.environ['PHASE_FILE']) if os.environ.get('PHASE_FILE') else None
if phase_file:
    phase_file.write_text('setup')
settings_path = Path(os.environ['SETTINGS'])
settings = json.loads(settings_path.read_text())
out = Path(os.environ['OUTPUT'])
out.mkdir(parents=True, exist_ok=False)
(out/'settings.json').write_text(settings_path.read_text())
params = settings['albis']
if phase_file:
    phase_file.write_text('generation')
started = time.perf_counter()
result = ab.simulate_3d_molecule_sphere_base(**params)
generation_seconds = time.perf_counter() - started
if phase_file:
    phase_file.write_text('validation')

# Read-only validation/measurement; do not alter the returned dictionary or arrays.
assert result['adata_cell_true'].shape == (settings['n_cells'], 556)
assert result['molecules']['full_xyz'].shape == (len(result['molecules']['full_gene']), 3)
assert len(result['molecules']['assigned_gene']) == 0
report = dict(status='ok', method='albis', n_cells=settings['n_cells'], seed=settings['seed'],
    n_genes=result['adata_cell_true'].n_vars, n_molecules=len(result['molecules']['full_gene']),
    generation_seconds=generation_seconds,
    process_peak_rss_bytes=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss*1024,
    output_policy='native returned dictionary and arrays unchanged; no simulated-data export',
    endpoint='simulate_3d_molecule_sphere_base return', albis_version=ab.__version__,
    script_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest())
(out/'measurement.json').write_text(json.dumps(report, indent=2)+'\n')
print(json.dumps(report), flush=True)
