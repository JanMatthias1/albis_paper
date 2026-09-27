#!/usr/bin/env bash
#SBATCH --job-name=figure5a_sccube600k
#SBATCH --partition=shared
#SBATCH --cpus-per-task=4
#SBATCH --mem=128G
#SBATCH --time=3-00:00:00
#SBATCH --output=/dcs04/hicks/data/Jan/sim_project/sim_paper/code/comparison_methods/overview/logs/sccube600k_%j.out
# Generate a fresh scCube realization using the current Figure 5A settings.
# Usage: sbatch submit_sccube600k.sh /absolute/path/to/new/output
set -euo pipefail
project=/dcs04/hicks/data/Jan/sim_project
code="$project/sim_paper/code/comparison_methods/overview"
reference="$project/sim_paper/data/comparison_methods/figure_5A_600k"
out="${1:?Supply a new output directory}"
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMBA_NUM_THREADS=1 PYTHONHASHSEED=20260921
"$project/comparison_methods/env/analysis/bin/python" - "$reference" "$out" <<'PY'
import json
import sys
from pathlib import Path

reference, out = map(lambda p: Path(p).resolve(), sys.argv[1:])
if out.exists():
    raise FileExistsError(f'Refusing to overwrite existing output: {out}')
settings = json.loads((reference / 'settings.json').read_text())
assert settings['n_cells'] == 600000 and settings['n_slices'] == 10
pool = (reference / 'splatter').resolve(strict=True)
assert json.loads((pool / 'manifest.json').read_text())['status'] == 'ok'
sources = {method: (reference / method).resolve(strict=True) for method in ['albis', 'spider']}
for source in sources.values():
    for modality in ['cell', 'bin', 'spot']:
        assert (source / f'{modality}.h5ad').is_file()
out.mkdir(parents=True)
(out / 'settings.json').write_text(json.dumps(settings, indent=2) + '\n')
(out / 'splatter').symlink_to(pool, target_is_directory=True)
for method, source in sources.items():
    (out / method).symlink_to(source, target_is_directory=True)
provenance = json.loads((reference / 'input_provenance.json').read_text())
provenance.pop('sccube_source', None)
provenance['sccube_regeneration'] = {
    'settings_source': str(reference / 'settings.json'),
    'splatter_source': str(pool),
    'reused_methods': {method: str(source) for method, source in sources.items()},
    'description': 'Fresh scCube VAE training and 600k-cell generation; reuse ALBIS and SPIDER.'
}
(out / 'input_provenance.json').write_text(json.dumps(provenance, indent=2) + '\n')
PY
export MPLCONFIGDIR="$out/cache/matplotlib" NUMBA_CACHE_DIR="$out/cache/numba" XDG_CACHE_HOME="$out/cache"
mkdir -p "$MPLCONFIGDIR" "$NUMBA_CACHE_DIR"
"$project/comparison_methods/env/sccube/bin/python" -u "$code/generate.py" --method sccube --settings "$out/settings.json" --out "$out" > "$out/sccube.log" 2>&1
bash "$code/render.sh" "$out" > "$out/plot.log" 2>&1
echo "Finished: $out/overview_combined.png"
