#!/usr/bin/env bash
# Regenerate Figure 5A scCube + SPIDER (600k cells) through the native-only
# generate.py (native 3D layout, native slicing, native capture/composition),
# reuse the current Splatter pool and ALBIS links, then render.
# Usage: bash submit_native600k.sh [/absolute/path/to/new/output]
# Submits three jobs: spider, sccube, render (runs after both succeed).
set -euo pipefail
project=/dcs04/hicks/data/Jan/sim_project
code="$project/sim_paper/code/comparison_methods/overview"
envs="$project/comparison_methods/env"
reference="$project/sim_paper/data/comparison_methods/figure_5A_600k"
out="${1:-$project/sim_paper/data/comparison_methods/figure_5A_600k_native_$(date +%Y%m%d)}"

"$envs/analysis/bin/python" - "$reference" "$out" <<'PY'
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
albis = (reference / 'albis').resolve(strict=True)
for modality in ['cell', 'bin', 'spot']:
    assert (albis / f'{modality}.h5ad').resolve(strict=True).is_file()
out.mkdir(parents=True)
(out / 'settings.json').write_text(json.dumps(settings, indent=2) + '\n')
(out / 'splatter').symlink_to(pool, target_is_directory=True)
(out / 'albis').symlink_to(albis, target_is_directory=True)
provenance = json.loads((reference / 'input_provenance.json').read_text())
provenance.pop('sccube_source', None)
provenance['native_regeneration'] = {
    'settings_source': str(reference / 'settings.json'),
    'splatter_source': str(pool),
    'reused_methods': {'albis': str(albis)},
    'description': ('Fresh scCube and SPIDER 600k runs with native slicing '
                    '(scCube is_split, spider.slice_anndata_by_z), native per-capture '
                    'composition (scCube calculate_spot_prop, SPIDER W) and '
                    'spider.make_transition_matrix. SPIDER circle-spot W uses the single '
                    'user-approved exception in comparison_methods/AGENTS.md.'),
}
(out / 'input_provenance.json').write_text(json.dumps(provenance, indent=2) + '\n')
PY

mkdir -p "$out/logs" "$out/cache/matplotlib" "$out/cache/numba"
common="set -euo pipefail; export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMBA_NUM_THREADS=1 PYTHONHASHSEED=20260921 MPLCONFIGDIR=$out/cache/matplotlib NUMBA_CACHE_DIR=$out/cache/numba XDG_CACHE_HOME=$out/cache"
slurm=(--parsable --partition=shared --cpus-per-task=4 --mem=128G)

spider=$(sbatch "${slurm[@]}" --job-name=fig5a_native_spider --time=1-00:00:00 \
  --output="$out/logs/spider_%j.out" \
  --wrap "$common; $envs/spider/bin/python -u $code/generate.py --method spider --settings $out/settings.json --out $out > $out/spider.log 2>&1")
# compute-111/095/116 hung scCube VAE training in earlier Figure 5 runs.
sccube=$(sbatch "${slurm[@]}" --job-name=fig5a_native_sccube --time=3-00:00:00 \
  --exclude=compute-111,compute-095,compute-116 --output="$out/logs/sccube_%j.out" \
  --wrap "$common; $envs/sccube/bin/python -u $code/generate.py --method sccube --settings $out/settings.json --out $out > $out/sccube.log 2>&1")
render=$(sbatch --parsable --partition=shared --cpus-per-task=2 --mem=64G --time=04:00:00 \
  --job-name=fig5a_native_render --dependency=afterok:$spider:$sccube \
  --output="$out/logs/render_%j.out" --wrap "bash $code/render.sh $out > $out/plot.log 2>&1")

echo "Output: $out"
echo "Jobs: spider=$spider sccube=$sccube render=$render (after both)"
