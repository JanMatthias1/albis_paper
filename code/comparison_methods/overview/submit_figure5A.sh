#!/usr/bin/env bash
# Figure 5A, fully regenerated in one folder (nothing linked from earlier runs):
#   1. splatter  556-gene Splatter expression pool (seed 20260921, 10k cells)
#   2. albis     ALBIS shared tissue (generate_strongmix_offsets.py defaults, seed 2025)
#   3. spider    SPIDER 600k, native 3D/slicing/capture, after 1
#   4. sccube    scCube 600k, native 3D/slicing/capture, after 1
#   5. finish    reproducibility check + render, after 2-4
# Usage: bash submit_figure5A.sh [/absolute/new/output]   (default: figure_5A_600k)
set -euo pipefail
project=/dcs04/hicks/data/Jan/sim_project
code="$project/sim_paper/code/comparison_methods/overview"
envs="$project/comparison_methods/env"
albis_python="$project/sim_paper/env/albis-tutorial/bin/python"
albis_generator="$project/sim_paper/code/applications_albis/cross_modality_alignment/strong_domain_mix/generate_strongmix_offsets.py"
out="${1:-$project/sim_paper/data/figure_5/figure_5A_600k}"
[[ -e "$out" ]] && { echo "Refusing to overwrite existing output: $out" >&2; exit 1; }

mkdir -p "$out/logs" "$out/albis" "$out/cache/matplotlib" "$out/cache/numba"
cp "$code/settings_figure5A.json" "$out/settings.json"
"$envs/analysis/bin/python" - "$out" <<'PY'
import json, sys
from pathlib import Path
out = Path(sys.argv[1]); s = json.loads((out / "settings.json").read_text())
(out / "splatter_contract.json").write_text(json.dumps(dict(
    n_cells=s["splatter_pool_cells"], n_genes=s["n_genes"], n_cell_types=len(s["proportions"]),
    cell_type_proportions=s["proportions"]), indent=2) + "\n")
PY
# ALBIS outputs are addressed as albis/{cell,bin,spot}.h5ad inside this folder.
ln -s ../albis_generation/data/cell/simulation_cell_z.h5ad "$out/albis/cell.h5ad"
ln -s ../albis_generation/data/bin16um/simulation_bin_z.h5ad "$out/albis/bin.h5ad"
ln -s ../albis_generation/data/spot/simulation_spot_z.h5ad "$out/albis/spot.h5ad"

common="set -euo pipefail; export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMBA_NUM_THREADS=1 PYTHONHASHSEED=20260921 PYTHONDONTWRITEBYTECODE=1 MPLCONFIGDIR=$out/cache/matplotlib NUMBA_CACHE_DIR=$out/cache/numba XDG_CACHE_HOME=$out/cache"
job() { sbatch --parsable --partition=shared "$@"; }

splatter=$(job --job-name=fig5a_splatter --cpus-per-task=1 --mem=16G --time=02:00:00 \
  --output="$out/logs/splatter_%j.out" --wrap "$common; \
  $envs/splatter/bin/Rscript $project/comparison_methods/code/workflows/splatter_expression.R \
  --contract $out/splatter_contract.json --seed 20260921 --out-dir $out/splatter > $out/splatter.log 2>&1; \
  test \$(wc -l < $out/splatter/genes.tsv) -eq 556")
albis=$(job --job-name=fig5a_albis --cpus-per-task=4 --mem=150G --time=06:00:00 \
  --output="$out/logs/albis_%j.out" --wrap "$common; \
  $albis_python -u $albis_generator --outdir $out/albis_generation > $out/albis.log 2>&1")
spider=$(job --job-name=fig5a_spider --cpus-per-task=4 --mem=128G --time=1-00:00:00 \
  --dependency=afterok:$splatter --output="$out/logs/spider_%j.out" --wrap "$common; \
  $envs/spider/bin/python -u $code/generate.py --method spider --settings $out/settings.json --out $out > $out/spider.log 2>&1")
# compute-111/095/116 hung scCube VAE training in earlier Figure 5 runs.
sccube=$(job --job-name=fig5a_sccube --cpus-per-task=4 --mem=128G --time=3-00:00:00 \
  --exclude=compute-111,compute-095,compute-116 --dependency=afterok:$splatter \
  --output="$out/logs/sccube_%j.out" --wrap "$common; \
  $envs/sccube/bin/python -u $code/generate.py --method sccube --settings $out/settings.json --out $out > $out/sccube.log 2>&1")
finish=$(job --job-name=fig5a_finish --cpus-per-task=2 --mem=128G --time=06:00:00 \
  --dependency=afterok:$albis:$spider:$sccube --output="$out/logs/finish_%j.out" --wrap "$common; \
  $envs/analysis/bin/python -u $code/verify_reproducibility.py --run $out > $out/reproducibility.log 2>&1; \
  bash $code/render.sh $out > $out/plot.log 2>&1")

printf 'stage\tjob\nsplatter\t%s\nalbis\t%s\nspider\t%s\nsccube\t%s\nfinish\t%s\n' \
  "$splatter" "$albis" "$spider" "$sccube" "$finish" > "$out/submitted_jobs.tsv"
echo "Output: $out"; cat "$out/submitted_jobs.tsv"
