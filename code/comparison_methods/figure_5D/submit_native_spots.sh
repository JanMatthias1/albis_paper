#!/usr/bin/env bash
# Figure 5D: each method's own spots on slice 5 (stored slice_id 4).
# ALBIS and scCube: figure_5A_600k. SPIDER: its own run with seed 20260922, because
# with Figure 5A's shared seed scCube and SPIDER draw identical uniform cell positions.
#   1. splatter  556-gene pool, seed 20260921; must be byte-identical to Figure 5A's
#   2. spider    SPIDER 600k cells + bins + spots, Figure 5A settings, seed 20260922, after 1
#   3. sccube    scCube's native slice memberships (sccube_slice_membership.py)
#   4. plot      native_spots.py, after 2 and 3
# Usage: bash submit_native_spots.sh [/absolute/new/output]
set -euo pipefail
project=/dcs04/hicks/data/Jan/sim_project
data="$project/sim_paper/data/figure_5"
overview="$project/sim_paper/code/comparison_methods/overview"
code="$project/sim_paper/code/comparison_methods/figure_5D"
envs="$project/comparison_methods/env"
spider_seed=20260922
out="${1:-$data/figure_5D_actual_coordinates/native_spots}"
[[ -e "$out" ]] && { echo "Refusing to overwrite existing output: $out" >&2; exit 1; }

mkdir -p "$out/logs" "$out/cache/matplotlib" "$out/cache/numba"
"$envs/analysis/bin/python" - "$overview/settings_figure5A.json" "$out" "$spider_seed" <<'PY'
import json, sys
from pathlib import Path
source, out, seed = Path(sys.argv[1]), Path(sys.argv[2]), int(sys.argv[3])
s = json.loads(source.read_text())
s.update(seed=seed)
(out / "settings.json").write_text(json.dumps(s, indent=2) + "\n")
(out / "splatter_contract.json").write_text(json.dumps(dict(
    n_cells=s["splatter_pool_cells"], n_genes=s["n_genes"], n_cell_types=len(s["proportions"]),
    cell_type_proportions=s["proportions"]), indent=2) + "\n")
PY

common="set -euo pipefail; export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMBA_NUM_THREADS=1 PYTHONHASHSEED=20260921 PYTHONDONTWRITEBYTECODE=1 MPLCONFIGDIR=$out/cache/matplotlib NUMBA_CACHE_DIR=$out/cache/numba XDG_CACHE_HOME=$out/cache"
job() { sbatch --parsable --partition=shared "$@"; }

splatter=$(job --job-name=fig5d_splatter --cpus-per-task=1 --mem=16G --time=02:00:00 \
  --output="$out/logs/splatter_%j.out" --wrap "$common; \
  $envs/splatter/bin/Rscript $project/comparison_methods/code/workflows/splatter_expression.R \
  --contract $out/splatter_contract.json --seed 20260921 --out-dir $out/splatter > $out/splatter.log 2>&1; \
  for f in counts.mtx genes.tsv cells.tsv; do cmp $out/splatter/\$f $data/figure_5A_600k/splatter/\$f; done")
spider=$(job --job-name=fig5d_spider --cpus-per-task=4 --mem=128G --time=1-00:00:00 \
  --dependency=afterok:$splatter --output="$out/logs/spider_%j.out" --wrap "$common; \
  $envs/spider/bin/python -u $overview/generate.py --method spider --settings $out/settings.json --out $out > $out/spider.log 2>&1")
sccube=$(job --job-name=fig5d_sccube --cpus-per-task=2 --mem=64G --time=04:00:00 \
  --exclude=compute-111,compute-095,compute-116 --output="$out/logs/sccube_%j.out" --wrap "$common; \
  $envs/sccube/bin/python -u $code/sccube_slice_membership.py --out $out > $out/sccube.log 2>&1")
plot=$(job --job-name=fig5d_plot --cpus-per-task=4 --mem=64G --time=02:00:00 \
  --dependency=afterok:$spider:$sccube --output="$out/logs/plot_%j.out" --wrap "$common; \
  $project/sim_paper/env/albis-tutorial/bin/python -u $code/native_spots.py --spider-base $out --out $out > $out/plot.log 2>&1")

printf 'stage\tjob\nsplatter\t%s\nspider\t%s\nsccube\t%s\nplot\t%s\n' "$splatter" "$spider" "$sccube" "$plot" > "$out/submitted_jobs.tsv"
echo "Output: $out"; cat "$out/submitted_jobs.tsv"
