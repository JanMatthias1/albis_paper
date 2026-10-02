#!/usr/bin/env bash
# Figure 5D variant with distinct simulation seeds for scCube and SPIDER.
# In figure_5A_600k both methods use seed 20260921 and both place cells with
# np.random.uniform in the same cube, so their 600k cell positions are identical
# (cell types/expression are independent). Here SPIDER alone is regenerated with
# seed 20260922 (cell level only, all other settings as Figure 5A); scCube and
# ALBIS stay the canonical figure_5A_600k data.
#   1. splatter  556-gene pool, seed 20260921; must be byte-identical to Figure 5A's
#   2. spider    SPIDER 600k cells, seed 20260922, after 1
#   3. plot      actual_capture.py --spider-base <out> --out <out>, after 2
# Usage: bash submit_distinct_seeds.sh [/absolute/new/output]
set -euo pipefail
project=/dcs04/hicks/data/Jan/sim_project
data="$project/sim_paper/data/figure_5"
overview="$project/sim_paper/code/comparison_methods/overview"
code="$project/sim_paper/code/comparison_methods/figure_5D"
envs="$project/comparison_methods/env"
spider_seed=20260922
out="${1:-$data/figure_5D_actual_coordinates/distinct_seeds}"
[[ -e "$out" ]] && { echo "Refusing to overwrite existing output: $out" >&2; exit 1; }

mkdir -p "$out/logs" "$out/cache/matplotlib" "$out/cache/numba"
"$envs/analysis/bin/python" - "$overview/settings_figure5A.json" "$out" "$spider_seed" <<'PY'
import json, sys
from pathlib import Path
source, out, seed = Path(sys.argv[1]), Path(sys.argv[2]), int(sys.argv[3])
s = json.loads(source.read_text())
s.update(seed=seed, cell_only=True)
(out / "settings.json").write_text(json.dumps(s, indent=2) + "\n")
(out / "splatter_contract.json").write_text(json.dumps(dict(
    n_cells=s["splatter_pool_cells"], n_genes=s["n_genes"], n_cell_types=len(s["proportions"]),
    cell_type_proportions=s["proportions"]), indent=2) + "\n")
PY

common="set -euo pipefail; export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMBA_NUM_THREADS=1 PYTHONHASHSEED=20260921 PYTHONDONTWRITEBYTECODE=1 MPLCONFIGDIR=$out/cache/matplotlib NUMBA_CACHE_DIR=$out/cache/numba XDG_CACHE_HOME=$out/cache"
job() { sbatch --parsable --partition=shared "$@"; }

splatter=$(job --job-name=fig5d_seeds_splatter --cpus-per-task=1 --mem=16G --time=02:00:00 \
  --output="$out/logs/splatter_%j.out" --wrap "$common; \
  $envs/splatter/bin/Rscript $project/sim_paper/code/comparison_methods/overview/splatter_expression.R \
  --contract $out/splatter_contract.json --seed 20260921 --out-dir $out/splatter > $out/splatter.log 2>&1; \
  for f in counts.mtx genes.tsv cells.tsv; do cmp $out/splatter/\$f $data/figure_5A_600k/splatter/\$f; done")
spider=$(job --job-name=fig5d_seeds_spider --cpus-per-task=4 --mem=128G --time=1-00:00:00 \
  --dependency=afterok:$splatter --output="$out/logs/spider_%j.out" --wrap "$common; \
  $envs/spider/bin/python -u $overview/generate.py --method spider --settings $out/settings.json --out $out > $out/spider.log 2>&1")
plot=$(job --job-name=fig5d_seeds_plot --cpus-per-task=4 --mem=64G --time=02:00:00 \
  --dependency=afterok:$spider --output="$out/logs/plot_%j.out" --wrap "set -euo pipefail; \
  export MPLCONFIGDIR=$out/cache/matplotlib PYTHONDONTWRITEBYTECODE=1; \
  $project/sim_paper/env/albis-tutorial/bin/python -u $code/actual_capture.py --spider-base $out --out $out > $out/plot.log 2>&1")

printf 'stage\tjob\nsplatter\t%s\nspider\t%s\nplot\t%s\n' "$splatter" "$spider" "$plot" > "$out/submitted_jobs.tsv"
echo "Output: $out"; cat "$out/submitted_jobs.tsv"
