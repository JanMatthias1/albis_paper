#!/usr/bin/env bash
#SBATCH --job-name=figure5a_spider600k
#SBATCH --partition=shared
#SBATCH --cpus-per-task=4
#SBATCH --mem=128G
#SBATCH --time=3-00:00:00
#SBATCH --output=/dcs04/hicks/data/Jan/sim_project/sim_paper/code/comparison_methods/overview/logs/spider600k_%j.out
set -euo pipefail
project=/dcs04/hicks/data/Jan/sim_project
code="$project/sim_paper/code/comparison_methods/overview"
out="${1:?Supply prepared overview directory}"
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMBA_NUM_THREADS=1 PYTHONHASHSEED=20260921
export MPLCONFIGDIR="$out/cache/matplotlib" NUMBA_CACHE_DIR="$out/cache/numba" XDG_CACHE_HOME="$out/cache"
mkdir -p "$MPLCONFIGDIR" "$NUMBA_CACHE_DIR"
"$project/comparison_methods/env/spider/bin/python" -u "$code/generate.py" --method spider --settings "$out/settings.json" --out "$out" > "$out/spider.log" 2>&1
"$project/comparison_methods/env/analysis/bin/python" "$code/plot_combined.py" --input "$out" > "$out/plot.log" 2>&1
echo "Finished: $out/overview_combined.png"
