#!/usr/bin/env bash
#SBATCH --job-name=sccube_100k
#SBATCH --partition=shared
#SBATCH --cpus-per-task=4
#SBATCH --mem=128G
#SBATCH --time=3-00:00:00
#SBATCH --output=/dcs04/hicks/data/Jan/sim_project/sim_paper/code/comparison_methods/overview/logs/sccube_100k_%j.out
set -euo pipefail
project=/dcs04/hicks/data/Jan/sim_project
code="$project/sim_paper/code/comparison_methods/overview"
output="$project/sim_paper/data/comparison_methods/overview_sccube_100k_${SLURM_JOB_ID}"
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMBA_NUM_THREADS=1 PYTHONHASHSEED=20260921
export MPLCONFIGDIR="$output/cache/matplotlib" NUMBA_CACHE_DIR="$output/cache/numba" XDG_CACHE_HOME="$output/cache"
"$project/comparison_methods/env/analysis/bin/python" "$code/prepare_sccube_100k.py" --out "$output"
mkdir -p "$MPLCONFIGDIR" "$NUMBA_CACHE_DIR"
"$project/comparison_methods/env/sccube/bin/python" -u "$code/generate.py" --method sccube --settings "$output/settings.json" --out "$output" > "$output/sccube.log" 2>&1
"$project/comparison_methods/env/analysis/bin/python" "$code/plot_combined.py" --input "$output" --slice-id 4 > "$output/plot.log" 2>&1
echo "Finished: $output/overview_combined.png"
