#!/usr/bin/env bash
#SBATCH --job-name=overview_fig2_600k
#SBATCH --partition=shared
#SBATCH --cpus-per-task=4
#SBATCH --mem=128G
#SBATCH --time=3-00:00:00
#SBATCH --exclude=compute-111,compute-095,compute-116
#SBATCH --output=/dcs04/hicks/data/Jan/sim_project/sim_paper/code/comparison_methods/overview/logs/figure2_%j.out
set -euo pipefail
project=/dcs04/hicks/data/Jan/sim_project
overview="$project/sim_paper/code/comparison_methods/overview"
output="${1:-$project/sim_paper/data/comparison_methods/overview_figure2_600k_${SLURM_JOB_ID}}"
cd "$project"
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMBA_NUM_THREADS=1 PYTHONHASHSEED=20260921
export MPLCONFIGDIR="$output/cache/matplotlib" NUMBA_CACHE_DIR="$output/cache/numba" XDG_CACHE_HOME="$output/cache"
echo "Preparing Figure 2 comparison: $output"
"$project/comparison_methods/env/analysis/bin/python" "$overview/prepare_figure2.py" --out "$output"
mkdir -p "$MPLCONFIGDIR" "$NUMBA_CACHE_DIR"
echo "Generating 600000 scCube cells; progress in $output/sccube.log"
"$project/comparison_methods/env/sccube/bin/python" -u "$overview/generate.py" --method sccube --settings "$output/settings.json" --out "$output" > "$output/sccube.log" 2>&1
echo "Plotting combined overview, stored slice_id=5"
"$project/comparison_methods/env/analysis/bin/python" "$overview/plot_combined.py" --input "$output" --slice-id 5 > "$output/plot.log" 2>&1
echo "Finished: $output/overview_combined.png"
