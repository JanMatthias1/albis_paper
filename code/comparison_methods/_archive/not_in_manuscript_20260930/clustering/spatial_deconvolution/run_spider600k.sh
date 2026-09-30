#!/usr/bin/env bash
#SBATCH --job-name=rctd_spider600k
#SBATCH --partition=shared
#SBATCH --cpus-per-task=8
#SBATCH --mem=96G
#SBATCH --time=1-00:00:00
#SBATCH --output=/dcs04/hicks/data/Jan/sim_project/sim_paper/code/comparison_methods/clustering/spatial_deconvolution/logs/spider600k_%A_%a.out
set -euo pipefail
project=/dcs04/hicks/data/Jan/sim_project
paper="$project/sim_paper"
code="$paper/code/comparison_methods/clustering/spatial_deconvolution"
seeds=(2025 101 202)
seed=${seeds[$((SLURM_ARRAY_TASK_ID / 2))]}
if (( SLURM_ARRAY_TASK_ID % 2 == 0 )); then genes=556; else genes=2000; fi
source_run="$paper/data/figure_5/clustering/cell_type/experiment_600k/runs/spider_genes${genes}_seed${seed}"
run="$paper/data/figure_5/spatial_deconvolution/spider600k_genes${genes}_seed${seed}"
mkdir -p "$run"
export OPENBLAS_NUM_THREADS=4 OMP_NUM_THREADS=4 MPLCONFIGDIR="$run/cache/mpl"
"$project/comparison_methods/env/spider/bin/python" -u "$code/prepare_spider600k.py" --run "$source_run" --out "$run/inputs" > "$run/preparation.log" 2>&1
rctd_env="$paper/env/rctd"
export LD_LIBRARY_PATH="$rctd_env/lib:${LD_LIBRARY_PATH:-}" R_HOME="$rctd_env/lib/R" RETICULATE_PYTHON="$rctd_env/bin/python"
export RCTD_CELL_H5AD="$run/inputs/cell.h5ad" RCTD_SPOT_H5AD="$run/inputs/spot.h5ad" RCTD_OUT_DIR="$run/results"
script="$paper/code/clustering/spatial_deconvolution/weak_domain/run_rctd_spot.R"
sha256sum "$script" > "$run/rctd_script.sha256"
"$rctd_env/bin/Rscript" -e 'set.seed(0); source(commandArgs(trailingOnly=TRUE)[1])' "$script" > "$run/rctd.log" 2>&1
