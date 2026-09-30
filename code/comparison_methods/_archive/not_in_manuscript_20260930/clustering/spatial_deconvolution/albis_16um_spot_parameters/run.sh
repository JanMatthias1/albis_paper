#!/usr/bin/env bash
#SBATCH --job-name=albis16um_rctd
#SBATCH --partition=shared
#SBATCH --cpus-per-task=8
#SBATCH --mem=160G
#SBATCH --time=1-00:00:00
set -euo pipefail
project=/dcs04/hicks/data/Jan/sim_project
code="$project/sim_paper/code/comparison_methods/clustering/spatial_deconvolution/albis_16um_spot_parameters"
base="$project/sim_paper/data/figure_5/spatial_deconvolution/albis_16um_spot_parameters"
seeds=(2025 101 202)
seed=${seeds[$SLURM_ARRAY_TASK_ID]}
run="$base/seed$seed"
export OPENBLAS_NUM_THREADS=8 OMP_NUM_THREADS=8 NUMBA_NUM_THREADS=8
export MPLCONFIGDIR="/tmp/albis16um_${SLURM_JOB_ID}" PYTHONDONTWRITEBYTECODE=1
"$project/sim_paper/env/albis-tutorial/bin/python" -u "$code/generate.py" --seed "$seed"
"$project/comparison_methods/env/analysis/bin/python" -u "$code/prepare_inputs.py" --seed "$seed"
envdir="$project/sim_paper/env/rctd"
export LD_LIBRARY_PATH="$envdir/lib:${LD_LIBRARY_PATH:-}" R_HOME="$envdir/lib/R" RETICULATE_PYTHON="$envdir/bin/python"
export RCTD_CELL_H5AD="$run/inputs/cell.h5ad" RCTD_SPOT_H5AD="$run/inputs/spot.h5ad" RCTD_OUT_DIR="$run/results"
script="$project/sim_paper/code/clustering/spatial_deconvolution/weak_domain/run_rctd_spot.R"
sha256sum "$script" > "$run/rctd_script.sha256"
"$envdir/bin/Rscript" -e 'set.seed(0); source(commandArgs(trailingOnly=TRUE)[1])' "$script" > "$run/rctd.log" 2>&1
