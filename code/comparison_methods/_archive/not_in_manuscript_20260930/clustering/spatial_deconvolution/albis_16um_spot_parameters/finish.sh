#!/usr/bin/env bash
#SBATCH --job-name=albis16um_figures
#SBATCH --partition=shared
#SBATCH --cpus-per-task=2
#SBATCH --mem=32G
#SBATCH --time=01:00:00
set -euo pipefail
project=/dcs04/hicks/data/Jan/sim_project
code="$project/sim_paper/code/comparison_methods/clustering/spatial_deconvolution/albis_16um_spot_parameters"
export OPENBLAS_NUM_THREADS=2 OMP_NUM_THREADS=2 MPLCONFIGDIR="/tmp/albis16um_figures_${SLURM_JOB_ID}" PYTHONDONTWRITEBYTECODE=1
"$project/comparison_methods/env/analysis/bin/python" -u "$code/summarize.py"
"$project/comparison_methods/env/analysis/bin/python" -u "$code/plot_deconvolution.py"
