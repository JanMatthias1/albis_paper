#!/usr/bin/env bash
#SBATCH --job-name=albis16um_cluster
#SBATCH --partition=shared
#SBATCH --cpus-per-task=8
#SBATCH --mem=192G
#SBATCH --time=12:00:00
set -euo pipefail
project=/dcs04/hicks/data/Jan/sim_project
code="$project/sim_paper/code/comparison_methods/clustering/spatial_deconvolution/albis_16um_spot_parameters"
export OPENBLAS_NUM_THREADS=8 OMP_NUM_THREADS=8 NUMBA_NUM_THREADS=8 PYTHONDONTWRITEBYTECODE=1
export MPLCONFIGDIR="/tmp/albis16um_cluster_${SLURM_JOB_ID}"
seeds=(2025 101 202)
"$project/sim_paper/env/albis-tutorial/bin/python" -u "$code/cluster_cells.py" --seed "${seeds[$SLURM_ARRAY_TASK_ID]}"
