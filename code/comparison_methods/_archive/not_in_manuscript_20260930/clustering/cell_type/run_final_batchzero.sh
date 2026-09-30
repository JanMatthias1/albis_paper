#!/usr/bin/env bash
#SBATCH --job-name=fig5c_albis_batchzero
#SBATCH --partition=shared
#SBATCH --cpus-per-task=8
#SBATCH --mem=192G
#SBATCH --time=1-00:00:00
set -euo pipefail
project=/dcs04/hicks/data/Jan/sim_project
code="$project/sim_paper/code/comparison_methods/clustering/cell_type"
export OPENBLAS_NUM_THREADS=8 OMP_NUM_THREADS=8 NUMBA_NUM_THREADS=8 PYTHONDONTWRITEBYTECODE=1
export MPLCONFIGDIR="/tmp/fig5c_albis_batchzero_${SLURM_JOB_ID}"
seeds=(2025 101 202)
"$project/sim_paper/env/albis-tutorial/bin/python" -u "$code/generate_final_batchzero.py" --seed "${seeds[$SLURM_ARRAY_TASK_ID]}"
"$project/sim_paper/env/albis-tutorial/bin/python" -u "$code/cluster_final_batchzero.py" --seed "${seeds[$SLURM_ARRAY_TASK_ID]}"
