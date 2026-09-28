#!/usr/bin/env bash
#SBATCH --job-name=batchzero_noharmony
#SBATCH --partition=shared
#SBATCH --cpus-per-task=8
#SBATCH --mem=192G
#SBATCH --time=12:00:00
#SBATCH --output=/dcs04/hicks/data/Jan/sim_project/sim_paper/code/comparison_methods/clustering/cell_type/logs/batchzero_%A_%a.out
set -euo pipefail
project=/dcs04/hicks/data/Jan/sim_project
export OPENBLAS_NUM_THREADS=8 OMP_NUM_THREADS=8 NUMBA_NUM_THREADS=8 MPLCONFIGDIR=/tmp/batchzero_${SLURM_JOB_ID}
if [[ "${1:-}" == aggregate ]]; then args=(--aggregate); else args=(--task "$SLURM_ARRAY_TASK_ID"); fi
"$project/sim_paper/env/albis-tutorial/bin/python" -u "$project/sim_paper/code/comparison_methods/clustering/cell_type/batchzero_trial.py" "${args[@]}"
