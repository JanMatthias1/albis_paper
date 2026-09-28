#!/usr/bin/env bash
#SBATCH --job-name=full_knn
#SBATCH --partition=shared
#SBATCH --cpus-per-task=8
#SBATCH --mem=192G
#SBATCH --time=12:00:00
#SBATCH --output=/dcs04/hicks/data/Jan/sim_project/sim_paper/code/comparison_methods/clustering/cell_type/logs/knn_%A_%a.out
set -euo pipefail
project=/dcs04/hicks/data/Jan/sim_project
export OPENBLAS_NUM_THREADS=8 OMP_NUM_THREADS=8 NUMBA_NUM_THREADS=8 MPLCONFIGDIR=/tmp/fullknn_${SLURM_JOB_ID}
if [[ "$1" == select ]]; then args=(--select); elif [[ "$1" == pilot ]]; then args=(--pilot "$SLURM_ARRAY_TASK_ID"); else args=(--task "$SLURM_ARRAY_TASK_ID"); fi
"$project/sim_paper/env/albis-tutorial/bin/python" -u "$project/sim_paper/code/comparison_methods/clustering/cell_type/knn_graph.py" "${args[@]}"
