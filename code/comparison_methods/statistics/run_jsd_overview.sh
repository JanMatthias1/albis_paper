#!/usr/bin/env bash
#SBATCH --job-name=figure5B_jsd
#SBATCH --partition=shared
#SBATCH --cpus-per-task=2
#SBATCH --mem=32G
#SBATCH --time=02:00:00
#SBATCH --output=/dcs04/hicks/data/Jan/sim_project/sim_paper/code/comparison_methods/statistics/logs/jsd_%A_%a.out
set -euo pipefail
project=/dcs04/hicks/data/Jan/sim_project
export MPLCONFIGDIR=/tmp/figure5B_jsd_${SLURM_JOB_ID}
export OPENBLAS_NUM_THREADS=2 OMP_NUM_THREADS=2
if [[ "${1:-}" == aggregate ]]; then args=(--aggregate); else args=(--slice-id "$SLURM_ARRAY_TASK_ID"); fi
"$project/comparison_methods/env/analysis/bin/python" -u "$project/sim_paper/code/comparison_methods/statistics/jsd_overview.py" "${args[@]}"
