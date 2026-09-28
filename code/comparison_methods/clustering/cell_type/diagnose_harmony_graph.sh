#!/usr/bin/env bash
#SBATCH --job-name=harmony_diagnostic
#SBATCH --partition=shared
#SBATCH --cpus-per-task=8
#SBATCH --mem=96G
#SBATCH --time=04:00:00
#SBATCH --output=/dcs04/hicks/data/Jan/sim_project/sim_paper/code/comparison_methods/clustering/cell_type/logs/harmony_diag_%A_%a.out
set -euo pipefail
project=/dcs04/hicks/data/Jan/sim_project
export OPENBLAS_NUM_THREADS=8 OMP_NUM_THREADS=8 NUMBA_NUM_THREADS=8 MPLCONFIGDIR=/tmp/harmony_diag_${SLURM_JOB_ID}
if [[ "$SLURM_ARRAY_TASK_ID" == 0 ]]; then genes=556; else genes=2000; fi
"$project/sim_paper/env/albis-tutorial/bin/python" -u "$project/sim_paper/code/comparison_methods/clustering/cell_type/diagnose_harmony_graph.py" --genes "$genes"
