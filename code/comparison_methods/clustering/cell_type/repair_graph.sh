#!/usr/bin/env bash
#SBATCH --job-name=repair_graph
#SBATCH --partition=shared
#SBATCH --cpus-per-task=4
#SBATCH --mem=96G
#SBATCH --time=04:00:00
#SBATCH --output=/dcs04/hicks/data/Jan/sim_project/sim_paper/code/comparison_methods/clustering/cell_type/logs/repair_%j.out
set -euo pipefail
project=/dcs04/hicks/data/Jan/sim_project
export OPENBLAS_NUM_THREADS=4 OMP_NUM_THREADS=4 NUMBA_NUM_THREADS=4 MPLCONFIGDIR=/tmp/repair_${SLURM_JOB_ID}
"$project/sim_paper/env/albis-tutorial/bin/python" -u "$project/sim_paper/code/comparison_methods/clustering/cell_type/repair_graph.py" --run "$1" "${@:2}"
