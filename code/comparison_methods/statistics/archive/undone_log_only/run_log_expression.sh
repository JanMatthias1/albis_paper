#!/usr/bin/env bash
#SBATCH --job-name=figure5B_logonly
#SBATCH --partition=shared
#SBATCH --cpus-per-task=2
#SBATCH --mem=32G
#SBATCH --time=01:00:00
#SBATCH --output=/dcs04/hicks/data/Jan/sim_project/sim_paper/code/comparison_methods/statistics/logs/logonly_%j.out
set -euo pipefail
project=/dcs04/hicks/data/Jan/sim_project
export MPLCONFIGDIR=/tmp/figure5B_logonly_${SLURM_JOB_ID}
export OPENBLAS_NUM_THREADS=2 OMP_NUM_THREADS=2
"$project/comparison_methods/env/analysis/bin/python" -u "$project/sim_paper/code/comparison_methods/statistics/plot_log_expression.py"
