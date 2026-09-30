#!/usr/bin/env bash
#SBATCH --job-name=celltype_summary
#SBATCH --partition=shared
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=00:10:00
#SBATCH --output=/dcs04/hicks/data/Jan/sim_project/sim_paper/code/comparison_methods/clustering/cell_type/logs/summary_%j.out
set -euo pipefail
project=/dcs04/hicks/data/Jan/sim_project
export MPLCONFIGDIR=/tmp/celltype_matplotlib_${SLURM_JOB_ID}
"$project/sim_paper/env/albis-tutorial/bin/python" "$project/sim_paper/code/comparison_methods/clustering/cell_type/aggregate.py"
