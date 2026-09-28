#!/usr/bin/env bash
#SBATCH --job-name=celltype600k
#SBATCH --partition=shared
#SBATCH --cpus-per-task=4
#SBATCH --mem=128G
#SBATCH --time=3-00:00:00
#SBATCH --output=/dcs04/hicks/data/Jan/sim_project/sim_paper/code/comparison_methods/clustering/cell_type/logs/%x_%A_%a.out
set -euo pipefail
project=/dcs04/hicks/data/Jan/sim_project
"$project/sim_paper/env/albis-tutorial/bin/python" -u "$project/sim_paper/code/comparison_methods/clustering/cell_type/run_task.py" "$1" "$SLURM_ARRAY_TASK_ID"
