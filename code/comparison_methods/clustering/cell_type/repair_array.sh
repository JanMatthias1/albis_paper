#!/usr/bin/env bash
#SBATCH --job-name=unique_profiles
#SBATCH --partition=shared
#SBATCH --cpus-per-task=4
#SBATCH --mem=96G
#SBATCH --time=04:00:00
#SBATCH --output=/dcs04/hicks/data/Jan/sim_project/sim_paper/code/comparison_methods/clustering/cell_type/logs/unique_%A_%a.out
set -euo pipefail
project=/dcs04/hicks/data/Jan/sim_project
code="$project/sim_paper/code/comparison_methods/clustering/cell_type"
run=$("$project/sim_paper/env/albis-tutorial/bin/python" -c 'import json,sys; print(json.load(open(sys.argv[1]))[int(sys.argv[2])]["directory"])' "$project/sim_paper/data/comparison_methods/clustering/cell_type/experiment_600k/runs.json" "$SLURM_ARRAY_TASK_ID")
if [[ -f "$run/clustering_unique_profiles/metrics.json" ]]; then exit 0; fi
bash "$code/repair_graph.sh" "$run" --group-across-slices
