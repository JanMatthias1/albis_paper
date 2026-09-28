#!/usr/bin/env bash
#SBATCH --job-name=comparison_overview
#SBATCH --partition=shared
#SBATCH --cpus-per-task=2
#SBATCH --mem=32G
#SBATCH --time=12:00:00
#SBATCH --exclude=compute-111,compute-095,compute-116
#SBATCH --output=/dcs04/hicks/data/Jan/sim_project/sim_paper/code/comparison_methods/overview/logs/overview_%j.out

set -euo pipefail
project=/dcs04/hicks/data/Jan/sim_project
overview="$project/sim_paper/code/comparison_methods/overview"
output="$project/sim_paper/data/comparison_methods/overview_full_${SLURM_JOB_ID}"
cd "$project"
echo "Job ${SLURM_JOB_ID}; host $(hostname); output $output"
"$project/comparison_methods/env/analysis/bin/python" "$overview/generate.py" \
    --settings "$overview/settings.json" --spider-spots circle --out "$output"
