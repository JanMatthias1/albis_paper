#!/bin/bash
#SBATCH --job-name=sim_figure2_final
#SBATCH --output=/dcs04/hicks/data/Jan/sim_project/sim_paper/code/count_distribution/logs/sim_figure2_final_%j.out
#SBATCH --time=04:00:00
#SBATCH --mem=256G
#SBATCH --cpus-per-task=4
#SBATCH --partition=shared

set -euo pipefail

mkdir -p /dcs04/hicks/data/Jan/sim_project/sim_paper/code/count_distribution/logs

source /dcs04/hicks/data/Jan/sim_project/sim_paper/code/count_distribution/_env.sh

cd /dcs04/hicks/data/Jan/sim_project

echo "Job started: $(date)"
echo "Host: $(hostname)"
echo "Python: ${PYTHON_BIN}"

echo "Command: ${PYTHON_BIN} sim_paper/code/count_distribution/generate_figure2_final.py"
"${PYTHON_BIN}" sim_paper/code/count_distribution/generate_figure2_final.py

echo "Job finished: $(date)"
