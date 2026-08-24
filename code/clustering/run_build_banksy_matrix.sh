#!/bin/bash
#SBATCH --job-name=sim_build_banksy_matrix
#SBATCH --output=/dcs04/hicks/data/Jan/sim_project/sim_paper/code/clustering/logs/sim_build_banksy_matrix_%A_%a.out
#SBATCH --time=08:00:00
#SBATCH --mem=256G
#SBATCH --cpus-per-task=4
#SBATCH --partition=shared

set -euo pipefail

mkdir -p /dcs04/hicks/data/Jan/sim_project/sim_paper/code/clustering/logs

source /dcs04/hicks/data/Jan/sim_project/sim_paper/code/clustering/_env_banksy.sh

cd /dcs04/hicks/data/Jan/sim_project

echo "Job started: $(date)"
echo "Host: $(hostname)"
echo "Python: ${PYTHON_BIN}"

echo "Command: ${PYTHON_BIN} sim_paper/code/clustering/01_build_banksy_matrix.py $*"
"${PYTHON_BIN}" sim_paper/code/clustering/01_build_banksy_matrix.py "$@"

echo "Job finished: $(date)"
