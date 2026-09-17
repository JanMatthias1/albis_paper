#!/bin/bash
#SBATCH --job-name=pre_harmony_ari
#SBATCH --output=/dcs04/hicks/data/Jan/sim_project/sim_paper/code/clustering/logs/pre_harmony_ari_%A_%a.out
#SBATCH --array=0-2%3
#SBATCH --time=04:00:00
#SBATCH --mem=250G
#SBATCH --cpus-per-task=8
#SBATCH --partition=shared
set -euo pipefail

mkdir -p /dcs04/hicks/data/Jan/sim_project/sim_paper/code/clustering/logs
source /dcs04/hicks/data/Jan/sim_project/sim_paper/code/clustering/_env.sh
cd /dcs04/hicks/data/Jan/sim_project

MODS=(cell bin bin16um)
MOD="${MODS[${SLURM_ARRAY_TASK_ID:-0}]}"

echo "[task] modality=${MOD}"
"${PYTHON_BIN}" sim_paper/code/clustering/misc/check_pre_harmony_ari.py --modality "${MOD}"
