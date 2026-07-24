#!/bin/bash
#SBATCH --job-name=sim_louvain_umap
#SBATCH --output=/dcs04/hicks/data/Jan/sim_project/sim_paper/code/clustering/logs/sim_louvain_umap_%A_%a.out
#SBATCH --array=0-2
#SBATCH --time=08:00:00
#SBATCH --mem=256G
#SBATCH --cpus-per-task=4
#SBATCH --partition=shared      # uncomment / edit to match your cluster's partition name

set -euo pipefail

mkdir -p /dcs04/hicks/data/Jan/sim_project/sim_paper/code/clustering/logs

source /dcs04/hicks/data/Jan/sim_project/sim_paper/code/clustering/_env.sh

cd /dcs04/hicks/data/Jan/sim_project

echo "Job started: $(date)"
echo "Host: $(hostname)"
echo "Python: ${PYTHON_BIN}"

MODALITIES=(spot bin cell)
MODALITY="${MODALITIES[${SLURM_ARRAY_TASK_ID:-0}]}"

if [[ "$#" -gt 0 ]]; then
    echo "Modality: ${MODALITY}"
    echo "Command: ${PYTHON_BIN} sim_paper/code/clustering/leiden_umap.py --modality ${MODALITY} --algorithm louvain $*"
    "${PYTHON_BIN}" sim_paper/code/clustering/leiden_umap.py --modality "${MODALITY}" --algorithm louvain "$@"
else
    RESOLUTIONS=(0.5 0.1 0.2)
    echo "Modality: ${MODALITY}"
    for resolution in "${RESOLUTIONS[@]}"; do
        echo "Command: ${PYTHON_BIN} sim_paper/code/clustering/leiden_umap.py --modality ${MODALITY} --algorithm louvain --resolution ${resolution}"
        "${PYTHON_BIN}" sim_paper/code/clustering/leiden_umap.py --modality "${MODALITY}" --algorithm louvain --resolution "${resolution}"
    done
fi

echo "Job finished: $(date)"
