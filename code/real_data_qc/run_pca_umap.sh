#!/bin/bash
#SBATCH --job-name=real_data_pca_umap
#SBATCH --output=/dcs04/hicks/data/Jan/sim_project/sim_paper/code/real_data_qc/logs/real_data_pca_umap_%A_%a.out
#SBATCH --array=0-3
#SBATCH --time=08:00:00
#SBATCH --mem=200G
#SBATCH --cpus-per-task=4
#SBATCH --partition=shared

set -euo pipefail

mkdir -p /dcs04/hicks/data/Jan/sim_project/sim_paper/code/real_data_qc/logs

source /dcs04/hicks/data/Jan/sim_project/sim_paper/code/real_data_qc/_env.sh

cd /dcs04/hicks/data/Jan/sim_project

echo "Job started: $(date)"
echo "Host: $(hostname)"
echo "Python: ${PYTHON_BIN}"

DATASETS=(human_pancreas_visium_hd breast_cancer_visium_hd non_diseased_lung lung_cancer)
DATASET="${DATASETS[${SLURM_ARRAY_TASK_ID:-0}]}"

echo "Dataset: ${DATASET}"
echo "Command: ${PYTHON_BIN} sim_paper/code/real_data_qc/pca_umap.py --dataset ${DATASET} $*"
"${PYTHON_BIN}" sim_paper/code/real_data_qc/pca_umap.py --dataset "${DATASET}" "$@"

echo "Job finished: $(date)"
