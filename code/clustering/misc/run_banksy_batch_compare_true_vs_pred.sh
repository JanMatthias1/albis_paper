#!/bin/bash
#SBATCH --job-name=banksy_batch_tvp
#SBATCH --output=/dcs04/hicks/data/Jan/sim_project/sim_paper/code/clustering/logs/banksy_batch_tvp_%A_%a.out
#SBATCH --array=0-5%6
#SBATCH --time=04:00:00
#SBATCH --mem=120G
#SBATCH --cpus-per-task=8
#SBATCH --partition=shared
#
# After-Harmony ground-truth-vs-predicted UMAP + contingency heatmap for the
# completed BANKSY batch-compare job (35608415), one task per (modality,
# batch) combo -- see plot_banksy_batch_compare_true_vs_pred.py.
set -euo pipefail

mkdir -p /dcs04/hicks/data/Jan/sim_project/sim_paper/code/clustering/logs
source /dcs04/hicks/data/Jan/sim_project/sim_paper/code/clustering/_env.sh
cd /dcs04/hicks/data/Jan/sim_project

MODS=(cell cell bin bin spot spot)
BATCHES=(prebatch tuned prebatch tuned prebatch tuned)

i="${SLURM_ARRAY_TASK_ID:-0}"
MOD="${MODS[$i]}"
BATCH="${BATCHES[$i]}"

echo "[task ${i}] modality=${MOD} batch=${BATCH}"
"${PYTHON_BIN}" sim_paper/code/clustering/misc/plot_banksy_batch_compare_true_vs_pred.py \
    --modality "${MOD}" --batch "${BATCH}"
echo "[done] task ${i}"
