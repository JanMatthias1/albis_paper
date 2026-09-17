#!/bin/bash
#SBATCH --job-name=domain_check_diag2
#SBATCH --output=/dcs04/hicks/data/Jan/sim_project/sim_paper/code/clustering/logs/domain_check_diag2_%A_%a.out
#SBATCH --array=0-5%6
#SBATCH --time=04:00:00
#SBATCH --mem=150G
#SBATCH --cpus-per-task=8
#SBATCH --partition=shared
#
# Same as domain_check_diagnostics.sh but for the 3 pca_harmony_domain_check
# combos that were still PENDING in the queue as of 2026-09-14 (job 35699479
# tasks 3-5: bin16um/tuned, spot/prebatch, spot/tuned). Submitted with
# --dependency=afterany:35699479 so it fires automatically once that whole
# array finishes (queue was congested behind an unrelated project's job
# array on the shared partition).
set -euo pipefail

mkdir -p /dcs04/hicks/data/Jan/sim_project/sim_paper/code/clustering/logs
source /dcs04/hicks/data/Jan/sim_project/sim_paper/code/clustering/_env.sh
cd /dcs04/hicks/data/Jan/sim_project

ROOT="sim_paper/data/figure_3/pca_harmony_domain_figure2_data"

MODS=(bin spot spot)
SUBDIRS=(bin16um spot spot)
LABELS=(tuned prebatch tuned)

i="${SLURM_ARRAY_TASK_ID:-0}"
j=$(( i % 3 ))
STEP=$(( i / 3 ))   # 0 = umap, 1 = classification

MOD="${MODS[$j]}"
SUBDIR="${SUBDIRS[$j]}"
LABEL="${LABELS[$j]}"

BANKSY_H5AD="${ROOT}/${SUBDIR}/${LABEL}/banksy_matrix/simulation_${MOD}_z_banksy_pca_harmony_qc.h5ad"
ARI_H5AD="${ROOT}/${SUBDIR}/${LABEL}/ari/simulation_${MOD}_z_ari_recovery.h5ad"

if [[ "${STEP}" == "0" ]]; then
    echo "[task ${i}] umap  modality=${MOD} subdir=${SUBDIR} label=${LABEL}"
    "${PYTHON_BIN}" sim_paper/code/clustering/misc/banksy_batch_compare/plot_banksy_batch_compare_umap.py \
        --modality "${MOD}" --input "${BANKSY_H5AD}" --label "${SUBDIR}_${LABEL}"
else
    echo "[task ${i}] classification  modality=${MOD} subdir=${SUBDIR} label=${LABEL}"
    "${PYTHON_BIN}" sim_paper/code/clustering/misc/banksy_batch_compare/plot_banksy_batch_compare_true_vs_pred.py \
        --modality "${MOD}" --input "${ARI_H5AD}" --label "${SUBDIR}_${LABEL}"
fi
echo "[done] task ${i}"
