#!/bin/bash
#SBATCH --job-name=domain_check_diag
#SBATCH --output=/dcs04/hicks/data/Jan/sim_project/sim_paper/code/clustering/logs/domain_check_diag_%A_%a.out
#SBATCH --array=0-5%6
#SBATCH --time=04:00:00
#SBATCH --mem=150G
#SBATCH --cpus-per-task=8
#SBATCH --partition=shared
#
# Before/after-Harmony UMAP + ground-truth-vs-predicted classification UMAP
# for the 3 ALREADY-COMPLETE pca_harmony_domain_check combos (job 35699479
# tasks 0-2: cell/prebatch, cell/tuned, bin16um/prebatch) -- the real,
# non-strengthened (manuscript-baseline) domain-mix data. See
# figure3_banksy_domain_sweep memory, 2026-09-14 "user's ask" section for
# background on this dataset. Task index maps 1:1 onto COMBOS below (0,1,2 =
# umap for the 3 combos; 3,4,5 = classification for the same 3 combos).
set -euo pipefail

mkdir -p /dcs04/hicks/data/Jan/sim_project/sim_paper/code/clustering/logs
source /dcs04/hicks/data/Jan/sim_project/sim_paper/code/clustering/_env.sh
cd /dcs04/hicks/data/Jan/sim_project

ROOT="sim_paper/data/figure_3/pca_harmony_domain_figure2_data"

# (banksy_modality, subdir, label)
MODS=(cell cell bin)
SUBDIRS=(cell cell bin16um)
LABELS=(prebatch tuned prebatch)

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
