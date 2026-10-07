#!/bin/bash
#SBATCH --job-name=weakmix_bs0_celltype_panel
#SBATCH --output=/dcs04/hicks/data/Jan/sim_project/albis_paper/data/figure_3/weak_domain_mix/logs/bs0_celltype_panel_%A_%a.out
#SBATCH --array=0-2
#SBATCH --time=06:00:00
#SBATCH --mem=96G
#SBATCH --cpus-per-task=4
#SBATCH --partition=shared
#
# Plain PCA->Harmony->Leiden on the weak domain-mix data at batch_sigma 0 (from
# ../generate/generate_weakmix_bs0.sh): resolution-matched ARI vs cell_type_true
# and domain_true, then UMAP + contingency plots. The batch_sigma-0 counterpart
# of the {cell,bin16um,spot}_celltype_panel.sh runs next to it (those use the
# Figure 2 data at its usual batch_sigma), and the weak-mix counterpart of
# ../../strong_domain_mix/pca_harmony/strongmix_celltype_panel.sh.
#
# PCA/Harmony/ARI are skipped when already present; this script adds the plots.
#
# Outputs: data/figure_3/weak_domain_mix/pca_harmony/<modality>/bs0/
set -euo pipefail

source /dcs04/hicks/data/Jan/sim_project/albis_paper/code/clustering/_env.sh
cd /dcs04/hicks/data/Jan/sim_project

UMAP_FLAGS=(--no-umap-sample)   # full-data UMAP; bin16um overrides below
case "${SLURM_ARRAY_TASK_ID:?run as an array task}" in
    0) NAME=cell; MOD=cell ;;
    1) NAME=bin16um; MOD=bin
       # full UMAP on ~366k bins took >5h and timed out job 35904914_1;
       # subsample it (plots only -- PCA/Harmony/Leiden/ARI use every bin)
       UMAP_FLAGS=() ;;
    2) NAME=spot; MOD=spot ;;
esac
DATA_DIR="albis_paper/data/figure_3/weak_domain_mix/data/${NAME}"
QC_FILES=("${DATA_DIR}"/*_bsigma0_qc.h5ad)
if [[ ${#QC_FILES[@]} -ne 1 || ! -f "${QC_FILES[0]}" ]]; then
    echo "[error] expected exactly one batch_sigma-0 QC h5ad in ${DATA_DIR}, found: ${QC_FILES[*]}" >&2
    exit 1
fi
SIM_QC="${QC_FILES[0]}"
CLUSTER_ROOT="albis_paper/data/figure_3/weak_domain_mix/pca_harmony/${NAME}/bs0"
PCA_H5AD="${CLUSTER_ROOT}/pca_harmony_qc/simulation_${MOD}_z_pca_harmony_qc.h5ad"
ARI_JSON="${CLUSTER_ROOT}/ari_recovery_qc/ari_summary_${MOD}.json"
echo "[input] ${SIM_QC} -> ${CLUSTER_ROOT}"

if [[ ! -f "${PCA_H5AD}" ]]; then
    echo "[pca_harmony] running"
    "${PYTHON_BIN}" albis_paper/code/clustering/step01_pca_harmony.py \
        --modality "${MOD}" --input "${SIM_QC}" "${UMAP_FLAGS[@]}" --output "${PCA_H5AD}"
fi

if [[ ! -f "${ARI_JSON}" ]]; then
    echo "[ari] resolution-matched ARI recovery"
    "${PYTHON_BIN}" albis_paper/code/clustering/step02_leiden_resolution_sweep.py \
        --modality "${MOD}" --input "${PCA_H5AD}" --output-dir "${CLUSTER_ROOT}/ari_recovery_qc"
fi

# Same cell_type_true-matched resolution as the other panels (see
# cell_celltype_panel.sh for why).
RESOLUTION=$("${PYTHON_BIN}" -c "
import json
with open('${ARI_JSON}') as f:
    summary = json.load(f)
print(next(r['resolution'] for r in summary if r['ground_truth'] == 'cell_type_true'))
")
echo "[leiden] cell_type_true-matched-resolution qualitative plots (resolution=${RESOLUTION})"
"${PYTHON_BIN}" albis_paper/code/clustering/step03_cluster_and_plot.py \
    --modality "${MOD}" \
    --input "${PCA_H5AD}" \
    --output-dir "${CLUSTER_ROOT}/leiden_pca_qc_celltype_matched" \
    --pipeline pca_harmony --resolution "${RESOLUTION}"

echo "[done] cell_type_true ARI -> ${ARI_JSON}"
