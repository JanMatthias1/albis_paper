#!/bin/bash
#SBATCH --job-name=bin16um_celltype_panel
#SBATCH --output=/dcs04/hicks/data/Jan/sim_project/albis_paper/code/clustering/logs/bin16um_celltype_panel_%j.out
#SBATCH --time=12:00:00
#SBATCH --mem=250G
#SBATCH --cpus-per-task=4
#SBATCH --partition=shared
#
# Figure 3 cell-type recovery, BIN 16 µm (Visium HD 16 µm bins).
# Uses the Figure 2 bin16um dataset (weak domain mix, batch_sigma 0.7);
# generates and QCs it with the Figure 2 settings if it is missing. Pipeline:
# plain PCA -> Harmony -> Leiden at the true cell-type count, ARI against
# cell_type_true, then UMAP (true vs predicted) and contingency plots.
# BANKSY is not used for bin cell types (see bin_celltype_panel.sh).
set -euo pipefail

mkdir -p /dcs04/hicks/data/Jan/sim_project/albis_paper/code/clustering/logs
source /dcs04/hicks/data/Jan/sim_project/albis_paper/code/clustering/_env.sh
cd /dcs04/hicks/data/Jan/sim_project

SIM_TAG="packing_pf0p04_bin16um_log_mu_-2.5_jitter0.6_bsigma07"
MODALITY="bin"
SIM_RAW="albis_paper/data/figure_2/smaller_sphere/data/${SIM_TAG}/simulation_${MODALITY}_z.h5ad"
SIM_QC="albis_paper/data/figure_2/smaller_sphere/data/${SIM_TAG}/simulation_${MODALITY}_z_qc.h5ad"
NOISY_DIR="albis_paper/data/noisy/${SIM_TAG}"
FIG2_DIR="albis_paper/data/figure_2/smaller_sphere/data/${SIM_TAG}"
# Output folder "bin16um" keeps this apart from the 8 µm panel; the tools still
# take --modality bin.
CLUSTER_ROOT="albis_paper/data/figure_3/weak_domain_mix/pca_harmony/bin16um/bs0.7"

if [[ ! -f "${SIM_QC}" ]]; then
    if [[ ! -f "${SIM_RAW}" ]]; then
        echo "[generate] ${SIM_TAG} not found under data/figure_2/, generating (same config as Figure 2)"
        "${PYTHON_BIN}" albis_paper/code/data/generate_simulation_noisy.py \
            --modality "${MODALITY}" --sphere-r-um 2050 --bin-size-um 16 \
            --base-gene-lognormal -2.5 0.7 --theta 2.0 --theta-jitter 0.6 --batch-sigma 0.7 \
            --sync-unaligned-seed \
            --out-tag "${SIM_TAG}"
        echo "[qc] ${SIM_TAG}"
        "${PYTHON_BIN}" albis_paper/code/clustering/step00_qc_filter.py --modality "${MODALITY}" --packing-tag "${SIM_TAG}"
        echo "[move] ${NOISY_DIR} -> ${FIG2_DIR}"
        mkdir -p albis_paper/data/figure_2/smaller_sphere/data
        mv "${NOISY_DIR}" "${FIG2_DIR}"
    else
        echo "[qc] ${SIM_QC} not found but raw exists, running step00_qc_filter.py directly against figure_2/"
        "${PYTHON_BIN}" albis_paper/code/clustering/step00_qc_filter.py --modality "${MODALITY}" \
            --input "${SIM_RAW}" --output "${SIM_QC}"
    fi
fi

if [[ ! -f "${CLUSTER_ROOT}/pca_harmony_qc/simulation_${MODALITY}_z_pca_harmony_qc.h5ad" ]]; then
    echo "[pca_harmony] running"
    "${PYTHON_BIN}" albis_paper/code/clustering/step01_pca_harmony.py \
        --modality "${MODALITY}" --input "${SIM_QC}" \
        --no-umap-sample \
        --output "${CLUSTER_ROOT}/pca_harmony_qc/simulation_${MODALITY}_z_pca_harmony_qc.h5ad"
fi

echo "[ari] resolution-matched ARI recovery"
"${PYTHON_BIN}" albis_paper/code/clustering/step02_leiden_resolution_sweep.py \
    --modality "${MODALITY}" \
    --packing-tag "${SIM_TAG}" \
    --input "${CLUSTER_ROOT}/pca_harmony_qc/simulation_${MODALITY}_z_pca_harmony_qc.h5ad" \
    --output-dir "${CLUSTER_ROOT}/ari_recovery_qc"

# Plot at the resolution step02 found for cell_type_true (true category count).
RESOLUTION=$("${PYTHON_BIN}" -c "
import json
with open('${CLUSTER_ROOT}/ari_recovery_qc/ari_summary_${MODALITY}.json') as f:
    summary = json.load(f)
print(next(r['resolution'] for r in summary if r['ground_truth'] == 'cell_type_true'))
")
echo "[leiden] cell_type_true-matched-resolution qualitative plots (resolution=${RESOLUTION}, UMAP true-vs-predicted, contingency heatmap)"
"${PYTHON_BIN}" albis_paper/code/clustering/step03_cluster_and_plot.py \
    --modality "${MODALITY}" \
    --input "${CLUSTER_ROOT}/pca_harmony_qc/simulation_${MODALITY}_z_pca_harmony_qc.h5ad" \
    --output-dir "${CLUSTER_ROOT}/leiden_pca_qc_celltype_matched" \
    --pipeline pca_harmony --resolution "${RESOLUTION}"

echo "[done] cell_type_true ARI -> ${CLUSTER_ROOT}/ari_recovery_qc/ari_summary_${MODALITY}.json"
