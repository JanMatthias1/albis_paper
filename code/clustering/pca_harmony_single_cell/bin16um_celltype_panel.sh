#!/bin/bash
#SBATCH --job-name=bin16um_celltype_panel
#SBATCH --output=/dcs04/hicks/data/Jan/sim_project/sim_paper/code/clustering/logs/bin16um_celltype_panel_%j.out
#SBATCH --time=12:00:00
#SBATCH --mem=250G
#SBATCH --cpus-per-task=4
#SBATCH --partition=shared
#
# Figure 3, cell-type recovery panel, BIN modality at the 16um Visium HD
# bin resolution (new 2026-08-26, alongside the existing 8um
# bin_celltype_panel.sh -- both resolutions are canonical, not a
# replacement, matching Figure 2's bin-at-both-resolutions decision).
# Fully self-contained: generates the sim data if not already present,
# QC-filters, runs plain PCA->Harmony->Leiden, reports ARI against
# cell_type_true, and produces the ground-truth-vs-predicted UMAP +
# contingency heatmap plots.
#
# Uses the EXACT SAME dataset as Figure 2's
# figure_2/bin_vs_breast_cancer_visium_hd_16um comparison
# (packing_pf0p04_bin16um_log_mu_-2.5_bsigma07) -- per user decision, this
# panel reads Figure 2's data directly (SIM_RAW/SIM_QC point at
# data/figure_2/<tag>/, not a separate data/noisy/<tag>/ copy). If that file
# isn't there yet, the fallback below generates+QCs it with the identical
# config Figure 2 uses, then moves it into data/figure_2/ itself (same
# pattern as
# code/count_distribution/bin16um_vs_breast_cancer_visium_hd_with_batch_tuned.sh).
# batch_sigma=0.5, same as the 8um bin panel (bin's Figure-3-decided value,
# applied at both resolutions per 2026-08-26 user decision) -- see
# bin_celltype_panel.sh for the batch_sigma=0.5 rationale (small-footprint
# aggregation needs a stronger shift than the manuscript default to produce
# a visible pre/post-Harmony correction story); not yet separately verified
# whether 0.5 is still the right value at 16um's larger aggregation
# footprint (16um pools ~4x the area of 8um), left as-is for now.
#
# Uses the WEAK (manuscript-baseline) domain_type_mix, same as the 8um
# celltype panel -- separate dataset from any future 16um domain panel.
# Pipeline: plain PCA+Harmony, NOT BANKSY -- same rationale as the 8um bin
# panel (BANKSY collapses bin's cell_type_true ARI; see bin_celltype_panel.sh).

set -euo pipefail

mkdir -p /dcs04/hicks/data/Jan/sim_project/sim_paper/code/clustering/logs
source /dcs04/hicks/data/Jan/sim_project/sim_paper/code/clustering/_env.sh
cd /dcs04/hicks/data/Jan/sim_project

SIM_TAG="packing_pf0p04_bin16um_log_mu_-2.5_bsigma07"
MODALITY="bin"
SIM_RAW="sim_paper/data/figure_2/smaller_sphere/data/${SIM_TAG}/simulation_${MODALITY}_z.h5ad"
SIM_QC="sim_paper/data/figure_2/smaller_sphere/data/${SIM_TAG}/simulation_${MODALITY}_z_qc.h5ad"
NOISY_DIR="sim_paper/data/noisy/${SIM_TAG}"
FIG2_DIR="sim_paper/data/figure_2/smaller_sphere/data/${SIM_TAG}"
# "bin16um" (not "bin") keeps this panel's output separate from the 8um
# bin_celltype_panel.sh output -- the underlying --modality passed to the
# python tools below is still "bin" (bin/spot/cell are the only valid
# modality values), only the output directory name is resolution-qualified.
CLUSTER_ROOT="sim_paper/data/figure_3/pca_harmony_single_cell/bin16um"

if [[ ! -f "${SIM_QC}" ]]; then
    if [[ ! -f "${SIM_RAW}" ]]; then
        echo "[generate] ${SIM_TAG} not found under data/figure_2/, generating (same config as Figure 2)"
        "${PYTHON_BIN}" sim_paper/code/data/generate_simulation_noisy.py \
            --modality "${MODALITY}" --sphere-r-um 2050 --bin-size-um 16 \
            --base-gene-lognormal -2.5 0.7 --batch-sigma 0.7 \
            --out-tag "${SIM_TAG}"
        echo "[qc] ${SIM_TAG}"
        "${PYTHON_BIN}" sim_paper/code/clustering/00_qc_filter.py --modality "${MODALITY}" --packing-tag "${SIM_TAG}"
        echo "[move] ${NOISY_DIR} -> ${FIG2_DIR}"
        mkdir -p sim_paper/data/figure_2/smaller_sphere/data
        mv "${NOISY_DIR}" "${FIG2_DIR}"
    else
        echo "[qc] ${SIM_QC} not found but raw exists, running 00_qc_filter.py directly against figure_2/"
        "${PYTHON_BIN}" sim_paper/code/clustering/00_qc_filter.py --modality "${MODALITY}" \
            --input "${SIM_RAW}" --output "${SIM_QC}"
    fi
fi

if [[ ! -f "${CLUSTER_ROOT}/pca_harmony_qc/simulation_${MODALITY}_z_pca_harmony_qc.h5ad" ]]; then
    echo "[pca_harmony] running"
    "${PYTHON_BIN}" sim_paper/code/clustering/01.2_pca_harmony.py \
        --modality "${MODALITY}" --input "${SIM_QC}" \
        --no-umap-sample \
        --output "${CLUSTER_ROOT}/pca_harmony_qc/simulation_${MODALITY}_z_pca_harmony_qc.h5ad"
fi

echo "[ari] resolution-matched ARI recovery"
"${PYTHON_BIN}" sim_paper/code/clustering/02_leiden_resolution_sweep.py \
    --modality "${MODALITY}" \
    --packing-tag "${SIM_TAG}" \
    --input "${CLUSTER_ROOT}/pca_harmony_qc/simulation_${MODALITY}_z_pca_harmony_qc.h5ad" \
    --output-dir "${CLUSTER_ROOT}/ari_recovery_qc"

# Use the SAME resolution 02_leiden_resolution_sweep.py's binary search already found
# for cell_type_true (achieves the true category count exactly), rather than a
# fixed guess -- see bin_celltype_panel.sh for why.
RESOLUTION=$("${PYTHON_BIN}" -c "
import json
with open('${CLUSTER_ROOT}/ari_recovery_qc/ari_summary_${MODALITY}.json') as f:
    summary = json.load(f)
print(next(r['resolution'] for r in summary if r['ground_truth'] == 'cell_type_true'))
")
echo "[leiden] cell_type_true-matched-resolution qualitative plots (resolution=${RESOLUTION}, UMAP true-vs-predicted, contingency heatmap)"
"${PYTHON_BIN}" sim_paper/code/clustering/step03_cluster_and_plot.py \
    --modality "${MODALITY}" \
    --input "${CLUSTER_ROOT}/pca_harmony_qc/simulation_${MODALITY}_z_pca_harmony_qc.h5ad" \
    --output-dir "${CLUSTER_ROOT}/leiden_pca_qc_celltype_matched" \
    --pipeline pca_harmony --resolution "${RESOLUTION}"

echo "[done] cell_type_true ARI -> ${CLUSTER_ROOT}/ari_recovery_qc/ari_summary_${MODALITY}.json"
