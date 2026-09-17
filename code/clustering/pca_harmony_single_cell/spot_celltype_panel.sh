#!/bin/bash
#SBATCH --job-name=spot_celltype_panel
#SBATCH --output=/dcs04/hicks/data/Jan/sim_project/sim_paper/code/clustering/logs/spot_celltype_panel_%j.out
#SBATCH --time=08:00:00
#SBATCH --mem=250G
#SBATCH --cpus-per-task=4
#SBATCH --partition=shared
#
# Figure 3, cell-type recovery panel, SPOT modality (Visium-like, 100um
# spacing / 27.5um capture radius).
# Fully self-contained: generates the sim data if not already present,
# QC-filters, runs plain PCA->Harmony->Leiden, reports ARI against
# cell_type_true, and produces the ground-truth-vs-predicted UMAP +
# contingency heatmap plots.
#
# Uses the WEAK (manuscript-baseline) domain_type_mix -- separate dataset
# from the domain panel (spot_domain_panel.sh), see cell_celltype_panel.sh
# for why the two panels don't share one dataset.
#
# 2026-08-26: switched to the EXACT SAME dataset as Figure 2's spot
# count-distribution comparison, instead of this panel's own
# previously-separate packing_pf0p04_bsigma015 tag. Per user decision, this
# panel reads Figure 2's data directly (SIM_RAW/SIM_QC point at
# data/figure_2/<tag>/, not a separate data/noisy/<tag>/ copy). If that file
# isn't there yet, the fallback below generates+QCs it with the identical
# config Figure 2 uses, then moves it into data/figure_2/ itself.
#
# 2026-09-05: spot's Figure 2 config was retuned onto the two 18k CytAssist
# probe references (lymph_node_visium + tonsil_visium; breast_cancer_visium
# dropped 2026-09-04). SIM_TAG here follows Figure 2:
#   packing_pf0p04_log_mu_-2.5_bsigma03 -> packing_pf0p04_log_mu_-2.0_theta_0.25_bsigma03
# (log_mu -2.5 -> -2.0, theta 2.0 -> 0.25; sphere_r_um 2050 + batch_sigma 0.3
# unchanged -- batch_sigma stays fixed, it is set on this panel's Harmony
# demo). Winner of code/data/misc/sweep_spot_*_probe.sh.
#
# 2026-09-06: --theta-jitter added, SIM_TAG bumped to
#   packing_pf0p04_log_mu_-2.0_theta_0.25_jitter0.10_bsigma03
# The shared spot config generated with --theta 0.25 and NO --theta-jitter,
# so generate_simulation_noisy.py fell back to its default THETA_JITTER=1.0 --
# per-gene NB dispersion ~N(0.25, 1.0), ~40% of genes floored to 1e-3, which
# shows up as a disjoint upper cloud in Figure 2's mean_variance /
# mean_dropout panels. Same artifact fixed for `cell` on 2026-08-25
# (--theta-jitter 0.15). --theta-jitter 0.10 (largest of the swept
# 0.10/0.15/0.25, code/data/misc/sweep_spot_jitter_probe.sh job 35536033)
# fully collapses the two clouds; theta unchanged at 0.25. This resubmit
# recomputes Figure 3's spot cell-typing on the fixed shared dataset; the
# pre-jitter-fix outputs are archived at
# data/figure_3/pca_harmony_single_cell/spot_pre_jitter_fix_20260906/
# (the 2026-09-05 pre-probe-retune outputs at spot_pre_probe_retune_20260905/).
# NOTE: the Figure 4/5 application scripts (code/applications_albis/**,
# code/clustering/domain_celltype_composition.py, the _strongmix derivative)
# still point at the old packing_pf0p04_log_mu_-2.5_bsigma03 tag, which stays
# on disk -- migrating those is a separate pass.
#
# batch_sigma=0.3 (per-modality value for the shared Fig 2 / Fig 3 spot
# dataset, finalized 2026-08-27): spot's large capture radius pools far more
# cells per observation than bin/cell, diluting real biological signal
# relative to the same fixed-size batch shift, so spot needs a larger
# per-slice shift than bin/cell to show a pre-Harmony separation that Harmony
# then fully resolves. Held fixed through the 2026-09-05 probe retune (only
# log_mu + theta moved). See figure.md, 2026-08-24/25/27, spot batch_sigma
# diagnosis.
# Pipeline: plain PCA+Harmony (this panel; domain panel uses BANKSY instead,
# see spot_domain_panel.sh -- spot is the one modality where the two panels'
# best pipelines actually differ).

set -euo pipefail

mkdir -p /dcs04/hicks/data/Jan/sim_project/sim_paper/code/clustering/logs
source /dcs04/hicks/data/Jan/sim_project/sim_paper/code/clustering/_env.sh
cd /dcs04/hicks/data/Jan/sim_project

SIM_TAG="packing_pf0p04_log_mu_-2.0_theta_0.25_jitter0.10_bsigma03"
MODALITY="spot"
SIM_RAW="sim_paper/data/figure_2/smaller_sphere/data/${SIM_TAG}/simulation_${MODALITY}_z.h5ad"
SIM_QC="sim_paper/data/figure_2/smaller_sphere/data/${SIM_TAG}/simulation_${MODALITY}_z_qc.h5ad"
NOISY_DIR="sim_paper/data/noisy/${SIM_TAG}"
FIG2_DIR="sim_paper/data/figure_2/smaller_sphere/data/${SIM_TAG}"
# Final-config output lives under figure_3/ (2026-08-25 reorg, same convention
# as the figure_2/ move) -- not the generic clustering_<tag>/ sweep location.
CLUSTER_ROOT="sim_paper/data/figure_3/pca_harmony_single_cell/${MODALITY}"

if [[ ! -f "${SIM_QC}" ]]; then
    if [[ ! -f "${SIM_RAW}" ]]; then
        echo "[generate] ${SIM_TAG} not found under data/figure_2/, generating (same config as Figure 2)"
        "${PYTHON_BIN}" sim_paper/code/data/generate_simulation_noisy.py \
            --modality "${MODALITY}" --sphere-r-um 2050 \
            --base-gene-lognormal -2.0 0.7 --theta 0.25 --theta-jitter 0.10 --batch-sigma 0.3 \
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
# fixed guess -- otherwise this qualitative plot's predicted-cluster count can
# drift from the true count and look like a mismatch that isn't really there
# (found 2026-08-25 on spot: fixed res=0.5 landed on 7 clusters vs. 8 true
# types, while the matched res=0.524 hits 8/8).
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
