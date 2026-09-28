#!/bin/bash
#SBATCH --job-name=fig2_cell_vs_lung_cancer
#SBATCH --output=/dcs04/hicks/data/Jan/sim_project/sim_paper/code/count_distribution/smaller_sphere/logs/fig2_cell_vs_lung_cancer_%j.out
#SBATCH --time=06:00:00
#SBATCH --mem=250G
#SBATCH --cpus-per-task=4
#SBATCH --partition=shared
#
# Figure 2 final comparison: cell vs. Xenium lung_cancer.
# Fully self-contained: generates the sim data (if not already present),
# QC-filters it, then runs all 4 count_distribution.py modes.
#
# Sim dataset: simulated single CELLS (no spatial aggregation), compared
# against Xenium lung cancer tissue (392-gene targeted panel).
# Current sim config (2026-09-17): log_mu=-2.5, theta=0.40, theta_jitter=0.15,
# batch_sigma=1.5 -- see cell_vs_non_diseased_lung.sh for the full rationale
# (same sim config, different real reference).

# 2026-08-27: SIM_TAG now carries the per-modality batch_sigma finalized for
# the SHARED Figure 2 / Figure 3 dataset (cell 1.5, bin8 0.8, bin16 0.7,
# spot 0.3), tuned on the Figure 3 pre/post-Harmony demo then confirmed here
# to still match the real count distribution. count_distribution.py now
# defaults to --slice-id 5 and post-batch counts, so those flags are no
# longer passed per-call below. The matching clustering panel is
# code/clustering/<modality>_celltype_panel.sh (same SIM_TAG).
#
# 2026-08-31: realwindow is now the default. generate_simulation_noisy.py no
# longer scales the capture window with --sphere-r-um; with --capture-window-um
# unset, cell crops to the real Xenium window (12 x 24 mm) via
# xenium_capture_window_um. cell has no off-tissue "empty" observations (cells
# are only placed within the tissue sphere), so this is effectively a no-op
# for cell beyond the window crop -- all 4 panels stay valid here. Sim data is
# (re)generated under data/noisy/<tag>/ then moved into data/figure_2/<tag>/;
# pre-realwindow data is archived at data/figure_2_oldwindow_20260831/.
#
# 2026-09-15: --sphere-r-um/--n-cells now explicit (2050 / 24207) instead of
# the script default (6000 / 600000), matching cross_modality_alignment/
# cross_tech_stair.py's cell dataset -- that script needs cell's disc the
# SAME physical size as bin16um/spot (all three are "the same tissue slice"
# captured at different technologies), n_cells scaled down to preserve cell's
# original ~0.16% 3D packing fraction (600000 * (2050/6000)**3) rather than
# forced up to bin/spot's ~4%. This tag is shared by both consumers -- until
# 2026-09-16 the flags lived only in a one-off manual rerun, never in this
# script, so a from-scratch regen here would have silently reverted to the
# old 6000um/600k config (which is exactly what happened -- see figure.md).
#
# 2026-09-17: log_mu -2.3 -> -2.5, dispersion re-tuned via an actual joint
# log_mu x theta grid search against both Xenium refs (code/misc/data/misc/
# sweep_cell_logmu_theta_joint_xenium.sh, job 35759008), mirroring spot's
# probe-sweep methodology -- cell had never gotten this treatment before,
# just one manual 2026-08-27 adjustment (see FIGURE2_METHODOLOGY.md's "Known
# open gaps": total_counts overshoot "accepted... revisit if priority").
# Composite score (sum of |ln(sim/real ratio)| over theta_hat/
# total_counts_median/matrix_zero_frac, both refs): -2.5/0.40 = 2.733 vs the
# old -2.3/0.40 = 2.996 -- meaningfully better, mainly from total_counts
# fidelity (ratio 1.04-1.52 vs 1.38-2.03). -2.5/0.50 scored marginally better
# still (2.702) but -2.5/0.40 wins on total_counts specifically, the
# previously-flagged soft spot -- picked over the outright-lowest point for
# that reason. Old data archived to figure_2/misc/archive/
# log_mu_-2.3_theta_0.40_jitter0.15_bsigma15_pre_dispersion_retune_20260917/.
# This retune cascades to every other cell-based figure (Fig3 celltype panel,
# Fig3/4B strong-mix, Fig4B STAGATE, Fig4C/E alignment) -- staged separately,
# deliberately sequenced AFTER this Figure 2 rerun lands and is checked.
set -euo pipefail

mkdir -p /dcs04/hicks/data/Jan/sim_project/sim_paper/code/count_distribution/smaller_sphere/logs
source /dcs04/hicks/data/Jan/sim_project/sim_paper/code/count_distribution/_env.sh
cd /dcs04/hicks/data/Jan/sim_project

SIM_TAG="log_mu_-2.5_theta_0.40_jitter0.15_bsigma15"
MODALITY="cell"
REAL_LABEL="lung_cancer"
SIM_RAW="sim_paper/data/figure_2/smaller_sphere/data/${SIM_TAG}/simulation_${MODALITY}_z.h5ad"
SIM_QC="sim_paper/data/figure_2/smaller_sphere/data/${SIM_TAG}/simulation_${MODALITY}_z_qc.h5ad"
NOISY_DIR="sim_paper/data/noisy/${SIM_TAG}"
REAL_INPUT="sim_paper/data/real_data_qc/${REAL_LABEL}/${REAL_LABEL}_qc.h5ad"
OUT_ROOT="sim_paper/data/figure_2/smaller_sphere/plots/${MODALITY}_vs_${REAL_LABEL}"

# generate_simulation_noisy.py only writes under data/noisy/<out-tag>/, so
# (re)generate + QC there and then move the directory into data/figure_2/.
if [[ ! -f "${SIM_RAW}" ]]; then
    echo "[generate] ${SIM_TAG} not present under figure_2/, building"
    if [[ ! -f "${NOISY_DIR}/simulation_${MODALITY}_z.h5ad" ]]; then
        "${PYTHON_BIN}" sim_paper/code/data/generate_simulation_noisy.py \
            --modality "${MODALITY}" \
            --sphere-r-um 2050 \
            --n-cells 24207 \
            --base-gene-lognormal -2.5 0.7 \
            --theta 0.40 \
            --theta-jitter 0.15 \
            --batch-sigma 1.5 \
            --sync-unaligned-seed \
            --out-tag "${SIM_TAG}"
    fi
    if [[ ! -f "${NOISY_DIR}/simulation_${MODALITY}_z_qc.h5ad" ]]; then
        "${PYTHON_BIN}" sim_paper/code/clustering/00_qc_filter.py \
            --modality "${MODALITY}" --packing-tag "${SIM_TAG}"
    fi
    mkdir -p sim_paper/data/figure_2/smaller_sphere/data
    mv "${NOISY_DIR}" "sim_paper/data/figure_2/smaller_sphere/data/${SIM_TAG}"
fi

if [[ ! -f "${SIM_QC}" ]]; then
    echo "[qc] ${SIM_QC} not found, running 00_qc_filter.py"
    "${PYTHON_BIN}" sim_paper/code/clustering/00_qc_filter.py \
        --modality "${MODALITY}" --input "${SIM_RAW}" --output "${SIM_QC}"
fi

"${PYTHON_BIN}" sim_paper/code/count_distribution/count_distribution.py \
    --modality "${MODALITY}" --input "${SIM_RAW}" \
    --compare-input "${REAL_INPUT}" --compare-label "${REAL_LABEL}" \
    --output-dir "${OUT_ROOT}/full_panel"

"${PYTHON_BIN}" sim_paper/code/count_distribution/count_distribution.py \
    --modality "${MODALITY}" --input "${SIM_RAW}" \
    --compare-input "${REAL_INPUT}" --compare-label "${REAL_LABEL}" \
    --match-panel-size \
    --output-dir "${OUT_ROOT}/hvg_matched"

"${PYTHON_BIN}" sim_paper/code/count_distribution/count_distribution.py \
    --modality "${MODALITY}" --input "${SIM_QC}" \
    --compare-input "${REAL_INPUT}" --compare-label "${REAL_LABEL}" \
    --output-dir "${OUT_ROOT}/qc_filtered"

"${PYTHON_BIN}" sim_paper/code/count_distribution/count_distribution.py \
    --modality "${MODALITY}" --input "${SIM_QC}" \
    --compare-input "${REAL_INPUT}" --compare-label "${REAL_LABEL}" \
    --match-panel-size \
    --output-dir "${OUT_ROOT}/qc_and_hvg_matched"

echo "[done] All 4 modes written under ${OUT_ROOT}"
