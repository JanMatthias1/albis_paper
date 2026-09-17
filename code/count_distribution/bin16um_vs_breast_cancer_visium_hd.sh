#!/bin/bash
#SBATCH --job-name=fig2_bin16um_vs_breast_cancer_visium_hd
#SBATCH --output=/dcs04/hicks/data/Jan/sim_project/sim_paper/code/count_distribution/logs/fig2_bin16um_vs_breast_cancer_visium_hd_%j.out
#SBATCH --time=08:00:00
#SBATCH --mem=250G
#SBATCH --cpus-per-task=4
#SBATCH --partition=shared
#
# Figure 2 final comparison: bin vs. Visium HD breast_cancer_visium_hd, at the
# 16um bin resolution (vs. the manuscript-baseline 8um used by
# bin_vs_breast_cancer_visium_hd.sh). Fully self-contained: generates the sim
# data (if not already present), QC-filters it, then runs all 4
# count_distribution.py modes.
#
# Sim dataset: simulated 16um BINS (real Visium HD ships 2um/8um/16um
# together, so this is a legitimate alternate real configuration, not just a
# tuning knob -- see generate_simulation_noisy.py --bin-size-um help and
# figure.md 2026-08-25 "bin's near-binary empty/non-empty bin sampling").
# Compared against real Visium HD breast cancer resampled to 16um
# (code/real_data_qc/misc/run_visium_hd_qc_16um.sh).
#
# Current winning sim config (2026-08-26, bin16um_logmu_sweep -- see
# figure.md / DATA_VERSIONS.md): packing_pf0p04 (sphere_r_um=2050, ~4% 3D
# packing, same as 8um -- packing is a 3D property independent of bin size),
# log_mu=-2.5 (deliberately more negative than 8um's winning log_mu=0.0,
# since a 16um bin covers ~4x the area of an 8um bin at fixed molecular
# density). theta/theta_jitter/noise_scale/marker_foldchange left at
# generate_simulation_noisy.py defaults (2.0/1.0/1.3/3.5), matching 8um --
# the 16um sweep only varied log_mu (and, in later rounds, theta/jitter,
# neither of which improved on the default).

# 2026-08-27: SIM_TAG now carries the per-modality batch_sigma finalized for
# the SHARED Figure 2 / Figure 3 dataset (cell 1.5, bin8 0.8, bin16 0.7,
# spot 0.3), tuned on the Figure 3 pre/post-Harmony demo then confirmed here
# to still match the real count distribution. count_distribution.py now
# defaults to --slice-id 5 and post-batch counts, so those flags are no
# longer passed per-call below. The matching clustering panel is
# code/clustering/<modality>_celltype_panel.sh (same SIM_TAG).
#
# 2026-08-31: realwindow is now the default (the old sphere-scaled tight
# capture window is retired). generate_simulation_noisy.py no longer scales
# the window with --sphere-r-um; with --capture-window-um unset it uses the
# real instrument window (6.5 x 6.5 mm for Visium / Visium HD). The ~2050 um
# tissue disc sits inside that window with a wide empty border: off-tissue
# bins get domain_true / cell_type_true = "unassigned", obs["is_empty"] = True,
# and are dropped by 00_qc_filter.py (~77% of rows). Only the qc_filtered /
# qc_and_hvg_matched panels are meaningful for the figure; full_panel /
# hvg_matched are now empty-swamped diagnostics. Sim data is (re)generated
# under data/noisy/<tag>/ then moved into data/figure_2/<tag>/; pre-realwindow
# data is archived at data/figure_2_oldwindow_20260831/.
set -euo pipefail

mkdir -p /dcs04/hicks/data/Jan/sim_project/sim_paper/code/count_distribution/logs
source /dcs04/hicks/data/Jan/sim_project/sim_paper/code/count_distribution/_env.sh
cd /dcs04/hicks/data/Jan/sim_project

SIM_TAG="packing_pf0p04_bin16um_log_mu_-2.5_bsigma07"
MODALITY="bin"
REAL_LABEL="breast_cancer_visium_hd_16um"
SIM_RAW="sim_paper/data/figure_2/${SIM_TAG}/simulation_${MODALITY}_z.h5ad"
SIM_QC="sim_paper/data/figure_2/${SIM_TAG}/simulation_${MODALITY}_z_qc.h5ad"
NOISY_DIR="sim_paper/data/noisy/${SIM_TAG}"
REAL_INPUT="sim_paper/data/real_data_qc/${REAL_LABEL}/breast_cancer_visium_hd_qc.h5ad"
OUT_ROOT="sim_paper/data/count_distribution/figure_2/${MODALITY}_vs_${REAL_LABEL}"

# generate_simulation_noisy.py only writes under data/noisy/<out-tag>/, so
# (re)generate + QC there and then move the directory into data/figure_2/.
if [[ ! -f "${SIM_RAW}" ]]; then
    echo "[generate] ${SIM_TAG} not present under figure_2/, building"
    if [[ ! -f "${NOISY_DIR}/simulation_${MODALITY}_z.h5ad" ]]; then
        "${PYTHON_BIN}" sim_paper/code/data/generate_simulation_noisy.py \
            --modality "${MODALITY}" \
            --sphere-r-um 2050 \
            --bin-size-um 16 \
            --base-gene-lognormal -2.5 0.7 \
            --batch-sigma 0.7 \
            --sync-unaligned-seed \
            --out-tag "${SIM_TAG}"
    fi
    if [[ ! -f "${NOISY_DIR}/simulation_${MODALITY}_z_qc.h5ad" ]]; then
        "${PYTHON_BIN}" sim_paper/code/clustering/00_qc_filter.py \
            --modality "${MODALITY}" --packing-tag "${SIM_TAG}"
    fi
    mkdir -p sim_paper/data/figure_2
    mv "${NOISY_DIR}" "sim_paper/data/figure_2/${SIM_TAG}"
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
