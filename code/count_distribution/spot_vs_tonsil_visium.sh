#!/bin/bash
#SBATCH --job-name=fig2_spot_vs_tonsil_visium
#SBATCH --output=/dcs04/hicks/data/Jan/sim_project/sim_paper/code/count_distribution/logs/fig2_spot_vs_tonsil_visium_%j.out
#SBATCH --time=08:00:00
#SBATCH --mem=250G
#SBATCH --cpus-per-task=4
#SBATCH --partition=shared
#
# Figure 2, "spot" modality, real reference 2 of 2: simulated Visium-like
# SPOTS (100um spacing, 27.5um capture radius) vs. CytAssist FFPE Protein
# Expression Human Tonsil AddOns (probe-based CytAssist, 18k-gene panel --
# see real_data_qc/run_visium_tonsil_qc.sh for provenance). Run alongside
# spot_vs_lymph_node_visium.sh, which uses the SAME SIM_TAG against the
# first probe reference (Visium V2 Human Lymph Node); spot's config is tuned
# jointly against both, not either one alone.
#
# 2026-09-05: spot retuned onto the two 18k CytAssist probe references.
# breast_cancer_visium (36k whole-transcriptome, fresh-frozen) was dropped
# as spot's tuning target on 2026-09-04.
# base_gene_lognormal log_mu -2.5 -> -2.0, theta (NB dispersion) 2.0 -> 0.25;
# sphere_r_um 2050 (packing_pf0p04) and batch_sigma 0.3 unchanged (batch_sigma
# stays fixed -- it is set on the Figure 3 Harmony demo, explicit user
# constraint). Winner of the probe sweeps (code/data/misc/sweep_spot_*_probe.sh,
# tabulated by summary_spot_logmu_theta_joint.sh): lowest composite
# = sum of |ln(sim/real ratio)| over {theta_hat, total_counts_median,
# matrix_zero_frac}, summed across BOTH probe refs. Tonsil is the LOOSER half
# of the joint fit -- its real HVG-matched theta_hat (~0.34) is well below
# lymph_node's (~0.76), so every config in the sweep overshoots tonsil
# dispersion by >=1.8x; the config is the joint optimum, not a tonsil-specific
# fit. Residual: genes_per_cell close to tonsil here -- pre-existing spot
# genes_per_cell gap, not in the composite, see figure.md Panel C.
#
# 2026-09-06: --theta-jitter added, SIM_TAG bumped to
#   packing_pf0p04_log_mu_-2.0_theta_0.25_jitter0.10_bsigma03
# The shared spot config generated with --theta 0.25 and NO --theta-jitter,
# so generate_simulation_noisy.py fell back to its default THETA_JITTER=1.0.
# Per-gene NB dispersion ~N(0.25, 1.0) -> ~40% of genes floored to 1e-3
# (extreme overdispersion), producing the disjoint upper cloud far above the
# NB fit in mean_variance_compare.png / mean_dropout_compare.png. Same
# artifact fixed for `cell` on 2026-08-25 (--theta-jitter 0.15); the fix was
# never carried to spot through the 2026-09-05 probe retune. --theta-jitter
# 0.10 is the largest of the swept values (0.10/0.15/0.25,
# code/data/misc/sweep_spot_jitter_probe.sh, job 35536033) that fully
# collapses the two clouds into one continuous locus (0.4% of genes floored,
# vs 40%); per-gene theta_hat median ~0.20. theta itself unchanged at 0.25 --
# a re-check now that jitter is sane is still open (figure.md 2026-09-06).
#
# The SIM_TAG data was staged from the jitter sweep build
# (data/noisy/spot_jitter_sweep_0.10/, byte-identical generate flags + the
# --out-tag, rng seeded) into data/figure_2/<tag>/. If absent, the block
# below regenerates it identically. Figure 3's
# code/clustering/spot_celltype_panel.sh reads the same SIM_TAG
# (byte-identical shared dataset).
#
# 2026-08-31 realwindow (still in force): generate_simulation_noisy.py uses
# the real 6.5x6.5 mm Visium window (no --capture-window-um); the ~2050 um
# tissue disc sits inside it with a wide empty border, so ~77% of spots are
# off-tissue before QC and get domain_true/cell_type_true = "unassigned",
# obs["is_empty"] = True, dropped by 00_qc_filter.py. Only qc_filtered /
# qc_and_hvg_matched are meaningful for the figure; full_panel / hvg_matched
# are empty-swamped diagnostics.
#
# Prereq: real_data_qc/run_visium_tonsil_qc.sh (writes
# data/real_data_qc/tonsil_visium/tonsil_visium_qc.h5ad).
set -euo pipefail

mkdir -p /dcs04/hicks/data/Jan/sim_project/sim_paper/code/count_distribution/logs
source /dcs04/hicks/data/Jan/sim_project/sim_paper/code/count_distribution/_env.sh
cd /dcs04/hicks/data/Jan/sim_project

SIM_TAG="packing_pf0p04_log_mu_-2.0_theta_0.25_jitter0.10_bsigma03"
MODALITY="spot"
REAL_LABEL="tonsil_visium"
SIM_RAW="sim_paper/data/figure_2/${SIM_TAG}/simulation_${MODALITY}_z.h5ad"
SIM_QC="sim_paper/data/figure_2/${SIM_TAG}/simulation_${MODALITY}_z_qc.h5ad"
NOISY_DIR="sim_paper/data/noisy/${SIM_TAG}"
REAL_INPUT="sim_paper/data/real_data_qc/${REAL_LABEL}/${REAL_LABEL}_qc.h5ad"
OUT_ROOT="sim_paper/data/count_distribution/figure_2/${MODALITY}_vs_${REAL_LABEL}"

# generate_simulation_noisy.py only writes under data/noisy/<out-tag>/, so
# (re)generate + QC there and then move the directory into data/figure_2/.
if [[ ! -f "${SIM_RAW}" ]]; then
    echo "[generate] ${SIM_TAG} not present under figure_2/, building"
    if [[ ! -f "${NOISY_DIR}/simulation_${MODALITY}_z.h5ad" ]]; then
        "${PYTHON_BIN}" sim_paper/code/data/generate_simulation_noisy.py \
            --modality "${MODALITY}" \
            --sphere-r-um 2050 \
            --base-gene-lognormal -2.0 0.7 \
            --theta 0.25 \
            --theta-jitter 0.10 \
            --batch-sigma 0.3 \
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

if [[ ! -f "${REAL_INPUT}" ]]; then
    echo "ERROR: ${REAL_INPUT} not found -- run real_data_qc/run_visium_tonsil_qc.sh first" >&2
    exit 1
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
