#!/bin/bash
#SBATCH --job-name=sweep_spot_jitter_probe
#SBATCH --output=/dcs04/hicks/data/Jan/sim_project/sim_paper/code/data/logs/sweep_spot_jitter_probe_%A_%a.out
#SBATCH --array=0-2
#SBATCH --time=08:00:00
#SBATCH --mem=250G
#SBATCH --cpus-per-task=4
#SBATCH --partition=shared
#
# spot theta_jitter sweep vs. the two 18k CytAssist probe-based real refs
# (lymph_node_visium, tonsil_visium).
#
# WHY: spot's mean_variance_compare.png shows two disjoint gene clouds -- a
# dense upper cloud sitting far above the NB fit, and the main cloud on it.
# Root-caused 2026-09-06: the shared spot config generates with --theta 0.25
# and NO --theta-jitter, so it falls back to generate_simulation_noisy.py's
# default NOISY_THETA_JITTER=1.0. Per-gene NB dispersion is then drawn ~
# N(0.25, 1.0) -- mostly <= 0, floored to 1e-3 -- so ~40% of genes end up
# pinned at the floor (extreme overdispersion -> upper cloud) and the rest
# sit at the nominal ~0.25 (main cloud). Nothing in between. This is the
# SAME artifact root-caused for the `cell` modality on 2026-08-25 and fixed
# there with --theta-jitter 0.15 (cell tag ..._jitter0.15_...); the fix was
# never propagated to spot when it was retuned onto the probe refs 2026-09-05.
#
# This sweep brackets --theta-jitter (0.10 / 0.15 / 0.25) at the current
# shared spot config otherwise held fixed: log_mu=-2.0, theta=0.25,
# batch_sigma=0.3, sphere_r_um=2050 (packing_pf0p04). Goal: the largest
# jitter that still collapses the two clouds into one continuous locus and
# keeps a realistic per-gene scatter (cell found 0.05 over-flattened), while
# theta_hat stays in range vs both probe refs. Whichever wins gets promoted
# to the shared SIM_TAG (..._theta_0.25_jitter<val>_bsigma03) in
# code/count_distribution/spot_vs_{lymph_node,tonsil}_visium.sh AND
# code/clustering/spot_celltype_panel.sh (byte-identical shared dataset),
# with theta re-checked afterward -- the 2026-09-05 joint theta sweep was run
# with jitter=1.0 and the floored-gene subpopulation drags the median
# theta_hat down, so theta=0.25 may not be the optimum once jitter is sane.
#
# Two comparison modes per (jitter, ref), matching the naming the other spot
# sweeps + summary_spot_logmu_theta_joint.sh use:
#   qc_filtered   -- QC'd sim, full 556-gene panel. This is the mode the
#                    figure reads spot's mean_variance / mean_dropout from
#                    (full_panel is empty-swamped post-realwindow), so it's
#                    where the two-cloud fix has to show up visually.
#   hvg_matched   -- QC'd sim + --match-panel-size (real refs HVG-subset to
#                    556 genes). Apples-to-apples theta_hat / zero_frac; the
#                    number the cross-sweep composite is scored on.
#
# Prereq: real_data_qc/run_visium_tonsil_qc.sh + run_visium_lymph_qc.sh
# already run (data/real_data_qc/{tonsil,lymph_node}_visium/*_qc.h5ad).
#
# Tabulate with: sim_paper/code/data/misc/summary_spot_jitter.sh
set -euo pipefail

mkdir -p /dcs04/hicks/data/Jan/sim_project/sim_paper/code/data/logs
source /dcs04/hicks/data/Jan/sim_project/sim_paper/code/count_distribution/_env.sh
cd /dcs04/hicks/data/Jan/sim_project

JITTERS=(0.10 0.15 0.25)
JITTER="${JITTERS[${SLURM_ARRAY_TASK_ID:-0}]}"
LOG_MU=-2.0
THETA=0.25
BATCH_SIGMA=0.3
SPHERE_R_UM=2050
TAG="spot_jitter_sweep_${JITTER}"
MODALITY="spot"

SIM_RAW="sim_paper/data/noisy/${TAG}/simulation_${MODALITY}_z.h5ad"
SIM_QC="sim_paper/data/noisy/${TAG}/simulation_${MODALITY}_z_qc.h5ad"
OUT_ROOT="sim_paper/data/count_distribution/sweeps/spot_jitter_sweep_probe"

echo "[config] jitter=${JITTER} theta=${THETA} log_mu=${LOG_MU} batch_sigma=${BATCH_SIGMA} tag=${TAG}"

if [[ ! -f "${SIM_RAW}" ]]; then
    echo "[generate] ${TAG} not found, generating"
    "${PYTHON_BIN}" sim_paper/code/data/generate_simulation_noisy.py \
        --modality "${MODALITY}" \
        --sphere-r-um "${SPHERE_R_UM}" \
        --base-gene-lognormal "${LOG_MU}" 0.7 \
        --theta "${THETA}" \
        --theta-jitter "${JITTER}" \
        --batch-sigma "${BATCH_SIGMA}" \
        --out-tag "${TAG}"
fi

if [[ ! -f "${SIM_QC}" ]]; then
    echo "[qc] ${SIM_QC} not found, running 00_qc_filter.py"
    "${PYTHON_BIN}" sim_paper/code/clustering/00_qc_filter.py --modality "${MODALITY}" --packing-tag "${TAG}"
fi

for real in tonsil_visium lymph_node_visium; do
    REAL_INPUT="sim_paper/data/real_data_qc/${real}/${real}_qc.h5ad"

    echo "[compare] ${TAG} vs ${real} (qc_filtered)"
    "${PYTHON_BIN}" sim_paper/code/count_distribution/count_distribution.py \
        --modality "${MODALITY}" --input "${SIM_QC}" \
        --compare-input "${REAL_INPUT}" --compare-label "${real}" \
        --output-dir "${OUT_ROOT}/${TAG}_qc_filtered_vs_${real}"

    echo "[compare] ${TAG} vs ${real} (hvg_matched)"
    "${PYTHON_BIN}" sim_paper/code/count_distribution/count_distribution.py \
        --modality "${MODALITY}" --input "${SIM_QC}" \
        --compare-input "${REAL_INPUT}" --compare-label "${real}" \
        --match-panel-size \
        --output-dir "${OUT_ROOT}/${TAG}_hvg_matched_vs_${real}"
done

echo "[done] ${TAG}"
