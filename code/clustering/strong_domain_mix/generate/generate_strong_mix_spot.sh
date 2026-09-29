#!/bin/bash
#SBATCH --job-name=strong_mix_spot
#SBATCH --output=/dcs04/hicks/data/Jan/sim_project/sim_paper/code/clustering/strong_domain_mix/generate/logs/strong_mix_spot_%A_%a.out
#SBATCH --array=0-15%8
#SBATCH --time=02:00:00
#SBATCH --mem=64G
#SBATCH --cpus-per-task=4
#SBATCH --partition=shared
#
# Spot's strong-domain-mix batch_sigma slide (no batch effect -> canonical
# 0.3), consolidated 2026-09-17 from code/misc/clustering/misc/
# spot_batch_sigma_slide_corrected.sh into data/strong_mix/ (was
# data/figure_3/cellbin_batch_sigma_slide/). Logic unchanged -- spot's config
# (BANKSY lambda=0.1/k_geom=8, dispersion log_mu=-2.0/theta=0.25/jitter=0.10)
# was already a single clean array script covering the full 16-point grid,
# just relocated.
#
# Full pipeline per point: generate (strong-domain-mix) -> QC ->
# BANKSY+Harmony -> step02_leiden_resolution_sweep.py -> composition_recovery.py
# (ARI + slice_id leakage).
# 2026-09-23: generate flags matched to the retuned Figure 2 config so the
# normal- and strong-mix datasets differ ONLY by --strong-domain-mix
# (bin16um: --theta 2.0 --theta-jitter 0.6; spot: --base-gene-lognormal
# -2.25 1.0, no per-domain depth factors; all: --sync-unaligned-seed, which only changes
# obsm['spatial_unaligned'] -- verified byte-identical counts/spatial/labels).
# Old figure_3/ and data/noisy/*_strong_mix_*/*_batch_slide_* inputs archived to
# data/figure_3_archive_20260923/ and data/noisy/_archive_figure3_20260923/,
# so every skip-if-exists step below regenerates from scratch.
# 2026-09-23 (later): ONE TREE PER POINT -- raw/QC h5ad are written straight
# into the point folder ${RUN} (next to banksy_matrix/ and ari/) via
# --output-dir, instead of data/noisy/<tag>/ + hand-made symlinks. Existing
# inputs were moved there from data/noisy/ (verified: Figure 2 config +
# --strong-domain-mix, check_matches_figure2.py).
# 2026-09-23 (latest): spot --domain-size-factors DROPPED by user decision
# (per-domain depth, spot-only, made domains partly identifiable from depth);
# spot now = Figure 2 packing_pf0p04_log_mu_-2.25_sigma1.0_theta_0.25_jitter0.10_bsigma03.
set -euo pipefail

mkdir -p /dcs04/hicks/data/Jan/sim_project/sim_paper/code/clustering/strong_domain_mix/generate/logs
source /dcs04/hicks/data/Jan/sim_project/sim_paper/code/clustering/_env_banksy.sh
BANKSY_PYTHON="${PYTHON_BIN}"
source /dcs04/hicks/data/Jan/sim_project/sim_paper/code/clustering/_env.sh
TUTORIAL_PYTHON="${PYTHON_BIN}"
cd /dcs04/hicks/data/Jan/sim_project

OUT_ROOT="sim_paper/data/figure_3/strong_domain_mix/batch_sigma_slide"
LAM=0.1
KG=8

BATCH_SIGMAS=(0.0 0.05 0.1 0.12 0.13 0.14 0.15 0.18 0.2 0.22 0.25 0.26 0.27 0.28 0.29 0.3)
BS="${BATCH_SIGMAS[${SLURM_ARRAY_TASK_ID:-0}]}"
TAG="spot_strong_mix_bs${BS}"

RUN="${OUT_ROOT}/spot/bs${BS}"
# The data files are named, Figure 2 style, by the parameters they were
# generated with: <Fig2 tag body>_strongmix_bsigma<batch_sigma, no dot>
# [_seed<seed>] (default seed 2025 not written), e.g.
# log_mu_-2.5_theta_0.40_jitter0.15_strongmix_bsigma05_seed101.h5ad / ..._qc.h5ad.
# The skip-if-exists checks use that name, so a parameter change regenerates.
DATA_TAG="packing_pf0p04_log_mu_-2.25_sigma1.0_theta_0.25_jitter0.10_strongmix_bsigma$(printf '%g' "${BS}" | tr -d .)"
SIM_RAW="${RUN}/${DATA_TAG}.h5ad"
SIM_QC="${RUN}/${DATA_TAG}_qc.h5ad"
mkdir -p "${RUN}"

echo "[task ${SLURM_ARRAY_TASK_ID:-0}] spot batch_sigma=${BS} lambda=${LAM} k_geom=${KG}"

if [[ ! -f "${SIM_RAW}" ]]; then
    "${BANKSY_PYTHON}" sim_paper/code/data/generate_simulation_noisy.py \
        --modality spot --sphere-r-um 2050 \
        --base-gene-lognormal -2.25 1.0 --theta 0.25 --theta-jitter 0.10 \
        --strong-domain-mix --batch-sigma "${BS}" --sync-unaligned-seed \
        --output-dir "${RUN}" --output-stem "${DATA_TAG}"
fi
if [[ ! -f "${SIM_QC}" ]]; then
    "${TUTORIAL_PYTHON}" sim_paper/code/clustering/step00_qc_filter.py \
        --modality spot --input "${SIM_RAW}" --output "${SIM_QC}"
fi

mkdir -p "${RUN}"
BANKSY_H5AD="${RUN}/banksy_matrix/simulation_spot_z_banksy_pca_harmony_qc.h5ad"
if [[ ! -f "${BANKSY_H5AD}" ]]; then
    "${BANKSY_PYTHON}" sim_paper/code/clustering/step01_build_banksy_matrix.py \
        --modality spot --input "${SIM_QC}" \
        --lambda "${LAM}" --k-geom "${KG}" --max-m 1 \
        --stagger-scale 5 --skip-umap \
        --output-dir "${RUN}/banksy_matrix"
fi

"${TUTORIAL_PYTHON}" sim_paper/code/clustering/step02_leiden_resolution_sweep.py \
    --modality spot --packing-tag "strong_mix_spot_bs${BS}" \
    --input "${BANKSY_H5AD}" --output-dir "${RUN}/ari"

"${TUTORIAL_PYTHON}" sim_paper/code/clustering/composition_recovery.py \
    --h5ad "${RUN}/ari/simulation_spot_z_ari_recovery.h5ad" --tag "spot_bs${BS}" \
    --out-dir "${OUT_ROOT}/scores"

# Render the manuscript panels from the completed embeddings and labels.
"${TUTORIAL_PYTHON}" sim_paper/code/clustering/plot_banksy_results.py \
    --input "${BANKSY_H5AD}" \
    --cluster-input "${RUN}/ari/simulation_spot_z_ari_recovery.h5ad"

echo "[done] task ${SLURM_ARRAY_TASK_ID:-0} -> ${RUN}"
