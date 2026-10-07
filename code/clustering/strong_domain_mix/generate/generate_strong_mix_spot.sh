#!/bin/bash
#SBATCH --job-name=strong_mix_spot
#SBATCH --output=/dcs04/hicks/data/Jan/sim_project/albis_paper/code/clustering/strong_domain_mix/generate/logs/strong_mix_spot_%A_%a.out
#SBATCH --array=0-15%8
#SBATCH --time=02:00:00
#SBATCH --mem=64G
#SBATCH --cpus-per-task=4
#SBATCH --partition=shared
#
# Spot strong-domain-mix batch-σ sweep, one array task per batch_sigma
# (16 points from 0 to the canonical 0.3), seed 2025.
# Data = Figure 2 spot settings (r = 2050 µm, log_mu -2.25 / sigma 1.0,
# theta 0.25, jitter 0.10) + --strong-domain-mix; nothing else differs.
# Per point: generate -> QC -> BANKSY (λ 0.5, k_geom 8) + Harmony ->
# resolution-matched Leiden ARI -> composition_recovery (ARI + slice leakage)
# -> plots. Raw/QC h5ad are written into the point folder; existing files are reused.
set -euo pipefail

mkdir -p /dcs04/hicks/data/Jan/sim_project/albis_paper/code/clustering/strong_domain_mix/generate/logs
source /dcs04/hicks/data/Jan/sim_project/albis_paper/code/clustering/_env_banksy.sh
BANKSY_PYTHON="${PYTHON_BIN}"
source /dcs04/hicks/data/Jan/sim_project/albis_paper/code/clustering/_env.sh
TUTORIAL_PYTHON="${PYTHON_BIN}"
cd /dcs04/hicks/data/Jan/sim_project

OUT_ROOT="albis_paper/data/figure_3/strong_domain_mix/batch_sigma_slide"
LAM=0.5
KG=8

BATCH_SIGMAS=(0.0 0.05 0.1 0.12 0.13 0.14 0.15 0.18 0.2 0.22 0.25 0.26 0.27 0.28 0.29 0.3)
BS="${BATCH_SIGMAS[${SLURM_ARRAY_TASK_ID:-0}]}"
TAG="spot_strong_mix_bs${BS}"

RUN="${OUT_ROOT}/spot/bs${BS}"
# Data files are named by their settings: <Figure 2 tag>_strongmix_bsigma<σ>[_seed<n>];
# existing files are reused, so a settings change regenerates.
DATA_TAG="packing_pf0p04_log_mu_-2.25_sigma1.0_theta_0.25_jitter0.10_strongmix_bsigma$(printf '%g' "${BS}" | tr -d .)"
SIM_RAW="${RUN}/${DATA_TAG}.h5ad"
SIM_QC="${RUN}/${DATA_TAG}_qc.h5ad"
mkdir -p "${RUN}"

echo "[task ${SLURM_ARRAY_TASK_ID:-0}] spot batch_sigma=${BS} lambda=${LAM} k_geom=${KG}"

if [[ ! -f "${SIM_RAW}" ]]; then
    "${BANKSY_PYTHON}" albis_paper/code/data/generate_simulation_noisy.py \
        --modality spot --sphere-r-um 2050 \
        --base-gene-lognormal -2.25 1.0 --theta 0.25 --theta-jitter 0.10 \
        --strong-domain-mix --batch-sigma "${BS}" --sync-unaligned-seed \
        --output-dir "${RUN}" --output-stem "${DATA_TAG}"
fi
if [[ ! -f "${SIM_QC}" ]]; then
    "${TUTORIAL_PYTHON}" albis_paper/code/clustering/step00_qc_filter.py \
        --modality spot --input "${SIM_RAW}" --output "${SIM_QC}"
fi

mkdir -p "${RUN}"
BANKSY_H5AD="${RUN}/banksy_matrix/simulation_spot_z_banksy_pca_harmony_qc.h5ad"
if [[ ! -f "${BANKSY_H5AD}" ]]; then
    "${BANKSY_PYTHON}" albis_paper/code/clustering/step01_build_banksy_matrix.py \
        --modality spot --input "${SIM_QC}" \
        --lambda "${LAM}" --k-geom "${KG}" --max-m 1 \
        --stagger-scale 5 --skip-umap \
        --output-dir "${RUN}/banksy_matrix"
fi

"${TUTORIAL_PYTHON}" albis_paper/code/clustering/step02_leiden_resolution_sweep.py \
    --modality spot --packing-tag "strong_mix_spot_bs${BS}" \
    --input "${BANKSY_H5AD}" --output-dir "${RUN}/ari"

"${TUTORIAL_PYTHON}" albis_paper/code/clustering/composition_recovery.py \
    --h5ad "${RUN}/ari/simulation_spot_z_ari_recovery.h5ad" --tag "spot_bs${BS}" \
    --out-dir "${OUT_ROOT}/scores"

# Render the manuscript panels from the completed embeddings and labels.
"${TUTORIAL_PYTHON}" albis_paper/code/clustering/plot_banksy_results.py \
    --input "${BANKSY_H5AD}" \
    --cluster-input "${RUN}/ari/simulation_spot_z_ari_recovery.h5ad"

echo "[done] task ${SLURM_ARRAY_TASK_ID:-0} -> ${RUN}"
