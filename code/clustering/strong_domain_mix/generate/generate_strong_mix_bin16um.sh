#!/bin/bash
#SBATCH --job-name=strong_mix_bin16um
#SBATCH --output=/dcs04/hicks/data/Jan/sim_project/sim_paper/code/clustering/strong_domain_mix/generate/logs/strong_mix_bin16um_%A_%a.out
#SBATCH --array=0-7%4
#SBATCH --time=18:00:00
#SBATCH --mem=320G
#SBATCH --cpus-per-task=8
#SBATCH --partition=shared
#
# Bin16um strong-domain-mix batch-σ sweep, one array task per batch_sigma
# (0, 0.25, 0.30, 0.35, 0.40, 0.45, 0.7, 0.05), seed 2025.
# Data = Figure 2 bin16um settings (r = 2050 µm, 16 µm bins, log_mu -2.5,
# theta 2.0, jitter 0.6) + --strong-domain-mix; nothing else differs.
# Per point: generate -> QC -> BANKSY (λ 0.5, k_geom 100) + Harmony ->
# resolution-matched Leiden ARI -> composition_recovery (ARI + slice leakage)
# -> plots. Raw/QC h5ad are written into the point folder; existing files are reused.
set -euo pipefail

mkdir -p /dcs04/hicks/data/Jan/sim_project/sim_paper/code/clustering/strong_domain_mix/generate/logs
source /dcs04/hicks/data/Jan/sim_project/sim_paper/code/clustering/_env_banksy.sh
BANKSY_PYTHON="${PYTHON_BIN}"
source /dcs04/hicks/data/Jan/sim_project/sim_paper/code/clustering/_env.sh
TUTORIAL_PYTHON="${PYTHON_BIN}"
cd /dcs04/hicks/data/Jan/sim_project

OUT_ROOT="sim_paper/data/figure_3/strong_domain_mix/batch_sigma_slide"
LAM=0.5
KG=100

BATCH_SIGMAS=(0 0.25 0.30 0.35 0.40 0.45 0.7 0.05)
BS="${BATCH_SIGMAS[${SLURM_ARRAY_TASK_ID:-0}]}"
TAG="bin16um_strong_mix_bs${BS}"

RUN="${OUT_ROOT}/bin16um/bs${BS}"
# Data files are named by their settings: <Figure 2 tag>_strongmix_bsigma<σ>[_seed<n>];
# existing files are reused, so a settings change regenerates.
DATA_TAG="packing_pf0p04_bin16um_log_mu_-2.5_jitter0.6_strongmix_bsigma$(printf '%g' "${BS}" | tr -d .)"
SIM_RAW="${RUN}/${DATA_TAG}.h5ad"
SIM_QC="${RUN}/${DATA_TAG}_qc.h5ad"
mkdir -p "${RUN}"

echo "[task ${SLURM_ARRAY_TASK_ID:-0}] bin16um batch_sigma=${BS} lambda=${LAM} k_geom=${KG}"

if [[ ! -f "${SIM_RAW}" ]]; then
    "${BANKSY_PYTHON}" sim_paper/code/data/generate_simulation_noisy.py \
        --modality bin --sphere-r-um 2050 --bin-size-um 16 \
        --base-gene-lognormal -2.5 0.7 --theta 2.0 --theta-jitter 0.6 \
        --strong-domain-mix --batch-sigma "${BS}" --sync-unaligned-seed \
        --output-dir "${RUN}" --output-stem "${DATA_TAG}"
fi
if [[ ! -f "${SIM_QC}" ]]; then
    "${TUTORIAL_PYTHON}" sim_paper/code/clustering/step00_qc_filter.py \
        --modality bin --input "${SIM_RAW}" --output "${SIM_QC}"
fi

mkdir -p "${RUN}"
BANKSY_H5AD="${RUN}/banksy_matrix/simulation_bin_z_banksy_pca_harmony_qc.h5ad"
if [[ ! -f "${BANKSY_H5AD}" ]]; then
    "${BANKSY_PYTHON}" sim_paper/code/clustering/step01_build_banksy_matrix.py \
        --modality bin --input "${SIM_QC}" \
        --lambda "${LAM}" --k-geom "${KG}" --max-m 1 \
        --stagger-scale 5 --skip-umap \
        --output-dir "${RUN}/banksy_matrix"
fi

"${TUTORIAL_PYTHON}" sim_paper/code/clustering/step02_leiden_resolution_sweep.py \
    --modality bin --packing-tag "strong_mix_bin16um_bs${BS}" \
    --input "${BANKSY_H5AD}" --output-dir "${RUN}/ari"

"${TUTORIAL_PYTHON}" sim_paper/code/clustering/composition_recovery.py \
    --h5ad "${RUN}/ari/simulation_bin_z_ari_recovery.h5ad" --tag "bin16um_bs${BS}" \
    --out-dir "${OUT_ROOT}/scores"

# Render the manuscript panels from the completed embeddings and labels.
"${TUTORIAL_PYTHON}" sim_paper/code/clustering/plot_banksy_results.py \
    --input "${BANKSY_H5AD}" \
    --cluster-input "${RUN}/ari/simulation_bin_z_ari_recovery.h5ad"

echo "[done] task ${SLURM_ARRAY_TASK_ID:-0} -> ${RUN}"
