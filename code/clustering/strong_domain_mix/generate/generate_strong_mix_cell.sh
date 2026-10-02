#!/bin/bash
#SBATCH --job-name=strong_mix_cell
#SBATCH --output=/dcs04/hicks/data/Jan/sim_project/sim_paper/code/clustering/strong_domain_mix/generate/logs/strong_mix_cell_%A_%a.out
#SBATCH --array=0-8%4
#SBATCH --time=06:00:00
#SBATCH --mem=96G
#SBATCH --cpus-per-task=8
#SBATCH --partition=shared
#
# Cell strong-domain-mix batch-σ sweep, one array task per batch_sigma
# (0, 0.05, 0.1, 0.2, 0.3, 0.4, 0.5, 1.0, 1.5), seed 2025.
# Data = Figure 2 cell settings (r = 2050 µm, 24,207 cells, log_mu -2.5,
# theta 0.40, jitter 0.15) + --strong-domain-mix; nothing else differs.
# Per point: generate -> QC -> BANKSY (λ 0.3, k_geom 60) + Harmony ->
# resolution-matched Leiden ARI -> composition_recovery (ARI + slice leakage)
# -> plots. Raw/QC h5ad are written into the point folder; existing files are
# reused. k_geom 60 suits this cell density (larger values smooth over a
# large part of a slice and collapse domain recovery).
# bs 1.5 and 0.05 are also Figure 4B (STAGATE) inputs.
set -euo pipefail

mkdir -p /dcs04/hicks/data/Jan/sim_project/sim_paper/code/clustering/strong_domain_mix/generate/logs
source /dcs04/hicks/data/Jan/sim_project/sim_paper/code/clustering/_env_banksy.sh
BANKSY_PYTHON="${PYTHON_BIN}"
source /dcs04/hicks/data/Jan/sim_project/sim_paper/code/clustering/_env.sh
TUTORIAL_PYTHON="${PYTHON_BIN}"
cd /dcs04/hicks/data/Jan/sim_project

OUT_ROOT="sim_paper/data/figure_3/strong_domain_mix/batch_sigma_slide"
LAM=0.3
KG=60

BATCH_SIGMAS=(0 0.05 0.1 0.2 0.3 0.4 0.5 1.0 1.5)
BS="${BATCH_SIGMAS[${SLURM_ARRAY_TASK_ID:-0}]}"
TAG="cell_strong_mix_bs${BS}"

RUN="${OUT_ROOT}/cell/bs${BS}"
# Data files are named by their settings: <Figure 2 tag>_strongmix_bsigma<σ>[_seed<n>];
# existing files are reused, so a settings change regenerates.
DATA_TAG="log_mu_-2.5_theta_0.40_jitter0.15_strongmix_bsigma$(printf '%g' "${BS}" | tr -d .)"
SIM_RAW="${RUN}/${DATA_TAG}.h5ad"
SIM_QC="${RUN}/${DATA_TAG}_qc.h5ad"
mkdir -p "${RUN}"

echo "[task ${SLURM_ARRAY_TASK_ID:-0}] cell batch_sigma=${BS} lambda=${LAM} k_geom=${KG}"

if [[ ! -f "${SIM_RAW}" ]]; then
    "${BANKSY_PYTHON}" sim_paper/code/data/generate_simulation_noisy.py \
        --modality cell --sphere-r-um 2050 --n-cells 24207 \
        --base-gene-lognormal -2.5 0.7 --theta 0.40 --theta-jitter 0.15 \
        --strong-domain-mix --batch-sigma "${BS}" --sync-unaligned-seed \
        --output-dir "${RUN}" --output-stem "${DATA_TAG}"
fi
if [[ ! -f "${SIM_QC}" ]]; then
    "${TUTORIAL_PYTHON}" sim_paper/code/clustering/step00_qc_filter.py \
        --modality cell --input "${SIM_RAW}" --output "${SIM_QC}"
fi

mkdir -p "${RUN}"
BANKSY_H5AD="${RUN}/banksy_matrix/simulation_cell_z_banksy_pca_harmony_qc.h5ad"
if [[ ! -f "${BANKSY_H5AD}" ]]; then
    "${BANKSY_PYTHON}" sim_paper/code/clustering/step01_build_banksy_matrix.py \
        --modality cell --input "${SIM_QC}" \
        --lambda "${LAM}" --k-geom "${KG}" --max-m 1 \
        --stagger-scale 5 --skip-umap \
        --output-dir "${RUN}/banksy_matrix"
fi

"${TUTORIAL_PYTHON}" sim_paper/code/clustering/step02_leiden_resolution_sweep.py \
    --modality cell --packing-tag "strong_mix_cell_bs${BS}" \
    --input "${BANKSY_H5AD}" --output-dir "${RUN}/ari"

"${TUTORIAL_PYTHON}" sim_paper/code/clustering/composition_recovery.py \
    --h5ad "${RUN}/ari/simulation_cell_z_ari_recovery.h5ad" --tag "cell_bs${BS}" \
    --out-dir "${OUT_ROOT}/scores"

# Render the manuscript panels from the completed embeddings and labels.
"${TUTORIAL_PYTHON}" sim_paper/code/clustering/plot_banksy_results.py \
    --input "${BANKSY_H5AD}" \
    --cluster-input "${RUN}/ari/simulation_cell_z_ari_recovery.h5ad"

echo "[done] task ${SLURM_ARRAY_TASK_ID:-0} -> ${RUN}"
