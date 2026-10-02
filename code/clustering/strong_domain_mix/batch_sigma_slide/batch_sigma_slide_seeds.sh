#!/bin/bash
#SBATCH --job-name=batch_slide_seeds
#SBATCH --output=/dcs04/hicks/data/Jan/sim_project/sim_paper/code/clustering/logs/batch_slide_seeds_%A_%a.out
#SBATCH --array=0-71%8
#SBATCH --time=18:00:00
#SBATCH --mem=320G
#SBATCH --cpus-per-task=8
#SBATCH --partition=shared
#
# Seed replication of the batch-σ sweep: seeds 101 and 202 for every point of
# all three modalities (plus extra spot points), so the final plot shows
# mean ± spread. One array task per row of batch_sigma_slide_seed_tasks.tsv
# (written by gen_batch_sigma_slide_seed_tasks.py). Same pipeline as the
# seed-2025 generators: generate -> QC -> BANKSY + Harmony -> ARI ->
# composition_recovery. Each task simulates fresh data (the batch shift is part
# of the counts); settings = Figure 2 + --strong-domain-mix.
set -euo pipefail

mkdir -p /dcs04/hicks/data/Jan/sim_project/sim_paper/code/clustering/logs
source /dcs04/hicks/data/Jan/sim_project/sim_paper/code/clustering/_env_banksy.sh
BANKSY_PYTHON="${PYTHON_BIN}"
source /dcs04/hicks/data/Jan/sim_project/sim_paper/code/clustering/_env.sh
TUTORIAL_PYTHON="${PYTHON_BIN}"
cd /dcs04/hicks/data/Jan/sim_project

OUT_ROOT="sim_paper/data/figure_3/strong_domain_mix/batch_sigma_slide"
TASKS_TSV="sim_paper/code/clustering/strong_domain_mix/batch_sigma_slide/batch_sigma_slide_seed_tasks.tsv"

i="${SLURM_ARRAY_TASK_ID:-0}"
row=$((i + 2))  # +1 for 1-indexing, +1 for the header line
LINE="$(sed -n "${row}p" "${TASKS_TSV}")"
MOD="$(echo "${LINE}" | cut -f1)"
BS="$(echo "${LINE}" | cut -f2)"
SEED="$(echo "${LINE}" | cut -f3)"

FOLDER="${MOD}"
case "${MOD}" in
  cell)
    SPHERE_R_UM=2050
    LAM=0.3; KG=60
    TAG_BODY="log_mu_-2.5_theta_0.40_jitter0.15"
    GEN_FLAGS=(--n-cells 24207 --base-gene-lognormal -2.5 0.7 --theta 0.40 --theta-jitter 0.15 --strong-domain-mix --sync-unaligned-seed)
    ;;
  bin)
    SPHERE_R_UM=2050
    LAM=0.3; KG=100
    TAG_BODY="packing_pf0p04_bin16um_log_mu_-2.5_jitter0.6"
    GEN_FLAGS=(--bin-size-um 16 --base-gene-lognormal -2.5 0.7 --theta 2.0 --theta-jitter 0.6 --strong-domain-mix --sync-unaligned-seed)
    # on-disk folder is bin16um; the tools take --modality bin
    FOLDER="bin16um"
    ;;
  spot)
    SPHERE_R_UM=2050
    LAM=0.5; KG=8
    TAG_BODY="packing_pf0p04_log_mu_-2.25_sigma1.0_theta_0.25_jitter0.10"
    GEN_FLAGS=(--base-gene-lognormal -2.25 1.0 --theta 0.25 --theta-jitter 0.10 --strong-domain-mix --sync-unaligned-seed)
    ;;
esac

if [[ -n "${SEED}" ]]; then
    TAG="${MOD}_batch_slide_bs${BS}_seed${SEED}"
    RUN="${OUT_ROOT}/${FOLDER}/bs${BS}_seed${SEED}"
    SEED_FLAG=(--seed "${SEED}")
else
    TAG="${MOD}_batch_slide_bs${BS}"
    RUN="${OUT_ROOT}/${FOLDER}/bs${BS}"
    SEED_FLAG=()
fi
# Data files are named by their settings: <Figure 2 tag>_strongmix_bsigma<σ>[_seed<n>];
# existing files are reused, so a settings change regenerates.
DATA_TAG="${TAG_BODY}_strongmix_bsigma$(printf '%g' "${BS}" | tr -d .)${SEED:+_seed${SEED}}"
SIM_RAW="${RUN}/${DATA_TAG}.h5ad"
SIM_QC="${RUN}/${DATA_TAG}_qc.h5ad"
mkdir -p "${RUN}"

echo "[task ${i}] modality=${MOD} batch_sigma=${BS} seed=${SEED:-<default>} lambda=${LAM} k_geom=${KG} -> ${RUN}"

if [[ ! -f "${SIM_RAW}" ]]; then
    echo "[generate] ${TAG}"
    "${BANKSY_PYTHON}" sim_paper/code/data/generate_simulation_noisy.py \
        --modality "${MOD}" --sphere-r-um "${SPHERE_R_UM}" \
        "${GEN_FLAGS[@]}" \
        --batch-sigma "${BS}" \
        "${SEED_FLAG[@]}" \
        --output-dir "${RUN}" --output-stem "${DATA_TAG}"
fi
if [[ ! -f "${SIM_QC}" ]]; then
    echo "[qc] ${TAG}"
    "${TUTORIAL_PYTHON}" sim_paper/code/clustering/step00_qc_filter.py \
        --modality "${MOD}" --input "${SIM_RAW}" --output "${SIM_QC}"
fi

mkdir -p "${RUN}"
BANKSY_H5AD="${RUN}/banksy_matrix/simulation_${MOD}_z_banksy_pca_harmony_qc.h5ad"
if [[ ! -f "${BANKSY_H5AD}" ]]; then
    "${BANKSY_PYTHON}" sim_paper/code/clustering/step01_build_banksy_matrix.py \
        --modality "${MOD}" --input "${SIM_QC}" \
        --lambda "${LAM}" --k-geom "${KG}" --max-m 1 \
        --stagger-scale 5 --skip-umap \
        --output-dir "${RUN}/banksy_matrix"
fi

"${TUTORIAL_PYTHON}" sim_paper/code/clustering/step02_leiden_resolution_sweep.py \
    --modality "${MOD}" --packing-tag "cellbin_batch_slide_${TAG}" \
    --input "${BANKSY_H5AD}" --output-dir "${RUN}/ari"

"${TUTORIAL_PYTHON}" sim_paper/code/clustering/composition_recovery.py \
    --h5ad "${RUN}/ari/simulation_${MOD}_z_ari_recovery.h5ad" --tag "${FOLDER}_$(basename "${RUN}")" \
    --out-dir "${OUT_ROOT}/scores"

# Render the manuscript panels from the completed embeddings and labels.
"${TUTORIAL_PYTHON}" sim_paper/code/clustering/plot_banksy_results.py \
    --input "${BANKSY_H5AD}" \
    --cluster-input "${RUN}/ari/simulation_${MOD}_z_ari_recovery.h5ad"

echo "[done] task ${i}  ->  ${RUN}"
