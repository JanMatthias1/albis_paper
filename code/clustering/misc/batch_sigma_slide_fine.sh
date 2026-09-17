#!/bin/bash
#SBATCH --job-name=batch_slide_fine
#SBATCH --output=/dcs04/hicks/data/Jan/sim_project/sim_paper/code/clustering/logs/batch_slide_fine_%A_%a.out
#SBATCH --array=0-12%8
#SBATCH --time=18:00:00
#SBATCH --mem=320G
#SBATCH --cpus-per-task=8
#SBATCH --partition=shared
#
# Finer batch_sigma grid to build a proper curve for all 3 modalities, at
# each modality's own Phase-1 winning BANKSY config:
#   cell:    lambda=0.5/k_geom=200 -- fills the 0-0.5 gap (cliff somewhere in there)
#   bin16um: lambda=0.5/k_geom=100 -- fills the 0.25-0.45 gap (cliff somewhere in there)
#   spot:    lambda=0.1/k_geom=8   -- brackets the anomalous bs=0.2 point densely
# Requires fresh simulations (batch shift baked in at generation time).
set -euo pipefail

mkdir -p /dcs04/hicks/data/Jan/sim_project/sim_paper/code/clustering/logs
source /dcs04/hicks/data/Jan/sim_project/sim_paper/code/clustering/_env_banksy.sh
BANKSY_PYTHON="${PYTHON_BIN}"
source /dcs04/hicks/data/Jan/sim_project/sim_paper/code/clustering/_env.sh
TUTORIAL_PYTHON="${PYTHON_BIN}"
cd /dcs04/hicks/data/Jan/sim_project

OUT_ROOT="sim_paper/data/figure_3/cellbin_batch_sigma_slide"

# 4 cell + 3 bin + 6 spot = 13 tasks
MODS=(cell cell cell cell   bin bin bin   spot spot spot spot spot spot)
BATCH_SIGMAS=(0.1 0.2 0.3 0.4   0.30 0.35 0.40   0.05 0.12 0.15 0.18 0.22 0.25)

i="${SLURM_ARRAY_TASK_ID:-0}"
MOD="${MODS[$i]}"
BS="${BATCH_SIGMAS[$i]}"

case "${MOD}" in
  cell)
    SPHERE_R_UM=6000
    LAM=0.5; KG=200
    GEN_FLAGS=(--base-gene-lognormal -2.3 0.7 --theta 0.40 --theta-jitter 0.15 --strong-domain-mix)
    ;;
  bin)
    SPHERE_R_UM=2050
    LAM=0.5; KG=100
    GEN_FLAGS=(--bin-size-um 16 --base-gene-lognormal -2.5 0.7 --strong-domain-mix)
    ;;
  spot)
    SPHERE_R_UM=2050
    LAM=0.1; KG=8
    GEN_FLAGS=(--base-gene-lognormal -2.5 0.7 --theta 2.0 --theta-jitter 1.0 --strong-domain-mix)
    ;;
esac

TAG="${MOD}_batch_slide_bs${BS}"
SIM_RAW="sim_paper/data/noisy/${TAG}/simulation_${MOD}_z.h5ad"
SIM_QC="sim_paper/data/noisy/${TAG}/simulation_${MOD}_z_qc.h5ad"
RUN="${OUT_ROOT}/${MOD}/bs${BS}"

echo "[task ${i}] modality=${MOD} batch_sigma=${BS} lambda=${LAM} k_geom=${KG}"

if [[ ! -f "${SIM_RAW}" ]]; then
    "${BANKSY_PYTHON}" sim_paper/code/data/generate_simulation_noisy.py \
        --modality "${MOD}" --sphere-r-um "${SPHERE_R_UM}" \
        "${GEN_FLAGS[@]}" \
        --batch-sigma "${BS}" \
        --out-tag "${TAG}"
fi
if [[ ! -f "${SIM_QC}" ]]; then
    "${TUTORIAL_PYTHON}" sim_paper/code/clustering/00_qc_filter.py \
        --modality "${MOD}" --input "${SIM_RAW}" --output "${SIM_QC}"
fi

mkdir -p "${RUN}"
BANKSY_H5AD="${RUN}/banksy_matrix/simulation_${MOD}_z_banksy_pca_harmony_qc.h5ad"
if [[ ! -f "${BANKSY_H5AD}" ]]; then
    "${BANKSY_PYTHON}" sim_paper/code/clustering/01_build_banksy_matrix.py \
        --modality "${MOD}" --input "${SIM_QC}" \
        --lambda "${LAM}" --k-geom "${KG}" --max-m 1 \
        --stagger-scale 5 --skip-umap \
        --output-dir "${RUN}/banksy_matrix"
fi

"${TUTORIAL_PYTHON}" sim_paper/code/clustering/02_leiden_resolution_sweep.py \
    --modality "${MOD}" --packing-tag "cellbin_batch_slide_${MOD}_bs${BS}" \
    --input "${BANKSY_H5AD}" --output-dir "${RUN}/ari"

"${TUTORIAL_PYTHON}" sim_paper/code/clustering/composition_recovery.py \
    --h5ad "${RUN}/ari/simulation_${MOD}_z_ari_recovery.h5ad" --tag "${MOD}_bs${BS}" \
    --out-dir "${OUT_ROOT}/scores"

echo "[done] task ${i}  ->  ${RUN}"
