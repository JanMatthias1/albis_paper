#!/bin/bash
#SBATCH --job-name=cellbin_batch_slide
#SBATCH --output=/dcs04/hicks/data/Jan/sim_project/sim_paper/code/clustering/logs/cellbin_batch_slide_%A_%a.out
#SBATCH --array=0-3%4
#SBATCH --time=18:00:00
#SBATCH --mem=320G
#SBATCH --cpus-per-task=8
#SBATCH --partition=shared
#
# Phase 2: slide batch_sigma at each modality's own Phase-1 ceiling config
# (cell: lambda=0.5/k_geom=200, unchanged, already near-ceiling; bin16um:
# lambda=0.5/k_geom=100, the Phase-1 winner, 0.806 vs the old k_geom=200's
# 0.784) to find the recoverability threshold -- does domain ARI degrade
# smoothly toward the canonical batch_sigma, or fall off a cliff somewhere
# below it? Requires FRESH simulations (the batch shift is baked into counts
# at generation time, not reprocessable from existing h5ads) at 2
# intermediate batch_sigma points per modality:
#   cell:    canonical batch_sigma=1.5 -> test 0.5, 1.0   (already have 0, 1.5)
#   bin16um: canonical batch_sigma=0.7 -> test 0.25, 0.45 (already have 0, 0.7)
# Full pipeline per task: generate (strong-domain-mix, same base counts
# config as the existing strongmix tags) -> QC -> BANKSY+Harmony (at the
# config above) -> 02_leiden_resolution_sweep.py -> composition_recovery.py (leak).
set -euo pipefail

mkdir -p /dcs04/hicks/data/Jan/sim_project/sim_paper/code/clustering/logs
source /dcs04/hicks/data/Jan/sim_project/sim_paper/code/clustering/_env_banksy.sh
BANKSY_PYTHON="${PYTHON_BIN}"
source /dcs04/hicks/data/Jan/sim_project/sim_paper/code/clustering/_env.sh
TUTORIAL_PYTHON="${PYTHON_BIN}"
cd /dcs04/hicks/data/Jan/sim_project

OUT_ROOT="sim_paper/data/figure_3/cellbin_batch_sigma_slide"

# task -> (modality, batch_sigma, lambda, k_geom)
MODS=(cell cell bin bin)
BATCH_SIGMAS=(0.5 1.0 0.25 0.45)
LAMBDAS=(0.5 0.5 0.5 0.5)
KGEOMS=(200 200 100 100)

i="${SLURM_ARRAY_TASK_ID:-0}"
MOD="${MODS[$i]}"
BS="${BATCH_SIGMAS[$i]}"
LAM="${LAMBDAS[$i]}"
KG="${KGEOMS[$i]}"

if [[ "${MOD}" == "cell" ]]; then
    SPHERE_R_UM=6000
    GEN_FLAGS=(--base-gene-lognormal -2.3 0.7 --theta 0.40 --theta-jitter 0.15 --strong-domain-mix)
else
    SPHERE_R_UM=2050
    GEN_FLAGS=(--bin-size-um 16 --base-gene-lognormal -2.5 0.7 --strong-domain-mix)
fi

TAG="${MOD}_batch_slide_bs${BS}"
SIM_RAW="sim_paper/data/noisy/${TAG}/simulation_${MOD}_z.h5ad"
SIM_QC="sim_paper/data/noisy/${TAG}/simulation_${MOD}_z_qc.h5ad"
RUN="${OUT_ROOT}/${MOD}/bs${BS}"

echo "[task ${i}] modality=${MOD} batch_sigma=${BS} lambda=${LAM} k_geom=${KG}"

if [[ ! -f "${SIM_RAW}" ]]; then
    echo "[generate] ${TAG}"
    "${BANKSY_PYTHON}" sim_paper/code/data/generate_simulation_noisy.py \
        --modality "${MOD}" --sphere-r-um "${SPHERE_R_UM}" \
        "${GEN_FLAGS[@]}" \
        --batch-sigma "${BS}" \
        --out-tag "${TAG}"
fi
if [[ ! -f "${SIM_QC}" ]]; then
    echo "[qc] ${TAG}"
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
