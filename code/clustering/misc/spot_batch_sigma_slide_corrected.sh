#!/bin/bash
#SBATCH --job-name=spot_slide_corrected
#SBATCH --output=/dcs04/hicks/data/Jan/sim_project/sim_paper/code/clustering/logs/spot_slide_corrected_%A_%a.out
#SBATCH --array=0-15%8
#SBATCH --time=02:00:00
#SBATCH --mem=64G
#SBATCH --cpus-per-task=4
#SBATCH --partition=shared
#
# Re-does the entire spot batch_sigma slide (all 16 levels tested so far:
# the original coarse+fine grid plus the dense dip-bracket) using spot's
# ACTUAL current calibrated dispersion instead of generate_strongmix.sh's
# stale defaults. Found 2026-09-14: every prior spot domain-sweep dataset
# (including the one that found lambda=0.1/k_geom=8 in the first place, and
# the whole dip/cliff curve) was generated with log_mu=-2.5, theta=2.0,
# theta-jitter=1.0 (generate_simulation_noisy.py's bare defaults) because
# generate_strongmix.sh only ever added explicit --theta/--theta-jitter for
# the cell branch, not spot. Figure 2's spot config was retuned 2026-09-05/06
# to log_mu=-2.0, theta=0.25, jitter=0.10 (fixing a "two clouds" artifact) and
# generate_strongmix.sh was never updated to match. batch_sigma itself (0.3
# canonical) is unaffected -- this only fixes the dispersion/mean it sits on.
#
# Same lambda=0.1/k_geom=8 BANKSY config as before (not re-tuned here --
# that's a separate follow-up if this curve's shape/optimum looks off).
# Output is a SEPARATE directory from the original (stale-param) curve, not
# a replacement, so the two can be compared directly.
set -euo pipefail

mkdir -p /dcs04/hicks/data/Jan/sim_project/sim_paper/code/clustering/logs
source /dcs04/hicks/data/Jan/sim_project/sim_paper/code/clustering/_env_banksy.sh
BANKSY_PYTHON="${PYTHON_BIN}"
source /dcs04/hicks/data/Jan/sim_project/sim_paper/code/clustering/_env.sh
TUTORIAL_PYTHON="${PYTHON_BIN}"
cd /dcs04/hicks/data/Jan/sim_project

OUT_ROOT="sim_paper/data/figure_3/cellbin_batch_sigma_slide"
LAM=0.1
KG=8

# Same 16 batch_sigma values as the stale-param curve (original + fine +
# dense bracket combined), for a direct point-by-point comparison.
BATCH_SIGMAS=(0.0 0.05 0.1 0.12 0.13 0.14 0.15 0.18 0.2 0.22 0.25 0.26 0.27 0.28 0.29 0.3)
BS="${BATCH_SIGMAS[${SLURM_ARRAY_TASK_ID:-0}]}"
TAG="spot_batch_slide_corrected_bs${BS}"

SIM_RAW="sim_paper/data/noisy/${TAG}/simulation_spot_z.h5ad"
SIM_QC="sim_paper/data/noisy/${TAG}/simulation_spot_z_qc.h5ad"
RUN="${OUT_ROOT}/spot/bs${BS}"

echo "[task ${SLURM_ARRAY_TASK_ID:-0}] spot(corrected) batch_sigma=${BS} lambda=${LAM} k_geom=${KG}"

if [[ ! -f "${SIM_RAW}" ]]; then
    "${BANKSY_PYTHON}" sim_paper/code/data/generate_simulation_noisy.py \
        --modality spot --sphere-r-um 2050 \
        --base-gene-lognormal -2.0 0.7 --theta 0.25 --theta-jitter 0.10 \
        --strong-domain-mix --batch-sigma "${BS}" \
        --out-tag "${TAG}"
fi
if [[ ! -f "${SIM_QC}" ]]; then
    "${TUTORIAL_PYTHON}" sim_paper/code/clustering/00_qc_filter.py \
        --modality spot --input "${SIM_RAW}" --output "${SIM_QC}"
fi

mkdir -p "${RUN}"
BANKSY_H5AD="${RUN}/banksy_matrix/simulation_spot_z_banksy_pca_harmony_qc.h5ad"
if [[ ! -f "${BANKSY_H5AD}" ]]; then
    "${BANKSY_PYTHON}" sim_paper/code/clustering/01_build_banksy_matrix.py \
        --modality spot --input "${SIM_QC}" \
        --lambda "${LAM}" --k-geom "${KG}" --max-m 1 \
        --stagger-scale 5 --skip-umap \
        --output-dir "${RUN}/banksy_matrix"
fi

"${TUTORIAL_PYTHON}" sim_paper/code/clustering/02_leiden_resolution_sweep.py \
    --modality spot --packing-tag "cellbin_batch_slide_spot_corrected_bs${BS}" \
    --input "${BANKSY_H5AD}" --output-dir "${RUN}/ari"

"${TUTORIAL_PYTHON}" sim_paper/code/clustering/composition_recovery.py \
    --h5ad "${RUN}/ari/simulation_spot_z_ari_recovery.h5ad" --tag "spot_bs${BS}" \
    --out-dir "${OUT_ROOT}/scores"

echo "[done] bs${BS}  ->  ${RUN}"
