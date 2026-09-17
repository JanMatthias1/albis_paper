#!/bin/bash
#SBATCH --job-name=spot_slide_weakmix
#SBATCH --output=/dcs04/hicks/data/Jan/sim_project/sim_paper/code/clustering/logs/spot_slide_weakmix_%A_%a.out
#SBATCH --array=0-15%8
#SBATCH --time=02:00:00
#SBATCH --mem=64G
#SBATCH --cpus-per-task=4
#SBATCH --partition=shared
#
# Third spot batch_sigma curve: same 16 batch_sigma levels and corrected
# dispersion (log_mu=-2.0, theta=0.25, jitter=0.10) as
# spot_batch_sigma_slide_corrected.sh, but WITHOUT --strong-domain-mix --
# i.e. the real, weak/manuscript-baseline domain_type_mix that
# spot_celltype_panel.sh's actual panel data uses, not the artificially
# strengthened testbed the rest of this investigation is built on.
# Tests whether the dip/cliff shape found on strongmix data (both stale- and
# corrected-dispersion versions) also shows up on real domain data, or is
# itself partly an artifact of strengthening the domain signal.
# Same lambda=0.1/k_geom=8 (not re-tuned). Output kept in its own sibling
# tree, spot_weakmix/, alongside spot/ (corrected-dispersion strongmix,
# renamed from spot_corrected/ 2026-09-15 once adopted as the official
# curve -- the old stale-dispersion spot/ was archived to
# _archive_20260915/cellbin_batch_sigma_slide_spot_staledispersion/).
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

BATCH_SIGMAS=(0.0 0.05 0.1 0.12 0.13 0.14 0.15 0.18 0.2 0.22 0.25 0.26 0.27 0.28 0.29 0.3)
BS="${BATCH_SIGMAS[${SLURM_ARRAY_TASK_ID:-0}]}"
TAG="spot_batch_slide_weakmix_bs${BS}"

SIM_RAW="sim_paper/data/noisy/${TAG}/simulation_spot_z.h5ad"
SIM_QC="sim_paper/data/noisy/${TAG}/simulation_spot_z_qc.h5ad"
RUN="${OUT_ROOT}/spot_weakmix/bs${BS}"

echo "[task ${SLURM_ARRAY_TASK_ID:-0}] spot(weakmix) batch_sigma=${BS} lambda=${LAM} k_geom=${KG}"

if [[ ! -f "${SIM_RAW}" ]]; then
    "${BANKSY_PYTHON}" sim_paper/code/data/generate_simulation_noisy.py \
        --modality spot --sphere-r-um 2050 \
        --base-gene-lognormal -2.0 0.7 --theta 0.25 --theta-jitter 0.10 \
        --batch-sigma "${BS}" \
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
    --modality spot --packing-tag "cellbin_batch_slide_spot_weakmix_bs${BS}" \
    --input "${BANKSY_H5AD}" --output-dir "${RUN}/ari"

"${TUTORIAL_PYTHON}" sim_paper/code/clustering/misc/composition_recovery.py \
    --h5ad "${RUN}/ari/simulation_spot_z_ari_recovery.h5ad" --tag "spot_weakmix_bs${BS}" \
    --out-dir "${OUT_ROOT}/scores"

echo "[done] bs${BS}  ->  ${RUN}"
