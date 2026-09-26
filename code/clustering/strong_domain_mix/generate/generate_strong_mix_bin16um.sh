#!/bin/bash
#SBATCH --job-name=strong_mix_bin16um
#SBATCH --output=/dcs04/hicks/data/Jan/sim_project/sim_paper/code/clustering/strong_domain_mix/generate/logs/strong_mix_bin16um_%A_%a.out
#SBATCH --array=0-7%4
#SBATCH --time=18:00:00
#SBATCH --mem=320G
#SBATCH --cpus-per-task=8
#SBATCH --partition=shared
#
# Bin16um's strong-domain-mix batch_sigma slide, consolidated 2026-09-17
# from the scattered points built across code/clustering/misc/
# bin16um_kgeom100_endpoints.sh (bs0/bs0.7) + cellbin_batch_sigma_slide.sh
# (bs0.25/0.45) + batch_sigma_slide_fine.sh (bs0.30/0.35/0.40) into one
# self-contained array script, output moved from
# data/figure_3/cellbin_batch_sigma_slide/bin/ to data/strong_mix/bin16um/.
# Config unchanged: BANKSY lambda=0.5/k_geom=100 (bin16um's own Phase-1
# winner, beats the older k_geom=200 -- 0.806 vs 0.784 domain ARI at
# bs=0), dispersion log_mu=-2.5 (defaults for theta/jitter), batch_sigma
# canonical=0.7, sphere_r_um=2050 (smaller_sphere disc), bin_size_um=16.
#
# Full pipeline per point: generate (strong-domain-mix) -> QC ->
# BANKSY+Harmony -> 02_leiden_resolution_sweep.py -> composition_recovery.py
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

# 0.05 appended 2026-09-23 (index 7): its baseline used to exist only as a
# hand-built folder; the seed script already expects bs0.05 + seeds 101/202.
BATCH_SIGMAS=(0 0.25 0.30 0.35 0.40 0.45 0.7 0.05)
BS="${BATCH_SIGMAS[${SLURM_ARRAY_TASK_ID:-0}]}"
TAG="bin16um_strong_mix_bs${BS}"

RUN="${OUT_ROOT}/bin16um/bs${BS}"
# The data files are named, Figure 2 style, by the parameters they were
# generated with: <Fig2 tag body>_strongmix_bsigma<batch_sigma, no dot>
# [_seed<seed>] (default seed 2025 not written), e.g.
# log_mu_-2.5_theta_0.40_jitter0.15_strongmix_bsigma05_seed101.h5ad / ..._qc.h5ad.
# The skip-if-exists checks use that name, so a parameter change regenerates.
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
    "${TUTORIAL_PYTHON}" sim_paper/code/clustering/00_qc_filter.py \
        --modality bin --input "${SIM_RAW}" --output "${SIM_QC}"
fi

mkdir -p "${RUN}"
BANKSY_H5AD="${RUN}/banksy_matrix/simulation_bin_z_banksy_pca_harmony_qc.h5ad"
if [[ ! -f "${BANKSY_H5AD}" ]]; then
    "${BANKSY_PYTHON}" sim_paper/code/clustering/01_build_banksy_matrix.py \
        --modality bin --input "${SIM_QC}" \
        --lambda "${LAM}" --k-geom "${KG}" --max-m 1 \
        --stagger-scale 5 --skip-umap \
        --output-dir "${RUN}/banksy_matrix"
fi

"${TUTORIAL_PYTHON}" sim_paper/code/clustering/02_leiden_resolution_sweep.py \
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
