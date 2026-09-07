#!/bin/bash
#SBATCH --job-name=batch_sigma_sweep
#SBATCH --output=/dcs04/hicks/data/Jan/sim_project/sim_paper/code/clustering/logs/batch_sigma_sweep_%A_%a.out
#SBATCH --array=0-9%5
#SBATCH --time=08:00:00
#SBATCH --mem=250G
#SBATCH --cpus-per-task=4
#SBATCH --partition=shared
#
# Figure 3 batch_sigma sweep (2026-08-27).
#
# The Figure-2 batch-tuned datasets (data/figure_2/<tag>/, adopted for Figure 3
# on 2026-08-26) are ~2-25x shallower than the datasets the previous batch_sigma
# values were tuned on, so on the new data there is no visible pre/post-Harmony
# slice separation for spot (0.15) or cell (0.22) -- the Harmony-correction demo
# panel shows nothing to correct. bin 8um at 0.5 still shows the "flower petal"
# separation; bin 16um at 0.5 is untested for the demo.
#
# This sweep regenerates each modality with the SAME biology as its Figure 2 tag
# (log_mu / theta / theta_jitter / packing / bin size all unchanged) and ONLY a
# stronger --batch-sigma, into throwaway sweep tags -- Figure 2 is left untouched.
# Once a value per modality is picked from the before/after UMAPs, propagate it
# back to the Figure 2 tags + the *_celltype_panel.sh scripts and re-run both
# figures (user decision 2026-08-27: "do a sweep here for figure 3, once the
# batch is confirmed we can then go back and re-run figure 2").
#
# Per point: generate_simulation_noisy.py -> 00_qc_filter.py -> pca_harmony.py.
# No Leiden/ARI -- the only deliverable is
#   data/figure_3/batch_sigma_sweep/<tag>/pca_harmony_qc/plots/umap_pca_harmony_before_after_by_slice_id.png
#
# Sweep matrix (array index -> modality / batch_sigma):
#   0 spot 0.30   1 spot 0.45   2 spot 0.60
#   3 cell 0.6    4 cell 1.0    5 cell 1.5
#   6 bin8 0.65   7 bin8 0.80
#   8 bin16 0.50  9 bin16 0.70

set -euo pipefail

mkdir -p /dcs04/hicks/data/Jan/sim_project/sim_paper/code/clustering/logs
source /dcs04/hicks/data/Jan/sim_project/sim_paper/code/clustering/_env.sh
cd /dcs04/hicks/data/Jan/sim_project

IDX="${SLURM_ARRAY_TASK_ID:-0}"

# label:   short name used in the sweep tags (bin8/bin16 disambiguate the two
#          bin resolutions; the actual --modality passed to the tools is in MOD)
# MOD:     spot|bin|cell (the only valid --modality values)
# SIGMA:   swept --batch-sigma
# GEN_ARGS: everything else generate_simulation_noisy.py needs, copied verbatim
#          from the matching *_celltype_panel.sh so the biology matches Figure 2.
case "${IDX}" in
  0) LABEL=spot;  MOD=spot; SIGMA=0.30; GEN_ARGS=(--sphere-r-um 2050 --base-gene-lognormal -2.5 0.7) ;;
  1) LABEL=spot;  MOD=spot; SIGMA=0.45; GEN_ARGS=(--sphere-r-um 2050 --base-gene-lognormal -2.5 0.7) ;;
  2) LABEL=spot;  MOD=spot; SIGMA=0.60; GEN_ARGS=(--sphere-r-um 2050 --base-gene-lognormal -2.5 0.7) ;;
  3) LABEL=cell;  MOD=cell; SIGMA=0.6;  GEN_ARGS=(--base-gene-lognormal -2.5 0.7 --theta 0.25 --theta-jitter 0.15) ;;
  4) LABEL=cell;  MOD=cell; SIGMA=1.0;  GEN_ARGS=(--base-gene-lognormal -2.5 0.7 --theta 0.25 --theta-jitter 0.15) ;;
  5) LABEL=cell;  MOD=cell; SIGMA=1.5;  GEN_ARGS=(--base-gene-lognormal -2.5 0.7 --theta 0.25 --theta-jitter 0.15) ;;
  6) LABEL=bin8;  MOD=bin;  SIGMA=0.65; GEN_ARGS=(--sphere-r-um 2050 --base-gene-lognormal 0.0 0.7) ;;
  7) LABEL=bin8;  MOD=bin;  SIGMA=0.80; GEN_ARGS=(--sphere-r-um 2050 --base-gene-lognormal 0.0 0.7) ;;
  8) LABEL=bin16; MOD=bin;  SIGMA=0.50; GEN_ARGS=(--sphere-r-um 2050 --bin-size-um 16 --base-gene-lognormal -2.5 0.7) ;;
  9) LABEL=bin16; MOD=bin;  SIGMA=0.70; GEN_ARGS=(--sphere-r-um 2050 --bin-size-um 16 --base-gene-lognormal -2.5 0.7) ;;
  *) echo "bad array index ${IDX}" >&2; exit 1 ;;
esac

SIGTAG="${SIGMA/./p}"                        # 0.30 -> 0p30, for filesystem-safe tags
TAG="batch_sweep_${LABEL}_bsigma${SIGTAG}"
NOISY_TAG="${TAG}"
NOISY_DIR="sim_paper/data/noisy/${NOISY_TAG}"
SIM_RAW="${NOISY_DIR}/simulation_${MOD}_z.h5ad"
SIM_QC="${NOISY_DIR}/simulation_${MOD}_z_qc.h5ad"
OUT_ROOT="sim_paper/data/figure_3/batch_sigma_sweep/${TAG}"
PCA_OUT="${OUT_ROOT}/pca_harmony_qc/simulation_${MOD}_z_pca_harmony_qc.h5ad"

echo "Job started: $(date)"
echo "Host: $(hostname)   idx=${IDX}"
echo "sweep point: label=${LABEL} modality=${MOD} batch_sigma=${SIGMA}"
echo "tag: ${TAG}"
echo "gen args: ${GEN_ARGS[*]} --batch-sigma ${SIGMA} --out-tag ${NOISY_TAG}"

if [[ ! -f "${SIM_RAW}" ]]; then
    echo "[generate] ${NOISY_TAG}"
    "${PYTHON_BIN}" sim_paper/code/data/generate_simulation_noisy.py \
        --modality "${MOD}" "${GEN_ARGS[@]}" --batch-sigma "${SIGMA}" \
        --out-tag "${NOISY_TAG}"
else
    echo "[generate] ${SIM_RAW} already exists, skipping"
fi

if [[ ! -f "${SIM_QC}" ]]; then
    echo "[qc] ${NOISY_TAG}"
    "${PYTHON_BIN}" sim_paper/code/clustering/00_qc_filter.py \
        --modality "${MOD}" --packing-tag "${NOISY_TAG}"
else
    echo "[qc] ${SIM_QC} already exists, skipping"
fi

echo "[pca_harmony] -> ${PCA_OUT}"
"${PYTHON_BIN}" sim_paper/code/clustering/pca_harmony.py \
    --modality "${MOD}" --input "${SIM_QC}" --output "${PCA_OUT}"

echo "[done] before/after UMAP: ${OUT_ROOT}/pca_harmony_qc/plots/umap_pca_harmony_before_after_by_slice_id.png"
echo "Job finished: $(date)"
