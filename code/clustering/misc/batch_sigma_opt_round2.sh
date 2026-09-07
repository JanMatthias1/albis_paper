#!/bin/bash
#SBATCH --job-name=batch_sigma_opt_r2
#SBATCH --output=/dcs04/hicks/data/Jan/sim_project/sim_paper/code/clustering/logs/batch_sigma_opt_r2_%A_%a.out
#SBATCH --array=0-9%5
#SBATCH --time=10:00:00
#SBATCH --mem=250G
#SBATCH --cpus-per-task=4
#SBATCH --partition=shared
#
# Figure 3 batch_sigma round 2 (2026-08-27) -- follows batch_sigma_sweep.sh.
#
# Round 1 picks (before/after-Harmony UMAP only):
#   spot   0.30   clean petals -> full merge
#   bin16  0.70   clear petals -> clean merge, no residual
#   bin8   0.50   already worked; higher only adds a residual knot -> UNCHANGED, not swept here
#   cell   1.5    only value that shows anything, BUT batch at that strength craters
#          cell's Figure 2 dispersion match (theta_hat 0.056 @ sigma 0.22 already
#          vs real lung ~0.16; 1.5 makes it far worse).
#
# This round:
#   idx 0-1  spot 0.30 / bin16 0.70 -- regenerate at the round-1 pick AND run
#            count_distribution vs the real ref, to confirm Figure 2 still matches
#            (spot's batch effect is load-bearing for its theta_hat match:
#            6.9 no-batch -> 1.47 @ 0.15; 0.30 will move it, need to see how far).
#   idx 2-9  cell optimization grid: log_mu {-2.5,-2.3} x theta {0.25,0.40} x
#            batch_sigma {1.0,1.5}, theta_jitter fixed 0.15.
#            - log_mu -2.3 splits the two real cell refs (non_diseased_lung ~73
#              counts / lung_cancer ~107) instead of sitting on non_diseased_lung
#              and undershooting lung_cancer; also lifts genes/cell 33 -> ~52
#              (real 49/63) and marginally helps demo visibility.
#            - theta 0.40 pre-compensates for the batch effect's theta_hat drop so
#              the WITH-batch theta_hat lands back near real's 0.16 rather than 0.056.
#            Each cell point runs count_distribution vs BOTH non_diseased_lung and
#            lung_cancer (full_panel only for this pass; qc/hvg-matched deferred to
#            the final run).
#
# Per point: generate_simulation_noisy.py -> 00_qc_filter.py -> pca_harmony.py
#            -> count_distribution.py (1x for spot/bin16, 2x for cell).
# Deliverables per point:
#   data/figure_3/batch_sigma_sweep/<tag>/pca_harmony_qc/plots/umap_pca_harmony_before_after_by_slice_id.png
#   data/figure_3/batch_sigma_sweep/<tag>/countdist_<reallabel>/comparison_summary.json

set -euo pipefail

mkdir -p /dcs04/hicks/data/Jan/sim_project/sim_paper/code/clustering/logs
source /dcs04/hicks/data/Jan/sim_project/sim_paper/code/clustering/_env.sh
cd /dcs04/hicks/data/Jan/sim_project

IDX="${SLURM_ARRAY_TASK_ID:-0}"
REAL_QC_ROOT="sim_paper/data/real_data_qc"

# LABEL: sweep-tag prefix   MOD: real --modality   GEN_ARGS: generator flags (no --batch-sigma)
# SIGMA: swept --batch-sigma
# COMPARE: space-separated "reallabel:path/to/real_qc.h5ad" entries for count_distribution
case "${IDX}" in
  0) LABEL="spot_bs0.30"; MOD=spot; SIGMA=0.30
     GEN_ARGS=(--sphere-r-um 2050 --base-gene-lognormal -2.5 0.7)
     COMPARE=("breast_cancer_visium:${REAL_QC_ROOT}/breast_cancer_visium/breast_cancer_visium_qc.h5ad") ;;
  1) LABEL="bin16_bs0.70"; MOD=bin; SIGMA=0.70
     GEN_ARGS=(--sphere-r-um 2050 --bin-size-um 16 --base-gene-lognormal -2.5 0.7)
     COMPARE=("breast_cancer_visium_hd_16um:${REAL_QC_ROOT}/breast_cancer_visium_hd_16um/breast_cancer_visium_hd_qc.h5ad") ;;
  2) LABEL="cell_lm-2.5_th0.25_bs1.0"; MOD=cell; SIGMA=1.0
     GEN_ARGS=(--base-gene-lognormal -2.5 0.7 --theta 0.25 --theta-jitter 0.15) ;;
  3) LABEL="cell_lm-2.5_th0.25_bs1.5"; MOD=cell; SIGMA=1.5
     GEN_ARGS=(--base-gene-lognormal -2.5 0.7 --theta 0.25 --theta-jitter 0.15) ;;
  4) LABEL="cell_lm-2.5_th0.40_bs1.0"; MOD=cell; SIGMA=1.0
     GEN_ARGS=(--base-gene-lognormal -2.5 0.7 --theta 0.40 --theta-jitter 0.15) ;;
  5) LABEL="cell_lm-2.5_th0.40_bs1.5"; MOD=cell; SIGMA=1.5
     GEN_ARGS=(--base-gene-lognormal -2.5 0.7 --theta 0.40 --theta-jitter 0.15) ;;
  6) LABEL="cell_lm-2.3_th0.25_bs1.0"; MOD=cell; SIGMA=1.0
     GEN_ARGS=(--base-gene-lognormal -2.3 0.7 --theta 0.25 --theta-jitter 0.15) ;;
  7) LABEL="cell_lm-2.3_th0.25_bs1.5"; MOD=cell; SIGMA=1.5
     GEN_ARGS=(--base-gene-lognormal -2.3 0.7 --theta 0.25 --theta-jitter 0.15) ;;
  8) LABEL="cell_lm-2.3_th0.40_bs1.0"; MOD=cell; SIGMA=1.0
     GEN_ARGS=(--base-gene-lognormal -2.3 0.7 --theta 0.40 --theta-jitter 0.15) ;;
  9) LABEL="cell_lm-2.3_th0.40_bs1.5"; MOD=cell; SIGMA=1.5
     GEN_ARGS=(--base-gene-lognormal -2.3 0.7 --theta 0.40 --theta-jitter 0.15) ;;
  *) echo "bad array index ${IDX}" >&2; exit 1 ;;
esac
COMPARE=("${COMPARE[@]:-}")
[ "${MOD}" = cell ] && COMPARE=(
  "non_diseased_lung:${REAL_QC_ROOT}/non_diseased_lung/non_diseased_lung_qc.h5ad"
  "lung_cancer:${REAL_QC_ROOT}/lung_cancer/lung_cancer_qc.h5ad"
)

TAG="batch_opt2_${LABEL}"
NOISY_DIR="sim_paper/data/noisy/${TAG}"
SIM_RAW="${NOISY_DIR}/simulation_${MOD}_z.h5ad"
SIM_QC="${NOISY_DIR}/simulation_${MOD}_z_qc.h5ad"
OUT_ROOT="sim_paper/data/figure_3/batch_sigma_sweep/${TAG}"
PCA_OUT="${OUT_ROOT}/pca_harmony_qc/simulation_${MOD}_z_pca_harmony_qc.h5ad"

echo "Job started: $(date)   host=$(hostname)   idx=${IDX}"
echo "point: ${LABEL}   modality=${MOD}   batch_sigma=${SIGMA}"
echo "gen: ${GEN_ARGS[*]} --batch-sigma ${SIGMA} --out-tag ${TAG}"

if [[ ! -f "${SIM_RAW}" ]]; then
    echo "[generate] ${TAG}"
    "${PYTHON_BIN}" sim_paper/code/data/generate_simulation_noisy.py \
        --modality "${MOD}" "${GEN_ARGS[@]}" --batch-sigma "${SIGMA}" --out-tag "${TAG}"
else
    echo "[generate] ${SIM_RAW} exists, skipping"
fi

if [[ ! -f "${SIM_QC}" ]]; then
    echo "[qc] ${TAG}"
    "${PYTHON_BIN}" sim_paper/code/clustering/00_qc_filter.py --modality "${MOD}" --packing-tag "${TAG}"
else
    echo "[qc] ${SIM_QC} exists, skipping"
fi

echo "[pca_harmony] -> ${PCA_OUT}"
"${PYTHON_BIN}" sim_paper/code/clustering/pca_harmony.py \
    --modality "${MOD}" --input "${SIM_QC}" --output "${PCA_OUT}"

for entry in "${COMPARE[@]}"; do
    [ -z "${entry}" ] && continue
    rlabel="${entry%%:*}"; rpath="${entry#*:}"
    echo "[count_distribution] vs ${rlabel} (full_panel, slice 5, --use-batch-effect)"
    "${PYTHON_BIN}" sim_paper/code/count_distribution/count_distribution.py \
        --modality "${MOD}" --input "${SIM_RAW}" --slice-id 5 --use-batch-effect \
        --compare-input "${rpath}" --compare-label "${rlabel}" \
        --output-dir "${OUT_ROOT}/countdist_${rlabel}"
done

echo "[done] UMAP: ${OUT_ROOT}/pca_harmony_qc/plots/umap_pca_harmony_before_after_by_slice_id.png"
echo "Job finished: $(date)"
