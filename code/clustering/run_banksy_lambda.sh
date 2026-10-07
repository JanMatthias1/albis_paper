#!/bin/bash
#SBATCH --time=12:00:00
#SBATCH --cpus-per-task=8
#SBATCH --partition=shared
#
# One BANKSY run at a fixed lambda = one array task. Used by the Figure 3
# strong_domain_mix/banksy_harmony_batch_zero/ and weak_domain_mix/banksy_harmony_batch_zero/ panels; submit through
# submit_banksy_lambda.sh, which sets --array/--mem/--output per modality.
#
#   sbatch --array=<ids> run_banksy_lambda.sh TASKS_TSV
#
# TASKS_TSV (tab-separated, header row; paths relative to sim_project/):
#   task  modality  k_geom  lambda  qc_h5ad  out_dir  [seed]
# (seed = the simulation seed of qc_h5ad; only the plots group by it)
#
# Each task: step01_build_banksy_matrix.py (BANKSY -> PCA -> Harmony, seeded) ->
# step02_leiden_resolution_sweep.py (Leiden resolution bisected to the true
# category count; scores BOTH domain_true and cell_type_true) -> labels.csv.gz.
# --max-m 1 / --stagger-scale 5 match strong_domain_mix/generate/generate_strong_mix_*.sh, so
# lambda is the only BANKSY setting that differs from the batch-sigma plot.
#
# MODE=sweep (default): ARI only; bin h5ads are deleted afterwards (GBs each).
# MODE=final: keep every h5ad and draw UMAP + spatial plots (plot_banksy_results.py).
# Steps whose outputs exist are skipped; delete out_dir to rerun a task.
set -euo pipefail

TASKS="${1:?usage: run_banksy_lambda.sh TASKS_TSV}"
MODE="${MODE:-sweep}"
CODE=/dcs04/hicks/data/Jan/sim_project/albis_paper/code/clustering
source "${CODE}/_env_banksy.sh"
BANKSY_PYTHON="${PYTHON_BIN}"
source "${CODE}/_env.sh"
TUTORIAL_PYTHON="${PYTHON_BIN}"
export MPLCONFIGDIR="/tmp/fig3-mpl-${USER}-${SLURM_JOB_ID:-local}"
cd /dcs04/hicks/data/Jan/sim_project

TASK="${SLURM_ARRAY_TASK_ID:?run as an array task (see submit_banksy_lambda.sh)}"
LINE="$(awk -F'\t' -v t="${TASK}" 'NR > 1 && $1 == t' "${TASKS}")"
[[ -n "${LINE}" ]] || { echo "no row for task ${TASK} in ${TASKS}" >&2; exit 1; }
IFS=$'\t' read -r _ MOD KG LAM QC RUN _ <<< "${LINE}"

# bin16um is --modality bin to the pipeline scripts
case "${MOD}" in bin8um|bin16um) PMOD=bin ;; *) PMOD="${MOD}" ;; esac

BANKSY_H5AD="${RUN}/banksy_matrix/simulation_${PMOD}_z_banksy_pca_harmony_qc.h5ad"
ARI_H5AD="${RUN}/ari/simulation_${PMOD}_z_ari_recovery.h5ad"
ARI_JSON="${RUN}/ari/ari_summary_${PMOD}.json"
[[ -f "${QC}" ]] || { echo "input missing: ${QC}" >&2; exit 1; }
mkdir -p "${RUN}"

echo "[task ${TASK}] ${MODE} ${MOD} lambda=${LAM} k_geom=${KG}"
echo "[input] ${QC}"
echo "[output] ${RUN}"

if [[ ! -f "${BANKSY_H5AD}" && ( "${MODE}" == final || ! -f "${ARI_JSON}" ) ]]; then
    "${BANKSY_PYTHON}" albis_paper/code/clustering/step01_build_banksy_matrix.py \
        --modality "${PMOD}" --input "${QC}" \
        --lambda "${LAM}" --k-geom "${KG}" --max-m 1 \
        --stagger-scale 5 --skip-umap \
        --output-dir "${RUN}/banksy_matrix"
fi

if [[ ! -f "${ARI_JSON}" || ( "${MODE}" == final && ! -f "${ARI_H5AD}" ) ]]; then
    "${TUTORIAL_PYTHON}" albis_paper/code/clustering/step02_leiden_resolution_sweep.py \
        --modality "${PMOD}" --input "${BANKSY_H5AD}" --output-dir "${RUN}/ari"

    "${TUTORIAL_PYTHON}" - "${ARI_H5AD}" "${RUN}/labels.csv.gz" <<'EOF'
import sys
import anndata as ad
a = ad.read_h5ad(sys.argv[1], backed="r")
cols = [c for c in ("slice_id", "cell_type_true", "domain_true",
                    "leiden_cell_type_true", "leiden_domain_true") if c in a.obs]
df = a.obs[cols].copy()
xy = a.obsm["spatial"][:, :2]
df["x"], df["y"] = xy[:, 0], xy[:, 1]
df.to_csv(sys.argv[2], compression="gzip")
print(f"[save] {sys.argv[2]} ({len(df)} rows)")
EOF
fi

if [[ "${MODE}" == final ]]; then
    "${TUTORIAL_PYTHON}" albis_paper/code/clustering/plot_banksy_results.py \
        --input "${BANKSY_H5AD}" --cluster-input "${ARI_H5AD}"
elif [[ "${PMOD}" == bin && -z "${KEEP_H5AD:-}" ]]; then
    rm -f "${BANKSY_H5AD}" "${ARI_H5AD}"
    echo "[cleanup] removed bin h5ads (labels.csv.gz + ari_summary kept)"
fi

echo "[done] ${RUN}"
