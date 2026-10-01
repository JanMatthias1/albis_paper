#!/usr/bin/env bash
# Figure 4B replot -- regenerate domains_3d_true_vs_stagate.png / umap_stagate3d*.png
# for every existing STAGATE run, then the fig4b summary figure.
# Pure CPU: plot_stagate_figures.py / plot_stagate_fig4b.py only read the
# already-written adata_results/{h5ad,metrics.json} -- no GPU, no retraining.
# Use this after a plotting-code change (as opposed to run_stagate_3D.sh, which
# retrains).
#
# Discovers run dirs from what's on disk under STAGATE_DIR (matching
# <mod>_strongmix_bs<value>_seed<seed>, the current stagate_inputs.py grid)
# rather than hardcoding, so it stays in sync with whatever 3D_stagate.py has
# actually produced.
#
# Optional $1: a job id to gate every submission on (--dependency=afterok:<id>).
#   "none" / unset = no dependency.
set -euo pipefail

DEP_JOB="${1:-none}"
DEP_ARG=()
[[ "${DEP_JOB}" != "none" ]] && DEP_ARG=(--dependency="afterok:${DEP_JOB}")

CODE_DIR="/dcs04/hicks/data/Jan/sim_project/sim_paper/code/applications_albis/spatial_clustering"
STAGATE_DIR="/dcs04/hicks/data/Jan/sim_project/sim_paper/data/figure_4/spatial_clustering/STAGATE"
CONDA_SH="/jhpce/shared/jhpce/core/anaconda3/2023.03/etc/profile.d/conda.sh"
# plot_stagate_figures.py only needs scanpy/anndata/matplotlib/sklearn (no
# torch/STAGATE_pyG/R), but those all live in the same stagate-pyg env used
# for training, so just reuse it rather than standing up a second env.
STAGATE_ENV="/dcs04/hicks/data/Jan/sim_project/sim_paper/env/stagate-pyg"
PYTHON_BIN="${STAGATE_ENV}/bin/python"
LOGDIR="${CODE_DIR}/logs_stagate_plots_3D"
mkdir -p "${LOGDIR}"

PARTITION="${PARTITION:-shared}"

# ---- discover run dirs -> (dataset, use_pre_batch) pairs -------------------
RUN_TAGS=()
for d in "${STAGATE_DIR}"/*/; do
  tag="$(basename "${d}")"
  [[ "${tag}" =~ ^(cell|bin16um|spot)_strongmix_bs[0-9.]+_seed(2025|101|202)$ ]] || continue
  [[ -f "${d}adata_results/metrics.json" ]] || continue
  RUN_TAGS+=("${tag}")
done
if [[ ${#RUN_TAGS[@]} -eq 0 ]]; then
  echo "no STAGATE run dirs with adata_results/metrics.json found under ${STAGATE_DIR}" >&2
  exit 1
fi

echo "replotting ${#RUN_TAGS[@]} run(s): ${RUN_TAGS[*]}"

JOB_IDS=()
for tag in "${RUN_TAGS[@]}"; do
  dataset="${tag}"

  jid=$(sbatch --parsable \
    "${DEP_ARG[@]}" \
    --job-name="stagate_plot" \
    --partition="${PARTITION}" \
    --mem=64G \
    --cpus-per-task=8 \
    --time=02:00:00 \
    --output="${LOGDIR}/plot_${tag}_%j.out" \
    --wrap="bash -c 'source ${CONDA_SH} && conda activate ${STAGATE_ENV} && ${PYTHON_BIN} ${CODE_DIR}/plot_stagate_figures.py --dataset ${dataset}'")
  echo "  ${tag} -> job ${jid}"
  JOB_IDS+=("${jid}")
done

# ---- final job: regenerate the two cross-run summary figures once every ---
# per-dataset replot above has finished (afterany: they should all succeed,
# but one flaky node shouldn't block the summary from covering the rest).
DEP_LIST=$(IFS=:; echo "${JOB_IDS[*]}")
sbatch \
  --dependency="afterany:${DEP_LIST}" \
  --job-name="stagate_plot_summary" \
  --partition="${PARTITION}" \
  --mem=16G \
  --cpus-per-task=2 \
  --time=00:30:00 \
  --output="${LOGDIR}/summary_%j.out" \
  --wrap="bash -c 'source ${CONDA_SH} && conda activate ${STAGATE_ENV} && ${PYTHON_BIN} ${CODE_DIR}/plot_stagate_fig4b.py'"

echo "submitted ${#JOB_IDS[@]} replot jobs + 1 summary job (partition=${PARTITION})"
