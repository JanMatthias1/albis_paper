#!/usr/bin/env bash
# Figure 4B -- run STAGATE 3D spatial-domain identification on the three ALBIS
# z-stacks (bin16um / spot / cell), same Figure 2 tags as Figure 4C.
# Optional $1: a job id to gate all three submissions on
# (--dependency=afterok:<id>), e.g. the QC-regen job. "none" / unset = no dep.
set -euo pipefail

DEP_JOB="${1:-none}"
DEP_ARG=()
[[ "${DEP_JOB}" != "none" ]] && DEP_ARG=(--dependency="afterok:${DEP_JOB}")

SCRIPT="/dcs04/hicks/data/Jan/sim_project/sim_paper/code/applications_albis/spatial_clustering/3D_stagate.py"
CONDA_SH="/jhpce/shared/jhpce/core/anaconda3/2023.03/etc/profile.d/conda.sh"
STAGATE_ENV="/dcs04/hicks/data/Jan/sim_project/sim_paper/env/stagate-pyg"
PYTHON_BIN="${STAGATE_ENV}/bin/python"
OUTDIR="/dcs04/hicks/data/Jan/sim_project/sim_paper/code/applications_albis/spatial_clustering"
LOGDIR="${OUTDIR}/logs_stagate_3D"
mkdir -p "${LOGDIR}"

DATASETS=(
  "bin16um"
  "spot"
  "cell"
)

for d in "${DATASETS[@]}"; do
  sbatch \
  "${DEP_ARG[@]}" \
  --job-name="3D_stagate" \
  --partition=gpu \
  --gres=gpu:l40s:1 \
  --mem=150G \
  --cpus-per-task=4 \
  --time=48:00:00 \
  --output="${LOGDIR}/stagate_${d}_%j.out" \
  --wrap="bash -c 'source ${CONDA_SH} && conda activate ${STAGATE_ENV} && ${PYTHON_BIN} ${SCRIPT} --dataset ${d}'"
done
