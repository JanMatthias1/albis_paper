#!/usr/bin/env bash
# Figure 4C -- run STAIR 3D alignment on the three ALBIS z-stacks.
# Optional $1: a job id to gate all three submissions on
# (--dependency=afterok:<id>), e.g. the QC-regen job. "none" / unset = no dep.
set -euo pipefail

DEP_JOB="${1:-none}"
DEP_ARG=()
[[ "${DEP_JOB}" != "none" ]] && DEP_ARG=(--dependency="afterok:${DEP_JOB}")

SCRIPT="/dcs04/hicks/data/Jan/sim_project/sim_paper/code/applications_albis/alignment/3D_stair.py"
CONDA_SH="/jhpce/shared/jhpce/core/anaconda3/2023.03/etc/profile.d/conda.sh"
STAIR_ENV="/dcs04/hicks/data/multi-sample-alignment-benchmark/envs/STAIR"
PYTHON_BIN="${STAIR_ENV}/bin/python"
OUTDIR="/dcs04/hicks/data/Jan/sim_project/sim_paper/code/applications_albis/alignment"
LOGDIR="${OUTDIR}/logs_stair_3D"
mkdir -p "${LOGDIR}"

DATASETS=(
  "bin16um"
  "spot"
  "cell"
)

for d in "${DATASETS[@]}"; do
  sbatch \
  "${DEP_ARG[@]}" \
  --job-name="3D_stair" \
  --partition=gpu \
  --gres=gpu:l40s:1 \
  --mem=150G \
  --cpus-per-task=4 \
  --time=48:00:00 \
  --output="${LOGDIR}/stair_${d}_%j.out" \
  --wrap="bash -c 'source ${CONDA_SH} && conda activate ${STAIR_ENV} && ${PYTHON_BIN} ${SCRIPT} --dataset ${d}'"
done
