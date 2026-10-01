#!/usr/bin/env bash
# Figure 4C -- run STAIR 3D alignment on the four ALBIS z-stacks.
# Optional $1: a job id to gate all submissions on
# (--dependency=afterok:<id>), e.g. the QC-regen job. "none" / unset = no dep.
# Optional $2+: subset of datasets to submit (default: all of DATASETS).
set -euo pipefail

DEP_JOB="${1:-none}"
DEP_ARG=()
[[ "${DEP_JOB}" != "none" ]] && DEP_ARG=(--dependency="afterok:${DEP_JOB}")
shift || true

SCRIPT="/dcs04/hicks/data/Jan/sim_project/sim_paper/code/applications_albis/alignment/3D_stair.py"
CONDA_SH="/jhpce/shared/jhpce/core/anaconda3/2023.03/etc/profile.d/conda.sh"
STAIR_ENV="/dcs04/hicks/data/multi-sample-alignment-benchmark/envs/STAIR"
PYTHON_BIN="${STAIR_ENV}/bin/python"
OUTDIR="/dcs04/hicks/data/Jan/sim_project/sim_paper/code/applications_albis/alignment"
LOGDIR="${OUTDIR}/logs_stair_3D"
mkdir -p "${LOGDIR}"

# Figure 4A strong-mix data (see 3D_stair.py DATASETS).
DATASETS=(
  "bin16um"
  "spot"
  "cell"
)
[[ $# -gt 0 ]] && DATASETS=("$@")

# Memory/time per dataset; larger stacks (bin8um, 600k+ cells) get more memory.
for d in "${DATASETS[@]}"; do
  MEM="150G"
  TIME="48:00:00"
  [[ "${d}" == "bin8um" ]] && MEM="250G" && TIME="24:00:00"
  [[ "${d}" == "cell_r6000" ]] && MEM="200G" && TIME="24:00:00"
  [[ "${d}" == "cell" ]] && MEM="200G" && TIME="24:00:00"

  sbatch \
  "${DEP_ARG[@]}" \
  --job-name="3D_stair" \
  --partition=gpu \
  --gres=gpu:l40s:1 \
  --mem="${MEM}" \
  --cpus-per-task=4 \
  --time="${TIME}" \
  --output="${LOGDIR}/stair_${d}_%j.out" \
  --wrap="bash -c 'source ${CONDA_SH} && conda activate ${STAIR_ENV} && ${PYTHON_BIN} ${SCRIPT} --dataset ${d}'"
done
