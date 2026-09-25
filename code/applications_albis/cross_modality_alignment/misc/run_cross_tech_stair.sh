#!/usr/bin/env bash
# Figure 4E -- cross-technology STAIR alignment: bin16um + spot + cell of the
# SAME ALBIS slice_id, aligned as if they were three "slices" of one z-stack.
# Usage: ./run_cross_tech_stair.sh [slice_id ...]   (default: slice 5)
set -euo pipefail

SCRIPT="/dcs04/hicks/data/Jan/sim_project/sim_paper/code/applications_albis/cross_modality_alignment/cross_tech_stair.py"
CONDA_SH="/jhpce/shared/jhpce/core/anaconda3/2023.03/etc/profile.d/conda.sh"
STAIR_ENV="/dcs04/hicks/data/multi-sample-alignment-benchmark/envs/STAIR"
PYTHON_BIN="${STAIR_ENV}/bin/python"
OUTDIR="/dcs04/hicks/data/Jan/sim_project/sim_paper/code/applications_albis/cross_modality_alignment"
LOGDIR="${OUTDIR}/logs_cross_tech_stair"
mkdir -p "${LOGDIR}"

SLICES=("$@")
[[ ${#SLICES[@]} -eq 0 ]] && SLICES=(5)

for s in "${SLICES[@]}"; do
  sbatch \
  --job-name="cross_tech_stair" \
  --partition=gpu \
  --gres=gpu:l40s:1 \
  --mem=100G \
  --cpus-per-task=4 \
  --time=08:00:00 \
  --output="${LOGDIR}/cross_tech_slice${s}_%j.out" \
  --wrap="bash -c 'source ${CONDA_SH} && conda activate ${STAIR_ENV} && ${PYTHON_BIN} ${SCRIPT} --slice ${s}'"
done
