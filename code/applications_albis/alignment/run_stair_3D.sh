#!/usr/bin/env bash
# Figure 4C -- run STAIR 3D alignment on the three ALBIS z-stacks.
# Usage:  ./run_stair_3D.sh            # submit all three
#         ./run_stair_3D.sh spot cell  # submit a subset
set -euo pipefail

SCRIPT_DIR="/dcs04/hicks/data/Jan/sim_project/sim_paper/code/applications_albis/alignment"
SCRIPT="${SCRIPT_DIR}/3D_stair.py"
LOGDIR="${SCRIPT_DIR}/logs_stair_3D"
mkdir -p "${LOGDIR}"

DATASETS=("${@:-bin16um spot cell}")
# allow a single quoted default arg to expand into words
read -r -a DATASETS <<< "${DATASETS[*]}"

for d in "${DATASETS[@]}"; do
  case "$d" in
    bin16um) MEM=64G  ; TIME=12:00:00 ;;
    spot)    MEM=32G  ; TIME=04:00:00 ;;
    cell)    MEM=150G ; TIME=24:00:00 ;;
    *) echo "unknown dataset: $d" >&2; exit 1 ;;
  esac

  sbatch \
    --job-name="f4c_stair_${d}" \
    --partition=gpu \
    --gres=gpu:l40s:1 \
    --mem="${MEM}" \
    --cpus-per-task=4 \
    --time="${TIME}" \
    --output="${LOGDIR}/stair_${d}_%j.out" \
    --wrap="bash -lc 'source ~/.bashrc && conda activate STAIR && python ${SCRIPT} --dataset ${d}'"
done
