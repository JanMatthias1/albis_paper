#!/usr/bin/env bash
# Figure 4B -- run STAGATE 3D spatial-domain identification on the ALBIS
# z-stacks: the Figure 3 strong-mix grid shared with BANKSY.
# Optional $1: a job id to gate every submission on (--dependency=afterok:<id>),
#   e.g. a QC-regen job. "none" / unset = no dependency.
# Optional $2..$N: explicit dataset names to run instead of the full set.
set -euo pipefail

DEP_JOB="${1:-none}"
DEP_ARG=()
[[ "${DEP_JOB}" != "none" ]] && DEP_ARG=(--dependency="afterok:${DEP_JOB}")

# Full-batch bin16um exceeds 80GB A100 memory.
# Use the 96GB H100 for every bin16um condition/seed, A100 for cell/spot.
# GPU_GRES remains an explicit override for all submitted datasets.
GPU_OVERRIDE="${GPU_GRES:-}"

SCRIPT="/dcs04/hicks/data/Jan/sim_project/sim_paper/code/applications_albis/spatial_clustering/3D_stagate.py"
CONDA_SH="/jhpce/shared/jhpce/core/anaconda3/2023.03/etc/profile.d/conda.sh"
STAGATE_ENV="/dcs04/hicks/data/Jan/sim_project/sim_paper/env/stagate-pyg"
PYTHON_BIN="${STAGATE_ENV}/bin/python"
OUTDIR="/dcs04/hicks/data/Jan/sim_project/sim_paper/code/applications_albis/spatial_clustering"
LOGDIR="${OUTDIR}/logs_stagate_3D"
mkdir -p "${LOGDIR}"

# Default: all 27 Figure 3 strong-mix inputs. Use the grid submitter for
# per-seed generation dependencies when inputs are not ready yet.
DATASETS=("${@:2}")
if [[ ${#DATASETS[@]} -eq 0 ]]; then
  mapfile -t DATASETS < <("${PYTHON_BIN}" -c "import sys; sys.path.insert(0, '${OUTDIR}'); from stagate_inputs import DATASETS; print('\n'.join(DATASETS))")
fi

for d in "${DATASETS[@]}"; do
  GPU_GRES="${GPU_OVERRIDE:-gpu:tesa100:1}"
  if [[ -z "${GPU_OVERRIDE}" && "${d}" == bin16um_* ]]; then
    GPU_GRES="gpu:tesh100:1"
  fi
  sbatch \
  "${DEP_ARG[@]}" \
  --job-name="3D_stagate" \
  --partition=gpu \
  --gres="${GPU_GRES}" \
  --mem=150G \
  --cpus-per-task=4 \
  --time=06:00:00 \
  --output="${LOGDIR}/stagate_${d}_%j.out" \
  --wrap="bash -c 'source ${CONDA_SH} && conda activate ${STAGATE_ENV} && export LD_LIBRARY_PATH=${STAGATE_ENV}/lib:\${LD_LIBRARY_PATH:-} R_HOME=${STAGATE_ENV}/lib/R PYTORCH_CUDA_ALLOC_CONF=expandable_segments:True && ${PYTHON_BIN} ${SCRIPT} --dataset ${d}'"
done
