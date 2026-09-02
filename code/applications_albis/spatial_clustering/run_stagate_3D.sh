#!/usr/bin/env bash
# Figure 4B -- run STAGATE 3D spatial-domain identification on the ALBIS
# z-stacks: the weak-mix Figure 2 tags (bin16um / spot / cell) AND the
# strong-domain-mix tags (bin8 / spot / cell) used by the BANKSY domain panels.
# Optional $1: a job id to gate every submission on (--dependency=afterok:<id>),
#   e.g. a QC-regen job. "none" / unset = no dependency.
# Optional $2..$N: explicit dataset names to run instead of the full set.
set -euo pipefail

DEP_JOB="${1:-none}"
DEP_ARG=()
[[ "${DEP_JOB}" != "none" ]] && DEP_ARG=(--dependency="afterok:${DEP_JOB}")

# STAGATE_pyG trains full-batch: the whole 3D SNN graph lives on the GPU at once.
# bin16um (~366k nodes / 7.8M edges) and cell (~600k nodes) OOM on the 46GB L40S
# (peak ~48-50GB); only spot fits. Default to the 80GB A100 (7 on the gpu
# partition, ~30GB headroom); override for the 96GB H100 (only 2, longer queue):
#   GPU_GRES=gpu:tesh100:1 bash run_stagate_3D.sh ...
GPU_GRES="${GPU_GRES:-gpu:tesa100:1}"

SCRIPT="/dcs04/hicks/data/Jan/sim_project/sim_paper/code/applications_albis/spatial_clustering/3D_stagate.py"
CONDA_SH="/jhpce/shared/jhpce/core/anaconda3/2023.03/etc/profile.d/conda.sh"
STAGATE_ENV="/dcs04/hicks/data/Jan/sim_project/sim_paper/env/stagate-pyg"
PYTHON_BIN="${STAGATE_ENV}/bin/python"
OUTDIR="/dcs04/hicks/data/Jan/sim_project/sim_paper/code/applications_albis/spatial_clustering"
LOGDIR="${OUTDIR}/logs_stagate_3D"
mkdir -p "${LOGDIR}"

# weak-mix Figure 2 tags + strong-domain-mix tags (see 3D_stagate.py DATASETS).
# Override by passing dataset names as $2..$N, e.g.:
#   bash run_stagate_3D.sh none spot_strongmix cell_strongmix
DATASETS=("${@:2}")
if [[ ${#DATASETS[@]} -eq 0 ]]; then
  DATASETS=(
    # family 1: weak mix, tuned batch (canonical Figure 2, reused in place)
    "bin16um" "spot" "cell"
    # family 2: strong domain mix, tuned batch
    "bin16um_strongmix" "spot_strongmix" "cell_strongmix"
    # family 3: weak mix, very low batch (0.05)
    "bin16um_lowbatch" "spot_lowbatch" "cell_lowbatch"
    # family 4: strong mix, very low batch (0.05)
    "bin16um_strongmix_lowbatch" "spot_strongmix_lowbatch" "cell_strongmix_lowbatch"
  )
fi

for d in "${DATASETS[@]}"; do
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
