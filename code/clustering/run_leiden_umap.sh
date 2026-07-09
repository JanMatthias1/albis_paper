#!/bin/bash
#SBATCH --job-name=sim_leiden_umap
#SBATCH --output=/dcs04/hicks/data/Jan/sim_project/sim_paper/code/clustering/logs/sim_leiden_umap_%A_%a.out
#SBATCH --array=0-2
#SBATCH --time=02:00:00
#SBATCH --mem=32G
#SBATCH --cpus-per-task=4
#SBATCH --partition=shared      # uncomment / edit to match your cluster's partition name

set -euo pipefail

mkdir -p /dcs04/hicks/data/Jan/sim_project/sim_paper/code/clustering/logs

ENV_PREFIX="/dcs04/hicks/data/Jan/sim_project/sim_app/env/sim-app-tutorial"
PYTHON_BIN="${ENV_PREFIX}/bin/python"

if ! command -v conda >/dev/null 2>&1; then
    module load conda 2>/dev/null || true
fi

if command -v conda >/dev/null 2>&1; then
    CONDA_BASE="$(conda info --base)"
    # shellcheck source=/dev/null
    source "${CONDA_BASE}/etc/profile.d/conda.sh"
    conda activate "${ENV_PREFIX}"
    PYTHON_BIN="$(command -v python)"
else
    echo "WARNING: conda was not found; using ${PYTHON_BIN} directly." >&2
    if [[ ! -x "${PYTHON_BIN}" ]]; then
        echo "ERROR: Python not found at ${PYTHON_BIN}" >&2
        echo "Create the environment with /dcs04/hicks/data/Jan/sim_project/sim_app/env/create_tutorial_env.sh first." >&2
        exit 1
    fi
    export PATH="${ENV_PREFIX}/bin:${PATH}"
fi

cd /dcs04/hicks/data/Jan/sim_project

echo "Job started: $(date)"
echo "Host: $(hostname)"
echo "Python: ${PYTHON_BIN}"

MODALITIES=(spot bin cell)
MODALITY="${MODALITIES[${SLURM_ARRAY_TASK_ID:-0}]}"

if [[ "$#" -gt 0 ]]; then
    echo "Modality: ${MODALITY}"
    echo "Command: ${PYTHON_BIN} sim_paper/code/clustering/leiden_umap.py --modality ${MODALITY} $*"
    "${PYTHON_BIN}" sim_paper/code/clustering/leiden_umap.py --modality "${MODALITY}" "$@"
else
    RESOLUTIONS=(0.5 0.1 0.2)
    echo "Modality: ${MODALITY}"
    for resolution in "${RESOLUTIONS[@]}"; do
        echo "Command: ${PYTHON_BIN} sim_paper/code/clustering/leiden_umap.py --modality ${MODALITY} --resolution ${resolution}"
        "${PYTHON_BIN}" sim_paper/code/clustering/leiden_umap.py --modality "${MODALITY}" --resolution "${resolution}"
    done
fi

echo "Job finished: $(date)"
