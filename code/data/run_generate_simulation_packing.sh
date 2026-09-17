#!/bin/bash
#SBATCH --job-name=albis_generate_packing
#SBATCH --output=/dcs04/hicks/data/Jan/sim_project/sim_paper/code/data/logs/albis_generate_packing_%A_%a.out
#SBATCH --array=0-2
#SBATCH --time=04:00:00
#SBATCH --mem=250G
#SBATCH --cpus-per-task=4
##SBATCH --partition=shared      # uncomment / edit to match your cluster's partition name


set -euo pipefail

# Make sure the log directory exists (SLURM needs it to exist before the job starts writing)
mkdir -p /dcs04/hicks/data/Jan/sim_project/sim_paper/code/data/logs

ENV_PREFIX="/dcs04/hicks/data/Jan/sim_project/albis/env/albis-tutorial"
PYTHON_BIN="${ENV_PREFIX}/bin/python"

# Activate the conda environment used by albis. Batch shells often do not
# source .bashrc, so conda activate is unavailable until conda.sh is loaded.
if ! command -v conda >/dev/null 2>&1; then
    module load conda 2>/dev/null || true
fi

if command -v conda >/dev/null 2>&1; then
    CONDA_BASE="$(conda info --base)"
    # shellcheck source=/dev/null
    source "${CONDA_BASE}/etc/profile.d/conda.sh"
    conda activate "${ENV_PREFIX}"
    # Don't overwrite PYTHON_BIN with `command -v python` here -- on a busy
    # node this can silently resolve to a DIFFERENT env's interpreter still
    # ahead on PATH, even though `conda activate` itself "succeeded"
    # (2026-08-25). The absolute path set above is unambiguous.
else
    echo "WARNING: conda was not found; using ${PYTHON_BIN} directly." >&2
    if [[ ! -x "${PYTHON_BIN}" ]]; then
        echo "ERROR: Python not found at ${PYTHON_BIN}" >&2
        echo "Create the environment with /dcs04/hicks/data/Jan/sim_project/albis/env/create_tutorial_env.sh first." >&2
        exit 1
    fi
    export PATH="${ENV_PREFIX}/bin:${PATH}"
fi

cd /dcs04/hicks/data/Jan/sim_project/sim_paper/code/data

echo "Job started: $(date)"
echo "Host: $(hostname)"
echo "Python: ${PYTHON_BIN}"

if [[ "${1:-}" == "--check-env" ]]; then
    "${PYTHON_BIN}" -c "from importlib.metadata import version; print('sim-app installed:', version('sim-app'))"
    echo "Environment check passed."
    exit 0
fi

MODALITIES=(spot bin cell)

# --modality is not accepted positionally here (unlike batch_effect_high's
# wrapper) because every other flag (--sphere-r-um, --out-tag, ...) is a
# real --flag forwarded via "$@" -- a bare positional modality would be
# ambiguous with those. Pick the modality via SLURM array index instead
# (sbatch --array=1 ... for bin only, etc.), same as run_generate_simulation_noisy.sh.
MODALITY="${MODALITIES[${SLURM_ARRAY_TASK_ID:-0}]}"

echo "Modality: ${MODALITY}"
echo "Command: ${PYTHON_BIN} generate_simulation_packing.py --modality ${MODALITY} $*"
"${PYTHON_BIN}" generate_simulation_packing.py --modality "${MODALITY}" "$@"

echo "Job finished: $(date)"
