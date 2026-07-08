#!/bin/bash
#SBATCH --job-name=sim_app_generate
#SBATCH --output=/dcs04/hicks/data/Jan/sim_project/sim_paper/code/data/logs/sim_app_generate_%j.out
#SBATCH --time=04:00:00
#SBATCH --mem=32G
#SBATCH --cpus-per-task=4
##SBATCH --partition=shared      # uncomment / edit to match your cluster's partition name


set -euo pipefail

# Make sure the log directory exists (SLURM needs it to exist before the job starts writing)
mkdir -p /dcs04/hicks/data/Jan/sim_project/sim_paper/code/data/logs

ENV_PREFIX="/dcs04/hicks/data/Jan/sim_project/sim_app/env/sim-app-tutorial"
PYTHON_BIN="${ENV_PREFIX}/bin/python"

# Activate the conda environment used by sim_app. Batch shells often do not
# source .bashrc, so conda activate is unavailable until conda.sh is loaded.
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

cd /dcs04/hicks/data/Jan/sim_project/sim_paper/code/data

echo "Job started: $(date)"
echo "Host: $(hostname)"
echo "Python: ${PYTHON_BIN}"

if [[ "${1:-}" == "--check-env" ]]; then
    "${PYTHON_BIN}" -c "from importlib.metadata import version; print('sim-app installed:', version('sim-app'))"
    echo "Environment check passed."
    exit 0
fi

"${PYTHON_BIN}" generate_simulation.py

echo "Job finished: $(date)"
