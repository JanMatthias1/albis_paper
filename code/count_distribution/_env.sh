# Shared conda-activation logic for the count_distribution SLURM wrapper.
# Usage: source "$(dirname "${BASH_SOURCE[0]}")/_env.sh"
# Sets PYTHON_BIN to the albis tutorial env's python.

ENV_PREFIX="/dcs04/hicks/data/Jan/sim_project/sim_paper/env/albis-tutorial"
PYTHON_BIN="${ENV_PREFIX}/bin/python"

if ! command -v conda >/dev/null 2>&1; then
    module load conda 2>/dev/null || true
fi

if command -v conda >/dev/null 2>&1; then
    CONDA_BASE="$(conda info --base)"
    # shellcheck source=/dev/null
    source "${CONDA_BASE}/etc/profile.d/conda.sh"
    conda activate "${ENV_PREFIX}"
    # Keep the absolute PYTHON_BIN: `command -v python` can resolve to another
    # env's interpreter still earlier on PATH, even after conda activate.
else
    echo "WARNING: conda was not found; using ${PYTHON_BIN} directly." >&2
    if [[ ! -x "${PYTHON_BIN}" ]]; then
        echo "ERROR: Python not found at ${PYTHON_BIN}" >&2
        echo "Create the environment with /dcs04/hicks/data/Jan/sim_project/sim_paper/env/create_tutorial_env.sh first." >&2
        exit 1
    fi
    export PATH="${ENV_PREFIX}/bin:${PATH}"
fi
