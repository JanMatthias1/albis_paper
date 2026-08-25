# Shared conda-activation logic for the count_distribution SLURM wrapper.
# Usage: source "$(dirname "${BASH_SOURCE[0]}")/_env.sh"
# Sets PYTHON_BIN to the sim_app tutorial env's python.

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
    # Don't overwrite PYTHON_BIN with `command -v python` here -- on a busy
    # node this can silently resolve to a DIFFERENT env's interpreter still
    # ahead on PATH, even though `conda activate` itself "succeeded"
    # (2026-08-25). The absolute path set above is unambiguous.
else
    echo "WARNING: conda was not found; using ${PYTHON_BIN} directly." >&2
    if [[ ! -x "${PYTHON_BIN}" ]]; then
        echo "ERROR: Python not found at ${PYTHON_BIN}" >&2
        echo "Create the environment with /dcs04/hicks/data/Jan/sim_project/sim_app/env/create_tutorial_env.sh first." >&2
        exit 1
    fi
    export PATH="${ENV_PREFIX}/bin:${PATH}"
fi
