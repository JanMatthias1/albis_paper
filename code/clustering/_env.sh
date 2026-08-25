# Shared conda-activation logic for the clustering SLURM wrappers.
# Usage: source "$(dirname "${BASH_SOURCE[0]}")/_env.sh"
# Sets PYTHON_BIN to the sim_app tutorial env's python.

ENV_PREFIX="/dcs04/hicks/data/Jan/sim_project/sim_app/env/sim-app-tutorial"
PYTHON_BIN="${ENV_PREFIX}/bin/python"
export PYTHONUNBUFFERED=1

if ! command -v conda >/dev/null 2>&1; then
    module load conda 2>/dev/null || true
fi

if command -v conda >/dev/null 2>&1; then
    CONDA_BASE="$(conda info --base)"
    # shellcheck source=/dev/null
    source "${CONDA_BASE}/etc/profile.d/conda.sh"
    conda activate "${ENV_PREFIX}"
    # Don't overwrite PYTHON_BIN with `command -v python` here: if the job's
    # inherited PATH already has another env's bin/ (e.g. sim-app-banksy)
    # ahead of what `conda activate` prepends -- seen in practice on a busy
    # node, 2026-08-25 -- this silently resolves to the WRONG interpreter
    # even though `conda activate` itself "succeeded". The absolute path
    # set above is unambiguous regardless of PATH ordering.
else
    echo "WARNING: conda was not found; using ${PYTHON_BIN} directly." >&2
    if [[ ! -x "${PYTHON_BIN}" ]]; then
        echo "ERROR: Python not found at ${PYTHON_BIN}" >&2
        echo "Create the environment with /dcs04/hicks/data/Jan/sim_project/sim_app/env/create_tutorial_env.sh first." >&2
        exit 1
    fi
    export PATH="${ENV_PREFIX}/bin:${PATH}"
fi
