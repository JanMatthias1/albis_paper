# Shared conda-activation logic for the BANKSY-matrix SLURM wrapper.
# Usage: source "$(dirname "${BASH_SOURCE[0]}")/_env_banksy.sh"
# Sets PYTHON_BIN to the isolated albis banksy env's python -- see
# sim_paper/env/create_banksy_env.sh for why this is a separate env from the
# rest of the clustering pipeline's albis-tutorial env.

ENV_PREFIX="/dcs04/hicks/data/Jan/sim_project/sim_paper/env/sim-app-banksy"
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
    # Don't overwrite PYTHON_BIN with `command -v python` here -- see
    # _env.sh for why (2026-08-25: this silently resolved to the wrong
    # env's interpreter on a busy node despite `conda activate` succeeding).
else
    echo "WARNING: conda was not found; using ${PYTHON_BIN} directly." >&2
    if [[ ! -x "${PYTHON_BIN}" ]]; then
        echo "ERROR: Python not found at ${PYTHON_BIN}" >&2
        echo "Create the environment with /dcs04/hicks/data/Jan/sim_project/sim_paper/env/create_banksy_env.sh first." >&2
        exit 1
    fi
    export PATH="${ENV_PREFIX}/bin:${PATH}"
fi
