#!/bin/bash
#SBATCH --job-name=sphere_intact
#SBATCH --output=/dcs04/hicks/data/Jan/sim_project/albis_paper/code/sphere_figure/logs/sphere_intact_%j.out
#SBATCH --time=04:00:00
#SBATCH --mem=250G
#SBATCH --cpus-per-task=4
##SBATCH --partition=shared      # uncomment / edit to match your cluster's partition name

set -euo pipefail

SCRIPT_DIR="/dcs04/hicks/data/Jan/sim_project/albis_paper/code/sphere_figure"
mkdir -p "$SCRIPT_DIR/logs"
mkdir -p "$SCRIPT_DIR/.mplconfig"
export MPLCONFIGDIR="$SCRIPT_DIR/.mplconfig"

ENV_PREFIX="/dcs04/hicks/data/Jan/sim_project/albis/env/albis-tutorial"
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
    # ahead on PATH, even though `conda activate` itself "succeeded".
    # The absolute path set above is unambiguous.
else
    echo "WARNING: conda was not found; using ${PYTHON_BIN} directly." >&2
    if [[ ! -x "${PYTHON_BIN}" ]]; then
        echo "ERROR: Python not found at ${PYTHON_BIN}" >&2
        echo "Create the environment with /dcs04/hicks/data/Jan/sim_project/albis/env/create_tutorial_env.sh first." >&2
        exit 1
    fi
    export PATH="${ENV_PREFIX}/bin:${PATH}"
fi

cd "$SCRIPT_DIR"

echo "Job started: $(date)"
echo "Host: $(hostname)"
echo "Python: ${PYTHON_BIN}"

"${PYTHON_BIN}" "$SCRIPT_DIR/01_plot_intact_sphere.py" "$@"

echo "Job finished: $(date)"
