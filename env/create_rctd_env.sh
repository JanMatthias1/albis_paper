#!/usr/bin/env bash
set -euo pipefail

# =============================================================================
# create_rctd_env.sh
#
# Creates an ISOLATED conda environment for RCTD (Robust Cell Type
# Decomposition, https://github.com/dmcable/spacexr) -- spatial deconvolution
# of Visium spots (albis_paper/code/clustering/spatial_deconvolution/).
#
# RCTD ("spacexr") is GitHub-only, not on CRAN/Bioconductor/conda-forge, so it
# is installed via remotes::install_github inside a conda-forge R. The h5ad
# inputs are read directly in R via the "anndata" R package (CRAN, wraps
# Python's anndata through reticulate) -- a Python + anndata are installed
# into this SAME env so reticulate has a colocated interpreter, rather than
# depending on a second, separately-activated conda env (same
# self-contained-single-env approach as create_stagate_env.sh).
#
# Usage:
#   bash create_rctd_env.sh
#
# Override defaults with env vars:
#   CONDA_ENV_NAME=my-env bash create_rctd_env.sh
# =============================================================================

if [[ -n "${SLURM_SUBMIT_DIR:-}" ]]; then
    REPO_ROOT="${SLURM_SUBMIT_DIR}"
else
    REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
fi

CONDA_ENV_NAME="${CONDA_ENV_NAME:-rctd}"
ENV_PREFIX="${REPO_ROOT}/env/${CONDA_ENV_NAME}"

if ! command -v conda >/dev/null 2>&1; then
    echo "ERROR: conda not found. Run: module load conda (or your cluster's conda module)." >&2
    exit 1
fi

CONDA_BASE="$(conda info --base)"
# shellcheck disable=SC1091
source "${CONDA_BASE}/etc/profile.d/conda.sh"

if [[ -d "${ENV_PREFIX}" ]]; then
    echo "[setup] Environment already exists at ${ENV_PREFIX} — reusing (installs are idempotent)."
else
    echo "[setup] Creating conda environment at ${ENV_PREFIX}"
    conda create -y --prefix "${ENV_PREFIX}" python=3.10 pip
fi

# Same warning as every other _env.sh / create_*_env.sh in this project: never
# rely on `conda activate` + bare `R`/`python`/`pip` on a busy login node --
# drive every install through the absolute env paths.
PY="${ENV_PREFIX}/bin/python"

echo "[setup] Installing R 4.3 + core deps (Matrix, doParallel, pbapply, dplyr,"
echo "        ggplot2, remotes for install_github) + anndata (R + Python) from conda-forge"
conda install -y --prefix "${ENV_PREFIX}" -c conda-forge \
    "r-base=4.3*" r-matrix r-doparallel r-pbapply r-dplyr r-ggplot2 \
    r-remotes r-reticulate r-anndata r-jsonlite anndata

export LD_LIBRARY_PATH="${ENV_PREFIX}/lib:${LD_LIBRARY_PATH:-}"
export R_HOME="${ENV_PREFIX}/lib/R"
export PATH="${ENV_PREFIX}/bin:${PATH}"
R_BIN="${ENV_PREFIX}/bin/R"
RSCRIPT_BIN="${ENV_PREFIX}/bin/Rscript"

# Point the R "anndata" package's reticulate at this same env's python (which
# has Python anndata installed above) rather than letting reticulate search
# and possibly find/create a different one.
export RETICULATE_PYTHON="${PY}"

echo "[setup] Installing RCTD (spacexr) from GitHub via remotes::install_github"
"${RSCRIPT_BIN}" -e '
Sys.setenv(RETICULATE_PYTHON = Sys.getenv("RETICULATE_PYTHON"))
remotes::install_github("dmcable/spacexr", upgrade = "never", quiet = FALSE)
'

echo "[setup] Verifying imports"
"${RSCRIPT_BIN}" -e '
library(spacexr)
library(anndata)
library(Matrix)
cat("spacexr", as.character(packageVersion("spacexr")), "OK\n")
cat("anndata (R)", as.character(packageVersion("anndata")), "OK\n")
'

echo ""
echo "[setup] Done. Environment created at: ${ENV_PREFIX}"
echo ""
echo "Run RCTD spot deconvolution with:"
echo "  sbatch albis_paper/code/clustering/spatial_deconvolution/run_rctd_spot.sh"
