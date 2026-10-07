#!/usr/bin/env bash
set -euo pipefail

# =============================================================================
# create_banksy_env.sh
#
# Creates an ISOLATED conda environment for the BANKSY clustering step
# (albis_paper/code/clustering/01_build_banksy_matrix.py).
#
# banksy_py pins a much older scanpy/numpy/anndata/pandas/scikit-learn/scipy/
# matplotlib stack than the rest of this project uses (the albis-tutorial
# env) -- installing it there would downgrade every other script in this
# repo. This env exists to keep that pinned stack fully isolated.
#
# Nothing downstream of the BANKSY step needs banksy_py itself -- only the
# resulting .h5ad file, which PCA/Harmony/Leiden (and everything else) reads
# back in the normal albis-tutorial env. So the split is: run
# 01_build_banksy_matrix.py in THIS env, then switch back to
# albis-tutorial for every step after.
#
# Usage:
#   bash create_banksy_env.sh
#
# Override defaults with env vars:
#   CONDA_ENV_NAME=my-env PYTHON_VERSION=3.10 bash create_banksy_env.sh
# =============================================================================

if [[ -n "${SLURM_SUBMIT_DIR:-}" ]]; then
    REPO_ROOT="${SLURM_SUBMIT_DIR}"
else
    REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
fi

CONDA_ENV_NAME="${CONDA_ENV_NAME:-albis-banksy}"
PYTHON_VERSION="${PYTHON_VERSION:-3.10}"
ENV_PREFIX="${REPO_ROOT}/env/${CONDA_ENV_NAME}"

if ! command -v conda >/dev/null 2>&1; then
    echo "ERROR: conda not found. Run: module load conda (or your cluster's conda module)." >&2
    exit 1
fi

CONDA_BASE="$(conda info --base)"
# shellcheck disable=SC1091
source "${CONDA_BASE}/etc/profile.d/conda.sh"

if [[ -d "${ENV_PREFIX}" ]]; then
    echo "[setup] Environment already exists at ${ENV_PREFIX} — skipping creation."
else
    echo "[setup] Creating conda environment at ${ENV_PREFIX} (python=${PYTHON_VERSION})"
    conda create -y --prefix "${ENV_PREFIX}" "python=${PYTHON_VERSION}" pip
fi

conda activate "${ENV_PREFIX}"

echo "[setup] Upgrading pip/setuptools/wheel"
python -m pip install --upgrade pip setuptools wheel

echo "[setup] Installing banksy_py (pins scanpy==1.9.1, numpy==1.22.4, anndata==0.8.0,"
echo "        pandas==1.5.2, scikit-learn==1.2.0, scipy==1.10.0, matplotlib==3.6.2 --"
echo "        isolated in this env on purpose, see header comment)"
python -m pip install banksy_py

echo "[setup] Installing harmonypy for batch correction (matches"
echo "        albis_paper/code/clustering/01.2_pca_harmony.py's Harmony implementation, so"
echo "        results are consistent with the rest of the clustering pipeline)"
python -m pip install harmonypy

echo "[setup] Installing leidenalg (banksy.main imports it directly but banksy_py"
echo "        doesn't declare it as a dependency)"
python -m pip install leidenalg

echo "[setup] Verifying imports"
python -c "
import banksy
import banksy_utils
import harmonypy
import leidenalg
import scanpy as sc
import anndata as ad
import numpy as np
print('banksy OK')
print('banksy_utils OK')
print('harmonypy OK')
print('leidenalg OK')
print(f'scanpy {sc.__version__}')
print(f'anndata {ad.__version__}')
print(f'numpy {np.__version__}')
"

echo ""
echo "[setup] Done. Environment created at: ${ENV_PREFIX}"
echo ""
echo "Activate with:"
echo "  conda activate ${ENV_PREFIX}"
echo ""
echo "This env is for 01_build_banksy_matrix.py only. Switch back to"
echo "albis-tutorial for every step after that (PCA/Harmony/Leiden, etc.)."
