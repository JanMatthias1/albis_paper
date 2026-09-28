#!/usr/bin/env bash
set -euo pipefail

# =============================================================================
# create_stagate_env.sh
#
# Creates an ISOLATED conda environment for STAGATE 3D spatial-domain
# identification (Figure 4B, sim_paper/code/applications_albis/spatial_clustering/).
#
# We install the PyTorch-Geometric port (STAGATE_pyG), NOT the original
# TensorFlow-1.15 STAGATE:
#   * TF 1.15 caps at CUDA 10 / compute capability <= 7.5. The GPU partition
#     here is L40S (compute capability 8.9) -- TF 1.15 cannot use it.
#   * The Figure 2 z-stacks are large (bin16um ~366k bins, cell ~600k cells
#     after QC). STAGATE needs a GPU at that scale; CPU TF 1.15 is not viable.
#   * STAGATE_pyG exposes the same API used by the 3D tutorial
#     (Cal_Spatial_Net_3D / train_STAGATE / mclust_R). The only unsupported
#     piece is the cell-type-aware module (alpha>0), which the 3D tutorial
#     does not use (it calls train_STAGATE(alpha=0)).
#
# Torch stack is pinned to match the STAIR alignment env
# (torch 2.6.0+cu124) so both Figure 4 application envs behave the same on
# the L40S nodes.
#
# mclust clustering needs R + the mclust package, reached through rpy2. We
# install r-base + r-mclust from conda-forge into this env, then build rpy2
# against that R.
#
# Usage:
#   bash create_stagate_env.sh
#
# Override defaults with env vars:
#   CONDA_ENV_NAME=my-env PYTHON_VERSION=3.10 STAGATE_REF=main bash create_stagate_env.sh
# =============================================================================

if [[ -n "${SLURM_SUBMIT_DIR:-}" ]]; then
    REPO_ROOT="${SLURM_SUBMIT_DIR}"
else
    REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
fi

CONDA_ENV_NAME="${CONDA_ENV_NAME:-stagate-pyg}"
PYTHON_VERSION="${PYTHON_VERSION:-3.10}"
ENV_PREFIX="${REPO_ROOT}/env/${CONDA_ENV_NAME}"

TORCH_VERSION="${TORCH_VERSION:-2.6.0}"
TORCH_CUDA="${TORCH_CUDA:-cu124}"
STAGATE_REPO="${STAGATE_REPO:-https://github.com/QIFEIDKN/STAGATE_pyG.git}"
STAGATE_REF="${STAGATE_REF:-main}"
STAGATE_SRC="${REPO_ROOT}/env/STAGATE_pyG"

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
    echo "[setup] Creating conda environment at ${ENV_PREFIX} (python=${PYTHON_VERSION})"
    conda create -y --prefix "${ENV_PREFIX}" "python=${PYTHON_VERSION}" pip
fi

# IMPORTANT: never rely on `conda activate` + bare `python`/`pip` here. On a
# busy login node the previously-active env's bin dir can stay ahead of what
# `conda activate` prepends to PATH, so bare `pip` installs into the WRONG env
# (the project's own code/*/_env.sh files carry the same warning). Drive every
# install through the absolute env python and `conda ... --prefix`.
PY="${ENV_PREFIX}/bin/python"

# R + mclust + rpy2 ALL from conda-forge, in one solve: conda-forge ships rpy2
# prebuilt against its own r-base, so the compiled API extension matches the
# runtime R. (pip-installing rpy2 here instead builds its cffi extension
# against whatever R is on PATH at build time -- which is none, since we never
# `conda activate` -- producing a broken `undefined symbol: R_ClosureEnv` .so.)
echo "[setup] Installing R + mclust + rpy2 from conda-forge (needed by STAGATE.mclust_R)"
"${PY}" -m pip uninstall -y rpy2 rpy2-rinterface rpy2-robjects 2>/dev/null || true
conda install -y --prefix "${ENV_PREFIX}" -c conda-forge r-mclust rpy2

# rpy2 dlopens libR.so, which pulls the conda libicuuc -> needs the conda
# libstdc++ (GLIBCXX_3.4.30), not the older system /lib64 one. Put the env's
# lib dir first, and point rpy2 straight at the conda R so it doesn't need
# `R` on PATH.
export LD_LIBRARY_PATH="${ENV_PREFIX}/lib:${LD_LIBRARY_PATH:-}"
export R_HOME="${ENV_PREFIX}/lib/R"
export PATH="${ENV_PREFIX}/bin:${PATH}"

echo "[setup] Upgrading pip/setuptools/wheel"
"${PY}" -m pip install --upgrade pip setuptools wheel

echo "[setup] Installing PyTorch ${TORCH_VERSION}+${TORCH_CUDA}"
"${PY}" -m pip install "torch==${TORCH_VERSION}" \
    --index-url "https://download.pytorch.org/whl/${TORCH_CUDA}"

echo "[setup] Installing torch-geometric (+ prebuilt scatter/sparse/pyg-lib wheels for this torch build)"
"${PY}" -m pip install torch-geometric
"${PY}" -m pip install pyg-lib torch-scatter torch-sparse \
    -f "https://data.pyg.org/whl/torch-${TORCH_VERSION}+${TORCH_CUDA}.html" || \
    echo "[warn] optional pyg extension wheels not installed — full-batch STAGATE still works, NeighborLoader batching will not"

echo "[setup] Installing scientific stack (scanpy / anndata / sklearn); rpy2 already came from conda"
"${PY}" -m pip install "scanpy>=1.9" anndata scikit-learn scikit-misc leidenalg python-igraph

echo "[setup] Fetching STAGATE_pyG (${STAGATE_REPO} @ ${STAGATE_REF})"
if [[ -d "${STAGATE_SRC}/.git" ]]; then
    git -C "${STAGATE_SRC}" fetch --depth 1 origin "${STAGATE_REF}"
    git -C "${STAGATE_SRC}" checkout -q "${STAGATE_REF}"
    git -C "${STAGATE_SRC}" reset --hard -q "origin/${STAGATE_REF}" || true
else
    git clone --depth 1 --branch "${STAGATE_REF}" "${STAGATE_REPO}" "${STAGATE_SRC}"
fi

echo "[setup] Installing STAGATE_pyG"
"${PY}" -m pip install "${STAGATE_SRC}"

echo "[setup] Verifying imports + required 3D API"
"${PY}" -c "
import torch, torch_geometric, scanpy as sc, anndata as ad, numpy as np
import STAGATE_pyG as ST
print('torch', torch.__version__, '| CUDA build', torch.version.cuda)
print('torch_geometric', torch_geometric.__version__)
print('scanpy', sc.__version__, '| anndata', ad.__version__, '| numpy', np.__version__)
for fn in ('Cal_Spatial_Net', 'Cal_Spatial_Net_3D', 'train_STAGATE', 'mclust_R'):
    assert hasattr(ST, fn), f'STAGATE_pyG is missing {fn} — 3D tutorial cannot run'
    print('STAGATE_pyG.' + fn, 'OK')
import rpy2.robjects as ro
ro.r('library(mclust)')
print('rpy2 + R mclust OK')
"

echo ""
echo "[setup] Done. Environment created at: ${ENV_PREFIX}"
echo "STAGATE_pyG source checked out at: ${STAGATE_SRC}"
echo ""
echo "Activate with:"
echo "  conda activate ${ENV_PREFIX}"
echo ""
echo "Run Figure 4B with:"
echo "  bash sim_paper/code/applications_albis/spatial_clustering/run_stagate_3D.sh"
