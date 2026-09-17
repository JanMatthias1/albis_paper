#!/bin/bash
#SBATCH --job-name=kgeom_no_harmony
#SBATCH --output=/dcs04/hicks/data/Jan/sim_project/sim_paper/code/clustering/logs/kgeom_no_harmony_%A_%a.out
#SBATCH --array=0-27%8
#SBATCH --time=06:00:00
#SBATCH --mem=250G
#SBATCH --cpus-per-task=8
#SBATCH --partition=shared
#
# BANKSY -> PCA -> Leiden, NO Harmony, for every (modality, lambda, k_geom)
# combo already built by spot_kgeom_sweep_tuned.sh (12 combos) and
# cell_bin_kgeom_sweep_tuned.sh (16 combos, cell + bin16um) -- reuses the
# already-built banksy_matrix h5ads (X_pca_pre_harmony always stored), no
# rebuild. Submitted with --dependency so it waits for both sweeps to finish
# before running.
set -euo pipefail

mkdir -p /dcs04/hicks/data/Jan/sim_project/sim_paper/code/clustering/logs
source /dcs04/hicks/data/Jan/sim_project/sim_paper/code/clustering/_env.sh
cd /dcs04/hicks/data/Jan/sim_project

SPOT_ROOT="sim_paper/data/figure_3/spot_kgeom_sweep_tuned"
CELLBIN_ROOT="sim_paper/data/figure_3/cellbin_kgeom_sweep_tuned"

# 0-11: spot (lambda x k_geom in {0.2,0.5,0.8} x {6,15,30,60})
SPOT_LAMBDAS=(0.2 0.2 0.2 0.2 0.5 0.5 0.5 0.5 0.8 0.8 0.8 0.8)
SPOT_KGEOMS=(6 15 30 60 6 15 30 60 6 15 30 60)
# 12-19: cell, 20-27: bin (lambda x k_geom in {0.2,0.5} x {4,6,10,15})
CB_LAMBDAS=(0.2 0.2 0.2 0.2 0.5 0.5 0.5 0.5)
CB_KGEOMS=(4 6 10 15 4 6 10 15)

i="${SLURM_ARRAY_TASK_ID:-0}"

if [[ "${i}" -lt 12 ]]; then
    MOD="spot"
    LAM="${SPOT_LAMBDAS[$i]}"
    KG="${SPOT_KGEOMS[$i]}"
    INPUT="${SPOT_ROOT}/lambda${LAM}_kgeom${KG}/banksy_matrix/simulation_spot_z_banksy_pca_harmony_qc.h5ad"
    OUT_DIR="${SPOT_ROOT}/lambda${LAM}_kgeom${KG}/ari_no_harmony"
elif [[ "${i}" -lt 20 ]]; then
    j=$(( i - 12 ))
    MOD="cell"
    LAM="${CB_LAMBDAS[$j]}"
    KG="${CB_KGEOMS[$j]}"
    INPUT="${CELLBIN_ROOT}/cell/lambda${LAM}_kgeom${KG}/banksy_matrix/simulation_cell_z_banksy_pca_harmony_qc.h5ad"
    OUT_DIR="${CELLBIN_ROOT}/cell/lambda${LAM}_kgeom${KG}/ari_no_harmony"
else
    j=$(( i - 20 ))
    MOD="bin"
    LAM="${CB_LAMBDAS[$j]}"
    KG="${CB_KGEOMS[$j]}"
    INPUT="${CELLBIN_ROOT}/bin/lambda${LAM}_kgeom${KG}/banksy_matrix/simulation_bin_z_banksy_pca_harmony_qc.h5ad"
    OUT_DIR="${CELLBIN_ROOT}/bin/lambda${LAM}_kgeom${KG}/ari_no_harmony"
fi

echo "[task ${i}] modality=${MOD} lambda=${LAM} k_geom=${KG}"
if [[ ! -f "${INPUT}" ]]; then
    echo "ERROR: input not found: ${INPUT} -- upstream sweep task must have failed/not run"
    exit 1
fi

"${PYTHON_BIN}" sim_paper/code/clustering/misc/check_ari_no_harmony_generic.py \
    --input "${INPUT}" --output-dir "${OUT_DIR}" --modality "${MOD}"

echo "[done] task ${i}"
