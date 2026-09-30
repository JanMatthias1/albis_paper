#!/usr/bin/env bash
#SBATCH --job-name=native556_rctd
#SBATCH --partition=shared
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=3-00:00:00
set -euo pipefail
project=/dcs04/hicks/data/Jan/sim_project
out="${1:?Supply overview directory}"
method="${2:?Supply spider or sccube}"
case "$method" in spider|sccube) ;; *) exit 2 ;; esac
seed="${3:-20260921}"
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMBA_NUM_THREADS=1 PYTHONHASHSEED="$seed"
export MPLCONFIGDIR="$out/cache/$method/matplotlib" NUMBA_CACHE_DIR="$out/cache/$method/numba" XDG_CACHE_HOME="$out/cache/$method"
mkdir -p "$MPLCONFIGDIR" "$NUMBA_CACHE_DIR"
"$project/comparison_methods/env/$method/bin/python" -u "$project/sim_paper/code/comparison_methods/overview/generate.py" \
 --method "$method" --settings "$out/settings.json" --out "$out" > "$out/$method.log" 2>&1
"$project/comparison_methods/env/analysis/bin/python" - "$out" "$method" <<'PY'
import anndata as ad,sys
from pathlib import Path
p=Path(sys.argv[1])/sys.argv[2]
for modality in ['cell','bin','spot']:
 a=ad.read_h5ad(p/(modality+'.h5ad'),backed='r')
 assert a.n_vars==556 and a.var_names.is_unique
 if modality=='cell': assert a.n_obs==10000
 a.file.close()
print('Verified regenerated native 556-gene outputs')
PY
unset RCTD_N_GENES
export RCTD_SOURCE="$out" RCTD_RUN_PREFIX="${NATIVE556_RUN_PREFIX:-native556}"
bash "$project/sim_paper/code/comparison_methods/clustering/spatial_deconvolution/run_pilot.sh" "$method"
