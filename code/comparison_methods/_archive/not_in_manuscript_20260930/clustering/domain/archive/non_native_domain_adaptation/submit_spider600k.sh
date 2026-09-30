#!/usr/bin/env bash
#SBATCH --job-name=spider_gyrus600k
#SBATCH --partition=shared
#SBATCH --cpus-per-task=4
#SBATCH --mem=128G
#SBATCH --time=3-00:00:00
#SBATCH --output=/dcs04/hicks/data/Jan/sim_project/sim_paper/code/comparison_methods/clustering/domain/logs/gyrus600k_%j.out
set -euo pipefail
project=/dcs04/hicks/data/Jan/sim_project
paper="$project/sim_paper"
code="$paper/code/comparison_methods/clustering/domain"
run="${1:-$paper/data/comparison_methods/clustering/domain/spider/gyrus600k_seed20260921}"
iterations="${2:-50000}"
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMBA_NUM_THREADS=1
export MPLCONFIGDIR="$run/cache/matplotlib" NUMBA_CACHE_DIR="$run/cache/numba" XDG_CACHE_HOME="$run/cache"
mkdir -p "$MPLCONFIGDIR" "$NUMBA_CACHE_DIR"
"$project/comparison_methods/env/spider/bin/python" -u "$code/generate_spider.py" \
    --n-cells 600000 --n-slices 10 --seed 20260921 --iterations "$iterations" \
    --pool "$paper/data/comparison_methods/overview_native556_20260926/splatter" \
    --out "$run/spider" > "$run/generation.log" 2>&1
"$project/comparison_methods/env/analysis/bin/python" - "$run/spider/cell.h5ad" <<'PY'
import sys
import anndata as ad
a = ad.read_h5ad(sys.argv[1], backed='r')
assert a.n_obs == 600000 and a.n_vars == 556
assert a.obs.domain_true.nunique() == 4
assert a.obs.cell_type_true.nunique() == 8
assert set(a.obs.slice_id.astype(int)) == set(range(10))
assert not a.obs[['domain_true','cell_type_true','slice_id']].isna().any().any()
print('Verified: 600,000 cells, 556 genes, four domains, eight cell types, ten slices.')
a.file.close()
PY
printf 'Shared domain/cell-type clustering input: %s\n' "$run/spider/cell.h5ad"
if [[ "$iterations" != 50000 ]]; then
    "$project/comparison_methods/env/analysis/bin/python" "$code/compare_spider_iterations.py" \
        --baseline "$paper/data/comparison_methods/clustering/domain/spider/gyrus600k_seed20260921/spider" \
        --current "$run/spider"
fi
