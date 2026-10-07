#!/usr/bin/env bash
# Redraw the current Figure 5A from saved data, without resimulation.
set -euo pipefail
code="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
project="$(cd -- "$code/../../../.." && pwd)"
out="${1:-$project/albis_paper/data/figure_5/figure_5A_600k}"
export MPLCONFIGDIR="$out/cache/matplotlib" XDG_CACHE_HOME="$out/cache"
export OPENBLAS_NUM_THREADS=1
mkdir -p "$MPLCONFIGDIR"
exec "$project/comparison_methods/env/analysis/bin/python" "$code/plot_combined.py" --input "$out"
