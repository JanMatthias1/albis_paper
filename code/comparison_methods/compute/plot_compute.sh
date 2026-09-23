#!/bin/bash
set -euo pipefail
SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
PROJECT="$(cd "$SCRIPT_DIR/../../../.." && pwd)"
export MPLCONFIGDIR="${MPLCONFIGDIR:-/tmp/albis-compute-mpl-${USER}}"
export XDG_CACHE_HOME="${XDG_CACHE_HOME:-/tmp/albis-compute-cache-${USER}}"
exec "$PROJECT/comparison_methods/env/analysis/bin/python" "$SCRIPT_DIR/plot_compute.py" "$@"
