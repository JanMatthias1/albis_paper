#!/usr/bin/env bash
# Run after the per-modality/seed analysis jobs finish.
set -euo pipefail
SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
source "$SCRIPT_DIR/_env.sh"
export MPLCONFIGDIR="${MPLCONFIGDIR:-/tmp/fig3-mpl-${USER}}"
"$PYTHON_BIN" "$SCRIPT_DIR/ari_recovery_summary/plot_ari_recovery.py"
"$PYTHON_BIN" "$SCRIPT_DIR/ari_recovery_summary/plot_domain_vs_celltype.py"
"$PYTHON_BIN" "$SCRIPT_DIR/strong_domain_mix/batch_sigma_slide/plot_batch_sigma_slide_final.py"
