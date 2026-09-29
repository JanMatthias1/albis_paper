#!/usr/bin/env bash
# Run after the per-modality/seed analysis jobs finish. The batch-zero bars need
# the no_harmony_domain/ and no_harmony_cell_type/ experiments to be complete
# (their own run.sh; not submitted by run_figure3_global.sh).
set -euo pipefail
SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
source "$SCRIPT_DIR/_env.sh"
export MPLCONFIGDIR="${MPLCONFIGDIR:-/tmp/fig3-mpl-${USER}}"
"$PYTHON_BIN" "$SCRIPT_DIR/ari_recovery_summary/plot_banksy_vs_pca_recovery.py"
"$PYTHON_BIN" "$SCRIPT_DIR/ari_recovery_summary/plot_banksy_vs_pca_condensed.py"
"$PYTHON_BIN" "$SCRIPT_DIR/strong_domain_mix/batch_sigma_slide/plot_batch_sigma_slide_final.py"
