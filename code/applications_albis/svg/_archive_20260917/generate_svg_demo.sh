#!/bin/bash
#SBATCH --job-name=svg_demo
#SBATCH --output=/dcs04/hicks/data/Jan/sim_project/sim_paper/code/applications_albis/svg/logs_nnsvg/svg_demo_%j.out
#SBATCH --time=10:00:00
#SBATCH --mem=250G
#SBATCH --cpus-per-task=8
#SBATCH --partition=shared
#
# Figure 4A *illustrative* dataset -- NOT a benchmark input.
#
# The AUROC/ROC panel (plot_figure4a_nnsvg.py) stays on the real shared
# datasets and honestly reports "SVG signal is detectable but modest"
# (marker genes reach only ~1.3x domain enrichment at marker_foldchange=3.5 /
# ~49% domain purity -> nnSVG AUROC ~0.65-0.75, nothing that "lights up" on a
# map). This script builds a SEPARATE, deliberately strong-signal sim used
# only for the qualitative "here is how you would run an SVG study with the
# package" panel (plot_figure4a_nnsvg_examples.py) -- same benchmark-vs-
# capability split as Figure 5a's representative outputs.
#
# Config: the bin16um strong-domain-mix family, with the marker contrast
# cranked and batch noise dropped so marker genes visibly track their
# domain:
#   --base-gene-lognormal -1.0 0.7   (higher baseline: bins well-populated)
#   --strong-domain-mix              (concentrate cell types in domains)
#   --marker-foldchange 10           (vs 3.5 default)
#   --shared-marker-foldchange 5     (vs 2.5 default)
#   --batch-sigma 0.3                (vs 0.7; cleaner, it is a demo)
#   --sphere-r-um 2050, --bin-size-um 16   (as the bin16um family)
# theta / theta-jitter left at generate defaults (2.0 / 1.0), as the rest of
# the bin family.
#
# Output: data/figure_4/svg/demo/<TAG>/simulation_bin_z_qc.h5ad  and, chained,
# the nnSVG ranking under data/figure_4/svg/nnsvg/bin16um_svgdemo/.
#
# Usage:  sbatch sim_paper/code/applications_albis/svg/generate_svg_demo.sh

set -euo pipefail

mkdir -p /dcs04/hicks/data/Jan/sim_project/sim_paper/code/applications_albis/svg/logs_nnsvg
source /dcs04/hicks/data/Jan/sim_project/sim_paper/code/count_distribution/_env.sh   # sets PYTHON_BIN (albis env)
cd /dcs04/hicks/data/Jan/sim_project

TAG="svgdemo_bin16um_mfc10_strongmix"
MODALITY="bin"
NOISY_DIR="sim_paper/data/noisy/${TAG}"
DEST_DIR="sim_paper/data/figure_4/svg/demo/${TAG}"
QC_H5AD="${DEST_DIR}/simulation_${MODALITY}_z_qc.h5ad"

STAGATE_PY="/dcs04/hicks/data/Jan/sim_project/sim_paper/env/stagate-pyg/bin/python"
NNSVG_OUT="sim_paper/data/figure_4/svg/nnsvg/bin16um_svgdemo"

echo "[config] tag=${TAG} modality=${MODALITY}"

# ---- 1. generate + QC (albis env) --------------------------------------------
if [[ ! -f "${QC_H5AD}" ]]; then
    if [[ ! -f "${NOISY_DIR}/simulation_${MODALITY}_z.h5ad" ]]; then
        echo "[generate] ${TAG}"
        "${PYTHON_BIN}" sim_paper/code/data/generate_simulation_noisy.py \
            --modality "${MODALITY}" \
            --base-gene-lognormal -1.0 0.7 \
            --strong-domain-mix \
            --marker-foldchange 10 \
            --shared-marker-foldchange 5 \
            --batch-sigma 0.3 \
            --sphere-r-um 2050 \
            --bin-size-um 16 \
            --out-tag "${TAG}"
    fi
    if [[ ! -f "${NOISY_DIR}/simulation_${MODALITY}_z_qc.h5ad" ]]; then
        echo "[qc] ${TAG}"
        "${PYTHON_BIN}" sim_paper/code/clustering/00_qc_filter.py --modality "${MODALITY}" --packing-tag "${TAG}"
    fi
    mkdir -p "sim_paper/data/figure_4/svg/demo"
    echo "[move] ${NOISY_DIR} -> ${DEST_DIR}"
    mv "${NOISY_DIR}" "${DEST_DIR}"
fi

# ---- 2. nnSVG ranking on slice 5 (stagate-pyg env) --------------------------
echo "[nnsvg] ${TAG}"
"${STAGATE_PY}" sim_paper/code/applications_albis/svg/3D_nnsvg.py \
    --dataset bin16um \
    --input "${QC_H5AD}" \
    --run-tag bin16um_svgdemo \
    --output-dir "${NNSVG_OUT}" \
    --n-threads 8

echo "[done] ${TAG}"
echo "  sim:   ${QC_H5AD}"
echo "  nnsvg: ${NNSVG_OUT}/gene_results_bin16um_svgdemo.csv"
echo "  then:  python sim_paper/code/applications_albis/svg/plot_figure4a_nnsvg_examples.py \\"
echo "           --dataset bin16um --h5ad ${QC_H5AD} \\"
echo "           --nnsvg-results ${NNSVG_OUT}/gene_results_bin16um_svgdemo.csv \\"
echo "           --output sim_paper/data/figure_4/svg/nnsvg/figure_4a_nnsvg_examples_demo.png"
