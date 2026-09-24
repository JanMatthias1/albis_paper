#!/bin/bash
# Figure 3, end to end: submits every job with the right dependencies.
#   bash sim_paper/code/clustering/run_figure3_global.sh
#
# Options (environment variables):
#   STAGES="celltype rctd banksy summary"   which stages to submit (default: all)
#   EXCLUDE=compute-158
#   DRY_RUN=1                                print the sbatch commands, submit nothing
#
# Stages and where their data comes from:
#   celltype  pca_harmony_single_cell/{cell,bin,bin16um,spot}_celltype_panel.sh
#             reads the Figure 2 datasets in data/figure_2/smaller_sphere/data/
#             -> data/figure_3/pca_harmony_single_cell/<modality>/
#   rctd      spatial_deconvolution/: canonical RCTD (Figure 2 cell reference +
#             Figure 2 spot query), 2 independent-seed cell references (generated
#             into data/figure_2/smaller_sphere/data/*_rctdref_seed*/ if absent),
#             RCTD per seed, then the 3-reference average (canonical + 2 seeds) (spatial_deconvolution/
#             weak_domain/) -> RCTD/weak_mix{,_seed*,_seed_avg}/; plus, when
#             banksy is also submitted, RCTD on the strong-mix spot bs0.3 query
#             with the strong-mix cell bs1.5 references at seeds 2025/101/202
#             (strong_domain/), then their average
#             -> RCTD/strong_mix{,_seed101,_seed202,_seed_avg}/
#             (all under data/figure_3/spatial_deconvolution/RCTD/)
#   banksy    strong_mix/generate_strong_mix_{cell,bin16um,spot}.sh (baseline
#             seed 2025) + cellbin_batch_sigma_slide/batch_sigma_slide_seeds.sh
#             (seeds 101/202; 72 tasks). Each task GENERATES its own data (Figure 2
#             parameters + --strong-domain-mix + its batch_sigma/seed), named by
#             those parameters, then QC -> seeded BANKSY -> Leiden ARI -> scores
#             -> data/figure_3/cellbin_batch_sigma_slide/<modality>/bs<s>[_seed<n>]/
#   summary   after everything above: strong_mix/check_matches_figure2.py (every
#             generated dataset = Figure 2 config + strong mix), then
#             plot_figure3_summaries.sh (plot_ari_recovery, plot_domain_vs_celltype,
#             plot_batch_sigma_slide_final)
#
# Data generation, QC, PCA/Harmony and BANKSY skip themselves when their output
# file already exists (so a re-run after a partial failure does not regenerate
# data); Leiden/ARI, scoring, plots and RCTD always rerun. To force a full
# rebuild, move data/figure_3/<stage folder> aside first.
# The summary stage uses afterany: it still runs if a task failed, so check
# logs/figure3_summary_*.out (the checker and the n= labels on the batch-sigma
# plot show what is missing).
set -euo pipefail

CODE=/dcs04/hicks/data/Jan/sim_project/sim_paper/code/clustering
LOGS="${CODE}/logs"
STAGES="${STAGES:-celltype rctd banksy summary}"
EXCLUDE="${EXCLUDE-compute-158}"
mkdir -p "${LOGS}"

sb() {  # sbatch --parsable, or print-only under DRY_RUN=1 (returns a fake job ID)
    if [[ -n "${DRY_RUN:-}" ]]; then
        echo "  sbatch $*" >&2; echo "DRY${RANDOM}"
    else
        sbatch --parsable "$@"
    fi
}
SB=(sb)
[[ -n "${EXCLUDE}" ]] && SB+=(--exclude="${EXCLUDE}")
has() { [[ " ${STAGES} " == *" $1 "* ]]; }
ALL=()   # job IDs the summary stage waits on

if has celltype; then
    for m in cell bin bin16um spot; do
        j=$("${SB[@]}" "${CODE}/pca_harmony_single_cell/${m}_celltype_panel.sh")
        echo "[celltype] ${m}: ${j}"; ALL+=("${j}")
    done
fi

if has rctd; then
    j=$("${SB[@]}" "${CODE}/spatial_deconvolution/weak_domain/run_rctd_spot.sh")
    echo "[rctd] canonical: ${j}"; ALL+=("${j}")
    ref=$("${SB[@]}" "${CODE}/spatial_deconvolution/weak_domain/generate_rctd_reference_independent_seed.sh")
    echo "[rctd] reference seeds: ${ref}"
    seeds=$("${SB[@]}" --dependency=afterok:"${ref}" "${CODE}/spatial_deconvolution/weak_domain/run_rctd_spot_independent_seed.sh")
    echo "[rctd] per-seed RCTD: ${seeds} (after ${ref})"
    canon="${j}"
    avg=$("${SB[@]}" --dependency=afterok:"${canon}":"${seeds}" --job-name=rctd_seed_avg \
        --output="${CODE}/spatial_deconvolution/logs/rctd_seed_avg_%j.out" \
        --time=01:00:00 --mem=16G --cpus-per-task=1 --partition=shared \
        --wrap="source ${CODE}/_env.sh && export MPLCONFIGDIR=/tmp/fig3-mpl-\${USER} && \"\${PYTHON_BIN}\" ${CODE}/spatial_deconvolution/average_rctd_seeds.py --config weak_mix && \"\${PYTHON_BIN}\" ${CODE}/spatial_deconvolution/plot_rctd_results.py /dcs04/hicks/data/Jan/sim_project/sim_paper/data/figure_3/spatial_deconvolution/RCTD/weak_mix_seed_avg")
    echo "[rctd] weak mix 3-reference average: ${avg} (after ${canon}, ${seeds})"; ALL+=("${avg}")
fi

if has banksy; then
    declare -A MIX
    for m in cell bin16um spot; do
        j=$("${SB[@]}" "${CODE}/strong_mix/generate_strong_mix_${m}.sh")
        echo "[banksy] strong-mix ${m}: ${j}"; ALL+=("${j}"); MIX[${m}]="${j}"
    done
    j=$("${SB[@]}" "${CODE}/cellbin_batch_sigma_slide/batch_sigma_slide_seeds.sh")
    echo "[banksy] seeds (72 tasks): ${j}"; ALL+=("${j}"); MIX[seeds]="${j}"
    if has rctd; then
        # afterany: the RCTD inputs are single tasks of those arrays, so one
        # unrelated failed task shouldn't block it; each RCTD task exits with
        # an error itself if its inputs are missing (and then the average
        # below, afterok, doesn't run).
        s=$("${SB[@]}" --dependency=afterany:"${MIX[cell]}":"${MIX[spot]}":"${MIX[seeds]}" \
            "${CODE}/spatial_deconvolution/strong_domain/run_rctd_spot_strong_domain.sh")
        echo "[rctd] strong mix, 3 references: ${s} (after ${MIX[cell]}, ${MIX[spot]}, ${MIX[seeds]})"
        # ${s} is one 3-task array (task 0 = canonical seed 2025, tasks 1-2 =
        # seeds 101/202), so afterok on it already waits for the canonical run
        # too -- unlike weak mix, where canonical is a separate job.
        avg=$("${SB[@]}" --dependency=afterok:"${s}" --job-name=rctd_strong_avg \
            --output="${CODE}/spatial_deconvolution/logs/rctd_strong_avg_%j.out" \
            --time=01:00:00 --mem=16G --cpus-per-task=1 --partition=shared \
            --wrap="source ${CODE}/_env.sh && export MPLCONFIGDIR=/tmp/fig3-mpl-\${USER} && \"\${PYTHON_BIN}\" ${CODE}/spatial_deconvolution/average_rctd_seeds.py --config strong_mix && \"\${PYTHON_BIN}\" ${CODE}/spatial_deconvolution/plot_rctd_results.py /dcs04/hicks/data/Jan/sim_project/sim_paper/data/figure_3/spatial_deconvolution/RCTD/strong_mix_seed_avg")
        echo "[rctd] strong mix 3-reference average: ${avg} (after ${s})"; ALL+=("${avg}")
    fi
fi

if has summary; then
    DEP=()
    if [[ ${#ALL[@]} -gt 0 ]]; then
        DEP=(--dependency=afterany:$(IFS=:; echo "${ALL[*]}"))
    fi
    j=$("${SB[@]}" "${DEP[@]}" --job-name=figure3_summary \
        --output="${LOGS}/figure3_summary_%j.out" \
        --time=02:00:00 --mem=64G --cpus-per-task=2 --partition=shared \
        --wrap="source ${CODE}/_env.sh && \"\${PYTHON_BIN}\" ${CODE}/strong_mix/check_matches_figure2.py && bash ${CODE}/plot_figure3_summaries.sh")
    echo "[summary] checker + plots: ${j}${DEP:+ (after all of the above)}"
fi
