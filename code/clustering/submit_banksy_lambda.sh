#!/bin/bash
# Submit every row of a BANKSY lambda task table: one array per modality
# (memory differs), all running run_banksy_lambda.sh.
#
#   bash submit_banksy_lambda.sh TASKS_TSV LOG_DIR [sweep|final]
#
# Called by the run_*.sh wrappers in strong_domain_mix/banksy_harmony_batch_zero/ and
# weak_domain_mix/banksy_harmony_batch_zero/. Rows already done are not submitted (sweep: ARI
# summary exists; final: ARI summary and ari/plots/ exist). DRY_RUN=1 prints
# the sbatch commands only. DEPENDENCY=<sbatch dependency, e.g. afterany:123>
# gates every submitted array (the runner fails loudly if an input is missing).
# EXCLUDE=<nodes> avoids nodes (default: the known-bad compute-158).
set -euo pipefail

TASKS="$(realpath "${1:?usage: submit_banksy_lambda.sh TASKS_TSV LOG_DIR [sweep|final]}")"
LOG_DIR="${2:?LOG_DIR required}"
MODE="${3:-sweep}"
[[ "${MODE}" == sweep || "${MODE}" == final ]] || { echo "mode must be sweep or final" >&2; exit 2; }
CODE=/dcs04/hicks/data/Jan/sim_project/sim_paper/code/clustering
EXCLUDE="${EXCLUDE-compute-158}"
mkdir -p "${LOG_DIR}"

# final runs also draw UMAPs from the full bin16um h5ad -> more memory
mem_for() {
    case "$1" in
        bin16um) [[ "${MODE}" == final ]] && echo 120G || echo 96G ;;
        *) echo 32G ;;
    esac
}

# leave out rows that are already done
todo() {
    local m="$1" task mod lam out pmod run
    while IFS=$'\t' read -r task mod _ lam _ out _; do
        [[ "${mod}" == "${m}" ]] || continue
        pmod="${mod}"; [[ "${mod}" == bin16um ]] && pmod=bin
        run="/dcs04/hicks/data/Jan/sim_project/${out}"
        if [[ -f "${run}/ari/ari_summary_${pmod}.json" && ( "${MODE}" == sweep || -d "${run}/ari/plots" ) ]]; then
            continue
        fi
        echo "${task}"
    done < <(tail -n +2 "${TASKS}") | paste -sd,
}

for m in cell bin16um spot; do
    idx="$(todo "${m}")"
    [[ -n "${idx}" ]] || { echo "${m}: nothing to run"; continue; }
    args=(--parsable --array="${idx}" --mem="$(mem_for "${m}")" --export=ALL,MODE="${MODE}"
          --job-name="banksy_${MODE}_${m}" --output="${LOG_DIR}/banksy_${MODE}_${m}_%A_%a.out")
    [[ -n "${EXCLUDE}" ]] && args+=(--exclude="${EXCLUDE}")
    [[ -n "${DEPENDENCY:-}" ]] && args+=(--dependency="${DEPENDENCY}")
    if [[ -n "${DRY_RUN:-}" ]]; then
        echo "sbatch ${args[*]} ${CODE}/run_banksy_lambda.sh ${TASKS}"
    else
        echo "${m} tasks ${idx}: $(sbatch "${args[@]}" "${CODE}/run_banksy_lambda.sh" "${TASKS}")"
    fi
done
