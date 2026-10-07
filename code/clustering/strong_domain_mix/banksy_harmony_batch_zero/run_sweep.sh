#!/bin/bash
# Strong domain mix, batch_sigma 0: BANKSY lambda sweep (sweep_tasks.tsv)
# -> data/figure_3/strong_domain_mix/banksy_harmony_batch_zero/sweep/<modality>/lam<lambda>/
# Rows already run are skipped. DRY_RUN=1 prints the sbatch commands only.
set -euo pipefail
HERE=/dcs04/hicks/data/Jan/sim_project/albis_paper/code/clustering/strong_domain_mix/banksy_harmony_batch_zero
bash "${HERE}/../../submit_banksy_lambda.sh" "${HERE}/sweep_tasks.tsv" \
    /dcs04/hicks/data/Jan/sim_project/albis_paper/data/figure_3/strong_domain_mix/banksy_harmony_batch_zero/logs sweep
