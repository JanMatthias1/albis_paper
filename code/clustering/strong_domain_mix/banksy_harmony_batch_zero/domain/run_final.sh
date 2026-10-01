#!/bin/bash
# Strong domain mix, batch_sigma 0: final BANKSY runs for domain recovery at the
# chosen λ (final_tasks.tsv; chosen on the strong mix).
# Keeps every h5ad and draws UMAP + spatial plots (MODE=final)
# -> data/figure_3/strong_domain_mix/banksy_harmony_batch_zero/domain/<modality>/lam<λ>/
set -euo pipefail
HERE=/dcs04/hicks/data/Jan/sim_project/sim_paper/code/clustering/strong_domain_mix/banksy_harmony_batch_zero/domain
bash "${HERE}/../../../submit_banksy_lambda.sh" "${HERE}/final_tasks.tsv" \
    /dcs04/hicks/data/Jan/sim_project/sim_paper/data/figure_3/strong_domain_mix/banksy_harmony_batch_zero/logs final
