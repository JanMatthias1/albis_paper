#!/bin/bash
# Weak domain mix, batch_sigma 0: final BANKSY runs for domain recovery at the
# chosen λ (final_tasks.tsv; chosen on the strong mix, see ../../../strong_domain_mix/banksy_harmony_batch_zero/README.md).
# Keeps every h5ad and draws UMAP + spatial plots (MODE=final)
# -> data/figure_3/weak_domain_mix/banksy_harmony_batch_zero/domain/<modality>/lam<λ>/
set -euo pipefail
HERE=/dcs04/hicks/data/Jan/sim_project/sim_paper/code/clustering/weak_domain_mix/banksy_harmony_batch_zero/domain
bash "${HERE}/../../../submit_banksy_lambda.sh" "${HERE}/final_tasks.tsv" \
    /dcs04/hicks/data/Jan/sim_project/sim_paper/data/figure_3/weak_domain_mix/banksy_harmony_batch_zero/logs final
