"""Figure 3 strong-mix inputs shared with BANKSY; simulation seeds, fixed training seed."""
from pathlib import Path

PROJECT = Path('/dcs04/hicks/data/Jan/sim_project')
SLIDE = PROJECT / 'sim_paper/data/figure_3/strong_domain_mix/batch_sigma_slide'
# Since 2026-09-23 the strong-mix generators name each QC file by its generation
# parameters: <Figure 2 tag body>_strongmix_bsigma<batch, no dot>[_seed<seed>]_qc.h5ad
# (see code/clustering/strong_domain_mix/generate/generate_strong_mix_*.sh).
TAG_BODY = {
    'cell': 'log_mu_-2.5_theta_0.40_jitter0.15',
    'bin16um': 'packing_pf0p04_bin16um_log_mu_-2.5_jitter0.6',
    'spot': 'packing_pf0p04_log_mu_-2.25_sigma1.0_theta_0.25_jitter0.10',
}
DATASETS = {}
for technology, canonical in [('cell', '1.5'), ('bin16um', '0.7'), ('spot', '0.3')]:
    for batch in ['0', '0.05', canonical]:
        for seed in [2025, 101, 202]:
            folder_batch = '0.0' if technology == 'spot' and batch == '0' and seed == 2025 else batch
            suffix = '' if seed == 2025 else f'_seed{seed}'
            name = f'{technology}_strongmix_bs{batch}_seed{seed}'
            data_tag = f"{TAG_BODY[technology]}_strongmix_bsigma{f'{float(batch):g}'.replace('.', '')}{suffix}"
            DATASETS[name] = dict(
                h5ad=str(SLIDE / technology / f'bs{folder_batch}{suffix}' / f'{data_tag}_qc.h5ad'),
                technology=technology, batch_sigma=float(batch), simulation_seed=seed,
            )
            if technology == 'cell':
                DATASETS[name].update(graph_model='knn', k_2d=6, k_z=3)

