"""Figure 3 strong-mix inputs shared with BANKSY; simulation seeds, fixed training seed."""
from pathlib import Path

PROJECT = Path('/dcs04/hicks/data/Jan/sim_project')
SLIDE = PROJECT / 'sim_paper/data/figure_3/cellbin_batch_sigma_slide'
DATASETS = {}
for technology, canonical in [('cell', '1.5'), ('bin16um', '0.7'), ('spot', '0.3')]:
    modality = 'bin' if technology == 'bin16um' else technology
    for batch in ['0', '0.05', canonical]:
        for seed in [2025, 101, 202]:
            folder_batch = '0.0' if technology == 'spot' and batch == '0' and seed == 2025 else batch
            suffix = '' if seed == 2025 else f'_seed{seed}'
            name = f'{technology}_strongmix_bs{batch}_seed{seed}'
            DATASETS[name] = dict(
                h5ad=str(SLIDE / technology / f'bs{folder_batch}{suffix}' / f'simulation_{modality}_z_qc.h5ad'),
                technology=technology, batch_sigma=float(batch), simulation_seed=seed,
            )
            if technology == 'cell':
                DATASETS[name].update(graph_model='knn', k_2d=6, k_z=3)

