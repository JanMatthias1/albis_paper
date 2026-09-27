"""Reproduce Figure 5A ALBIS parameters with native ALBIS simulation calls."""
import argparse
import hashlib
import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[6]
SOURCE = ROOT / 'sim_paper/data/figure_4/cross_modality_alignment/strong_domain_mix_shift3x'
OUT = ROOT / 'sim_paper/data/comparison_methods/spatial_deconvolution/albis_16um_spot_parameters'
SEEDS = [2025, 101, 202]


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--seed', type=int, choices=SEEDS, required=True)
    args = parser.parse_args()
    run = OUT / f'seed{args.seed}'
    cfg = json.loads((run / 'config.json').read_text())
    assert cfg['seed'] == args.seed and cfg['n_cells'] == 600000
    dest = run / 'data'
    dest.mkdir(exist_ok=False)
    sys.path.insert(0, str(ROOT / 'albis'))
    import albis
    source_file = ROOT / 'albis/albis/simulation_sphere.py'
    source_hash = hashlib.sha256(source_file.read_bytes()).hexdigest()
    expected = json.loads((SOURCE / 'data/spot/generation_manifest.json').read_text())['source_sha256']['albis/albis/simulation_sphere.py']
    assert source_hash == expected, 'ALBIS source changed from the original Figure 4 run'
    (run / 'generation_status.json').write_text(json.dumps({'status': 'running', 'source_sha256': source_hash}, indent=2))
    result = albis.simulate_3d_molecule_sphere_multires(**cfg)
    records = {}
    for folder, modality, key in [('cell', 'cell', 'adata_cell_sectioned'), ('bin16um', 'bin', 'bin_adatas'), ('spot', 'spot', 'spot_adatas')]:
        a = result[key]['Z']
        assert a.n_vars == 556 and a.obs.slice_id.nunique() == 10
        if modality == 'cell':
            assert a.n_obs == 600000 and a.obs.cell_type_true.nunique() == 8
        target = dest / folder
        target.mkdir()
        a.uns['native_generation_parameters'] = cfg
        a.write_h5ad(target / f'simulation_{modality}_z.h5ad')
        records[modality] = {'n_obs': a.n_obs, 'n_genes': a.n_vars}
        print('Saved', modality, records[modality], flush=True)
    (run / 'generation_status.json').write_text(json.dumps({'status': 'complete', 'source_sha256': source_hash, 'outputs': records}, indent=2) + '\n')


if __name__ == '__main__':
    main()
