"""Check native ALBIS max_shift=3075 against the historical section-5 shift3x.

This replays native perturbation draws from existing truth; it does not claim
to validate a full tissue regeneration. A separate native generation run does that.
"""
import argparse
import json
from pathlib import Path
import sys

import h5py
import numpy as np
from anndata._io.specs import read_elem

ROOT = Path(__file__).resolve().parents[5]
sys.dont_write_bytecode = True
sys.path.insert(0, str(ROOT / 'albis'))
from albis.simulation_sphere import random_rigid_inplane_2d


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    base = ROOT / 'sim_paper/data/figure_4/cross_modality_alignment'
    parser.add_argument('--source', type=Path, default=base / 'independent_offsets')
    # The post-hoc shift3x this replays against was archived when the native
    # regen was promoted into independent_offsets_shift3x/ (2026-09-20).
    parser.add_argument('--comparison', type=Path,
                        default=base / 'independent_offsets_shift3x_before_native_20260920T162259708032Z')
    parser.add_argument('--output', type=Path, default=base / 'native_offsets_3075_seed12345/native_equivalence.json')
    parser.add_argument('--slice', type=int, default=5)
    args = parser.parse_args()
    rows = []
    for tech, modality in [('bin16um', 'bin'), ('spot', 'spot'), ('cell', 'cell')]:
        filename = f'simulation_{modality}_z_qc.h5ad'
        with h5py.File(args.source / 'data' / tech / filename) as fh:
            obs = read_elem(fh['obs'])
            xy = read_elem(fh['obsm/spatial'])
            params = read_elem(fh['uns/rigid_perturb_inplane/spatial_unaligned'])
        mask = obs.slice_id.astype(str).to_numpy() == str(args.slice)
        original = params[str(args.slice)]
        native, rotation, translation = random_rigid_inplane_2d(
            xy[mask], max_deg=270., max_shift=3075., rng=np.random.default_rng(int(original['seed'])))
        np.testing.assert_array_equal(rotation, original['R'])
        np.testing.assert_allclose(translation, np.asarray(original['t'])*3, atol=1e-10)
        with h5py.File(args.comparison / 'data' / tech / filename) as fh:
            comparison_obs = read_elem(fh['obs'])
            comparison_xy = read_elem(fh['obsm/spatial_unaligned'])
            truth = read_elem(fh['obsm/spatial'])
        np.testing.assert_array_equal(obs.index[mask], comparison_obs.index)
        np.testing.assert_array_equal(xy[mask], truth)
        max_error = float(np.max(np.abs(native.astype(float)-comparison_xy)))
        if max_error >= .001:
            raise ValueError(f'{tech}: native/shift3x mismatch {max_error} µm')
        rows.append({'modality': tech, 'slice_id': args.slice, 'native_seed': int(original['seed']),
                     'n_obs': int(mask.sum()), 'max_coordinate_difference_um': max_error,
                     'native_rotation': rotation.tolist(), 'native_translation_um': translation.tolist()})
    report = {'scope': 'Native coordinate RNG replay on existing truth; not full regeneration',
              'max_shift': 3075., 'max_deg': 270., 'base_seed_unaligned': 12345,
              'sync_unaligned_seed': False, 'tolerance_um': .001, 'results': rows}
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(report, indent=2) + '\n')
    print(json.dumps(report, indent=2))


if __name__ == '__main__':
    main()
