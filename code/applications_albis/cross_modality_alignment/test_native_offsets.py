"""Verify native perturbation magnitude and seed semantics independently of alignment."""
import unittest
from pathlib import Path
import sys
sys.dont_write_bytecode = True
sys.path.insert(0, str(Path(__file__).resolve().parents[4] / "albis"))
import numpy as np
import albis as ab
from albis.simulation_sphere import random_rigid_inplane_2d


class NativeOffsetsTest(unittest.TestCase):
    def test_same_seed_scales_translation_not_rotation(self):
        points = np.array([[0., 0.], [2000., -1000.], [-800., 300.], [200., 900.]])
        for seed in (22440, 32472, 42472):
            old, old_r, old_t = random_rigid_inplane_2d(points, 270., 1025., np.random.default_rng(seed))
            new, new_r, new_t = random_rigid_inplane_2d(points, 270., 3075., np.random.default_rng(seed))
            np.testing.assert_array_equal(old_r, new_r)
            np.testing.assert_allclose(new_t, 3*old_t, atol=1e-10)
            np.testing.assert_allclose(new, old.astype(float) + 2*old_t, atol=.001)

    def test_native_generation_separates_tissue_and_perturbation_seeds(self):
        kwargs = dict(output='cell', slice_axis='Z', n_cells=100, n_slices=3,
                      sphere_radius_um=250., seed=7, max_deg=270.)
        baseline = ab.generate_data(**kwargs, max_shift=1025., base_seed_unaligned=12345)
        stronger = ab.generate_data(**kwargs, max_shift=3075., base_seed_unaligned=12345)
        independent = ab.generate_data(**kwargs, max_shift=3075., base_seed_unaligned=54321)
        for other in (stronger, independent):
            self.assertEqual((baseline.X != other.X).nnz, 0)
            np.testing.assert_array_equal(baseline.obsm['spatial'], other.obsm['spatial'])
            np.testing.assert_array_equal(baseline.obs_names, other.obs_names)
        for sid, old in baseline.uns['rigid_perturb_inplane']['spatial_unaligned'].items():
            new = stronger.uns['rigid_perturb_inplane']['spatial_unaligned'][sid]
            np.testing.assert_array_equal(old['R'], new['R'])
            np.testing.assert_allclose(np.asarray(old['t'])*3, new['t'], atol=1e-10)
        self.assertFalse(np.allclose(stronger.obsm['spatial_unaligned'], independent.obsm['spatial_unaligned']))


if __name__ == '__main__':
    unittest.main()
