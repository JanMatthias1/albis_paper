"""Analytic checks for reference-only scoring; run with unittest."""
import unittest
import numpy as np
from reference_metrics import evaluate, rigid_fit, TECHS


class ReferenceMetricTests(unittest.TestCase):
    def setUp(self):
        self.points = np.array([[0., 0.], [3., 0.], [0., 2.], [2., 4.]])
        self.truth = np.tile(self.points, (3, 1))
        self.tech = np.repeat(TECHS, 4)

    def score(self, xy, truth=None, tech=None):
        return evaluate(self.truth if truth is None else truth, {"test": xy},
                        self.tech if tech is None else tech, "bin16um")[0]["stages"]["test"]

    def test_common_rigid_motion_is_removed(self):
        r = np.array([[0., -1.], [1., 0.]])
        row = self.score(self.truth @ r + [400., -700.])
        self.assertLess(row["balanced_target_rmse_um"], 1e-10)

    def test_target_error_survives_and_reference_is_excluded(self):
        xy = self.truth.copy()
        xy[self.tech == "spot"] += [3., 4.]
        xy[self.tech == "cell"] += [0., 12.]
        row = self.score(xy)
        self.assertAlmostEqual(row["per_modality"]["spot"]["rmse_um"], 5.)
        self.assertAlmostEqual(row["balanced_target_rmse_um"], np.sqrt((25+144)/2))

    def test_density_does_not_change_balanced_score(self):
        xy = self.truth.copy()
        xy[self.tech == "spot"] += [3., 4.]
        ix = np.r_[np.arange(12), np.tile(np.arange(4, 8), 20)]
        self.assertAlmostEqual(self.score(xy)["balanced_target_rmse_um"],
                               self.score(xy[ix], self.truth[ix], self.tech[ix])["balanced_target_rmse_um"])

    def test_reflection_and_scale_are_not_removed(self):
        for xy in (self.truth * [-1, 1], self.truth * 2):
            self.assertGreater(self.score(xy)["balanced_target_rmse_um"], .1)
        r, _ = rigid_fit(self.points * [-1, 1], self.points)
        self.assertAlmostEqual(np.linalg.det(r), 1.)

    def test_reference_changes_target_set(self):
        result, _ = evaluate(self.truth, {"test": self.truth}, self.tech, "spot")
        self.assertEqual(result["scored_modalities"], ["bin16um", "cell"])

    def test_invalid_inputs_fail(self):
        with self.assertRaises(ValueError):
            rigid_fit(np.zeros((4, 2)), self.points)
        xy = self.truth.copy(); xy[5, 0] = np.nan
        with self.assertRaises(ValueError):
            self.score(xy)
        with self.assertRaises(ValueError):
            evaluate(self.truth, {"test": self.truth}, self.tech, "unknown")


if __name__ == "__main__":
    unittest.main()
