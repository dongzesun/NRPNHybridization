"""Public radiation on ordered physical histories with unordered observer cuts."""

from concurrent.futures import ProcessPoolExecutor, ThreadPoolExecutor
from dataclasses import fields, replace
import multiprocessing
import unittest

import numpy as np

from eccpn import Accuracy, Binary, InitialState, compute_radiation, evolve
from eccpn.frame import BondiFrame
from eccpn.history import FiniteOrbitalPast, batched_log_weights, reuse_radiation_weights


def compare_fields(actual, expected):
    if isinstance(expected, dict):
        assert actual.keys() == expected.keys()
        for key in expected:
            compare_fields(actual[key], expected[key])
    elif isinstance(expected, BondiFrame):
        for item in fields(expected):
            compare_fields(getattr(actual, item.name), getattr(expected, item.name))
    else:
        np.testing.assert_array_equal(actual, expected)


class ParallelRadiation(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        binary = Binary(.71, .29, c=1.4)
        initial = InitialState((23., 0., 0.), (-.022, .262, .033))
        precursor = evolve(binary, initial, None, (-.4, 0.), Accuracy("newtonian"))
        initial = InitialState(*precursor.state(0.))
        spin = binary.mu * np.cross(initial.position, initial.velocity)
        past = FiniteOrbitalPast((-.4, -.2, 0.), precursor.state,
                                 {0: binary.mass, 2: -.0014, 4: .00005},
                                 {0: spin, 2: .018 * spin}, 1., "parallel-radiation")
        zero = (0., 0., 0.)
        cls.trajectory = evolve(binary, initial, past.initial_data(binary, zero, zero, zero),
                                (0., .36), Accuracy(max_step=.15))

    def test_cut_weights_respect_each_moment_grid(self):
        first = np.array([-2., -.3, .4])
        other = np.array([-2., -.7, .4])
        ends = [np.r_[grid[1:], 1.1] for grid in (first, other)]
        expected = [batched_log_weights(grid, end, .85, 2., 2) for grid, end in zip((first, other), ends)]
        with reuse_radiation_weights():
            a = batched_log_weights(first, ends[0], .85, 2., 2)
            self.assertIs(a, batched_log_weights(first.copy(), ends[0].copy(), .85, 2., 2))
            b = batched_log_weights(other, ends[1], .85, 2., 2)
        np.testing.assert_array_equal(a, expected[0])
        np.testing.assert_array_equal(b, expected[1])

    def test_processes_all_native_products_and_memory(self):
        times = np.array([[.35, .03], [.14, .35]])
        expected = compute_radiation(self.trajectory, times)
        with ProcessPoolExecutor(2, mp_context=multiprocessing.get_context("spawn")) as pool:
            actual = compute_radiation(self.trajectory, times, executor=pool)
            self.assertEqual(pool.submit(abs, -2).result(), 2)
        for item in fields(expected):
            with self.subTest(field=item.name):
                compare_fields(getattr(actual, item.name), getattr(expected, item.name))
        self.assertTrue(actual.retained_memory_strain)
        self.assertTrue(actual.energy_flux_parts)

    def test_threads_native_and_bondi_cuts_share_history(self):
        times = np.array([.32, .05, .32])
        accuracy = replace(self.trajectory.accuracy, waveform_pn=3., retain_memory_excess=False)
        frame = BondiFrame(ell_max=2, polar_points=12, azimuthal_points=23,
                           supertranslation={(0, 0): .01})
        expected = compute_radiation(self.trajectory, times, frame, accuracy)
        with ThreadPoolExecutor(2) as pool:
            actual = compute_radiation(self.trajectory, times, frame, accuracy, executor=pool)
        for item in fields(expected):
            with self.subTest(field=item.name):
                compare_fields(getattr(actual, item.name), getattr(expected, item.name))


if __name__ == "__main__":
    unittest.main()
