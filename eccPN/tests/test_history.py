"""Analytic and quadrature checks of the causal history kernel."""

import unittest

import numpy as np
from scipy.integrate import quad

from eccpn.history import integrate_log_panel


class HistoryKernel(unittest.TestCase):
    def test_constant_at_split_endpoint(self):
        duration = 3.4
        scale = 1.7
        expected = duration * (np.log(duration / scale) - 1)
        self.assertAlmostEqual(integrate_log_panel(10, 10 + duration, [1], 10 + duration, scale), expected)

    def test_quadrature_on_clipped_near_and_remote_panels(self):
        coefficients = np.array([1, 0.4, -0.22, 0.15])
        scale = 0.83
        for time in (-2, 0.63, 1, 1.07, 1.5, 3, 40, 1e7):
            with self.subTest(time=time):
                for power in (0, 1, 2):
                    actual = integrate_log_panel(0, 1, coefficients, time, scale, power)
                    expected = 0.0
                    if time > 0:
                        expected = quad(
                            lambda t: np.polynomial.polynomial.polyval(t, coefficients)
                            * np.log((time - t) / scale)**power,
                            0, min(time, 1), epsabs=3e-11, epsrel=1e-12,
                        )[0]
                    np.testing.assert_allclose(actual, expected, atol=4e-11, rtol=3e-12)

    def test_panel_split_and_tensor_components(self):
        source = [[1.0, -3.0], [2.0, 0.25]]
        time = 2.1
        scale = 0.7
        result = integrate_log_panel(0, 2, source, time, scale, 2)
        left = integrate_log_panel(0, 1, [[1, -3], [1, 0.125]], time, scale, 2)
        right = integrate_log_panel(1, 2, [[2, -2.875], [1, 0.125]], time, scale, 2)
        np.testing.assert_allclose(left + right, result, atol=3e-14, rtol=3e-14)
