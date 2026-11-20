"""Canonical gauge products and generic luminosity contractions."""

from fractions import Fraction
from math import comb, factorial
import unittest

import numpy as np
from scipy.integrate import quad

from eccpn import Accuracy, Binary, InitialState, compute_radiation, evolve
from eccpn.canonical import canonical_mass_quadrupole
from eccpn.gauge_moments import canonical_source_jets
from eccpn.history import PanelHistory, SourcePanel
from eccpn.luminosity import (
    canonical_local_luminosity, canonical_luminosity_from_jets, canonical_nonlocal_luminosity,
    causal_canonical_nonlocal_luminosity, nonlinear_luminosity,
)
from eccpn.newtonian_moments import newtonian_source_jets
from eccpn.tails import quadrupole_tail_sources
from eccpn.waveform_jets import canonical_waveform_jets


class CanonicalAssembly(unittest.TestCase):
    def test_published_nonlocal_power_differs_from_retarded_force_for_past_pulse(self):
        binary = Binary(0.62, 0.38, G=1.05, c=1.6)
        position, velocity = np.array([13., 0., 0.]), np.array([0., 0.29, 0.03])
        time, infrared = 0.7, 0.91
        current = newtonian_source_jets(binary, position, velocity, 9, 2)["M2"][0]
        shape = np.array([[0., 1., 0.2], [1., 0., 0.], [0.2, 0., 0.]])
        base = {9 + k: (-1)**k * comb(9, k) for k in range(10)}
        panels = []
        for piece in range(16):
            a, width = Fraction(piece, 16), Fraction(1, 16)
            coefficients = np.zeros((19, 10, 3, 3))
            for order in range(10):
                for k in range(19 - order):
                    exact = sum(
                        value * factorial(p) // factorial(p - order) * comb(p - order, k)
                        * a**(p - order - k) * width**k
                        for p, value in base.items() if p >= order + k
                    )
                    coefficients[k, order] = float(exact) * shape
            panels.append(SourcePanel(-2 + float(a), -2 + float(a + width), coefficients))
        # Local jets and initial endpoints are held fixed. Only the
        # hereditary response to the earlier smooth pulse varies.
        panels.append(SourcePanel(-1., 1., np.zeros((1, 10, 3, 3))))
        history = PanelHistory(-2., (10, 3, 3), "generic-past-pulse", tuple(panels))
        actual = causal_canonical_nonlocal_luminosity(
            binary, position, velocity, history, time, infrared, current,
        )
        printed = canonical_nonlocal_luminosity(
            binary, history, time, position, velocity, binary.mass, infrared, current,
        )
        # The endpoint-free variation in the dynamic term contains a
        # factor of two. At zero radial speed the moving-scale derivative
        # vanishes. Direct integration by parts avoids cancellation of
        # high-degree pulse derivatives in the independent evaluation.
        def convolution(order):
            return (-1)**(order + 1) * factorial(order - 1) * quad(
                lambda x: (x * (1 - x))**9 / (time + 2 - x)**order, 0, 1,
                epsabs=1e-18, epsrel=2e-12,
            )[0]

        position_term = np.outer(position, position) @ shape.T
        velocity_term = (3 * np.outer(velocity, position) + np.outer(position, velocity)) @ shape.T
        expected = -16 * binary.G**3 * binary.mass * binary.mu / (25 * binary.c**13) * np.sum(
            current[3] * (position_term * convolution(8) + velocity_term * convolution(7))
        )
        self.assertGreater(abs(expected), 1e-17)
        np.testing.assert_allclose(actual - printed, expected, rtol=1e-4, atol=4e-17)

    def test_all_local_terms_through_35pn_against_published_scalar(self):
        binary = Binary(0.65, 0.35, G=1.06, c=1.4)
        for speed in (-0.09, 0, 0.09):
            position, velocity = [15., 0., 0.], [speed, .31, .01]
            dipole = binary.mu * np.cross(position, velocity)
            jets = canonical_waveform_jets(binary, position, velocity, 0.7, 1.2, {0: dipole}, 7)
            independently_contracted = canonical_luminosity_from_jets(binary, jets, max_grade=7)
            published = canonical_local_luminosity(binary, position, velocity, 0.7, 1.2)
            for grade, flux in independently_contracted.items():
                np.testing.assert_allclose(flux, published[grade], rtol=3e-12, atol=2e-23)

    def test_general_harmonic_source_inventory_closes_canonical_product(self):
        binary = Binary(0.64, 0.36, G=1.2, c=1.4)
        position = np.array([14, -2, 3], dtype=float)
        velocity = np.array([-0.04, 0.22, -0.09])
        moments = canonical_source_jets(binary, position, velocity, 0.8, 0.76)
        independent = quadrupole_tail_sources(binary, np.concatenate((position, velocity)))
        np.testing.assert_allclose(moments["I2"][0][[4, 5, 7]], independent, rtol=3e-14, atol=4e-19)
        for order in range(4):
            canonical = canonical_mass_quadrupole(binary, moments, order, max_grade=7)
            self.assertEqual(set(canonical), {0, 2, 4, 5, 6, 7})
            np.testing.assert_allclose(canonical[2], moments["I2"][2][order], rtol=1e-14)
            np.testing.assert_allclose(
                np.trace(canonical[7]), 0, atol=3e-16,
            )
        equal = canonical_source_jets(Binary(0.5, 0.5), position, velocity, 1, 1)
        for name in ("I3", "J2", "W1", "Y1"):
            np.testing.assert_array_equal(equal[name][0], 0)

    def test_leibniz_and_explicit_grading(self):
        binary = Binary(0.67, 0.33, G=1.1, c=1.6)
        basis = np.diag([2.0, -1.0, -1.0])
        time = 1.1
        # Set I = Q t^2 and W = t^3. The grade-5 contact term is
        # 4G/c^5 * (6 t * Q t^2 - 3 t^2 * 2 Q t) and cancels exactly.
        # With W^(grade 2) = t^4, the residual is 4G/c^5 * 4 Q t^4.
        moments = {
            "W": {0: np.array([time**3, 3*time**2, 6*time, 6., 0., 0.]),
                  2: np.array([time**4, 4*time**3, 12*time**2, 24*time, 24.])},
            "I2": {0: np.array([time**2, 2 * time, 2., 0., 0., 0., 0.])[:, None, None] * basis},
            "I3": {0: np.zeros((6, 3, 3, 3))},
            "X": {0: np.zeros(7)},
        }
        for name in ("W1", "Y1", "Z1"):
            moments[name] = {0: np.zeros((8, 3))}
        for name in ("W2", "Y2", "J2"):
            moments[name] = {0: np.zeros((8, 3, 3))}
        for derivative in (0, 1, 2):
            actual = canonical_mass_quadrupole(binary, moments, derivative)
            np.testing.assert_allclose(actual[5], 0, atol=3e-14)
            derivative_of_fourth_power = (time**4, 4*time**3, 12*time**2)[derivative]
            np.testing.assert_allclose(
                actual[7], (16 * binary.G / binary.c**5) * derivative_of_fourth_power * basis,
                atol=2e-13, rtol=1e-14,
            )

    def test_generated_leading_luminosity_on_an_eccentric_trajectory(self):
        binary = Binary(0.8, 0.2, G=1.2, c=1.4)
        initial = InitialState((32., -1., 0.), (-0.04, 0.17, 0.06))
        trajectory = evolve(binary, initial, None, (0, 15), Accuracy("newtonian"))
        times = np.linspace(0, 15, 21)
        position, velocity = trajectory.state(times)
        flux = canonical_local_luminosity(binary, position, velocity, 0.7, 1.1)
        np.testing.assert_allclose(flux[0], compute_radiation(trajectory, times).energy_flux_infinity, rtol=3e-14)
        self.assertEqual(set(flux), {0, 2, 4, 5, 6, 7, 8})
        negative = canonical_local_luminosity(binary, position, -velocity, 0.7, 1.1)
        for grade, value in flux.items():
            np.testing.assert_allclose(negative[grade], (-1)**grade * value, rtol=2e-14)

    def test_power_keeps_cross_terms_and_discards_out_of_grade_terms(self):
        binary = Binary(0.7, 0.3)
        basis = np.diag([2.0, -1.0, -1.0])
        canonical = {("mass", 2): {0: basis}}
        corrections = {("mass", 2): {3: 0.5 * basis, 5: -0.1 * basis}}
        flux = nonlinear_luminosity(binary, canonical, corrections)
        self.assertEqual(set(flux), {3, 5, 6, 8})
        np.testing.assert_allclose([flux[g] for g in (3, 5, 6, 8)], [1.2, -0.24, 0.3, -0.12])
