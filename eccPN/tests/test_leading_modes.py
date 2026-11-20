"""Full leading multipole inventory and independent circular normalization."""

import unittest

import numpy as np

from eccpn import Accuracy, Binary, InitialState, compute_radiation, evolve
from eccpn.leading_modes import leading_mode_coefficients


class LeadingModes(unittest.TestCase):
    def test_quadrupole_and_odd_mass_circular_mode(self):
        binary = Binary(0.69, 0.31, G=1.3, c=1.4)
        radius = 24.0
        speed = np.sqrt(binary.G * binary.mass / radius)
        position, velocity = (radius, 0, 0), (0, speed, 0)
        trajectory = evolve(binary, InitialState(position, velocity), None, (0, 1), Accuracy("newtonian"))
        radiation = compute_radiation(trajectory, 0.0)
        strain, news = leading_mode_coefficients(binary, position, velocity)
        for m in range(-2, 3):
            np.testing.assert_allclose(strain[0, 2, m], radiation.strain[2, m], atol=2e-17, rtol=3e-14)
            np.testing.assert_allclose(news[0, 2, m], radiation.news[2, m], atol=2e-17, rtol=3e-14)
        delta = (binary.m1 - binary.m2) / binary.mass
        expected_h21 = -8j / 3 * np.sqrt(np.pi / 5) * binary.G * binary.mu * delta * speed**3 / binary.c**5
        np.testing.assert_allclose(strain[1, 2, 1], expected_h21, rtol=3e-14)
        self.assertIn((7, 9, 9), strain)
        self.assertIn((7, 8, 8), strain)

    def test_news_tracks_higher_leading_terms_on_an_eccentric_trajectory(self):
        binary = Binary(0.65, 0.35)
        initial = InitialState((15., -2., 0.), (0.07, 0.25, 0.))
        trajectory = evolve(binary, initial, None, (0, 3), Accuracy("newtonian"))
        delta = 1e-4
        earlier = leading_mode_coefficients(binary, *trajectory.state(1.0 - delta))[0]
        later = leading_mode_coefficients(binary, *trajectory.state(1.0 + delta))[0]
        news = leading_mode_coefficients(binary, *trajectory.state(1.0))[1]
        for mode in ((3, 5, 5), (5, 6, 5), (7, 9, 9), (7, 8, 7)):
            expected = (later[mode] - earlier[mode]) / (2 * delta)
            np.testing.assert_allclose(news[mode], expected, rtol=1e-7, atol=1e-15)
