"""Angular waveform flux and leading null memory without orbital averaging."""

import unittest

import numpy as np
from scipy.integrate import quad

from eccpn import Accuracy, Binary, InitialState, compute_radiation, evolve
from eccpn.memory import _memory_weights, angular_news_products, integrate_null_memory_panel, memory_boundary, memory_rate_polynomials, null_memory_news
from eccpn.mpm import formal_memory_rates
from eccpn.projection import project_stf
from eccpn.waveform_jets import canonical_waveform_jets


class NullMemory(unittest.TestCase):
    def test_logarithmic_finite_start_weights(self):
        for start, stop, origin, fraction in ((1.2, 1.7, 1.2, 0.37), (2.1, 3.7, 1.2, 0.73)):
            scale = 2.7
            weights = _memory_weights(start, stop, origin, scale, fraction)
            for n in range(5):
                for p in range(3):
                    width = stop - start
                    value, _ = quad(
                        lambda q: width * q**n * np.log((start - origin + width*q)/scale)**p,
                        0., fraction, epsabs=1e-13, epsrel=1e-13,
                    )
                    np.testing.assert_allclose(weights[n, p], value, atol=8e-15, rtol=2e-12)

    def test_formal_mpm_endpoint_matches_raw_on_generic_state(self):
        binary = Binary(0.69, 0.31, G=1.05, c=1.4)
        position, velocity = [13.0, -1.0, 2.0], [-0.03, 0.26, 0.1]
        jets = canonical_waveform_jets(binary, position, velocity, 1.2, 1.5, {0: [0, 0, 0.8]})
        news = {}
        for kind in ("mass", "current"):
            for rank in range(2, 5):
                offset = rank - 2 + (kind == "current")
                name = ("M" if kind == "mass" else "S") + str(rank)
                for grade, values in jets[name].items():
                    if grade + offset > 2:
                        continue
                    for (ell, m), mode in project_stf(binary, rank, values[rank + 1], kind).items():
                        news[grade + offset, ell, m] = news.get((grade + offset, ell, m), 0) + mode
        raw = null_memory_news(binary, news, max_grade=2)
        _, boundary_rate = memory_boundary(binary, jets)
        mpm = {}
        for rank in range(2, 7):
            for grade, tensor in formal_memory_rates(binary, f"u{rank}", jets).items():
                for (ell, m), mode in project_stf(binary, rank, tensor).items():
                    mpm[grade + rank - 2, ell, m] = mode
        for key in set(mpm) | set(raw):
            np.testing.assert_allclose(
                raw.get(key, 0) - boundary_rate.get(key, 0), mpm.get(key, 0),
                rtol=3e-11, atol=6e-20, err_msg=f"Memory rate {key}",
            )

    def test_complex_polynomial_memory_panel(self):
        binary = Binary(0.63, 0.37)
        news = {(0, 2, 2): np.array([3 + 2j, -1j, 0.3]), (0, 2, -2): np.array([1j, 2, 0])}
        duration, fraction = 2.1, 0.77
        integrated = integrate_null_memory_panel(binary, news, duration, fraction)
        nodes, weights = np.polynomial.legendre.leggauss(5)
        nodes = fraction * (1 + nodes) / 2
        sampled = {key: np.polynomial.polynomial.polyval(nodes, field) for key, field in news.items()}
        expected = null_memory_news(binary, sampled)
        for key in expected:
            np.testing.assert_allclose(
                integrated[key], duration * fraction / 2 * (expected[key] @ weights), rtol=3e-14, atol=6e-15,
            )

    def test_log_memory_polynomial_matches_sampled_angular_products(self):
        binary = Binary(0.63, 0.37)
        news = {
            (0, 2, 2): np.array([3 - 2j, 1, -1j]),
            (1, 3, 1): np.array([1, 2 - 1j, 1j]),
            (3, 2, -2): np.array([[1 - 2j, -1j], [2, 0.3], [-1, 0.7j]]),
            (4, 3, -1): np.array([[0.3j, 0.2], [-1, 0.5j], [1, -1j]]),
        }
        polynomials = memory_rate_polynomials(binary, news)
        for time in (0.19, 0.67, 0.99):
            log = np.log(time)
            samples = {
                key: (np.polynomial.polynomial.polyval(log, np.polynomial.polynomial.polyval(time, field))
                      if field.ndim == 2 else np.polynomial.polynomial.polyval(time, field))
                for key, field in news.items()
            }
            direct = null_memory_news(binary, samples)
            self.assertEqual(set(polynomials), set(direct))
            for key, series in polynomials.items():
                evaluated = np.polynomial.polynomial.polyval(
                    log, np.polynomial.polynomial.polyval(time, series),
                )
                np.testing.assert_allclose(evaluated, direct[key], rtol=8e-14, atol=5e-16)

    def test_circular_normalization_and_oscillatory_mode(self):
        binary = Binary(0.63, 0.37, G=1.13, c=1.47)
        radius = 24.0
        speed = np.sqrt(binary.G * binary.mass / radius)
        initial = InitialState((radius, 0, 0), (0, speed, 0))
        trajectory = evolve(binary, initial, None, (0, 1), Accuracy("newtonian"))
        radiation = compute_radiation(trajectory, [0.0, 0.5])
        news = {(0, ell, m): value for (ell, m), value in radiation.news.items()}
        rate = null_memory_news(binary, news)
        x = speed**2 / binary.c**2
        chirp = (64 / 5) * binary.c**3 / (binary.G * binary.mass) * binary.nu * x**5
        expected = (4 / 7) * np.sqrt(5 * np.pi / 6) * binary.G * binary.mu / binary.c**2 * chirp
        np.testing.assert_allclose(rate[5, 2, 0], expected, rtol=5e-13)

        z, weights = np.polynomial.legendre.leggauss(10)
        azimuth = np.arange(24) * (2 * np.pi / 24)
        phase = np.exp(2j * azimuth)
        normalization = np.sqrt(5 / (64 * np.pi))
        angular_news = normalization * (
            radiation.news[(2, 2)][0] * (1 + z[:, None])**2 * phase
            + radiation.news[(2, -2)][0] * (1 - z[:, None])**2 * phase.conj()
        )
        harmonic_44 = (3 / 16) * np.sqrt(35 / (2 * np.pi)) * (1 - z[:, None]**2)**2 * phase**2
        projected = np.sum(weights[:, None] * np.abs(angular_news)**2 * harmonic_44.conj()) * (2 * np.pi / 24)
        actual = angular_news_products(news)[0, 4, 4][0]
        self.assertGreater(abs(projected), 0.0)
        np.testing.assert_allclose(actual, projected, rtol=2e-14)

    def test_grade_cutoff(self):
        source = {(1, 3, 3): 1 - 1j, (7, 9, 9): 1 + 1j}
        products = angular_news_products(source)
        self.assertTrue(products)
        self.assertEqual({grade for grade, _, _ in products}, {2})
