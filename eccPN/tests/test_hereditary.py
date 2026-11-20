"""Committed history, checkpoint continuation and a finite-history tail fixture."""

from pathlib import Path
from tempfile import TemporaryDirectory
import unittest

import numpy as np
from scipy.integrate import quad

from eccpn import Accuracy, Binary, InitialState, compute_radiation, evolve
from eccpn.hereditary import AcceptedHistory, integrate_hereditary
from eccpn.history import PanelHistory, SourcePanel, TailInitialData
from eccpn.mpm import _LEVI_CIVITA
from eccpn.recoil import center_flux_rate, recoil_fluxes
from eccpn.tails import quadrupole_tail_sources
from eccpn.waveform_jets import canonical_waveform_jets


class HereditaryIntegration(unittest.TestCase):
    def test_pointwise_recoil_fluxes_against_independent_stf_products(self):
        binary = Binary(0.72, 0.28, G=1.09, c=1.37)
        position = np.array([13.2, -2.1, 1.3])
        for sign in (-1, 1):
            velocity = np.array([sign * 0.11, 0.29, -0.13])
            jets = canonical_waveform_jets(binary, position, velocity, 0.9, 1.4, {0: [0, 0, 1.]})
            m2, m3, m4 = (jets[name] for name in ("M2", "M3", "M4"))
            s2, s3 = (jets[name] for name in ("S2", "S3"))

            def leading(i, j):
                return (
                    2 / 63 * np.einsum("iab,ab->i", m3[i][4], m2[j][3])
                    + 16 / 45 * np.einsum("ijk,ja,ka->i", _LEVI_CIVITA, m2[j][3], s2[i][3])
                )

            prefix = binary.G / binary.c**7
            order_zero = prefix * leading(0, 0)
            order_two = prefix * (
                leading(0, 2) + leading(2, 0) + (
                    np.einsum("iabc,abc->i", m4[0][5], m3[0][4]) / 1134
                    + np.einsum("ijk,jab,kab->i", _LEVI_CIVITA, m3[0][4], s3[0][4]) / 126
                    + 4 / 63 * np.einsum("iab,ab->i", s3[0][4], s2[0][3])
                ) / binary.c**2
            )
            center = 2 / 21 * prefix * np.einsum("iab,ab->i", m3[0][3], m2[0][3])
            np.testing.assert_allclose(recoil_fluxes(binary, position, velocity), [order_zero, order_two, center],
                                       rtol=3e-12, atol=2e-23)
            step = 3e-4
            acceleration = -binary.G * binary.mass / np.linalg.norm(position)**3 * position
            previous = recoil_fluxes(binary, position - step * velocity, velocity - step * acceleration)[2]
            following = recoil_fluxes(binary, position + step * velocity, velocity + step * acceleration)[2]
            np.testing.assert_allclose(center_flux_rate(binary, position, velocity),
                                       (following - previous) / (2 * step), rtol=3e-8, atol=2e-18)

    def test_recoil_is_integrated_on_the_same_trajectory(self):
        binary = Binary(0.65, 0.35, c=1.9)
        initial = InitialState((22.0, 0, 0), (0.05, 0.21, 0))
        source = quadrupole_tail_sources(binary, np.array(initial.position + initial.velocity))
        prefix = PanelHistory(-1, (3, 3, 3), "finite-recoil-fixture", (
            SourcePanel(-1, 0, np.stack((0 * source, 0 * source, source))),
        ))
        zero = (0, 0, 0)
        prehistory = TailInitialData(prefix, zero, zero, zero)
        trajectory = evolve(
            binary, initial, prehistory, (0, 1),
            Accuracy("sh-through-4.5PN+recoil", max_step=0.125),
        )
        position, velocity = trajectory.state(0.5)
        expected = recoil_fluxes(binary, position, velocity)
        early = trajectory.recoil(0.5 - 1e-4)
        late = trajectory.recoil(0.5 + 1e-4)
        at_cut = trajectory.recoil(0.5)
        expected_rates = (expected[0], expected[1], at_cut[0] + at_cut[1] + expected[2])
        for before, after, rate in zip(early, late, expected_rates):
            np.testing.assert_allclose((after - before) / 2e-4, rate, rtol=1e-7, atol=1e-20)
        for flux in recoil_fluxes(Binary(0.5, 0.5), position, velocity):
            np.testing.assert_array_equal(flux, np.zeros(3))
        (first, first_rate), (second, second_rate) = trajectory.body_states(0.5)
        np.testing.assert_allclose(first - second, position, atol=2e-15)
        np.testing.assert_allclose(first_rate - second_rate, velocity, atol=2e-17)
        for side, rate in enumerate((first_rate, second_rate)):
            before = trajectory.body_states(0.5 - 1e-4)[side][0]
            after = trajectory.body_states(0.5 + 1e-4)[side][0]
            np.testing.assert_allclose((after - before) / 2e-4, rate, rtol=5e-7, atol=2e-11)

    def test_rejection_and_checkpoint_resume(self):
        def source(state):
            return state

        def rhs(time, state, history, trial):
            memory = history.integral(time, 1, trial=trial)
            return -10 * state + 0.4 * memory

        history = PanelHistory(0, (1,), "scalar-log-test")
        initial = AcceptedHistory(np.array([1.0]), history, 0.4)
        options = dict(rtol=2e-6, atol=2e-8, max_step=0.4)
        full = integrate_hereditary(rhs, source, initial, 0.4, **options)
        self.assertGreater(full.rejected_steps, 0)
        split = full.boundaries[4]
        first = integrate_hereditary(rhs, source, initial, split, **options)
        with TemporaryDirectory() as directory:
            path = Path(directory) / "checkpoint.json"
            first.checkpoint.save(path)
            with self.assertRaises(ValueError):
                AcceptedHistory.load(path, identity="incorrect")
            checkpoint = AcceptedHistory.load(path, identity=history.identity)
            continued = integrate_hereditary(rhs, source, checkpoint, 0.4, **options)
        np.testing.assert_allclose(
            continued(np.linspace(split, 0.4, 21)), full(np.linspace(split, 0.4, 21)),
            rtol=1e-5, atol=2e-8,
        )
        np.testing.assert_allclose(
            continued.checkpoint.history.integral(0.4, 1),
            full.checkpoint.history.integral(0.4, 1), rtol=1e-5,
        )
        self.assertEqual(history.panels, ())
        self.assertEqual(len(full.checkpoint.history.panels), len(full.boundaries) - 1)

    def test_shared_history_tail(self):
        binary = Binary(0.6, 0.4, c=1.8)
        initial = InitialState((18, 0, 0), (-0.035, 0.28, 0))
        state = np.array(initial.position + initial.velocity)
        joined = quadrupole_tail_sources(binary, state)
        # Artificial finite past with a smooth turn-on, matched in value at zero.
        prefix = PanelHistory(-7, (3, 3, 3), "finite-tail-fixture", (
            SourcePanel(-7, 0, np.stack((0 * joined, 0 * joined, joined))),
        ))
        trajectory = evolve(
            binary, initial, prefix, (0, 2),
            Accuracy("newtonian+2.5PN+4PN-tail", rtol=2e-8, max_step=0.5),
        )
        time = 1.75
        result = compute_radiation(trajectory, time)
        self.assertTrue(result.leading_quadrupole_tail)
        position, velocity = trajectory.state(time)
        radius = np.linalg.norm(position)
        acceleration = -binary.G * binary.mass * position / radius**3
        instant = 2 * binary.mu * (velocity[0] * velocity[1] + (
            position[0] * acceleration[1] + position[1] * acceleration[0]
        ) / 2)
        expected_tail = 0.0
        for panel in trajectory.history.panels:
            if panel.start >= time:
                break
            expected_tail += quad(
                lambda t: np.polynomial.polynomial.polyval(
                    (t - panel.start) / (panel.end - panel.start), panel.coefficients[:, 0, 0, 1],
                ) * (np.log((time - t) / 2) + 11 / 12),
                panel.start, min(panel.end, time), epsabs=1e-15,
            )[0]
        expected_tail *= 2 * binary.G * binary.mass / binary.c**3
        np.testing.assert_allclose(
            result.radiative_mass_quadrupole[0, 1] - instant, expected_tail,
            atol=2e-15, rtol=2e-8,
        )
