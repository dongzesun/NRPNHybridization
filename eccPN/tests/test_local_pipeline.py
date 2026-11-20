"""Analytic regressions through the two public calls."""

import unittest

import numpy as np

from eccpn import Accuracy, Binary, InitialState, compute_radiation, evolve
from eccpn.dynamics import acceleration


class LocalPipeline(unittest.TestCase):
    binary = Binary(0.7, 0.3, G=1.1, c=1.3)
    eccentricity = 0.55
    semimajor = 42.0
    gm = binary.G * binary.mass
    mean_motion = np.sqrt(gm / semimajor**3)
    period = 2 * np.pi / mean_motion
    initial = InitialState(
        (semimajor * (1 - eccentricity), 0, 0),
        (0, np.sqrt(gm / semimajor * (1 + eccentricity) / (1 - eccentricity)), 0),
    )

    def test_eccentric_orbit_and_waveform(self):
        trajectory = evolve(
            self.binary, self.initial, None, (0, self.period),
            Accuracy("newtonian", rtol=2e-12, atol=2e-13, max_step=self.period / 200),
        )
        eccentric_anomaly = np.linspace(0, 2 * np.pi, 13)
        eccentricity = self.eccentricity
        times = (eccentric_anomaly - eccentricity * np.sin(eccentric_anomaly)) / self.mean_motion
        times[-1] = self.period
        radius = self.semimajor * (1 - eccentricity * np.cos(eccentric_anomaly))
        x = self.semimajor * (np.cos(eccentric_anomaly) - eccentricity)
        y = self.semimajor * np.sqrt(1 - eccentricity**2) * np.sin(eccentric_anomaly)
        azimuth = np.arctan2(y, x)
        radial_speed = np.sqrt(self.gm * self.semimajor) * eccentricity * np.sin(eccentric_anomaly) / radius
        angular_momentum = np.sqrt(self.gm * self.semimajor * (1 - eccentricity**2))
        tangential_speed = angular_momentum / radius
        waveform = compute_radiation(trajectory, times)
        position, _ = trajectory.state(times)
        np.testing.assert_allclose(position[:, 0], x, atol=1e-8)
        np.testing.assert_allclose(position[:, 1], y, atol=1e-8)
        expected = (
            4 * self.binary.G * self.binary.mu / self.binary.c**4 * np.sqrt(np.pi / 5)
            * np.exp(-2j * azimuth)
            * ((radial_speed - 1j * tangential_speed)**2 - self.gm / radius)
        )
        np.testing.assert_allclose(waveform.h22, expected, rtol=2e-10)
        derotated = waveform.h22 * np.exp(2j * azimuth)
        np.testing.assert_allclose(derotated[2], derotated[-3].conjugate(), rtol=2e-10)
        self._check_news(trajectory)

    def test_harmonic_balance_and_news(self):
        # Tilt the ellipse to exercise the full native STF projection.
        tangential_speed = self.initial.velocity[1]
        initial = InitialState(self.initial.position, (0, 0.8 * tangential_speed, 0.6 * tangential_speed))
        trajectory = evolve(
            self.binary, initial, None, (0, self.period),
            Accuracy("newtonian+2.5PN", max_step=self.period / 100),
        )
        times = np.linspace(0, self.period, 101)
        position, velocity = trajectory.state(times)
        _, reaction = acceleration(self.binary, position, velocity, trajectory.accuracy.dynamics)
        mechanical_work = self.binary.mu * np.sum(velocity * reaction, axis=-1)
        mechanical_torque = self.binary.mu * np.cross(position, reaction)

        # Differentiate the 2.5PN Schott energy and angular momentum using
        # the Newtonian radial acceleration at each point on the orbit.
        radius = np.linalg.norm(position, axis=-1)
        speed_squared = np.sum(velocity**2, axis=-1)
        radial_speed = np.sum(position * velocity, axis=-1) / radius
        radial_acceleration = (speed_squared - radial_speed**2 - self.gm / radius) / radius
        factor = 8 * self.binary.mu * self.binary.nu * self.gm**2 / (5 * self.binary.c**5)
        schott_factor = factor * radial_speed / radius**2
        schott_factor_rate = factor * (radial_acceleration - 2 * radial_speed**2 / radius) / radius**2
        speed_squared_rate = -2 * self.gm * radial_speed / radius**2
        schott_work = schott_factor_rate * speed_squared + schott_factor * speed_squared_rate
        schott_torque = -schott_factor_rate[:, None] * np.cross(position, velocity)
        radiation = compute_radiation(trajectory, times)
        power = radiation.energy_flux_infinity
        torque = radiation.angular_momentum_flux_infinity
        self.assertTrue((power > 0).all())
        np.testing.assert_allclose(
            -(mechanical_work + schott_work), power, rtol=4e-14, atol=1e-23,
        )
        angular_residual = mechanical_torque + schott_torque + torque
        relative_angular_error = np.linalg.norm(angular_residual, axis=-1) / np.linalg.norm(torque, axis=-1)
        self.assertLess(relative_angular_error.max(), 4e-14)
        self._check_news(trajectory)

    def _check_news(self, trajectory):
        times = self.period * np.array([0.02, 0.27, 0.5, 0.73, 0.98])
        delta = 1e-6 * self.period
        earlier = compute_radiation(trajectory, times - delta)
        later = compute_radiation(trajectory, times + delta)
        central = compute_radiation(trajectory, times)
        for mode, derivative in central.news.items():
            numerical = (later.strain[mode] - earlier.strain[mode]) / (2 * delta)
            np.testing.assert_allclose(derivative, numerical, atol=2e-11, rtol=1e-7)
        if trajectory.accuracy.dynamics == "newtonian":
            news_squared = sum(np.abs(mode)**2 for mode in central.news.values())
            luminosity = self.binary.c**3 / (16 * np.pi * self.binary.G) * news_squared
            np.testing.assert_allclose(luminosity, central.energy_flux_infinity, rtol=2e-14)

    def test_full_request_and_domain_errors(self):
        with self.assertRaises(NotImplementedError):
            evolve(self.binary, self.initial, None, (0, 1))
        trajectory = evolve(
            self.binary, self.initial, None, (0, 1), Accuracy("newtonian"),
        )
        with self.assertRaises(ValueError):
            compute_radiation(trajectory, [1.1])
        with self.assertRaises(NotImplementedError):
            compute_radiation(trajectory, [0.5], frame="boosted")
