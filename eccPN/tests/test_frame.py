"""Bondi retiming, spin phase and source coverage."""

import unittest

import numpy as np
from scipy.integrate import quad

from eccpn import Accuracy, Binary, BondiFrame, InitialState, compute_radiation, evolve, incoming_bondi_data
from eccpn.frame import angular_grid, retime_and_project
from eccpn.incoming import IncomingOrbit, ParabolicOrbit


class FrameTransformation(unittest.TestCase):
    def test_incoming_rest_mass_aspect_and_early_shear(self):
        binary = Binary(0.76, 0.24, G=1.3, c=1.2)
        incoming = IncomingOrbit((0.17, 0.27, -0.06), (-3., 6., 12.), -1000.)
        data = incoming_bondi_data(binary, incoming, ell_max=8)
        fine = incoming_bondi_data(binary, incoming, ell_max=8, polar_points=32, azimuthal_points=64)
        gamma = 1 / np.sqrt(1 - np.sum(data.betas**2, axis=1))
        lengths = np.asarray(data.mass_lengths)
        rays, weights = angular_grid(24, 48)
        aspect = data.mass_aspect(rays)
        self.assertGreater(np.ptp(aspect), 0.01)
        np.testing.assert_allclose(aspect @ weights / (4 * np.pi), lengths @ gamma, rtol=2e-14)
        np.testing.assert_allclose(rays.T @ (weights * aspect) / (4 * np.pi), 0, atol=3e-16)
        np.testing.assert_allclose(binary.c * (data.betas[0] - data.betas[1]), incoming.velocity, rtol=3e-14)
        np.testing.assert_allclose(list(data.shear_strain.values()), list(fine.shear_strain.values()),
                                   atol=6e-16, rtol=3e-12)
        np.testing.assert_allclose(data.frame.rotation @ data.betas[0], [0, np.linalg.norm(data.betas[0]), 0],
                                   atol=2e-17)

        def inertial_strain(times):
            return (
                {key: np.full_like(times, value, dtype=complex) for key, value in data.shear_strain.items()},
                {key: np.zeros_like(times) for key in data.shear_strain},
            )

        strain, news = retime_and_project(data.frame, np.array(-700.), binary.c, inertial_strain, data.shear_strain)
        np.testing.assert_allclose(list(strain.values()), 0, atol=3e-17)
        np.testing.assert_array_equal(list(news.values()), 0)

        parabolic = incoming_bondi_data(binary, ParabolicOrbit(30., (1., 0., 0.), (0., 1., 0.), -1000.))
        np.testing.assert_array_equal(list(parabolic.shear_strain.values()), 0)
        np.testing.assert_allclose(parabolic.mass_aspect(rays), np.sum(lengths), rtol=2e-15)

    def test_native_time_shift_and_rotation(self):
        binary = Binary(0.63, 0.37, G=1.03, c=1.4)
        initial = InitialState((18.0, 0, 0), (0.025, 0.305, 0))
        trajectory = evolve(binary, initial, None, (0, 6), Accuracy("newtonian"))
        angle, delay = 0.34, 0.7
        rotation = np.array([[np.cos(angle), -np.sin(angle), 0],
                             [np.sin(angle), np.cos(angle), 0], [0, 0, 1]])
        frame = BondiFrame(rotation=rotation, supertranslation={(0, 0): delay * binary.c * np.sqrt(4 * np.pi)})
        times = np.array([1.0, 2.2])
        result = compute_radiation(trajectory, times, frame)
        direct = compute_radiation(trajectory, times + delay)
        phase = np.exp(-2j * angle)
        np.testing.assert_allclose(result.strain[2, 2], phase * direct.h22, rtol=2e-13)
        np.testing.assert_allclose(result.news[2, 2], phase * direct.news22, rtol=2e-13)
        np.testing.assert_allclose(result.h22, compute_radiation(trajectory, times).h22, rtol=2e-13)
        with self.assertRaises(ValueError):
            compute_radiation(trajectory, [5.8], frame)

    def test_electric_shear_offset(self):
        def zero_source(times):
            return ({(2, 2): np.zeros_like(times)}, {(2, 2): np.zeros_like(times)})

        alpha20 = 0.16
        frame = BondiFrame(supertranslation={(2, 0): alpha20})
        strain, news = retime_and_project(frame, np.array([0.4]), 1.2, zero_source, [(2, 2)])
        np.testing.assert_allclose(strain[2, 0], -np.sqrt(24) * alpha20, rtol=2e-14)
        np.testing.assert_array_equal(news[2, 0], [0])

    def test_axial_boost_at_independently_integrated_cut(self):
        beta = 0.16
        frame = BondiFrame(boost=(0, 0, beta), polar_points=20, azimuthal_points=30)
        time = 0.8

        def source(times):
            return ({(2, 2): times**2}, {(2, 2): 2 * times})

        result = retime_and_project(frame, np.array(time), 1.0, source, [(2, 2)])

        def integrand(z, derivative):
            conformal = (1 + beta * z) / np.sqrt(1 - beta**2)
            source_z = (z + beta) / (1 + beta * z)
            amplitude = time**2 if derivative == 0 else 2 * time
            return (5 / 32) * amplitude / conformal**3 * (1 + source_z)**2 * (1 + z)**2

        for derivative, modes in enumerate(result):
            np.testing.assert_allclose(modes[2, 2], quad(integrand, -1, 1, args=(derivative,))[0], rtol=3e-14)
