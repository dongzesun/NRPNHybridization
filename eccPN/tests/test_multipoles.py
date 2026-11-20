"""Tests for explicit PN grading and local SH source quadrupole derivatives."""

import unittest

import numpy as np

from eccpn.gauge_moments import canonical_source_jets
from eccpn.multipoles import known_local_quadrupole_grade_eight
from eccpn.radiative_quadrupole import _contract, _spin_contract

from eccpn import Accuracy, Binary, InitialState, compute_radiation, evolve
from eccpn.multipoles import quadrupole_infrared_commutator, source_mass_quadrupole


class SourceQuadrupole(unittest.TestCase):
    def test_generic_source_infrared_scale_coefficient(self):
        binary = Binary(0.62, 0.38, G=1.1, c=1.6)
        position = np.array([11., -0.7, 0.1])
        velocity = np.array([0.08, 0.3, -0.12])
        jets = canonical_source_jets(binary, position, velocity, 1.2, 0.9)["I2"]
        newtonian = jets[0]
        dipole = binary.mu * np.cross(position, velocity)
        mass_scale = binary.mu * (velocity @ velocity / 2 - binary.G * binary.mass / np.linalg.norm(position))
        mass_scale /= binary.c**2
        tail = (-214 / 105 * (binary.G * binary.mass)**2 / binary.c**6
                * (jets[2][2] + 2 * mass_scale / binary.mass * newtonian[2]))
        mixed = (
            20 / 7 * _contract(newtonian[0], newtonian[4])
            - 24 / 7 * _contract(newtonian[1], newtonian[3])
            + 12 / 7 * _contract(newtonian[2], newtonian[2])
            - 4 / 3 * _spin_contract(dipole, newtonian[3])
        ) * binary.G**2 * binary.mass / binary.c**8
        reference = known_local_quadrupole_grade_eight(binary, position, velocity, (1.2, 1.1), 0.9)[0]
        shifted = known_local_quadrupole_grade_eight(binary, position, velocity, (1.2, 1.1), 2.4)[0]
        np.testing.assert_allclose((shifted - reference) / np.log(2.4 / 0.9), tail + mixed,
                                   rtol=1e-12, atol=2e-16)

    def test_newtonian_jet(self):
        binary = Binary(0.68, 0.32, G=1.1, c=1.3)
        initial = InitialState((14.0, -2.0, 3.0), (0.03, 0.25, 0.16))
        trajectory = evolve(binary, initial, None, (0, 3), Accuracy("newtonian"))
        times = [0.0, 1.0, 3.0]
        position, velocity = trajectory.state(times)
        radiation = compute_radiation(trajectory, times)
        for derivative, expected in enumerate((
            radiation.radiative_mass_quadrupole, radiation.radiative_mass_quadrupole_derivative,
        ), start=2):
            result = source_mass_quadrupole(binary, position, velocity, derivative, 2.1, 1.4)
            self.assertEqual(set(result), {0, 2, 4, 5, 6, 7})
            np.testing.assert_allclose(result[0], expected, rtol=2e-14, atol=2e-17)
            for tensor in result.values():
                np.testing.assert_allclose(np.trace(tensor, axis1=-1, axis2=-2), 0, atol=1e-15)

    def test_first_pn_source_coefficient(self):
        binary = Binary(0.62, 0.38, G=1.07, c=1.31)
        position = np.array([27.0, 2.0, 6.0])
        velocity = np.array([0.06, 0.2, 0.15])
        radius = np.linalg.norm(position)
        p = binary.G * binary.mass / radius
        w = velocity @ velocity
        n = binary.nu
        radial_speed = position @ velocity / radius
        raw = (
            (w * (29 / 42 - 29 * n / 14) + p * (-5 / 7 + 8 * n / 7))
            * np.outer(position, position)
            + radius**2 * (11 / 21 - 11 * n / 7) * np.outer(velocity, velocity)
            + radius * radial_speed * (-2 / 7 + 6 * n / 7)
            * (np.outer(position, velocity) + np.outer(velocity, position))
        )
        expected = binary.mu / binary.c**2 * (raw - np.eye(3) * np.trace(raw) / 3)
        result = source_mass_quadrupole(binary, position, velocity, 0, 1.0, 1.0)
        np.testing.assert_allclose(result[2], expected, rtol=3e-14)

    def test_4pn_commutator_visible_only_in_generic_fixture(self):
        binary = Binary(0.66, 0.34)
        circular = [0., 0.25, 0.]
        np.testing.assert_array_equal(
            quadrupole_infrared_commutator(binary, [16., 0., 0.], circular), np.zeros((3, 3)),
        )
        radial = [0.1, 0.25, 0.]
        corrected = quadrupole_infrared_commutator(binary, [16., 0., 0.], radial)
        self.assertGreater(np.linalg.norm(corrected), 0)
