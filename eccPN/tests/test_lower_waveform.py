"""Development waveforms and continuation on a declared finite history."""

from dataclasses import replace
from math import factorial, prod
from pathlib import Path
from tempfile import TemporaryDirectory
import unittest
from unittest.mock import patch

import numpy as np
from scipy.optimize import brentq

from eccpn import Accuracy, Binary, InitialState, compute_radiation, evolve
from eccpn.gauge_moments import canonical_source_jets
from eccpn.hereditary import AcceptedHistory
from eccpn.history import FiniteOrbitalPast, PanelHistory, SourcePanel, TailInitialData
from eccpn.leading_modes import newtonian_position_series
from eccpn.luminosity import angular_momentum_flux, canonical_local_luminosity
from eccpn.memory import null_memory_news
from eccpn.multipoles import additional_local_quadrupole_grade_eight
from eccpn.projection import project_stf, stf_tensor
from eccpn.tails import quadrupole_tail_sources, radiative_tail
from eccpn.waveform_history import build_canonical_history, finite_quadrupole_grade_eight, finite_radiation_available


def independent_commutator_jet(binary, position, velocity):
    """Evaluate the QC2 commutator from a Cartesian Newtonian Taylor jet."""
    x = newtonian_position_series(binary,position,velocity,4)
    v = x[1:]*np.arange(1,5)[:,None]
    x = x[:4]

    def multiply(a,b):
        return np.convolve(a,b)[:4]

    def dot(a,b):
        return sum(multiply(first,second) for first,second in zip(a.T,b.T))

    r_squared = dot(x,x)
    change = r_squared/r_squared[0]
    change[0] = 0.

    def radial_power(half_power):
        powers = np.array([1.,0.,0.,0.])
        value = powers.copy()
        for order in range(1,4):
            powers = multiply(powers,change)
            value += prod(half_power-k for k in range(order))/factorial(order)*powers
        return value*r_squared[0]**half_power

    gm = binary.G*binary.mass
    xv = dot(x,v)
    scalar = dot(v,v) - 3*multiply(multiply(xv,xv),radial_power(-1)) - gm*radial_power(-0.5)
    amplitude = -4*binary.nu*binary.mu*gm**3/binary.c**8*multiply(radial_power(-1.5),scalar)
    tensor = np.stack([multiply(amplitude,multiply(first,second)) for first in x.T for second in x.T],axis=-1)
    return stf_tensor(tensor.reshape(4,3,3)*np.array([1,1,2,6])[:,None,None],2)


class LowerWaveform(unittest.TestCase):
    def test_public_4pn_quadrupole_and_total_at_eccentric_turning_points(self):
        binary = Binary(0.72,0.28,c=1.1)
        for radius,radial,tangential in ((22.,-0.0011,0.28),(53.,0.00011,0.085)):
            initial = InitialState((radius,0.,0.),(radial,tangential,0.))
            precursor = evolve(binary,initial,None,(-0.4,0),Accuracy("newtonian"))
            state = InitialState(*(tuple(row) for row in precursor.state(0)))
            current = binary.mu*np.cross(*precursor.state(0))
            past = FiniteOrbitalPast((-0.4,-0.2,0),precursor.state,{0: binary.mass,2: -0.001,4: -0.00002},
                                     {0: current,2: current*0.014},0.92,"turning-point-4pn")
            zero = (0.,0.,0.)
            model = Accuracy("sh-through-4.5PN+recoil",max_step=0.2,waveform_pn=4.,harmonic_scales=(0.91,0.91))
            path = evolve(binary,state,past.initial_data(binary,zero,zero,zero),(0.,1.4),model)
            def radial_numerator(t):
                return np.dot(*path.state(t))
            tp = brentq(radial_numerator,0.,1.4)
            times = np.array([tp/2,tp,(tp+1.4)/2])
            self.assertLess(radial_numerator(times[0])*radial_numerator(times[2]),0)
            actual = compute_radiation(path,times)
            lower_path = replace(path,accuracy=replace(model,waveform_pn=3.5))
            lower = compute_radiation(lower_path,times)
            self.assertEqual(actual.fixed_history_amplitude_pn,3.5)
            self.assertEqual(actual.fixed_history_quadrupole_pn,4.0)
            self.assertEqual(actual.fixed_history_luminosity_relative_pn,4.0)
            np.testing.assert_allclose(actual.energy_flux_infinity,actual.energy_flux_from_modes,rtol=3e-13)
            np.testing.assert_allclose(actual.h22,project_stf(binary,2,actual.radiative_mass_quadrupole)[2,2],
                                       rtol=3e-14)
            np.testing.assert_allclose(actual.news22,project_stf(binary,2,actual.radiative_mass_quadrupole_derivative)[2,2],
                                       rtol=3e-14)
            x,v = path.state(times)
            published = canonical_local_luminosity(binary,x,v,0.91,0.92)[8]
            np.testing.assert_allclose(actual.energy_flux_parts["canonical_local"][8],published,rtol=2e-14)
            self.assertGreater(np.max(np.abs(actual.energy_flux_infinity/lower.energy_flux_infinity-1)),1e-8)
            self.assertGreater(np.max(np.abs(actual.news22-lower.news22)),1e-12)
            self.assertGreater(np.max(np.abs(actual.energy_flux_parts["canonical_nonlocal_published"][8]
                                             -actual.energy_flux_parts["canonical_nonlocal_retarded"][8])),1e-17)

            def omit_commutator(binary,state_position,state_velocity,scale):
                return (additional_local_quadrupole_grade_eight(binary,state_position,state_velocity,scale)
                        - independent_commutator_jet(binary,state_position,state_velocity))

            with patch("eccpn.multipoles.additional_local_quadrupole_grade_eight",omit_commutator):
                removed = compute_radiation(path,tp)
            correction = independent_commutator_jet(binary,*path.state(tp))[2:4]
            expected_strain, expected_news = (project_stf(binary,2,field)[2,2] for field in correction)
            self.assertGreater(abs(expected_strain),1e-11)
            self.assertGreater(abs(expected_news),1e-13)
            np.testing.assert_allclose(actual.h22[1]-removed.h22,expected_strain,rtol=2e-7)
            np.testing.assert_allclose(actual.news22[1]-removed.news22,expected_news,rtol=2e-7)
            np.testing.assert_allclose(actual.energy_flux_infinity[1],removed.energy_flux_infinity)
            self.assertGreater(abs(removed.energy_flux_from_modes/removed.energy_flux_infinity-1),1e-7)

    def test_finite_35pn_all_modes_and_power(self):
        binary = Binary(0.71, 0.29)
        q, eccentricity = 26, 0.68
        peri_speed = np.sqrt(binary.G * binary.mass * (1 + eccentricity) / q)
        initial = InitialState((q, 0., 0.), (-0.022, peri_speed, 0.))
        past_path = evolve(binary, initial, None, (-0.7, 0), Accuracy("newtonian"))
        state = InitialState(*(tuple(row) for row in past_path.state(0)))
        dipole = binary.mu * np.cross(initial.position, initial.velocity)
        past = FiniteOrbitalPast(
            (-0.7, -0.35, 0), past_path.state, {0: binary.mass, 2: -0.0014, 4: 0.00005},
            {0: dipole, 2: dipole * 0.018}, 0.87, "complete-35pn-fixed-history",
        )
        zero = (0., 0., 0.)
        trajectory = evolve(binary, state, past.initial_data(binary, zero, zero, zero), (0, 0.5),
                            Accuracy("sh-through-4.5PN+recoil", max_step=0.21, waveform_pn=3.5))
        radiation = compute_radiation(trajectory, [0.13, 0.34, 0.46])
        self.assertEqual(radiation.fixed_history_amplitude_pn, 3.5)
        self.assertEqual(radiation.fixed_history_luminosity_relative_pn, 3.5)
        for mode in ((3, 2), (4, 3), (5, 5), (6, 5), (7, 7), (8, 7), (9, 9)):
            self.assertGreater(np.max(np.abs(radiation.strain[mode])), 0)
        np.testing.assert_allclose(radiation.energy_flux_infinity, radiation.energy_flux_from_modes, rtol=3e-12)
        history = trajectory._canonical_history
        panel = history.moments["M2"][0].panels[-1]
        cut = (panel.start + panel.end) / 2
        dt = (panel.end - panel.start) / 4096
        times = [cut - dt, cut, cut + dt]
        retained = compute_radiation(trajectory, times, accuracy=replace(trajectory.accuracy, retain_memory_excess=True))
        baseline = compute_radiation(trajectory, times)
        np.testing.assert_array_equal(retained.energy_flux_infinity, baseline.energy_flux_infinity)
        np.testing.assert_array_equal(retained.radiative_mass_quadrupole, baseline.radiative_mass_quadrupole)
        canonical_news = finite_radiation_available(binary, history, cut, *trajectory.state(cut), 1.0)[1]
        rates = null_memory_news(binary, canonical_news)
        self.assertEqual(max(ell for _, ell, _ in retained.retained_memory_news), 11)
        self.assertIn((8, 2, 0), retained.retained_memory_strain)
        for key, value in retained.retained_memory_strain.items():
            grade, _, _ = key
            self.assertGreater(grade, 7)
            expected = retained.retained_memory_news[key][1]
            np.testing.assert_allclose(expected, rates[key], rtol=3e-8, atol=3e-22)
            np.testing.assert_allclose((value[2] - value[0]) / (2 * dt), expected, rtol=3e-8, atol=3e-22)
        self.assertTrue(any(abs(values[1]) > 1e-20 for (grade, ell, m), values
                            in retained.retained_memory_news.items() if grade >= 9 and m != 0))

    def test_3pn_restart_replays_states_and_wave_history(self):
        binary = Binary(0.65, 0.35, c=1.52)
        initial = InitialState((20., 0., 0.), (0.08, 0.25, 0.0))
        prefix = evolve(binary, initial, None, (-0.4, 0), Accuracy("newtonian"))
        dipole = binary.mu * np.cross(initial.position, initial.velocity)
        past = FiniteOrbitalPast((-0.4, 0.), prefix.state, {0: binary.mass, 2: -0.002}, {0: dipole}, 1.0, "restart")
        state = InitialState(*(tuple(x) for x in prefix.state(0)))
        history = past.initial_data(binary, (0.,) * 3, (0.,) * 3, (0.,) * 3)
        accuracy = Accuracy("sh-through-4.5PN+recoil", max_step=0.08, waveform_pn=3)
        first = evolve(binary, state, history, (0., 0.46), accuracy)
        raw = first.checkpoint
        with TemporaryDirectory() as directory:
            path = Path(directory) / "finite-3pn.json"
            raw.save(path)
            restored = AcceptedHistory.load(path, identity="restart")
        state = InitialState(tuple(raw.state[:3]), tuple(raw.state[3:6]))
        resumed = evolve(binary, state, restored, (0.46, 0.78), accuracy)
        in_memory = evolve(binary, state, raw, (0.46, 0.78), accuracy)
        cuts = [0.52, 0.76]
        np.testing.assert_array_equal(resumed.state(cuts), in_memory.state(cuts))
        np.testing.assert_array_equal(resumed.recoil(cuts), in_memory.recoil(cuts))
        result, reference = compute_radiation(resumed, cuts), compute_radiation(in_memory, cuts)
        np.testing.assert_array_equal(result.h22, reference.h22)
        np.testing.assert_array_equal(result.news22, reference.news22)
        np.testing.assert_array_equal(result.energy_flux_infinity, reference.energy_flux_infinity)

    def test_grade_seven_assembly_preserves_lower_modes(self):
        binary = Binary(0.66, 0.34, G=1.05, c=1.4)
        prefix = evolve(
            binary, InitialState((17.0, 0.0, 4.0), (0.07, 0.23, -0.12)), None,
            (-1, 1), Accuracy("newtonian"),
        )
        position, velocity = prefix.state(0.53)
        dipole = binary.mu * np.cross(*prefix.state(-1))
        args = (binary, [-1, 0, 1], prefix.state, {0: binary.mass, 2: -0.005, 4: -0.00007},
                {0: dipole, 2: 0.01 * dipole}, 0.88, 0.97)
        lower = build_canonical_history(*args, maximum_absolute_grade=6)
        extended = build_canonical_history(*args, maximum_absolute_grade=7)
        reference = finite_radiation_available(binary, lower, 0.53, position, velocity, 1.17)
        result = finite_radiation_available(binary, extended, 0.53, position, velocity, 1.17)
        for before, after in zip(reference[:2], result[:2]):
            for key, coefficient in before.items():
                np.testing.assert_allclose(after[key], coefficient, atol=2e-16, rtol=2e-13)
            self.assertGreater(abs(after[7, 7, 1]), 0)
            self.assertGreater(abs(after[7, 8, 1]), 0)
            self.assertGreater(abs(after[7, 9, 1]), 0)
        torque = angular_momentum_flux(binary, result[2], 7)
        mode_z = {}
        for (a, ell, m), first in result[0].items():
            for b in range(8 - a):
                second = result[1].get((b, ell, m))
                if second is not None:
                    power = binary.c**3 / (16 * np.pi * binary.G) * (1j * m * np.conj(first) * second).real
                    mode_z[a + b] = mode_z.get(a + b, 0) + power
        for grade, component in torque.items():
            np.testing.assert_allclose(component[2], mode_z[grade], rtol=3e-12, atol=3e-22)
        parts = finite_quadrupole_grade_eight(binary, extended, 0.53, position, velocity, (0.88, 0.89), 1.17)
        self.assertEqual(set(parts), {"canonical_local", "canonical_nonlocal", "ordinary_tail", "tail_of_tail", "cubic"})
        for value in parts.values():
            self.assertGreater(np.linalg.norm(value), 0)
            np.testing.assert_allclose(value, value.swapaxes(-1, -2), atol=1e-16)
            np.testing.assert_allclose(np.trace(value, axis1=-1, axis2=-2), 0, atol=1e-16)

    def test_finite_3pn_waveform_power_recovers_independent_scalar(self):
        binary = Binary(0.61, 0.39, G=1.04, c=1.5)
        early = InitialState((19.0, 0.0, 0.0), (0.05, 0.25, 0.0))
        prefix = evolve(binary, early, None, (-2, 0), Accuracy("newtonian"))
        position, velocity = prefix.state(0)
        dipole = binary.mu * np.cross(early.position, early.velocity)
        past = FiniteOrbitalPast(
            (-2, -1, 0), prefix.state, {0: binary.mass, 2: -0.003}, {0: dipole},
            0.84, "shared-3pn-test-worldline",
        )
        zero = (0., 0., 0.)
        trajectory = evolve(
            binary, InitialState(tuple(position), tuple(velocity)), past.initial_data(binary, zero, zero, zero),
            (0, 1.3), Accuracy("sh-through-4.5PN+recoil", max_step=0.33, waveform_pn=3),
        )
        output = compute_radiation(trajectory, [0.11, 0.61, 1.17])
        self.assertEqual(output.fixed_history_amplitude_pn, 3)
        self.assertEqual(output.fixed_history_luminosity_relative_pn, 3)
        np.testing.assert_allclose(output.energy_flux_infinity, output.energy_flux_from_modes, rtol=2e-13)
        self.assertGreater(np.max(abs(output.strain[8, 8])), 0)

    def test_full_lower_order_modes_through_public_interface(self):
        binary = Binary(0.61, 0.39, G=1.04, c=1.5)
        initial = InitialState((19.0, 0.0, 0.0), (0.05, 0.25, 0.0))
        source = quadrupole_tail_sources(binary, np.array(initial.position + initial.velocity))
        prefix = PanelHistory(-1, (3, 3, 3), "1.5PN-finite-fixture", (
            SourcePanel(-1, 0, np.stack((0 * source, 0 * source, source))),
        ))
        zero = (0, 0, 0)
        history = TailInitialData(prefix, zero, zero, zero)
        accuracy = Accuracy("sh-through-4.5PN+recoil", max_step=0.2, waveform_pn=1.5)
        trajectory = evolve(binary, initial, history, (0, 0.9), accuracy)
        time = 0.71
        output = compute_radiation(trajectory, time)
        self.assertEqual(output.fixed_history_amplitude_pn, 1.5)
        self.assertEqual(set(output.strain), {(ell, m) for ell in range(2, 6) for m in range(-ell, ell + 1)})
        tensors = canonical_source_jets(binary, *trajectory.state(time), 1.0, 1.0)
        orbital_tail = radiative_tail(binary, time, trajectory.history, accuracy.radiative_scale)
        expect = tensors["I2"][0][2] + tensors["I2"][2][2] + orbital_tail[0]
        np.testing.assert_allclose(output.radiative_mass_quadrupole, expect, atol=2e-17)
        np.testing.assert_allclose(output.h22, project_stf(binary, 2, expect)[2, 2], rtol=5e-14)
        self.assertGreater(abs(output.strain[5, 5]), 0)
        for (ell, m), value in output.strain.items():
            np.testing.assert_allclose(output.strain[ell, -m], (-1)**ell * value.conj(), atol=2e-16)

    def test_circular_1pn_h22_coefficient(self):
        binary = Binary(0.64, 0.36, G=1.2, c=60)
        radius = 40
        potential = binary.G * binary.mass / radius
        gauge_radius = potential / binary.c**2
        speed = np.sqrt(potential * (1 + (-3 + binary.nu) * gauge_radius))
        phase_frequency = speed / radius
        invariant_frequency = (binary.G * binary.mass * phase_frequency / binary.c**3)**(2 / 3)
        source = canonical_source_jets(binary, [radius, 0, 0], [0, speed, 0], 1, 1)
        quadrupole = source["I2"][0][2] + source["I2"][2][2]
        waveform = project_stf(binary, 2, quadrupole)[2, 2]
        leading = -8 * np.sqrt(np.pi / 5) * binary.G * binary.mu / binary.c**2 * invariant_frequency
        measured = (waveform / leading - 1) / invariant_frequency
        np.testing.assert_allclose(measured, (-107 + 55 * binary.nu) / 42, atol=5e-5)
