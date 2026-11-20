"""Incoming quadrupole asymptote and analytic lag integrations."""

from dataclasses import replace
from pathlib import Path
from tempfile import TemporaryDirectory
import unittest
import warnings

import numpy as np
from scipy.integrate import IntegrationWarning, quad
from mpmath import mp

from eccpn import Accuracy, Binary, InitialState, compute_radiation, evolve
from eccpn.asymptotic_waveform import AsymptoticQuadrupole
from eccpn.hereditary import AcceptedHistory
from eccpn.incoming import IncomingOrbit, IncomingQuadrupole, ParabolicOrbit, _newtonian_hyperbolic, match_incoming_state, remote_log_moment
from eccpn.tails import quadrupole_tail_sources


class IncomingHistory(unittest.TestCase):
    def test_exact_newtonian_remainder_and_sampled_force_checkpoint(self):
        weak = Binary(0.62,0.38,c=1e6)
        incoming = IncomingOrbit((0.02,0.23,0.08),(29.,4.,-2.),-6000.,asymptotic_order=3)
        x,v = incoming.state_at(weak,incoming.cutoff)
        np.testing.assert_allclose(v@v/2-weak.G*weak.mass/np.linalg.norm(x),np.dot(incoming.velocity,incoming.velocity)/2,
                                   rtol=2e-13)
        np.testing.assert_allclose(np.cross(x,v),np.cross(incoming.impact,incoming.velocity),rtol=2e-12)
        x_exact,v_exact = _newtonian_hyperbolic(weak,np.asarray(incoming.velocity),np.asarray(incoming.impact),
                                               -incoming.cutoff,incoming.log_reference)
        np.testing.assert_allclose(x,x_exact,rtol=2e-14)
        np.testing.assert_allclose(v,v_exact,rtol=2e-14)
        binary = Binary(0.62,0.38,c=1.3)
        initial,data = incoming.prepare(binary)
        accuracy = Accuracy("sh-through-4.5PN+recoil",rtol=2e-8,max_step=1.7,
                            waveform_pn=4.0,incoming_waveform_nodes=12)
        first = evolve(binary,initial,data,(incoming.cutoff,incoming.cutoff+12),accuracy)
        with TemporaryDirectory() as directory:
            filename = Path(directory)/"sampled-prefix.json"
            first.checkpoint.save(filename)
            restored = AcceptedHistory.load(filename,identity=data.quadrupole.identity)
        for k in range(3):
            np.testing.assert_array_equal(
                first.checkpoint.history.integral(first.time_span[1],1.1,power=k),
                restored.history.integral(first.time_span[1],1.1,power=k),
            )
        self.assertGreater(np.linalg.norm(data.recoil_state()),0.)
        result = compute_radiation(first,incoming.cutoff+7.)
        self.assertEqual(result.incoming_asymptotic_order,3)
        self.assertEqual(result.fixed_history_quadrupole_pn,4)
        np.testing.assert_allclose(result.energy_flux_infinity,result.energy_flux_from_modes,rtol=3e-13)

    def test_fractional_remote_log_moments_near_join(self):
        cutoff, origin, source_scale, lag_scale = -24000., 0., 1.7, 2.3
        for t, alpha, p, k in ((cutoff + .73, 10 / 3, 0, 2),
                               (cutoff + 4e-8, 7 / 3, 2, 2),
                               (cutoff + .19, 13 / 3, 1, 1)):
            with mp.workdps(55):
                D = mp.mpf(origin) - cutoff
                d = (mp.mpf(t) - cutoff) / D
                def integrand(q):
                    return (q**(mp.mpf(alpha) - 2) * mp.log(D / (source_scale * q))**p
                            * mp.log(D * ((1 - q) + d * q) / (lag_scale * q))**k)

                exact = D**(1 - mp.mpf(alpha)) * mp.quad(integrand, [0, .5, 1 - d, 1])
            with warnings.catch_warnings():
                warnings.simplefilter("error", IntegrationWarning)
                result = remote_log_moment(cutoff, origin, t, source_scale, lag_scale, alpha, p, k)
            np.testing.assert_allclose(result, float(exact), rtol=2e-13, atol=0.)

    def test_parabolic_conservative_correction_and_waveform_charge(self):
        binary = Binary(0.6, 0.4, c=1.4)
        radius, peri, tangent = 9., (1., 0., 0.), (0., 1., 0.)
        for duration in (15000., 60000.):
            incoming = ParabolicOrbit(radius, peri, tangent, -duration)
            corrected = replace(incoming, asymptotic_order=2)

            def energy(data):
                position, velocity = data.state_at(binary, -1.1*duration)
                r = np.linalg.norm(position)
                u, w, s, nu = binary.G*binary.mass/r, velocity @ velocity, position @ velocity/r, binary.nu
                first = 3/8*(1-3*nu)*w*w + u*(nu*s*s+(3+nu)*w)/2 + u*u/2
                return w/2-u+first/binary.c**2

            self.assertLess(abs(energy(corrected)), abs(energy(incoming))/100)
            cut = -1.1*duration
            position, velocity = corrected.state_at(binary, cut)
            left, _ = corrected.state_at(binary, cut-0.1)
            right, _ = corrected.state_at(binary, cut+0.1)
            np.testing.assert_allclose((right-left)/0.2, velocity, rtol=2e-9, atol=1e-11)
            self.assertGreater(np.linalg.norm(position), 100.)
            initial, _ = incoming.prepare(binary)
            actual = quadrupole_tail_sources(binary, np.concatenate((initial.position, initial.velocity)))
            source = AsymptoticQuadrupole.parabolic_newtonian(binary, incoming)
            approximate = np.stack([source.value(-duration, order) for order in (4, 5, 7)])
            self.assertLess(np.linalg.norm(approximate[0]-actual[0])/np.linalg.norm(actual[0]), 8e-4)
        masses, angular = corrected.conserved_charges(binary)
        self.assertEqual(sum(masses.values()), binary.mass)
        np.testing.assert_allclose(angular[0], np.array([0., 0., binary.mu*np.sqrt(2*radius)]))

    def test_analytic_acceleration_of_both_incoming_asymptotes(self):
        binary = Binary(0.62,0.38,c=1.3)
        hyperbolic = IncomingOrbit((0.02,0.23,0.08),(29.,4.,-2.),-4000.)
        parabolic = ParabolicOrbit(21.,(0.6,0.,0.8),(0.,1.,0.),-4000.)
        cut,step = -7250.,2.4
        for template,orders in ((hyperbolic,(1,2,3)),(parabolic,(1,2))):
            for order in orders:
                with self.subTest(orbit=type(template).__name__,order=order):
                    orbit = replace(template,asymptotic_order=order)
                    position,rate,acceleration = orbit.state_at(binary,cut,include_acceleration=True)
                    np.testing.assert_array_equal(orbit.state_at(binary,cut),(position,rate))
                    samples = np.stack([orbit.state_at(binary,cut+k*step) for k in (-2,-1,1,2)])
                    derivatives = np.tensordot([1.,-8.,8.,-1.],samples,axes=1)/(12*step)
                    np.testing.assert_allclose(derivatives[0],rate,rtol=2e-10,atol=1e-13)
                    np.testing.assert_allclose(derivatives[1],acceleration,rtol=2e-9,atol=1e-15)

    def test_high_waveform_uses_remote_hyperbolic_mode_history(self):
        binary = Binary(0.63, 0.37, c=1.5)
        orbit = IncomingOrbit((0.23, -0.11, 0.07), (0., 31., 0.), -1600.)
        initial, data = orbit.prepare(binary)
        accuracy = Accuracy("sh-through-4.5PN+recoil", rtol=2e-9, max_step=10,
                            waveform_pn=3.5, incoming_waveform_nodes=24)
        path = evolve(binary, initial, data, (-1600., -1550.), accuracy)
        history = path._canonical_history
        refined = replace(path, _canonical_history=replace(history, incoming=replace(history.incoming, order=48)))
        cut = -1558.
        baseline = compute_radiation(path, cut)
        fine = compute_radiation(refined, cut)
        self.assertEqual(baseline.incoming_asymptotic_order, 2)
        self.assertEqual(baseline.incoming_quadrature_nodes, 24)
        self.assertEqual(fine.incoming_quadrature_nodes, 48)
        np.testing.assert_allclose(fine.news22, baseline.news22, rtol=1e-10)
        np.testing.assert_allclose(fine.energy_flux_infinity, fine.energy_flux_from_modes, rtol=3e-13)
        matched = history.covered_moments()["M2"][0]
        remote = matched.integrate_jet(cut, 2., 1, 4) - matched.finite.integrate_jet(cut, 2., 1, 4)
        self.assertGreater(np.linalg.norm(remote), np.linalg.norm(matched.finite.integrate_jet(cut, 2., 1, 4)))
        recorded = history.from_record(history.to_record())
        replay = compute_radiation(replace(path, _canonical_history=recorded), cut)
        np.testing.assert_array_equal(replay.h22, baseline.h22)
        np.testing.assert_array_equal(replay.news22, baseline.news22)

    def test_parabolic_branch_approaches_the_fractional_power_asymptote(self):
        binary = Binary(0.6, 0.4, c=1.7)
        peri, tangent = (1., 0., 0.), (0., 1., 0.)
        errors = []
        for cutoff in (-15000., -120000.):
            initial, past = ParabolicOrbit(9., peri, tangent, cutoff).prepare(binary)
            errors.append(np.max(past.asymptotic_join_residual))
            kinetic = np.sum(np.array(initial.velocity)**2) / 2
            potential = binary.G * binary.mass / np.linalg.norm(initial.position)
            np.testing.assert_allclose(kinetic, potential, rtol=5e-15)
        self.assertLess(errors[1], 0.7 * errors[0])
        path = evolve(
            binary, initial, past, (cutoff, cutoff + 1.0),
            Accuracy("sh-through-4.5PN+recoil", rtol=1e-7, max_step=0.3, waveform_pn=1.5),
        )
        self.assertTrue(np.isfinite(compute_radiation(path, cutoff + 0.6).news22))

    def test_incoming_charges_match_special_relativistic_energy_and_angular_momentum(self):
        from scipy.optimize import brentq

        orbit = IncomingOrbit((0.057, -0.033, 0.01), (1., 2., 3.), -20000.)
        binary = Binary(0.62, 0.38)
        w = np.linalg.norm(orbit.velocity) / binary.c
        mass, dipole = orbit.conserved_charges(binary)
        rest = np.array([binary.m1, binary.m2])
        p = brentq(lambda p: np.sum(p / np.hypot(rest, p)) - w, 0, binary.mass)
        energy = np.sum(np.hypot(rest, p))
        self.assertLess(abs(sum(mass.values()) - energy), 1e-8)
        exact = p / np.linalg.norm(orbit.velocity) * np.cross(orbit.impact, orbit.velocity)
        np.testing.assert_allclose(dipole[0] + dipole[2], exact, rtol=3e-6)

    def test_subleading_terms_control_quadrupole_join(self):
        binary = Binary(0.59, 0.41, c=1.4)
        velocity, impact = (0.32, -0.14, 0.08), (0., 12., 0.)
        errors = [
            [max(IncomingOrbit(velocity, impact, -duration, asymptotic_order=p)
                 .prepare(binary)[1].asymptotic_join_residual) for p in (1, 2)]
            for duration in (4800., 19200., 76800.)
        ]
        leading, subleading = np.asarray(errors).T
        self.assertLess(subleading[-1], leading[-1] / 100)
        self.assertLess(subleading[-1], subleading[-2] / 10)

    def test_forward_shooting_recovers_hyperbolic_target(self):
        binary = Binary(0.62, 0.38, c=1.6)
        actual = IncomingOrbit((0.31, -0.13, 0.06), (1., 14., -2.), -8000.0)
        guess = IncomingOrbit((0.309, -0.131, 0.061), (1.2, 13.8, -1.9), -8000.0)
        initial, past = actual.prepare(binary)
        accuracy = Accuracy("sh-through-4.5PN+recoil", max_step=6.0, rtol=2e-8)
        end = -7967.0
        forward = evolve(binary, initial, past, (actual.cutoff, end), accuracy)
        target = InitialState(*(tuple(row) for row in forward.state(end)))
        match = match_incoming_state(binary, guess, target, end, accuracy, target_tolerance=2e-9)
        self.assertLess(match.scaled_target_residual, 2e-9)
        np.testing.assert_allclose(match.incoming.velocity, actual.velocity, atol=2e-9)
        np.testing.assert_allclose(match.incoming.impact, actual.impact, atol=3e-5)
        self.assertLess(match.maximum_compactness, 0.001)
        self.assertLess(match.maximum_speed_ratio, 0.23)

    def test_near_cutoff_on_a_retarded_incoming_trajectory(self):
        binary = Binary(0.59, 0.41, c=1.4)
        velocity = (0.32, -0.14, 0.08)
        cut = -1200
        initial, prehistory = IncomingOrbit(velocity, (0.0, 12.0, 0.0), cut).prepare(binary)
        path = evolve(
            binary, initial, prehistory, (cut, cut + 8),
            Accuracy("sh-through-4.5PN+recoil", rtol=1e-7, max_step=2),
        )
        self.assertEqual(path.asymptotic_join_residual, prehistory.asymptotic_join_residual)
        self.assertIs(path.history.incoming, prehistory.quadrupole.incoming)
        output = compute_radiation(path, [cut + 3, cut + 6])
        self.assertTrue(np.isfinite(output.h22).all())
        momenta = path.recoil(cut)
        for actual, expected in zip(momenta, prehistory.recoil_state().reshape(3, 3)):
            np.testing.assert_array_equal(actual, expected)

    def test_leading_quadrupole_decay(self):
        binary = Binary(0.59, 0.41, G=1.02)
        velocity = np.array([0.25, -0.12, 0.08])
        relative = []
        for duration in (3e7, 6e7):
            prefix = IncomingQuadrupole.from_velocity(binary, -duration, velocity)
            position = -duration * velocity + np.array([0.7, 1.5, -2.0])
            actual = quadrupole_tail_sources(binary, np.concatenate((position, velocity)))
            error = np.linalg.norm(actual - prefix.at_cutoff(), axis=(-2, -1))
            relative.append(error / np.linalg.norm(actual, axis=(-2, -1)))
        np.testing.assert_allclose(relative[1] / relative[0], 0.5, rtol=2e-4)
        self.assertLess(relative[1].max(), 3e-6)

    def test_prefix_convolution_across_cutoff_and_at_late_cuts(self):
        binary = Binary(0.59, 0.41)
        prefix = IncomingQuadrupole.from_velocity(binary, -21.0, [0.25, -0.12, 0.08], time_origin=2.0)
        scale = 4.0
        duration = 23.0
        for time in (-21.0, -20.9999, -14.1, 2.0, 5.8, 1000.0):
            for power in (0, 1, 2):
                result = prefix.integral(time, scale, power)
                for index, (decay, factor) in enumerate(((3, 4), (4, 12), (6, 240))):
                    def integrand(lag):
                        return (duration + lag)**-decay * np.log((time + 21 + lag) / scale)**power

                    expected = quad(integrand, 0, np.inf, epsabs=2e-15, epsrel=2e-11)[0]
                    expected *= factor * prefix.tensor[0, 1]
                    np.testing.assert_allclose(result[index, 0, 1], expected, rtol=7e-8, atol=3e-15)
