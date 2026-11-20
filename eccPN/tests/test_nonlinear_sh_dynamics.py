"""Exact-grade and omission regressions at the staged public entry points."""

from dataclasses import replace
from pathlib import Path
from tempfile import TemporaryDirectory
import unittest

import numpy as np

from eccpn import Accuracy, Binary, InitialState, compute_radiation, evolve
from eccpn.conventions import NONLINEAR_5PN_SECTORS, SH_LOCAL_TAIL_MODEL, SH_RADIATIVE_5PN_MODEL
from eccpn.frame import BondiFrame
from eccpn.hereditary import AcceptedHistory
from eccpn.history import FiniteOrbitalPast
from eccpn.incoming import ParabolicOrbit
from eccpn.nonlinear_force import sh_nonlinear_5pn_terms


class NonlinearHarmonicDynamics(unittest.TestCase):
    def test_rotation_parity_and_fixed_grade(self):
        binary = Binary(.71,.29)
        x = np.array([[14.,3.,1.],[31.,5.,-2.]])
        v = np.array([[.14,.29,-.07],[-.07,.15,.04]])
        angle = .37
        sine,cosine = np.sin(angle),np.cos(angle)
        rotation = np.array([[cosine,-sine,0.],[sine,cosine,0.],[0.,0.,1.]])
        direct = sh_nonlinear_5pn_terms(binary,x,v)
        reverse = sh_nonlinear_5pn_terms(binary,x,-v)
        rotated = sh_nonlinear_5pn_terms(binary,x@rotation.T,v@rotation.T)
        weak_pn = sh_nonlinear_5pn_terms(replace(binary,c=2.),x,v)
        for name,actual in direct.items():
            with self.subTest(name=name):
                self.assertTrue(np.all(np.linalg.norm(actual,axis=1)>0))
                np.testing.assert_allclose(actual,reverse[name],rtol=5e-14,atol=1e-22)
                np.testing.assert_allclose(actual/2**10,weak_pn[name],rtol=5e-14,atol=1e-24)
                np.testing.assert_allclose(actual@rotation.T,rotated[name],rtol=1e-13,atol=1e-22)

    @staticmethod
    def declared_past(binary):
        initial = InitialState((17.,1.,0.),(.08,.26,0.))
        precursor = evolve(binary,initial,None,(-.2,0.),Accuracy("newtonian"))
        state = InitialState(*(tuple(value) for value in precursor.state(0.)))
        angular = binary.mu*np.cross(*precursor.state(0.))
        template = FiniteOrbitalPast((-.2,-.1,0.),precursor.state,{0:binary.mass,2:-.001,4:-.00002},
                                      {0:angular,2:.014*angular},1.,"nonlinear-force-test")
        zero = (0.,0.,0.)
        return state,template.initial_data(binary,zero,zero,zero)

    def test_default_uses_full_radiative_budget(self):
        binary = Binary(.71,.29)
        state,past = self.declared_past(binary)
        trajectory = evolve(binary,state,past,(0.,.02))
        self.assertEqual(trajectory.accuracy.dynamics,"full")
        self.assertEqual(tuple(trajectory.nonlinear_5pn_terms(.01)),NONLINEAR_5PN_SECTORS)
        self.assertEqual(trajectory.accuracy.causal_tail_split_time,480.)
        wave = compute_radiation(trajectory,[.005,.015])
        self.assertEqual((wave.fixed_history_amplitude_pn,wave.fixed_history_quadrupole_pn,
                          wave.fixed_history_luminosity_relative_pn),(3.5,4.,4.))
        self.assertGreater(len(wave.retained_memory_news),0)
        self.assertEqual(max(key[1] for key in wave.retained_memory_news),11)
        np.testing.assert_allclose(wave.energy_flux_infinity,wave.energy_flux_from_modes,rtol=1e-13)
        phase = .21
        rotation = np.array([[np.cos(phase),-np.sin(phase),0.],[np.sin(phase),np.cos(phase),0.],[0.,0.,1.]])
        frame = BondiFrame(rotation=rotation,ell_max=2,polar_points=9,azimuthal_points=18)
        with self.assertRaisesRegex(ValueError, "12 x 23"):
            compute_radiation(trajectory,[.005,.015],frame=frame)
        frame = BondiFrame(rotation=rotation,ell_max=2,polar_points=12,azimuthal_points=23)
        retimed = compute_radiation(trajectory,[.005,.015],frame=frame)
        np.testing.assert_array_equal(retimed.h22,wave.h22)
        np.testing.assert_array_equal(retimed.energy_flux_infinity,wave.energy_flux_infinity)
        self.assertGreater(np.max(abs(retimed.strain[2,2]-wave.h22)),1e-8)

    def test_default_checkpoint_matches_uninterrupted_physical_radiation(self):
        binary = Binary(.73, .27)
        orbit = ParabolicOrbit(31., (1., 0., 0.), (0., 1., 0.), -24000.)
        initial, past = orbit.prepare(binary)
        settings = Accuracy(max_step=.125, rtol=1e-9)
        start, stop = orbit.cutoff, orbit.cutoff + 1.0
        uninterrupted = evolve(binary, initial, past, (start, stop), settings)
        boundary = uninterrupted._solution.boundaries[4]
        first = evolve(binary, initial, past, (start, boundary), settings)
        with TemporaryDirectory() as directory:
            checkpoint = Path(directory) / "default-physical.json"
            first.checkpoint.save(checkpoint)
            saved = AcceptedHistory.load(checkpoint, identity=past.quadrupole.identity)
        state = InitialState(tuple(saved.state[:3]), tuple(saved.state[3:6]))
        resumed = evolve(binary, state, saved, (boundary, stop), settings)
        cuts = np.array([start + .73, start + .96])
        np.testing.assert_allclose(resumed.state(cuts), uninterrupted.state(cuts), rtol=3e-13, atol=1e-23)
        np.testing.assert_allclose(resumed.recoil(cuts), uninterrupted.recoil(cuts), rtol=3e-13, atol=1e-23)
        expected, actual = (compute_radiation(path, cuts) for path in (uninterrupted, resumed))
        for name in ("h22", "news22", "energy_flux_infinity"):
            np.testing.assert_allclose(getattr(actual, name), getattr(expected, name), rtol=5e-12, atol=1e-22)
        self.assertEqual(set(actual.retained_memory_news), set(expected.retained_memory_news))
        for name in ("retained_memory_strain", "retained_memory_news"):
            for key, sample in getattr(actual, name).items():
                np.testing.assert_allclose(sample, getattr(expected, name)[key], rtol=5e-12, atol=1e-24)

    def test_each_sector_enters_the_evolved_state_and_replay(self):
        binary = Binary(.71,.29)
        state,past = self.declared_past(binary)
        settings = Accuracy(SH_RADIATIVE_5PN_MODEL,rtol=1e-11,atol=1e-13,waveform_pn=4.)
        t = .025
        complete = evolve(binary,state,past,(0.,t),settings)
        terms = complete.nonlinear_5pn_terms(0.)
        self.assertEqual(tuple(terms),NONLINEAR_5PN_SECTORS)
        for omitted in NONLINEAR_5PN_SECTORS:
            with self.subTest(omitted=omitted):
                mask = tuple(name for name in NONLINEAR_5PN_SECTORS if name!=omitted)
                reduced = evolve(binary,state,past,(0.,t),replace(settings,nonlinear_5pn_sectors=mask))
                def difference(time):
                    return (complete.state(time)[1]-reduced.state(time)[1])/time
                observed = 2*difference(t/2)-difference(t)
                self.assertGreater(np.linalg.norm(observed),1e-11)
                np.testing.assert_allclose(observed,terms[omitted],rtol=4e-3,atol=1e-13)
        lower = evolve(binary,state,past,(0.,t),replace(settings,dynamics=SH_LOCAL_TAIL_MODEL))
        empty = evolve(binary,state,past,(0.,t),replace(settings,nonlinear_5pn_sectors=()))
        np.testing.assert_array_equal(lower.state(t),empty.state(t))
        native = compute_radiation(complete,[t/4,t/2,t])
        self.assertEqual(native.fixed_history_quadrupole_pn,4.)
        np.testing.assert_allclose(native.energy_flux_infinity,native.energy_flux_from_modes,rtol=8e-14)
        self.assertTrue(np.isfinite(native.h22).all())
        restart = InitialState(*(tuple(array) for array in complete.state(t)))
        with self.assertRaises(ValueError):
            evolve(binary,restart,complete.checkpoint,(t,2*t),replace(settings,nonlinear_5pn_sectors=()))
