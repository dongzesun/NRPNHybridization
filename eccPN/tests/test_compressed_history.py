"""Long-history compression check with independent direct panels and restart."""

from pathlib import Path
from tempfile import TemporaryDirectory
import unittest

import numpy as np
from scipy.integrate import quad_vec

from eccpn import Accuracy, Binary, InitialState, evolve
from eccpn.compressed_history import CompressedPanelHistory
from eccpn.hereditary import AcceptedHistory, integrate_hereditary
from eccpn.history import PanelHistory, SourcePanel
from eccpn.incoming import ParabolicOrbit
from eccpn.newtonian_moments import newtonian_source_jets
from eccpn.tail_derivative import fourfold_tail_jet
from eccpn.tail_parts import fourfold_tail_integral, inverse_four_panel
from eccpn.tails import PrimitiveIncomingTail, quadrupole_tail_sources, sixth_quadrupole


class CompressedHistory(unittest.TestCase):
    def test_irregular_long_prefix_clipping_and_provisional_panel(self):
        generator = np.random.default_rng(96813)
        widths = np.exp(generator.uniform(-6,4,1017))
        times = -92 + np.r_[0,np.cumsum(widths)]
        panels = tuple(SourcePanel(a,b,generator.normal(size=(3,3,2)))
                       for a,b in zip(times[:-1],times[1:]))
        direct = PanelHistory(times[0],(3,2),"compression-weights",panels)
        fast = CompressedPanelHistory.from_history(direct)
        trial = SourcePanel(times[-1],times[-1]+0.31,generator.normal(size=(3,3,2)))
        for time in (times[0],times[1],times[529]+widths[529]/2,times[-1],times[-1]+0.25):
            for scale in (0.008,2.,900.):
                for power in (0,1,2):
                    exact = direct.integral(time,scale,power,trial)
                    computed = fast.integral(time,scale,power,trial)
                    np.testing.assert_allclose(computed,exact,rtol=3e-13,atol=4e-12)
        for time in (times[0],times[1],times[300]+widths[300]/3,times[-1]):
            np.testing.assert_array_equal(fast.value(time),direct.value(time))

    def test_rejection_and_bitwise_checkpoint_continuation(self):
        def rhs(time,state,history,trial):
            return -4*state+0.4*history.integral(time,1.,trial=trial)

        history = CompressedPanelHistory.from_history(PanelHistory(0.,(1,),"compressed-resume"))
        initial = AcceptedHistory(np.array([1.]),history,0.36)
        options = dict(rtol=2e-7,atol=1e-10,max_step=0.36)
        full = integrate_hereditary(rhs,lambda x:x,initial,0.72,**options)
        self.assertGreater(full.rejected_steps,0)
        split = full.boundaries[8]
        first = integrate_hereditary(rhs,lambda x:x,initial,split,**options)
        with TemporaryDirectory() as directory:
            path = Path(directory)/"checkpoint.json"
            first.checkpoint.save(path)
            recovered = AcceptedHistory.load(path,identity="compressed-resume")
            self.assertIsInstance(recovered.history,CompressedPanelHistory)
            resumed = integrate_hereditary(rhs,lambda x:x,recovered,0.72,**options)
        np.testing.assert_array_equal(resumed.boundaries,full.boundaries[8:])
        np.testing.assert_array_equal(resumed.dense,full.dense[8:])
        np.testing.assert_array_equal(resumed.checkpoint.state,full.checkpoint.state)

    def test_remote_q6_primitive_and_retarded_restart(self):
        binary = Binary(0.72,0.28,G=1.13,c=1.41)
        orbit = ParabolicOrbit(28.,(1.,0.,0.),(0.,1.,0.),-260.)
        initial,past = orbit.prepare(binary)
        remote = PrimitiveIncomingTail(past.quadrupole.incoming,orbit,binary,nodes=128)
        exact = newtonian_source_jets(binary,*orbit.state_at(binary,orbit.cutoff),
                                      highest_derivative=7,highest_rank=2)["M2"][0]
        np.testing.assert_allclose(sixth_quadrupole(binary,np.r_[initial.position,initial.velocity]),
                                   exact[6],rtol=7e-14,atol=1e-24)
        for cut in (orbit.cutoff+2.3,orbit.cutoff+30.):
            for power in (0,1,2):
                def source(q):
                    duration = -orbit.cutoff*q/(1-q)
                    y = np.concatenate(orbit.state_at(binary,orbit.cutoff-duration))
                    return (-orbit.cutoff/(1-q)**2*np.log((cut-orbit.cutoff+duration)/3.)**power
                            *quadrupole_tail_sources(binary,y)[2])
                expected,_,info = quad_vec(source,0.,1.,epsabs=1e-21,epsrel=3e-11,full_output=True)
                self.assertTrue(info.success)
                np.testing.assert_allclose(remote.integral(cut,3.,power)[2],expected,
                                           rtol=2e-10,atol=2e-21)
        accuracy = Accuracy("sh-through-4.5PN+recoil",compress_history=True,
                            causal_tail_interpolation="primitive",max_step=0.125,rtol=1e-9)
        t0,t1 = orbit.cutoff,orbit.cutoff+3.
        full = evolve(binary,initial,past,(t0,t1),accuracy)
        q6_at_stop = sixth_quadrupole(binary,np.concatenate(full.state(t1)))
        np.testing.assert_allclose(full.history.integral(t1,1.,0)[2],q6_at_stop,rtol=2e-13,atol=8e-24)
        middle = full._solution.boundaries[24]
        before = evolve(binary,initial,past,(t0,middle),accuracy)
        with TemporaryDirectory() as directory:
            path = Path(directory)/"physical-checkpoint.json"
            before.checkpoint.save(path)
            checkpoint = AcceptedHistory.load(path,identity=past.quadrupole.identity)
            resumed = evolve(binary,InitialState(tuple(checkpoint.state[:3]),tuple(checkpoint.state[3:6])),
                             checkpoint,(middle,t1),accuracy)
        np.testing.assert_array_equal(resumed.checkpoint.state,full.checkpoint.state)

    def test_ratios_in_far_integral_with_clipped_quadratic_panels(self):
        generator = np.random.default_rng(8745)
        edges = np.cumsum(np.r_[0., np.exp(generator.uniform(-4,3,127))])
        data = tuple(SourcePanel(a,b,generator.normal(size=(3,5,3,3)))
                     for a,b in zip(edges[:-1],edges[1:]))
        finite = PanelHistory(edges[0],(5,3,3),"inverse-quartic",data)
        compressed = CompressedPanelHistory.from_history(finite)
        for time,stop in ((edges[-1]+10,edges[-1]-3.5),(3e6,edges[-1]),(edges[12]+2,edges[12]-0.01)):
            expected = np.zeros((3,3))
            for panel in data:
                if panel.start>=stop:
                    break
                width = panel.end-panel.start
                fraction = (min(stop,panel.end)-panel.start)/width
                integral,_,info = quad_vec(
                    lambda q:width*np.polynomial.polynomial.polyval(q,panel.coefficients[:,3])
                    /(time-panel.start-width*q)**4,0.,fraction,epsabs=1e-25,epsrel=3e-12,full_output=True,
                )
                self.assertTrue(info.success)
                expected += integral
            observation = compressed.inverse_four(time,stop,3,inverse_four_panel)
            np.testing.assert_allclose(observation,expected,rtol=3e-11,atol=3e-26)

    def test_fourfold_physical_restart(self):
        binary = Binary(0.73,0.27)
        orbit = ParabolicOrbit(31.,(1.,0.,0.),(0.,1.,0.),-220.)
        initial,past = orbit.prepare(binary)
        accuracy = Accuracy("sh-through-4.5PN+recoil",compress_history=True,
                            causal_tail_interpolation="integration_by_parts",causal_tail_split_time=90.,
                            max_step=0.125,rtol=1e-9)
        t0,t1 = orbit.cutoff,orbit.cutoff+3.
        full = evolve(binary,initial,past,(t0,t1),accuracy)
        middle = full._solution.boundaries[24]
        before = evolve(binary,initial,past,(t0,middle),accuracy)
        with TemporaryDirectory() as directory:
            path = Path(directory)/"fourfold.json"
            before.checkpoint.save(path)
            checkpoint = AcceptedHistory.load(path,identity=past.quadrupole.identity)
            resumed = evolve(binary,InitialState(tuple(checkpoint.state[:3]),tuple(checkpoint.state[3:6])),
                             checkpoint,(middle,t1),accuracy)
        np.testing.assert_array_equal(resumed.checkpoint.state,full.checkpoint.state)


    def test_fourfold_trial_cache_and_same_history_log_derivative(self):
        generator = np.random.default_rng(8124)
        cuts = (0.,0.7,1.5,2.6,4.1,4.8)
        panels = tuple(SourcePanel(a,b,generator.normal(size=(3,5,3,3)))
                       for a,b in zip(cuts[:-1],cuts[1:]))
        finite = PanelHistory(cuts[0],(5,3,3),"fourfold-log-jet",panels[:-1])
        logarithmic_rate = 0.17
        for committed in (finite,CompressedPanelHistory.from_history(finite)):
            extended = committed.append(panels[-1])
            for cut,lag in ((0.42,1.2),(2.28,1.0),(4.55,0.23),(4.52,7.),(4.63,1.6)):
                def value(t,scale):
                    provisional = panels[-1] if t>committed.end else None
                    return fourfold_tail_integral(committed,t,scale,lag,provisional)

                scale = lambda t: np.exp(logarithmic_rate*t)
                history = extended if cut>committed.end else committed
                observed,rate = fourfold_tail_jet(history,cut,scale(cut),logarithmic_rate,lag)
                for radial_scale in (0.03,scale(cut),500.):
                    expected = fourfold_tail_integral(history,cut,radial_scale,lag,None)
                    np.testing.assert_allclose(value(cut,radial_scale),expected,rtol=5e-14,atol=2e-12)
                np.testing.assert_allclose(observed,value(cut,scale(cut)),rtol=5e-14,atol=2e-12)
                step = 2e-5
                stencil = ((-2,1),(-1,-8),(1,8),(2,-1))
                expected = sum(coef*value(cut+k*step,scale(cut+k*step)) for k,coef in stencil)/(12*step)
                np.testing.assert_allclose(rate,expected,rtol=5e-8,atol=7e-8)

    def test_fourfold_remote_derivative_before_and_after_split_crosses_join(self):
        binary = Binary(0.73,0.27)
        orbit = ParabolicOrbit(31.,(1.,0.,0.),(0.,1.,0.),-220.)
        initial,past = orbit.prepare(binary)
        accuracy = Accuracy("sh-through-4.5PN+recoil",compress_history=True,
                            causal_tail_interpolation="integration_by_parts",causal_tail_split_time=1.5,
                            max_step=0.125,rtol=1e-9)
        t0 = orbit.cutoff
        result = evolve(binary,initial,past,(t0,t0+3.),accuracy)
        for cut in (t0+0.81,t0+2.57):
            scale = lambda time: np.exp((time-t0)*0.006)*58.
            integral = lambda time: fourfold_tail_integral(result.history,time,scale(time),1.5,None)
            first,derivative = fourfold_tail_jet(result.history,cut,scale(cut),0.006,1.5)
            np.testing.assert_allclose(first,integral(cut),rtol=1e-14,atol=0)
            step = 0.001
            reference = (integral(cut-2*step)-8*integral(cut-step)
                         +8*integral(cut+step)-integral(cut+2*step))/(12*step)
            self.assertLess(np.linalg.norm(reference-derivative)/np.linalg.norm(derivative),2e-5)


if __name__ == "__main__":
    unittest.main()
