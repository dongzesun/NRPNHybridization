"""Independent time-domain checks of paired logarithmic tails."""

import unittest

import numpy as np
from scipy.integrate import quad

from eccpn.consistent_tail import TailPrimitive
from eccpn.history import PanelHistory, SourcePanel
from eccpn.waveform_history import CanonicalHistory


class ConsistentTail(unittest.TestCase):
    def test_polynomial_causal_source(self):
        self.check_polynomial(1)

    def test_third_time_derivative(self):
        self.check_polynomial(3)

    def check_polynomial(self, derivatives):
        cuts = np.array([0., .125, .5, 1., 4., 16., 45., 64.])
        panels = []
        for a, b in zip(cuts[:-1],cuts[1:]):
            # Only endpoint jets matter; use their exact endpoint chord and
            # a zero quadratic coefficient to exercise the canonical reader.
            def jet(time):
                return np.array([time**3,3*time**2,6*time,6.,0.])
            left, right = jet(a), jet(b)
            panels.append(SourcePanel(a,b,np.stack((left,right-left,0*left))))
        past = PanelHistory(cuts[0],(5,),"polynomial",tuple(panels))
        history = CanonicalHistory({"M2":{0:past}},{0:1.},{},1.,1.)
        fitted = TailPrimitive.from_canonical(history,"M2",0,{0:1.},news_derivatives=derivatives)
        primitive = TailPrimitive(fitted.cuts,fitted.coefficients,(np.array([]),)*3)
        terms = ((.9,2.3,2),(-.4,3.7,1),(.2,1.,0))
        time = 53.

        def kernel(delay):
            return sum(weight*np.log(delay/scale)**power for weight,scale,power in terms)

        # Direct quadrature over the whole causal interval, without a split or
        # integration by parts. Cubic compact onset has P=P'=P''=0 at zero.
        expected = np.array([quad(lambda delay:6*(time-delay)*kernel(delay),0,time,
                                  epsabs=1e-10,epsrel=2e-13)[0],
                             quad(lambda delay:6*kernel(delay),0,time,
                                  epsabs=1e-10,epsrel=2e-13)[0]])
        for lag in (.25,1.,7.,28.,52.95):
            with self.subTest(lag=lag):
                actual = primitive.integrate(time,lag,terms)
                np.testing.assert_allclose(actual[:2],expected,rtol=2e-11,atol=1e-9)
                if derivatives==3:
                    # The causal cubic's third primitive derivative turns on
                    # at zero. Subsequent derivatives include that boundary.
                    kernel_rate = (.9*2*np.log(time/2.3)-.4)/time
                    np.testing.assert_allclose(actual[2:],[6*kernel(time),6*kernel_rate],
                                               rtol=2e-7,atol=1e-9)
                dt = 1/512
                evaluated = [primitive.integrate(time+i*dt,lag,terms)[0]
                             for i in (-2,-1,1,2)]
                stencil = np.array([1,-8,8,-1])@evaluated/(12*dt)
                np.testing.assert_allclose(stencil,actual[1],rtol=2e-8,atol=2e-7)

    def test_multiple_integrations_by_parts(self):
        from math import comb, factorial

        from eccpn.consistent_tail import SeparatedLogJet

        cuts = np.array([0., .2, 1., 2., 4., 8., 16., 32.])

        def polynomial(a,b):
            jets = np.zeros((10,9))
            for order in range(9):
                prefactor = factorial(9)/factorial(9-order)/16**9
                for power in range(10-order):
                    jets[power,order] = prefactor*comb(9-order,power)*a**(9-order-power)*(b-a)**power
            return jets

        jets = np.stack([polynomial(a,b) for a,b in zip(cuts[:-1],cuts[1:])])
        history = SeparatedLogJet(cuts,jets,np.array([]),np.array([]),np.empty((0,9)))
        time = 20.
        terms = ((.5,2.1,2),(-.3,1.7,1),(.6,1.,0))

        def kernel(delay):
            return sum(weight*np.log(delay/scale)**power for weight,scale,power in terms)

        exact = [quad(lambda lag:factorial(9)/factorial(9-order)*(time-lag)**(9-order)
                      /16**9*kernel(lag),0,time,epsabs=2e-14,epsrel=5e-13)[0]
                 for order in range(5,9)]
        for delay in (.4,1.7,7.,19.9):
            for integrations in (2,4):
                with self.subTest(delay=delay,integrations=integrations):
                    value = history.integrate(time,delay,terms,5,4,integrations)
                    np.testing.assert_allclose(value,exact,rtol=2e-9,atol=1e-11)

    def test_earlier_product_pulse_rate(self):
        from eccpn.consistent_tail import SeparatedLogRate

        cuts = np.array([0.,1.,3.,8.])
        polynomial = np.zeros((3,9))
        polynomial[0,4:] = [1.,-4.,6.,-4.,1.]
        history = SeparatedLogRate(cuts,polynomial,np.array([]),np.array([]),np.array([]))
        time = 5.
        terms = ((1.,2.,1),(.25,3.,2))
        target = quad(lambda t:t**4*(1-t)**4*(1+.5*np.log((time-t)/3))/(time-t),0,1,
                      epsabs=1e-15)[0]
        expected_value = quad(lambda t:t**4*(1-t)**4*sum(
            weight*np.log((time-t)/scale)**power for weight,scale,power in terms),0,1,
            epsabs=1e-15)[0]
        np.testing.assert_allclose(history.value(time,terms),expected_value,rtol=4e-11,atol=1e-14)
        for lag in (.25,2.,4.5,4.95):
            np.testing.assert_allclose(history.rate(time,lag,terms),target,rtol=4e-11,atol=1e-14)
