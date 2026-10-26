"""Numerical regressions for the EMR-supplemented flux and evolution."""
import json
from pathlib import Path
import unittest

import mpmath as mp
import numpy as np
import quaternionic
from PYPostNewtonian.Code import EMR22, PostNewtonian, PNEvolution
from PYPostNewtonian.Code.EMR22Coefficients import FLUX, ENERGY

ROOT = Path(__file__).resolve().parents[1]


class EMR22Tests(unittest.TestCase):
    def test_high_precision_flux_tail(self):
        data = json.loads((ROOT/'PYPostNewtonian/PNTerms/EMR22/coefficients_decimal.json').read_text())
        with mp.workdps(80):
            for v in [.05,.2,.22,.3,.4,.45,.5]:
                x,nu = mp.mpf(v),mp.mpf(3)/16
                for n in range(13,45):
                    exact = 32*nu**2*x**10/5 * mp.fsum(
                        mp.fsum(mp.mpf(c)*mp.log(x)**j for j,c in enumerate(data['flux'][k]))*x**k
                        for k in range(13,n+1))
                    actual = EMR22.flux_tail(v,float(nu),n)
                    self.assertAlmostEqual(actual/float(exact),1.,delta=2e-12)

    def test_binding_energy(self):
        with mp.workdps(60):
            for v in [.1,.2,.3,.4,.5]:
                exact = mp.fsum(-mp.mpf(3)/32*(2*n+2)*mp.mpf(str(ENERGY[n]))*mp.mpf(v)**(2*n+1) for n in range(7,23))
                np.testing.assert_allclose(EMR22.energy_derivative_tail(v,1.,3/16,44),float(exact),rtol=2e-15)
            # Through 22PN the Taylor series must recover the exact test-particle
            # binding energy in the weak field, including the leading -v^2/2.
            v=mp.mpf('.1')
            series=-v*v/2*mp.fsum(mp.mpf(str(e))*v**(2*n) for n,e in enumerate(ENERGY))
            exact=(1-2*v*v)/mp.sqrt(1-3*v*v)-1
            self.assertLess(abs(series-exact),mp.mpf('1e-20'))

    def test_preserved_api_and_flux(self):
        with np.load(ROOT/'tests/fixtures/flux_before22pn.npz') as z:
            y,orders,golden=z['y'],z['orders'],z['flux']
        for n,expected in zip(orders,golden):
            args=dict(return_energy_flux=True,y=y,PNEvolutionOrder=n,PNWaveformModeOrder=6)
            values=(3.,.01,[4e-6,0,1e-5],[0,4e-6,6e-6])
            np.testing.assert_array_equal(PostNewtonian.PNWaveform(*values,**args),expected)
        f6=PostNewtonian.PNWaveform(*values,**args)
        for n in [6.5,7,10,21.5,22]:
            args['PNEvolutionOrder']=n
            f=PostNewtonian.PNWaveform(*values,**args)
            tail=np.array([EMR22.flux_tail(v,3/16,int(2*n)) for v in y[0]])
            np.testing.assert_allclose(f,f6+tail,rtol=3e-15,atol=1e-22)

    def test_balance_and_evolution(self):
        args=(quaternionic.one,quaternionic.x,quaternionic.y,quaternionic.z,
              .75,.25,.22,0.*quaternionic.one,0.*quaternionic.one,np.array([1.,0,0,0]))
        cons,var=EMR22.initialize(*args[:9])
        for v in [.2,.3,.4]:
            y=np.array([v,0.,0.,0.,0.,1.,0.,0.,0.])
            dy=EMR22.rhs(cons,var,y,44)
            f,de,a=EMR22.balance(cons,var,44)
            self.assertLess(abs(de*dy[0]+f+a),1e-17)
            np.testing.assert_array_equal(dy[1:],PNEvolution.TaylorT1_6p0(cons,var)[1:])
        a=EMR22.Evolution(*args,t_PNStart=-100,t_PNEnd=100,StepsPerOrbit=None,dt=5,tol=1e-10)
        b=EMR22.Evolution(*args,t_PNStart=-100,t_PNEnd=100,StepsPerOrbit=None,dt=5,tol=1e-12)
        np.testing.assert_array_equal(a.t,b.t)
        np.testing.assert_allclose(a.y,b.y,rtol=1e-8,atol=2e-10)
        self.assertTrue(np.all(np.diff(a.y[0])>0))

    def test_invalid_orders(self):
        for n in [-1,6.1,22.5,np.nan]:
            with self.assertRaises(ValueError):EMR22.validate_order(n)

    def test_grid_beyond_newtonian_time(self):
        q=3.0
        v0=.15
        nu=q/(1+q)**2
        tnewton=5/(256*nu*v0**8)
        args=(quaternionic.one,quaternionic.x,quaternionic.y,quaternionic.z,
              q/(1+q),1/(1+q),v0,0.*quaternionic.one,0.*quaternionic.one,
              np.array([1.,0,0,0]))
        sol=EMR22.Evolution(*args,t_PNEnd=tnewton+1000,PNEvolutionOrder=22,StepsPerOrbit=4)
        self.assertLess(len(sol.t),20000)
        self.assertGreater(sol.t[-1],tnewton)
        self.assertTrue(np.all(np.diff(sol.t)>0))
        self.assertTrue(np.all(np.isfinite(sol.y)))


if __name__=='__main__':
    unittest.main(verbosity=2)
