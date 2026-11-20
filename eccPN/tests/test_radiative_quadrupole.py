"""Generic radiative quadrupole kernel and news checks on known source histories."""

from math import factorial
import unittest

import numpy as np
from scipy.integrate import quad, quad_vec

from eccpn import Binary
from eccpn.asymptotic_waveform import AsymptoticQuadrupole, MatchedQuadrupoleHistory
from eccpn.history import PanelHistory, SourcePanel
from eccpn.hereditary import _interpolate
from eccpn.memory import (
    null_memory_news, quadrupole_tail_memory_boundary, raw_quadrupole_overlap_grade_eight,
)
from eccpn.projection import project_stf
from eccpn.radiation_kernels import ordinary_tail
from eccpn.radiative_quadrupole import _nested_memory, mpm_quadrupole_interactions


def polynomial_panel(start, end, power, basis):
    grid = np.polynomial.Polynomial([start, end - start])
    coefficients = np.zeros((power + 1, 9, 3, 3))
    for derivative in range(power + 1):
        polynomial = grid**(power - derivative) * factorial(power) / factorial(power - derivative)
        coefficients[:len(polynomial.coef), derivative] = polynomial.coef[:, None, None] * basis
    return SourcePanel(start, end, coefficients)


class RadiativeQuadrupole(unittest.TestCase):
    def test_remote_and_mixed_ordered_domains_against_full_past(self):
        shape = np.array([[0.8, 0.2, 0.15], [0.2, -1.1, -0.33], [0.15, -0.33, 0.3]])
        product = shape @ shape
        product -= np.eye(3)*np.trace(product)/3
        scale, offset, time = 0.83, -1613/270, -2.0
        harmonic = sum(1/k for k in range(1, 6))
        constant = 1/5 - harmonic + offset
        expected_rate = 120**2/5 * (-time)**-11 * (np.log(-time/scale)+constant) * product
        expected_value = 120**2/50 * (-time)**-10 * (np.log(-time/scale)+constant+1/10) * product

        def moment(cut):
            return np.array([factorial(n+1)*(-cut)**(-n-2)*shape for n in range(10)])

        for join, n_panels in ((-3.1, 25), (-4.2, 65)):
            remote = AsymptoticQuadrupole.from_terms(join, 0.0, 1.1, {(2, 0): shape})
            times = np.linspace(join, time, n_panels+1)
            panels = tuple(SourcePanel(a, b, _interpolate((moment(a), moment((a+b)/2), moment(b))))
                           for a, b in zip(times[:-1], times[1:]))
            past = MatchedQuadrupoleHistory(PanelHistory(join, (10, 3, 3), "combined-ordered", panels), remote)
            actual = _nested_memory(past, time, scale, 2e-8, 2e-13)
            np.testing.assert_allclose(actual, [expected_value, expected_rate], rtol=7e-6, atol=2e-10)
            self.assertLess(np.linalg.norm(past.start_jump(4)),
                            4*np.finfo(float).eps*np.linalg.norm(moment(join)[4]))

    def test_four_pn_raw_memory_boundary_and_news(self):
        binary = Binary(0.64, 0.36, G=1.1, c=1.4)
        shape = np.array([[0.8, 0.2, 0.15], [0.2, -1.1, -0.33], [0.15, -0.33, 0.3]])
        history = PanelHistory(0, (9, 3, 3), "raw-overlap", (
            polynomial_panel(0, 0.31, 4, shape), polynomial_panel(0.31, 0.94, 4, shape),
        ))
        time, scale = 0.74, 0.88
        modes = [(2, m) for m in range(-2, 3)]

        def angular_projection(cut):
            newtonian = project_stf(binary, 2, history.value(cut)[3])
            tail = project_stf(binary, 2, ordinary_tail(binary, history, cut, 2, "mass", binary.mass, scale)[1])
            inventory = {(grade, *mode): value for grade, items in ((0, newtonian), (3, tail))
                         for mode, value in items.items()}
            raw = null_memory_news(binary, inventory, max_grade=3)
            return np.array([raw[(8, *mode)] for mode in modes])

        strain, rate = raw_quadrupole_overlap_grade_eight(binary, history, time, binary.mass, scale)
        strain_modes, news_modes = project_stf(binary, 2, strain), project_stf(binary, 2, rate)
        independent = quad_vec(angular_projection, 0, time, points=[0.31], epsabs=3e-10, epsrel=5e-10)[0]
        np.testing.assert_allclose([strain_modes[key] for key in modes], independent, rtol=2e-9, atol=1e-8)
        np.testing.assert_allclose([news_modes[key] for key in modes], angular_projection(time),
                                   rtol=2e-12, atol=2e-12)
        infrared = 1.07
        nested = _nested_memory(history, time, 2 * infrared / binary.c, 2e-11, 2e-12)
        nested *= (8 * binary.G**2 * binary.mass) / (7 * binary.c**8)
        boundary = quadrupole_tail_memory_boundary(binary, history, time, binary.mass, infrared, scale)
        np.testing.assert_allclose(nested + boundary, [strain, rate], rtol=4e-12, atol=5e-12)

    def test_ordered_nested_integral_under_panel_split(self):
        basis = np.diag([2.0, -1.0, -1.0])
        split = PanelHistory(0, (9, 3, 3), "ordered-test", (
            polynomial_panel(0, 0.4, 4, basis), polynomial_panel(0.4, 1.0, 4, basis),
        ))
        time, scale = 0.83, 0.74
        actual = _nested_memory(split, time, scale, 1e-11, 1e-12)
        product = basis @ basis
        product = product - np.eye(3) * np.trace(product) / 3
        # Each fourth derivative is 4! times the fixed STF basis.
        expected_value = 24**2 * quad(
            lambda t: t * (np.log(t / scale) - 1 - 1613 / 270), 0, time,
        )[0] * product
        expected_news = 24**2 * time * (np.log(time / scale) - 1 - 1613 / 270) * product
        np.testing.assert_allclose(actual, [expected_value, expected_news], rtol=4e-13)
        cache = {}
        for cut in (1., 0.83, 0.4, 0.17, 1.):
            cached = _nested_memory(split, cut, scale, 1e-11, 1e-12, prefix_cache=cache)
            fresh = _nested_memory(split, cut, scale, 1e-11, 1e-12)
            np.testing.assert_array_equal(cached, fresh)

    def test_analytic_news_including_current_dipole(self):
        binary = Binary(0.65, 0.35, G=1.05, c=1.4)
        first = np.array([[1.0, 0.2, 0.4], [0.2, -0.6, 0.1], [0.4, 0.1, -0.4]])
        second = np.array([[0.2, 0.8, 0.1], [0.8, -0.06, 0.3], [0.1, 0.3, -0.14]])
        intervals = ((0., 0.41), (0.41, 1.))
        panels = []
        for start, end in intervals:
            lower = polynomial_panel(start, end, 4, first).coefficients
            higher = polynomial_panel(start, end, 7, second).coefficients.copy()
            higher[:len(lower)] += lower
            panels.append(SourcePanel(start, end, higher))
        history = PanelHistory(0, (9, 3, 3), "polynomial-news", tuple(panels))

        def radiation(time):
            return mpm_quadrupole_interactions(
                binary, history, time, binary.mass, [0.4, 0.7, -0.2], 1.1, 2.0,
            )

        time, offset = 0.87, 2e-5
        middle = radiation(time)
        left, right = radiation(time - offset), radiation(time + offset)
        for grade, (_, derivative) in middle.items():
            numerical = (right[grade][0] - left[grade][0]) / (2 * offset)
            np.testing.assert_allclose(derivative, numerical, rtol=5e-6, atol=2e-5)
            np.testing.assert_allclose(np.trace(derivative), 0.0, atol=1e-8)

    def test_radiative_scale_is_time_translation_across_cubic_grades(self):
        binary = Binary(0.65,0.35,G=1.05,c=1.4)
        mass = binary.mass
        first = np.array([[1.0,0.2,0.4],[0.2,-0.6,0.1],[0.4,0.1,-0.4]])
        second = np.array([[0.2,0.8,0.1],[0.8,-0.06,0.3],[0.1,0.3,-0.14]])
        panels = []
        for a,b in ((0.,0.34),(0.34,0.66),(0.66,1.)):
            lower = polynomial_panel(a,b,6,first).coefficients
            higher = polynomial_panel(a,b,8,second).coefficients.copy()
            higher[:len(lower)] += lower
            panels.append(SourcePanel(a,b,higher))
        history = PanelHistory(0.,(9,3,3),"radiative-scale-translation",tuple(panels))
        scale = 1.2
        shifted_scale = 2.1
        shift = 2*binary.G*mass*np.log(shifted_scale/scale)/binary.c**3
        cache = {}
        spin = np.array([0.4,0.7,-0.2])
        def nonlinear(time,radiative_scale):
            return mpm_quadrupole_interactions(
                binary,history,time,mass,spin,1.4,radiative_scale,prefix_cache=cache,
            )

        for cut in (0.58,0.87):
            late = nonlinear(cut,shifted_scale)
            early = nonlinear(cut,scale)
            source = history.value(cut)
            tail = ordinary_tail(binary,history,cut,2,"mass",mass,scale)
            moved_tail = ordinary_tail(binary,history,cut,2,"mass",mass,shifted_scale)
            np.testing.assert_allclose(moved_tail-tail,-shift*source[3:5],rtol=2e-12,atol=2e-11)
            np.testing.assert_allclose(late[6][0]-early[6][0],-shift*tail[1]+shift**2*source[4]/2,
                                       rtol=1e-11,atol=1e-10)
            np.testing.assert_allclose(late[8][0]-early[8][0],-shift*early[5][1],rtol=1e-11,atol=1e-10)
            step = 2e-5
            def diff(function):
                return (function(cut-2*step)-8*function(cut-step)
                        +8*function(cut+step)-function(cut+2*step))/(12*step)

            tail_news_rate = diff(lambda t:ordinary_tail(binary,history,t,2,"mass",mass,scale)[1])
            quadratic_news_rate = diff(lambda t:nonlinear(t,scale)[5][1])
            np.testing.assert_allclose(late[6][1]-early[6][1],-shift*tail_news_rate+shift**2*source[5]/2,
                                       rtol=2e-9,atol=2e-8)
            np.testing.assert_allclose(late[8][1]-early[8][1],-shift*quadratic_news_rate,rtol=2e-9,atol=2e-8)
