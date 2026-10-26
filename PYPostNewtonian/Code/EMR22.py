"""Schwarzschild test-particle additions to the existing comparable-mass PN model.

N > 6 means EMR-supplemented flux and energy, not complete N-PN dynamics.
Known spin/precession and horizon absorption terms retain their original orders.
"""
import numpy as np
from numba import njit
from scipy.integrate import solve_ivp
from types import SimpleNamespace

from . import PNEvolution as base
from .EMR22Coefficients import FLUX, ENERGY
from .PNBalance6 import balance6


def validate_order(order):
    if not np.isfinite(order) or order < 0 or order > 22 or 2*order != int(2*order):
        raise ValueError('PNEvolutionOrder must be a half-integer between 0 and 22')


@njit(cache=True)
def flux_tail(v, nu, twice_order):
    """Only relative powers v^13 through v^(2N); no double counting."""
    logv = np.log(v)
    total = 0.0
    for k in range(twice_order, 12, -1):
        coefficient = 0.0
        for j in range(7, -1, -1):
            coefficient = coefficient*logv + FLUX[k, j]
        total = total*v + coefficient
    return (32.0/5.0)*nu*nu*v**23*total


@njit(cache=True)
def energy_derivative_tail(v, mass, nu, twice_order):
    total = 0.0
    for n in range(twice_order//2, 6, -1):
        total = total*v*v + (2*n+2)*ENERGY[n]
    return -0.5*mass*nu*v**15*total


def initialize(wHat, xHat, yHat, zHat, m1, m2, v0, s1, s2):
    z = np.array([0.0])
    cons = base.Cons(*([z]*(len(base.ConsSpec)-2)), True, True)
    var = base.Vars(*([z]*len(base.VarsSpec)))
    base.Initialization(cons, wHat, xHat, yHat, zHat, m1, m2, v0, s1, s2, 1., 1., 1., 1.)
    return cons, var


@njit(cache=True)
def balance(cons, var, order2):
    f, de, absorption = balance6(cons, var)
    v, nu, mass = var.v[0], cons.nu[0], cons.M[0]
    return f + flux_tail(v, nu, order2), de + energy_derivative_tail(v, mass, nu, order2), absorption


@njit(cache=True)
def rhs(cons, var, y, order2):
    base.Recalculate_6p0(cons, var, y)
    f, de, absorption = balance(cons, var, order2)
    dy = base.TaylorT1_6p0(cons, var)
    dy[0] = -(f + absorption)/de
    return dy


def energy_flux(y, wHat, xHat, yHat, zHat, m1, m2, v0, s1, s2, frame,
                omega_start=None, omega_end=None, t_PNStart=False, t_PNEnd=False,
                PNEvolutionOrder=22., TaylorTn=1, StepsPerOrbit=32, dt=None,
                tol=1e-8, MinStep=1e-7):
    validate_order(PNEvolutionOrder)
    y = np.asarray(y, dtype=float)
    if y.ndim != 2 or y.shape[0] != 9 or not np.all(np.isfinite(y)) or np.any(y[0] <= 0):
        raise ValueError('Flux states must be finite 9-by-N arrays with positive velocity')
    cons, var = initialize(wHat,xHat,yHat,zHat,m1,m2,v0,s1,s2)
    result = np.empty(y.shape[1])
    for i in range(y.shape[1]):
        base.Recalculate_6p0(cons, var, y[:,i])
        result[i] = balance(cons, var, int(2*PNEvolutionOrder))[0]
    return result


def Evolution(wHat, xHat, yHat, zHat, m1, m2, v0, s1, s2, frame,
              omega_start=None, omega_end=None, t_PNStart=False, t_PNEnd=False,
              PNEvolutionOrder=22., TaylorTn=1, StepsPerOrbit=32, dt=None,
              tol=1e-8, MinStep=1e-7):
    """DOP853 evolution using the original frame/spin equations and new balance.

    Integrate with v as the independent variable: dt/dv = -E'/F is regular at
    the energy turning point. Stop there, then sample an orbital time grid.
This also avoids stepping across the singularity of dv/dt near the ICO.
"""
    validate_order(PNEvolutionOrder)
    if TaylorTn != 1:
        raise ValueError('EMR supplementation currently supports TaylorT1 only')
    cons, var = initialize(wHat,xHat,yHat,zHat,m1,m2,v0,s1,s2)
    mass, nu = float(cons.M[0]), float(cons.nu[0])
    order2 = int(2*PNEvolutionOrder)
    y0 = np.array([v0,0.,0.,0.,0.,*np.asarray(frame, dtype=float)])
    base.Recalculate_6p0(cons,var,y0)
    f0, de0, a0 = balance(cons,var,order2)
    if de0 >= 0 or f0+a0 <= 0:
        raise ValueError('Initial state is beyond the inspiral branch of energy balance')
    tmerger = 5*mass/(256*nu*v0**8)
    tend = float(t_PNEnd) if t_PNEnd else tmerger
    tstart = float(t_PNStart) if t_PNStart else -3*tmerger
    vstart = (omega_start*mass)**(1/3) if omega_start else 0.01
    vend = (omega_end*mass)**(1/3) if omega_end else 0.7

    def state(v,z):
        return np.concatenate(([v],z[1:]))

    def integrand(v,z):
        y = state(v,z)
        base.Recalculate_6p0(cons,var,y)
        f,de,a = balance(cons,var,order2)
        if not np.isfinite(f+de+a) or f+a <= 0:
            raise ValueError('Nonfinite or nonpositive dissipative flux during evolution')
        dtdv = -de/(f+a)
        dy = base.TaylorT1_6p0(cons,var)
        return np.concatenate(([dtdv], dy[1:]*dtdv))

    def turning_point(v,z):
        base.Recalculate_6p0(cons,var,state(v,z))
        return -balance(cons,var,order2)[1]
    turning_point.terminal = True
    turning_point.direction = -1

    def forward_time(v,z): return z[0]-tend
    forward_time.terminal = True
    def backward_time(v,z): return z[0]-tstart
    backward_time.terminal = True

    initial = np.concatenate(([0.],y0[1:]))
    solutions = []
    directions = [(vend,forward_time)]
    if omega_start or t_PNStart:
        directions.insert(0,(vstart,backward_time))
    for limit,event in directions:
        sol = solve_ivp(integrand,(v0,limit),initial,method='DOP853',rtol=tol,atol=tol,
                        dense_output=True,events=(event,turning_point),max_step=0.002)
        if not sol.success:
            raise RuntimeError(f'High-order PN evolution failed: {sol.message}')
        solutions.append(sol)

    # Invert the strictly monotone t(v) separately on the two branches. Brent
    # solves against the dense ODE solution, without introducing spline errors.
    from scipy.optimize import brentq
    times, states = [], []
    for sol in solutions:
        endtime = float(sol.y[0,-1])
        sign = 1 if endtime > 0 else -1
        # The legacy Newtonian time-grid estimate vanishes at tmerger, which
        # need not coincide with this PN model's endpoint. Bound its step by
        # the shortest orbital step on the actual integrated branch.
        orbital_step = (2*mass/(max(v0,sol.t[-1])**3*StepsPerOrbit)
                        if dt is None else 0.0)
        samples = [0.]
        while sign*samples[-1] < sign*endtime:
            step = dt if dt is not None else max(orbital_step,
                2*mass*(256*nu*max(tmerger-samples[-1],0.)/(5*mass))**(3/8)/StepsPerOrbit)
            if step < MinStep: break
            samples.append(samples[-1]+sign*step)
        samples = np.array(samples[:-1])
        for t in samples:
            if t == 0:
                v = v0
            else:
                v = brentq(lambda u: sol.sol(u)[0]-t,min(v0,sol.t[-1]),max(v0,sol.t[-1]),xtol=1e-14)
            times.append(t)
            states.append(state(v,sol.sol(v)))
    t,index = np.unique(np.asarray(times),return_index=True)
    return SimpleNamespace(t=t,y=np.asarray(states)[index].T,success=True,
                           message='Completed on the physical inspiral branch')
