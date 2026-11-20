"""Leading radiative tail and 4PN dynamical tail on a supplied causal history."""

from dataclasses import dataclass
from functools import cached_property

import numpy as np


_IDENTITY = np.eye(3)
_IDENTITY.flags.writeable = False


def _stf(left, right):
    outer = np.outer(left, right)
    return (outer + outer.T) / 2 - _IDENTITY * np.trace(outer) / 3


def _tail_basis(binary, state):
    position, velocity = state[:3], state[3:]
    radius = np.linalg.norm(position)
    normal = position / radius
    radial_speed = normal @ velocity
    speed_squared = velocity @ velocity
    potential = binary.G * binary.mass / radius
    basis = np.stack((_stf(normal, normal), _stf(normal, velocity), _stf(velocity, velocity)))
    return radius, radial_speed, speed_squared, potential, basis


def _force_sources(binary, geometry):
    radius, s, w, u, basis = geometry
    coefficients = np.array([
        (-30 * s**2 + 2 * u + 6 * w, 36 * s, -8),
        (210 * s**3 - 90 * s * w, -240 * s**2 - 16 * u + 48 * w, 60 * s),
        (
            42 * s * (-4 * u**2 + 120 * u * s**2 - 48 * u * w
                      + 495 * s**4 - 450 * s**2 * w + 75 * w**2),
            -8 * (22 * u**2 + 840 * u * s**2 - 141 * u * w
                  + 2835 * s**4 - 1890 * s**2 * w + 135 * w**2),
            840 * s * (2 * u + 7 * s**2 - 3 * w),
        ),
    ])
    # Section 4's total-derivative recurrence, with G and c explicit.
    factors = binary.mu * (binary.G*binary.mass) / radius**np.array([3, 4, 6])
    return factors[:, None, None] * np.einsum("ab,bij->aij", coefficients, basis)


def quadrupole_tail_sources(binary, state):
    """Newtonian fourth, fifth and seventh derivatives on the instantaneous state."""
    return _force_sources(binary, _tail_basis(binary, state[:6]))


def _sixth_source(binary, geometry):
    radius,s,w,u,basis = geometry
    a = 2*(8*u**2-105*u*s**2+21*u*w-945*s**4+630*s**2*w-45*w**2)
    b = 60*s*(6*u+35*s**2-15*w)
    c = -4*(19*u+135*s**2-27*w)
    return binary.mu*u/radius**4*(a*basis[0]+b*basis[1]+c*basis[2])


def sixth_quadrupole(binary,state):
    """Newtonian sixth time derivative in Cartesian SH phase space."""
    return _sixth_source(binary, _tail_basis(binary, state[:6]))


def primitive_tail_panel(binary,source,states,width):
    """Use the cubic Q6 primitive with endpoint Newtonian Q7 slopes."""
    before,after = source(states[0]),source(states[-1])
    from .hereditary import _interpolate
    result = _interpolate((before,source(states[1]),after))
    if before.shape == (5, 3, 3):
        rate = (after[4] - before[4]) / width
    else:
        rate = (sixth_quadrupole(binary,states[-1])-sixth_quadrupole(binary,states[0]))/width
    result[1,2] = 6*rate-4*before[2]-2*after[2]
    result[2,2] = 3*(before[2]+after[2])-6*rate
    return result


@dataclass(frozen=True)
class PrimitiveIncomingTail:
    """Retarded prehistory for Q7 as the time derivative of the incoming Q6."""

    base: object
    orbit: object
    binary: object
    nodes: int=64

    @property
    def cutoff(self):
        return self.orbit.cutoff

    @cached_property
    def _quadrature(self):
        distance = self.orbit.time_origin-self.cutoff
        legendre, weights = np.polynomial.legendre.leggauss(self.nodes)
        legendre, weights = (1+legendre)/2,weights/2
        near = distance*legendre
        far = distance/(1-legendre)
        offsets = np.r_[near,far]
        weights = np.r_[distance*weights,distance*weights/(1-legendre)**2]

        def value(time):
            return sixth_quadrupole(self.binary,np.concatenate(self.orbit.state_at(self.binary,time)))

        first = value(self.cutoff)
        source = np.stack([value(self.cutoff-delay) for delay in offsets])
        source[:self.nodes] -= first
        return distance,offsets,weights,source,first

    def at_cutoff(self):
        """Expose the original Newtonian jet used to match the SH phase state."""
        return self.base.at_cutoff()

    def integral(self,time,scale,power=1):
        ordinary = self.base.integral(time,scale,power).copy()
        distance,offsets,weights,source,first = self._quadrature
        if power==0:
            ordinary[2] = first
            return ordinary
        lags = (time-self.cutoff)+offsets
        kernel = power*np.log(lags/scale)**(power-1)/lags
        boundary = np.log(((time-self.cutoff)+distance)/scale)**power
        ordinary[2] = boundary*first+np.tensordot(weights*kernel,source,axes=1)
        return ordinary


def dynamical_tail(binary, time, state, history, trial=None, *, split_time=None):
    """Relative acceleration of Section 4, on the covered finite history."""
    position, velocity = state[:3], state[3:]
    radius = np.linalg.norm(position)
    normal = position / radius
    s = normal @ velocity
    w = velocity @ velocity
    gm = binary.G * binary.mass
    u = gm / radius
    scale = 2*radius/binary.c
    if split_time is None:
        log_integral = history.integral(time,scale,trial=trial)[2]
    else:
        from .tail_parts import fourfold_tail_integral
        log_integral = fourfold_tail_integral(history,time,scale,split_time,trial)
    retarded = -8 * binary.G * gm / (5 * binary.c**8) * (log_integral @ position)
    local = 16 * gm**3 * binary.nu / (15 * binary.c**8 * radius**4) * (
        normal * (-6 * u**2 + 70 * u * s**2 - 36 * u * w + 608 * s**4 - 471 * s**2 * w + 30 * w**2)
        - 3 * s * velocity * (37 * u + 136 * s**2 - 63 * w)
    )
    return retarded + local


def radiative_tail(binary, time, history, radiative_scale):
    """Leading tail of the mass quadrupole and its news on the same history."""
    scale = 2 * radiative_scale
    tail = history.integral(time, scale)[:2]
    tail += (11 / 12) * history.integral(time, scale, power=0)[:2]
    # Differentiating the finite-history integral produces a start term.
    if time > history.start and history.panels:
        boundary = history.panels[0].coefficients[0, 0]
        if history.incoming is not None:
            boundary = boundary - history.incoming.at_cutoff()[0]
        tail[1] += boundary * (np.log((time - history.start) / scale) + 11 / 12)
    return (2 * binary.G * binary.mass / binary.c**3) * tail
