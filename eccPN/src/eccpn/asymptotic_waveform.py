"""Analytic remote terms for the Newtonian mass quadrupole and its products."""

from dataclasses import dataclass
from functools import cached_property

import numpy as np
from scipy.integrate import quad_vec

from .incoming import incoming_grid, remote_log_moment
from .tails import _stf


def _add(target, key, value):
    target[key] = target.get(key, 0) + value


def _product(left, right):
    product = left @ right
    return (product + product.T)/2 - np.trace(product)*np.eye(3)/3


@dataclass(frozen=True)
class AsymptoticQuadrupole:
    """Cartesian quadrupole represented as powers of T and log(T/T_ref)."""

    cutoff: float
    origin: float
    log_reference: float
    jets: tuple

    @property
    def jet_count(self):
        return len(self.jets)

    @classmethod
    def from_terms(cls, cutoff, origin, log_reference, terms, highest_derivative=9):
        if (not np.isfinite([cutoff, origin, log_reference]).all() or cutoff >= origin
                or log_reference <= 0 or highest_derivative < 0):
            raise ValueError("Supply a covered inverse-time interval and a positive log reference.")
        inventory = {key: np.asarray(value, dtype=float) for key, value in terms.items()}
        jets = [inventory]
        for _ in range(highest_derivative):
            derivative = {}
            for (decay, power), field in inventory.items():
                if decay:
                    _add(derivative, (decay+1, power), decay*field)
                if power:
                    _add(derivative, (decay+1, power-1), -power*field)
            inventory = derivative
            jets.append(inventory)
        return cls(cutoff, origin, log_reference, tuple(jets))

    @classmethod
    def hyperbolic_newtonian(cls, binary, incoming):
        """Coulomb-log expansion of the Newtonian hyperbola through T^0 in I."""
        velocity = np.asarray(incoming.velocity)
        impact = np.asarray(incoming.impact)
        speed = np.linalg.norm(velocity)
        if not speed > 0:
            raise ValueError("A hyperbolic leading quadrupole requires nonzero incoming speed.")
        direction = velocity/speed
        parallel = impact @ direction
        transverse = impact-parallel*direction
        g = binary.G*binary.mass
        drift = g/speed**3
        position_log = -g**2/speed**5 * direction
        position_inverse = (g/speed**3 * (2*parallel*direction-transverse) + 3*position_log)/2
        mass = binary.mu
        vv, vb = _stf(velocity, velocity), _stf(velocity, impact)
        terms = {
            (-2, 0): mass*vv,
            (-1, 0): -2*mass*vb,
            (-1, 1): 2*mass*drift*vv,
            (0, 0): mass*(_stf(impact, impact)-2*_stf(velocity, position_inverse)),
            (0, 1): -2*mass*(drift*vb+_stf(velocity, position_log)),
            (0, 2): mass*drift**2*vv,
        }
        return cls.from_terms(incoming.cutoff, incoming.time_origin, incoming.log_reference, terms)

    @classmethod
    def parabolic_newtonian(cls, binary, incoming):
        """Fractional powers of the Newtonian parabolic quadrupole."""
        radial, tangent = np.asarray((incoming.pericenter_direction, incoming.transverse_direction))
        far = (9*binary.G*binary.mass/2)**(1/3)
        peri = incoming.pericenter_distance
        rr, rt, tt = _stf(radial, radial), _stf(radial, tangent), _stf(tangent, tangent)
        terms = {
            (-4/3, 0): binary.mu*far**2*rr,
            (-1, 0): 4*binary.mu*np.sqrt(peri)*far**1.5*rt,
            (-2/3, 0): binary.mu*far*peri*(-6*rr+4*tt),
            (-1/3, 0): -16*binary.mu*np.sqrt(far)*peri**1.5*rt,
        }
        return cls.from_terms(incoming.cutoff, incoming.time_origin, 1., terms)

    def value(self, time, derivative):
        if not np.isfinite(time) or time > self.cutoff:
            raise ValueError("The asymptotic moments cover only the remote past.")
        duration = self.origin-time
        logarithm = np.log(duration/self.log_reference)
        return sum((duration**-p*logarithm**q*field for (p, q), field in self.jets[derivative].items()),
                   np.zeros((3, 3)))

    def _integral(self, time, scale, power, terms, cutoff):
        return sum((
            field*remote_log_moment(cutoff, self.origin, time, self.log_reference, scale, p, q, power)
            for (p, q), field in terms.items()
        ), np.zeros((3, 3)))

    def linear_integral(self, time, scale, power, derivative):
        """Integrate the derivative across the remote part of a retarded cut."""
        return self._integral(time, scale, power, self.jets[derivative], min(self.cutoff, time))

    def product_integral(self, time, scale, power, left, right):
        """Integrate a quadratic STF source before the finite-history join."""
        terms = {}
        for (p, q), first in self.jets[left].items():
            for (a, b), second in self.jets[right].items():
                _add(terms, (p+a, q+b), _product(first, second))
        return self._integral(time, scale, power, terms, min(self.cutoff, time))

    def nested_remote(self, scale, offset, *, rtol=1e-9, atol=1e-30):
        """Tail-of-memory sector with both ordered integration times remote."""
        duration = self.origin-self.cutoff

        def integrand(q):
            time = self.origin-duration/q
            source = self.value(time, 4)
            inner = self.linear_integral(time, scale, 1, 4)
            inner += offset*self.linear_integral(time, scale, 0, 4)
            return duration/q**2 * _product(source, inner)

        result, _, info = quad_vec(integrand, 0., 1., epsrel=rtol, epsabs=atol, full_output=True)
        if not info.success:
            raise RuntimeError(f"The remote ordered integral failed: {info.message}")
        return result


@dataclass(frozen=True)
class MatchedQuadrupoleHistory:
    """Combine a remote Newtonian quadrupole and its finite numerical panels."""

    finite: object
    remote: object

    def __post_init__(self):
        if (self.finite.start != self.remote.cutoff or self.finite.incoming is not None
                or self.finite.source_shape != (self.remote.jet_count, 3, 3)):
            raise ValueError("Remote and finite canonical jets must have the same order and join.")

    @property
    def start(self):
        return self.finite.start

    @property
    def end(self):
        return self.finite.end

    @property
    def source_shape(self):
        return self.finite.source_shape

    @property
    def panels(self):
        return self.finite.panels

    def value(self, time):
        if time >= self.start:
            return self.finite.value(time)
        return np.array([self.remote.value(time, k) for k in range(self.source_shape[0])])

    def integrate_jet(self, time, scale, power, orders):
        derivatives = range(self.source_shape[0])[orders]
        if isinstance(orders, slice):
            past = np.stack([self.remote.linear_integral(time, scale, power, k) for k in derivatives])
        else:
            past = self.remote.linear_integral(time, scale, power, derivatives)
        if time >= self.start:
            past += self.finite.integrate_jet(time, scale, power, orders)
        return past

    def start_jump(self, order):
        return self.finite.start_jump(order) - self.remote.value(self.start, order)


@dataclass(frozen=True)
class SampledWaveformQuadrupole:
    """Analytic leading kernels plus a Newtonian quadrupole remainder."""

    binary: object
    orbit: object
    asymptote: AsymptoticQuadrupole
    nodes: int

    @property
    def cutoff(self):
        return self.orbit.cutoff

    @property
    def jet_count(self):
        return self.asymptote.jet_count

    def _jet(self,time):
        from .newtonian_moments import newtonian_source_jets

        return newtonian_source_jets(self.binary,*self.orbit.state_at(self.binary,time),
                                     highest_derivative=self.jet_count-1,highest_rank=2)["M2"][0]

    def value(self,time,derivative):
        return self._jet(time)[derivative]

    def _quadrature(self,cut,start,end=None):
        times,weights = incoming_grid(self.nodes,cut,self.orbit.time_origin)
        jets = np.stack([self._jet(t)[start:end] for t in times])
        analytic = np.stack([np.stack([self.asymptote.value(t,k)
                                      for k in range(start,self.jet_count if end is None else end)]) for t in times])
        return times,weights,jets,analytic

    @cached_property
    def samples(self):
        return self._quadrature(self.cutoff,0)

    def linear_integral(self,time,scale,power,derivative):
        analytic = self.asymptote.linear_integral(time,scale,power,derivative)
        if time >= self.cutoff:
            times,weights,values,baseline = self.samples
            remainder = values[:,derivative]-baseline[:,derivative]
        else:
            times,weights,values,baseline = self._quadrature(time,derivative,derivative+1)
            remainder = values[:,0]-baseline[:,0]
        return analytic+np.tensordot(weights*np.log((time-times)/scale)**power,remainder,axes=1)

    def product_integral(self,time,scale,power,left,right):
        if time < self.cutoff:
            raise ValueError("The sampled product integral requires a cut after the numeric join.")
        analytic = self.asymptote.product_integral(time,scale,power,left,right)
        times,weights,jet,baseline = self.samples
        residual = sum(w*(_product(new[left],new[right])-_product(old[left],old[right]))
                       for w,new,old in zip(weights*np.log((time-times)/scale)**power,jet,baseline))
        return analytic+residual

    @cached_property
    def _ordered_integrals(self):
        return {}

    def nested_remote(self,scale,offset,*,rtol,atol):
        key = (scale,offset,rtol,atol)
        if key in self._ordered_integrals:
            return self._ordered_integrals[key].copy()
        # Keep the independently validated analytic leading integral. The
        # regular remainder converges on the common compactified grid.
        reference = self.asymptote.nested_remote(scale,offset,rtol=rtol,atol=atol)
        times,weights,jet,baseline = self.samples
        result = reference
        for t,weight,actual,leading in zip(times,weights,jet[:,4],baseline[:,4]):
            outer = self.asymptote.linear_integral(t,scale,1,4)+offset*self.asymptote.linear_integral(t,scale,0,4)
            inner_times,inner_weights,values,lower = self._quadrature(t,4,5)
            residual = np.tensordot(inner_weights*(np.log((t-inner_times)/scale)+offset),
                                    values[:,0]-lower[:,0],axes=1)
            result += weight*(_product(actual,outer+residual)-_product(leading,outer))
        self._ordered_integrals[key] = result.copy()
        return result
