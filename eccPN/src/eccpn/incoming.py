"""Incoming hyperbolic phase space and its quadrupole asymptote."""

from dataclasses import asdict, dataclass
from functools import cache, cached_property
from hashlib import sha256
import json
from math import comb, exp, factorial, fsum, log, log1p

import numpy as np
from scipy.integrate import quad
from scipy.signal import convolve2d
from scipy.special import polygamma

from .generated import sh_local as conservative


def incoming_parameters(binary, orbit):
    """JSON-safe physical inputs, independent of scalar and vector container type."""
    return {"type": type(orbit).__name__, "binary": asdict(binary), "orbit": asdict(orbit)}


def _normalize_orbit(orbit, vector_names):
    # Freeze numeric inputs at construction, including mutable NumPy vectors.
    for name in orbit.__dataclass_fields__:
        value = getattr(orbit, name)
        if name in vector_names:
            value = tuple(float(x) + 0. for x in value)
        elif name == "asymptotic_order":
            if int(value) != value:
                raise ValueError("The incoming asymptotic_order must be an integer.")
            value = int(value)
        else:
            value = float(value) + 0.
        object.__setattr__(orbit, name, value)


def incoming_identity(binary, orbit):
    payload = json.dumps(incoming_parameters(binary, orbit), sort_keys=True,
                         separators=(",", ":"), allow_nan=False).encode("ascii")
    return "eccpn.incoming.v1:" + sha256(payload).hexdigest()


def incoming_grid(order,cutoff,origin):
    """Positive physical-time weights on a compactified incoming interval."""
    nodes,weights = np.polynomial.legendre.leggauss(order)
    angle = (1+nodes)*np.pi/4
    q = np.sin(angle)**2
    duration = origin-cutoff
    times = origin-duration/q
    weights *= np.sin(2*angle)*np.pi/4*duration/q**2
    return times,weights


@cache
def _log_beta_series(a,source_power,kernel_power):
    """Taylor coefficients of log[B(a-u-v,1+v)/B(a,1)]."""
    result = np.zeros((source_power+1,kernel_power+1))
    for i in range(source_power+1):
        for j in range(kernel_power+1):
            order = i+j
            if not order:
                continue
            coefficient = (-1)**order*polygamma(order-1,a)/(factorial(i)*factorial(j))
            if i == 0:
                coefficient += polygamma(j-1,1.)/factorial(j)
            if j == 0:
                coefficient -= (-1)**i*polygamma(i-1,a+1)/factorial(i)
            result[i,j] = coefficient
    return result


def _join_log_integral(a, source_log, lag_log, source_power, kernel_power):
    """Mixed log moments at the matching cut via Euler's beta identity."""
    unit = _log_beta_series(a,source_power,kernel_power).copy()
    if source_power:
        unit[1,0] += source_log
    if kernel_power:
        unit[0,1] += lag_log
    power = np.zeros_like(unit)
    power[0,0] = 1.
    result = power.copy()
    for n in range(1,source_power+kernel_power+1):
        power = convolve2d(power,unit)[:source_power+1,:kernel_power+1]/n
        result += power
    return result[-1,-1]*factorial(source_power)*factorial(kernel_power)/a


def remote_log_moment(cutoff, origin, time, source_scale, lag_scale, decay, source_log_power, kernel_power):
    """Integrate one remote inverse-power/log monomial in the true retarded lag."""
    if (not np.isfinite([cutoff, origin, time, source_scale, lag_scale, decay]).all()
            or not cutoff < origin or time < cutoff or min(source_scale, lag_scale) <= 0 or decay <= 1
            or source_log_power not in range(5) or kernel_power not in range(3)):
        raise ValueError("Supply a covered remote cut, positive scales and supported inverse powers and logarithms.")
    duration = origin - cutoff
    ratio = (time - origin) / duration
    source_log = log(duration / source_scale)
    lag_log = log(duration / lag_scale)

    def primitive(index, extra_log):
        # After q=T_cut/T the weight is q**(decay-2+index).
        exponent = decay - 1 + index
        return fsum(
            coefficient * factorial(j + k) / exponent**(j + k + 1)
            for j in range(source_log_power + 1)
            for k in range(extra_log + 1)
            for coefficient in (comb(source_log_power, j) * source_log**(source_log_power - j)
                                * comb(extra_log, k) * lag_log**(extra_log - k),)
        )

    if kernel_power == 0 or ratio == 0:
        value = primitive(0, kernel_power)
    elif time == cutoff:
        value = _join_log_integral(decay-1,source_log,lag_log,source_log_power,kernel_power)
    elif kernel_power == 1 and source_log_power == 0 and float(decay).is_integer() and abs(ratio) > 0.65:
        degree = int(decay) - 2
        endpoint = log1p(ratio)
        reduced = endpoint / ratio
        for k in range(1, degree + 2):
            reduced = (1/k - reduced) / ratio
        varying_lag = (endpoint - ratio * reduced) / (degree + 1)
        value = primitive(0, 1) + varying_lag
    elif abs(ratio) <= 0.65:
        # The logarithm of the actual lag contributes log(1+ratio*q).
        # The log-square series has coefficients 2*(-1)**n*H[n-1]/n.
        harmonic = 0.
        terms = [primitive(0, kernel_power)]
        for n in range(1, 145):
            signed = (-ratio)**n / n
            term = -kernel_power * signed * primitive(n, kernel_power - 1)
            if kernel_power == 2 and n > 1:
                term += 2 * signed * harmonic * primitive(n, 0)
            terms.append(term)
            harmonic += 1/n
        value = fsum(terms)
    else:
        def integrand(q, varying_lag):
            logarithm = log(q)
            return (q**(decay - 2) * (source_log - logarithm)**source_log_power
                    * (lag_log - logarithm + varying_lag)**kernel_power)

        def direct(q):
            return integrand(q, log1p(ratio * q))

        if ratio < 0:
            distance = (time - cutoff) / duration
            left = quad(direct, 0., .5, epsabs=4e-13, epsrel=4e-13)[0]
            # Near the join the original q integral has a corner in a
            # layer of width distance. Resolve the layer using the exact
            # logarithm as the integration variable, including its Jacobian.
            span = 1 - distance
            def right(log_lag):
                y = exp(-log_lag)
                q = 1 - (y - distance) / span
                return y / span * integrand(q, -log_lag)

            upper = quad(right, -log((1 + distance) / 2), -log(distance),
                         epsabs=4e-13, epsrel=4e-13)[0]
            value = left + upper
        else:
            value = quad(direct, 0., 1., epsabs=4e-13, epsrel=4e-13)[0]
    return duration**(1 - decay) * value


def _is_stf(tensors):
    # Scale each tensor separately so the roundoff allowance has its physical units.
    tolerance = 32 * np.finfo(float).eps * np.max(np.abs(tensors), axis=(-2, -1))
    skew = np.max(np.abs(tensors - tensors.swapaxes(-1, -2)), axis=(-2, -1))
    return np.all(skew <= tolerance) and np.all(np.abs(np.trace(tensors, axis1=-1, axis2=-2)) <= tolerance)


@dataclass(frozen=True)
class IncomingQuadrupole:
    """Inverse-power and logarithmic asymptote of I^(4), I^(5) and I^(7)."""

    cutoff: float
    time_origin: float
    tensor: np.ndarray
    subleading: np.ndarray | None = None
    log_reference: float = 1.0
    orbit: object = None
    binary: object = None

    def __post_init__(self):
        tensor = np.array(self.tensor, dtype=float, copy=True)
        if (not np.isfinite([self.cutoff, self.time_origin]).all() or self.cutoff >= self.time_origin
                or tensor.shape != (3, 3) or not np.isfinite(tensor).all()):
            raise ValueError("An incoming asymptote requires a finite cutoff, origin and Cartesian STF tensor.")
        if not _is_stf(tensor):
            raise ValueError("The incoming quadrupole must be symmetric and trace-free.")
        tensor.flags.writeable = False
        object.__setattr__(self, "tensor", tensor)
        if not np.isfinite(self.log_reference) or self.log_reference <= 0:
            raise ValueError("The incoming log reference is a positive time.")
        if self.subleading is not None:
            correction = np.array(self.subleading, dtype=float, copy=True)
            if (correction.shape != (3, 2, 3, 3) or not np.isfinite(correction).all()
                    or not _is_stf(correction)):
                raise ValueError("Supply constant and logarithmic STF coefficients for each incoming derivative.")
            correction.flags.writeable = False
            object.__setattr__(self, "subleading", correction)

    @classmethod
    def from_velocity(cls, binary, cutoff, past_velocity, time_origin=0.0):
        """Construct the leading terms from the nonzero incoming relative velocity."""
        velocity = np.asarray(past_velocity, dtype=float)
        speed = np.linalg.norm(velocity)
        if velocity.shape != (3,) or not np.isfinite(speed) or speed <= 0:
            raise ValueError("A hyperbolic incoming state needs a finite nonzero relative velocity.")
        basis = np.outer(velocity, velocity) - np.eye(3) * speed**2 / 3
        return cls(cutoff, time_origin, binary.mu * binary.G * binary.mass / speed**3 * basis)

    def at_cutoff(self):
        duration = self.time_origin - self.cutoff
        value = np.array([4 / duration**3, 12 / duration**4, 240 / duration**6])[:, None, None] * self.tensor
        if self.subleading is not None:
            coefficients = self.subleading[:, 0] + log(duration / self.log_reference) * self.subleading[:, 1]
            value += coefficients * duration**-np.array([4, 5, 7], dtype=float)[:, None, None]
        return value

    def integral(self, time, scale, power):
        """Integrate the analytic prefix up to cutoff on a supported later cut."""
        if not np.isfinite([time, scale]).all() or time < self.cutoff or scale <= 0 or power not in (0, 1, 2):
            raise ValueError("The incoming history integral requires a later cut, positive scale and power 0–2.")
        def moment(decay, log_power):
            return remote_log_moment(
                self.cutoff, self.time_origin, time, self.log_reference, scale, decay, log_power, power,
            )

        results = [factor * moment(decay, 0) for decay, factor in ((3, 4), (4, 12), (6, 240))]
        result = np.asarray(results)[:, None, None] * self.tensor
        if self.subleading is not None:
            for index, decay in enumerate((4, 5, 7)):
                result[index] += moment(decay, 0) * self.subleading[index, 0]
                result[index] += moment(decay, 1) * self.subleading[index, 1]
        return result


def _incoming_pn_powers(binary):
    """Denominators for the Newtonian through 4PN conservative coefficients."""
    return binary.c ** (2 * np.arange(5))


def _incoming_phase_space(
    binary, velocity, impact, duration, log_reference, asymptotic_order, *, include_acceleration=False,
):
    speed = np.linalg.norm(velocity)
    direction = velocity / speed
    log_time = log(duration / log_reference)
    g = binary.G * binary.mass
    values = np.asarray(conservative.incoming_conservative(speed, binary.nu))
    values /= _incoming_pn_powers(binary)[:, None]
    radial, projected, derivative_u = values.T
    drift = g / speed**3 * sum(projected)
    position = -duration * velocity - drift * velocity * log_time + impact
    rate = velocity + drift * velocity / duration
    if asymptotic_order == 1:
        return (position,rate,drift*velocity/duration**2) if include_acceleration else (position,rate)

    # The correction to x(-T) is (p log(T/t*) + q)/T. Discard
    # conservative products with explicit powers beyond 4PN.
    parallel = impact @ direction
    transverse = impact - parallel * direction
    rate_derivative = sum(np.convolve(projected, 2 * np.arange(5) * projected)[:5])
    p = -g**2 / speed**5 * sum(np.convolve(projected, projected)[:5]) * direction
    constant = g / speed**2 * (
        (g / speed * derivative_u.sum() + g / speed**3 * rate_derivative
         + 2 * parallel / speed * projected.sum()) * direction - radial.sum() / speed * transverse
    )
    args = (0.0, -speed, speed**2, binary.nu)
    constant += binary.nu * g**2 / speed**3 * direction * (
        16 / 5 * speed**3 / binary.c**5
        + (-conservative.rr35n(*args) + speed * conservative.rr35v(*args)) / binary.c**7
        + (-conservative.rr45mattern(*args) + speed * conservative.rr45matterv(*args)) / binary.c**9
    )
    q = (constant + 3 * p) / 2
    polynomial = p * log_time + q
    result = position + polynomial/duration, rate + (polynomial-p)/duration**2
    if include_acceleration:
        acceleration = drift*velocity/duration**2+(2*polynomial-3*p)/duration**3
    if asymptotic_order == 3:
        newtonian = _newtonian_hyperbolic(binary,velocity,impact,duration,log_reference)
        # Sum the Newtonian remainder exactly while keeping the independently
        # matched relativistic inverse-time coefficients in SH.
        x_n = -duration*velocity - g/speed**3*velocity*log_time + impact
        v_n = velocity + g/speed**3*velocity/duration
        p_n = -g**2/speed**5*direction
        q_n = (g/speed**3*(2*parallel*direction-transverse)+3*p_n)/2
        polynomial_n = p_n*log_time+q_n
        x_n += polynomial_n/duration
        v_n += (polynomial_n-p_n)/duration**2
        result = newtonian[0]+(result[0]-x_n), newtonian[1]+(result[1]-v_n)
        if include_acceleration:
            acceleration += (-g*newtonian[0]/np.linalg.norm(newtonian[0])**3
                             -g/speed**3*velocity/duration**2-(2*polynomial_n-3*p_n)/duration**3)
    return (*result,acceleration) if include_acceleration else result


def _newtonian_hyperbolic(binary, velocity, impact, duration, log_reference):
    """Exact Cartesian Kepler hyperbola at fixed incoming Coulomb data."""
    from scipy.optimize import brentq

    speed = np.linalg.norm(velocity)
    direction = velocity/speed
    parallel = impact@direction
    transverse = impact-parallel*direction
    a = binary.G*binary.mass/speed**2
    beta = np.linalg.norm(transverse)/a
    eccentricity = np.hypot(1.,beta)
    axis = (direction+transverse/a)/eccentricity
    if beta:
        tangent = (beta*direction-transverse/(a*beta))/eccentricity
    else:
        tangent = np.zeros(3)
    offset = a/speed*(1-np.log(2*speed*log_reference/(a*eccentricity)))-parallel/speed
    mean = (-duration-offset)*speed/a
    if mean >= 0:
        raise ValueError("The incoming cut must precede the Newtonian periastron.")
    lower = np.arcsinh(mean/eccentricity)-2
    anomaly = brentq(lambda u: eccentricity*np.sinh(u)-u-mean,lower,0.,xtol=1e-13)
    cosine, sine = np.cosh(anomaly), np.sinh(anomaly)
    position = a*((eccentricity-cosine)*axis+beta*sine*tangent)
    rate = speed*(-sine*axis+beta*cosine*tangent)/(eccentricity*cosine-1)
    return position,rate


def _remote_correction(binary, velocity, impact):
    """Expand the order-reduced Newtonian jets on the chosen incoming SH worldline."""
    from .tails import _stf

    speed = np.linalg.norm(velocity)
    direction = velocity / speed
    parallel = direction @ impact
    axis = _stf(direction, direction)
    transverse = _stf(direction, impact - parallel * direction)
    g = binary.G * binary.mass
    projection = np.asarray(conservative.incoming_conservative(speed, binary.nu))[:, 1]
    drift = g / speed**3 * np.sum(projection / _incoming_pn_powers(binary))
    coefficients = []
    # At s=-v and u=0, the radial jet coefficient is B*v**(n-2).
    # The next two columns are its potential derivative and the response
    # of the tensor basis to the incoming transverse displacement.
    for n, b, potential_coefficient, transverse_coefficient in ((4, 4, 2, 12), (5, 12, 16, 48), (7, 240, 888, 1440)):
        base = b * speed**(n - 2)
        scale = binary.mu * g / speed**(n - 1)
        constant = (base * ((n - 2) * drift + (n - 1) * parallel / speed)
                    + g * potential_coefficient * speed**(n - 5))
        coefficients.append([
            scale * (constant * axis + transverse_coefficient * speed**(n - 3) * transverse),
            -scale * (n - 1) * base * drift * axis,
        ])
    return coefficients


@dataclass(frozen=True)
class IncomingOrbit:
    """Unbound incoming relative velocity and impact offset in SH source time."""

    velocity: tuple[float, float, float]
    impact: tuple[float, float, float]
    cutoff: float
    log_reference: float = 1.0
    time_origin: float = 0.0
    asymptotic_order: int = 2

    def __post_init__(self):
        _normalize_orbit(self, ('velocity', 'impact'))

    def state_at(self, binary, time, *, include_acceleration=False):
        """Evaluate the chosen SH asymptote in its finite-to-remote domain."""
        if not np.isfinite(time) or time > self.cutoff or not time < self.time_origin:
            raise ValueError("The incoming phase expansion covers finite cuts at and before its join.")
        return _incoming_phase_space(
            binary, np.asarray(self.velocity, dtype=float), np.asarray(self.impact, dtype=float), self.time_origin-time,
            self.log_reference, self.asymptotic_order, include_acceleration=include_acceleration,
        )

    def conserved_charges(self, binary):
        """Conserved mass through c^-4 and current dipole through c^-2 at infinity."""
        velocity, impact = np.asarray(self.velocity, dtype=float), np.asarray(self.impact, dtype=float)
        if velocity.shape != (3,) or impact.shape != (3,) or not np.isfinite([velocity, impact]).all():
            raise ValueError("Incoming velocity and impact must be finite Cartesian vectors.")
        w = velocity @ velocity / binary.c**2
        if not 0 < w < 1:
            raise ValueError("Use a nonzero subluminal incoming relative speed.")
        mass = {0: binary.mass, 2: binary.mu * w / 2,
                4: 3 * binary.mu * (1 - 3 * binary.nu) * w**2 / 8}
        dipole = binary.mu * np.cross(impact, velocity)
        return mass, {0: dipole, 2: (1 - 3 * binary.nu) * w * dipole / 2}

    def force_prefix(self, binary):
        if self.asymptotic_order == 3:
            return SampledIncomingTail(self, binary)
        velocity = np.asarray(self.velocity, dtype=float)
        impact = np.asarray(self.impact, dtype=float)
        leading = IncomingQuadrupole.from_velocity(binary, self.cutoff, velocity, self.time_origin)
        extra = _remote_correction(binary, velocity, impact) if self.asymptotic_order >= 2 else None
        return IncomingQuadrupole(self.cutoff, self.time_origin, leading.tensor, extra,
                                  self.log_reference, self, binary)

    def prepare(self, binary):
        """Construct the leading physical inner and outer matching data."""
        from .conventions import InitialState
        from .history import PanelHistory, TailInitialData
        from .tails import quadrupole_tail_sources

        velocity = np.asarray(self.velocity, dtype=float)
        impact = np.asarray(self.impact, dtype=float)
        speed = np.linalg.norm(velocity)
        duration = self.time_origin - self.cutoff
        if (velocity.shape != (3,) or impact.shape != (3,) or not np.isfinite(impact).all()
                or not np.isfinite([speed, duration, self.log_reference]).all()
                or not 0 < speed < binary.c or duration <= 0 or self.log_reference <= 0
                or self.asymptotic_order not in (1, 2, 3)):
            raise ValueError("Supply finite, subluminal incoming data with a positive remote duration.")
        g, c, nu = binary.G * binary.mass, binary.c, binary.nu
        position, evolved_velocity = self.state_at(binary, self.cutoff)
        remote = self.force_prefix(binary)
        phase_space = np.concatenate((position, evolved_velocity))
        at_cut = quadrupole_tail_sources(binary, phase_space)
        difference = np.linalg.norm(at_cut - remote.at_cutoff(), axis=(-2, -1))
        normalization = np.linalg.norm(remote.at_cutoff(), axis=(-2, -1))
        residual = tuple(difference / normalization)
        history = PanelHistory(self.cutoff, (3, 3, 3), incoming_identity(binary, self), incoming=remote)

        delta = (binary.m1 - binary.m2) / binary.mass
        prefactor = binary.mass * delta * nu**2 / c**7
        # The coefficients below multiply T^-4 (momentum flux) or T^-2
        # (centre flux) on the incoming straight inertial asymptote.
        leading_flux = -16 / 105 * prefactor * g**3 * velocity / speed**2
        next_flux = 8 * 12 / 1890 * (13 - 4 * nu) * prefactor * g**3 * velocity / c**2
        center_flux = -16 / 35 * prefactor * g**2 * speed * velocity
        recoil = leading_flux / (3 * duration**3)
        next_recoil = next_flux / (3 * duration**3)
        center = center_flux / duration + (leading_flux + next_flux) / (6 * duration**2)
        if isinstance(remote,SampledIncomingTail):
            recoil,next_recoil,center = remote.recoil()
        return InitialState(tuple(position), tuple(evolved_velocity)), TailInitialData(
            history, tuple(recoil), tuple(next_recoil), tuple(center),
            asymptotic_join_residual=residual, incoming_orbit=self,
        )


@dataclass(frozen=True)
class SampledIncomingTail:
    """Causal far-past force and recoil quadrature on the improved phase expansion."""

    orbit: IncomingOrbit
    binary: object
    nodes: int = 48

    @property
    def cutoff(self):
        return self.orbit.cutoff

    @cached_property
    def samples(self):
        from .tails import quadrupole_tail_sources

        times,weights = incoming_grid(self.nodes,self.cutoff,self.orbit.time_origin)
        values = np.stack([np.concatenate(self.orbit.state_at(self.binary,time)) for time in times])
        tensors = np.stack([quadrupole_tail_sources(self.binary,state) for state in values])
        return times,weights,values,tensors

    def at_cutoff(self):
        from .tails import quadrupole_tail_sources

        return quadrupole_tail_sources(self.binary,np.concatenate(self.orbit.state_at(self.binary,self.cutoff)))

    def integral(self,time,scale,power):
        if not np.isfinite([time,scale]).all() or time < self.cutoff or scale <= 0 or power not in range(3):
            raise ValueError("The incoming force integral requires a covered cut and positive scale.")
        times,weights,_,tensors = self.samples
        return np.tensordot(weights*np.log((time-times)/scale)**power,tensors,axes=1)

    def recoil(self):
        from .recoil import recoil_fluxes

        times,weights,values,_ = self.samples
        leading,next_order,center = recoil_fluxes(self.binary,values[:,:3],values[:,3:])
        return (weights @ leading, weights @ next_order,
                weights @ (center+(self.cutoff-times)[:,None]*(leading+next_order)))


@dataclass(frozen=True)
class ParabolicQuadrupole:
    """Leading Q^(4), Q^(5) and Q^(7) on a marginally bound incoming branch."""

    cutoff: float
    time_origin: float
    tensor: np.ndarray
    orbit: object = None
    binary: object = None

    def __post_init__(self):
        base = IncomingQuadrupole(self.cutoff, self.time_origin, self.tensor)
        object.__setattr__(self, "tensor", base.tensor)

    @staticmethod
    def _powers():
        from math import prod

        return [(n - 4 / 3, (-1)**n * prod(4 / 3 - i for i in range(n))) for n in (4, 5, 7)]

    def at_cutoff(self):
        duration = self.time_origin - self.cutoff
        scales = np.array([factor * duration**-p for p, factor in self._powers()])
        return scales[:, None, None] * self.tensor

    def integral(self, time, scale, power):
        if not np.isfinite([time, scale]).all() or time < self.cutoff or scale <= 0 or power not in (0, 1, 2):
            raise ValueError("The parabolic prefix requires a covered cut and a positive kernel scale.")
        values = [
            factor * remote_log_moment(self.cutoff, self.time_origin, time, 1, scale, p, 0, power)
            for p, factor in self._powers()
        ]
        return np.array(values)[:, None, None] * self.tensor


@dataclass(frozen=True)
class ParabolicOrbit:
    """Leading Newtonian incoming parabolic data, independent of hyperbolic speed."""

    pericenter_distance: float
    pericenter_direction: tuple[float, float, float]
    transverse_direction: tuple[float, float, float]
    cutoff: float
    time_origin: float = 0.0
    asymptotic_order: int = 1

    def __post_init__(self):
        _normalize_orbit(self, ('pericenter_direction', 'transverse_direction'))

    def conserved_charges(self, binary):
        """Use pericenter_distance to parameterize the zero-energy incoming angular momentum."""
        peri, tangent = np.asarray((self.pericenter_direction, self.transverse_direction))
        momentum = binary.mu * np.sqrt(2*binary.G*binary.mass*self.pericenter_distance)*np.cross(peri, tangent)
        return {0: binary.mass, 2: 0., 4: 0.}, {0: momentum, 2: np.zeros(3)}

    def state_at(self, binary, time, *, include_acceleration=False):
        """Evaluate the marginally bound SH asymptote and optional acceleration."""
        radius, duration = self.pericenter_distance, self.time_origin-time
        peri, tangent = np.asarray((self.pericenter_direction, self.transverse_direction), dtype=float)
        if (not np.isfinite([time, radius, duration]).all() or time > self.cutoff or min(radius, duration) <= 0
                or self.asymptotic_order not in (1, 2) or not np.allclose(
                    [peri @ peri, tangent @ tangent, peri @ tangent], [1., 1., 0.], atol=2e-12, rtol=0,
                )):
            raise ValueError("Supply a parabolic remote cut, positive radius, and orthonormal directions.")
        gm = binary.G*binary.mass
        tau = np.sqrt(2*radius**3/gm)
        parameter = 2*np.sinh(np.arcsinh(-3*duration/(2*tau))/3)
        position = radius*((1-parameter**2)*peri + 2*parameter*tangent)
        velocity = 2*radius/(tau*(1+parameter**2)) * (-parameter*peri+tangent)
        if self.asymptotic_order == 1:
            result = (position,velocity)
            return (*result,-gm*position/np.linalg.norm(position)**3) if include_acceleration else result
        return self._conservative_asymptote(
            binary,duration,radius,peri,tangent,position,velocity,include_acceleration=include_acceleration,
        )

    @staticmethod
    def _conservative_asymptote(
        binary, duration, peri_radius, peri, tangent, position, velocity, *, include_acceleration=False,
    ):
        # Expand the D1 SH E and J invariants at E=0 and fixed incoming J.
        gm, inv_c2, nu = binary.G*binary.mass, binary.c**-2, binary.nu
        J = np.sqrt(2*gm*peri_radius)
        angular2 = (J/gm)**2
        w2 = -5*(2-nu)*inv_c2
        first_energy = (6-7*nu)*w2/2 - nu*angular2/2
        second_energy = (68-119*nu+50*nu**2)/4
        w3 = -2*(inv_c2*first_energy + inv_c2**2*second_energy)
        angular1 = (4-2*nu)*inv_c2
        squared_radial2 = w2-angular2
        squared_radial3 = w3+2*angular2*angular1
        lag1 = -squared_radial2/4
        lag2 = 3*squared_radial2**2/32 - squared_radial3/4
        offset = -2*lag1*gm
        inverse = (lag1**2+2*lag2)*gm**2
        far = (9*gm*duration**2/2)**(1/3)
        far_rate = -2*far/(3*duration)
        r_newton = np.linalg.norm(position)
        rn_rate = position @ velocity/r_newton
        r = r_newton + offset+peri_radius + (inverse-peri_radius**2)/far
        r_rate = rn_rate - (inverse-peri_radius**2)*far_rate/far**2
        j_factor = J/np.sqrt(2*gm)
        theta_newton = np.arctan2(position @ tangent, position @ peri)
        theta = theta_newton + j_factor*(
            2*(r**-0.5-r_newton**-0.5) + 2*gm/3*((lag1-angular1)*r**-1.5-angular2/4*r_newton**-1.5)
        )
        theta_rate = J/r_newton**2 + j_factor*(
            -(r**-1.5+gm*(lag1-angular1)*r**-2.5)*r_rate
            +(r_newton**-1.5+gm*angular2/4*r_newton**-2.5)*rn_rate
        )
        normal = np.cos(theta)*peri + np.sin(theta)*tangent
        transverse = -np.sin(theta)*peri + np.cos(theta)*tangent
        result = r*normal, r_rate*normal+r*theta_rate*transverse
        if not include_acceleration:
            return result
        rn_acceleration = J**2/r_newton**3-gm/r_newton**2
        far_acceleration = -2*far/(9*duration**2)
        r_acceleration = rn_acceleration+(inverse-peri_radius**2)*(
            2*far_rate**2/far**3-far_acceleration/far**2
        )
        k = gm*(lag1-angular1)
        ell = gm*angular2/4
        d_radial = (3*r**-2.5+5*k*r**-3.5)*r_rate**2/2-(r**-1.5+k*r**-2.5)*r_acceleration
        d_newtonian = ((r_newton**-1.5+ell*r_newton**-2.5)*rn_acceleration
                       -(3*r_newton**-2.5+5*ell*r_newton**-3.5)*rn_rate**2/2)
        theta_acceleration = -2*J*rn_rate/r_newton**3+j_factor*(d_radial+d_newtonian)
        acceleration = ((r_acceleration-r*theta_rate**2)*normal
                        +(2*r_rate*theta_rate+r*theta_acceleration)*transverse)
        return (*result,acceleration)

    def force_prefix(self, binary):
        from .tails import _stf

        peri = np.asarray(self.pericenter_direction, dtype=float)
        tensor = binary.mu * (9 * binary.G * binary.mass / 2)**(2 / 3) * _stf(peri, peri)
        return ParabolicQuadrupole(self.cutoff, self.time_origin, tensor, self, binary)

    def prepare(self, binary):
        from scipy.integrate import quad_vec

        from .conventions import InitialState
        from .history import PanelHistory, TailInitialData
        from .recoil import recoil_fluxes
        from .tails import quadrupole_tail_sources

        peri, tangent = np.asarray((self.pericenter_direction, self.transverse_direction), dtype=float)
        radius = self.pericenter_distance
        duration = self.time_origin - self.cutoff
        if (peri.shape != (3,) or tangent.shape != (3,) or not np.isfinite([duration, radius]).all()
                or min(duration, radius) <= 0 or not np.allclose([peri @ peri, tangent @ tangent, peri @ tangent],
                                                                 [1., 1., 0.], rtol=0, atol=2e-12)):
            raise ValueError("Supply an earlier cutoff, positive pericenter distance and orthonormal directions.")
        g = binary.G * binary.mass
        tau = np.sqrt(2 * radius**3 / g)
        barker = 2 * np.sinh(np.arcsinh(-3 * duration / (2 * tau)) / 3)

        def state(parameter):
            position = radius * ((1 - parameter**2) * peri + 2 * parameter * tangent)
            velocity = 2 * radius / (tau * (1 + parameter**2)) * (-parameter * peri + tangent)
            return position, velocity

        position, velocity = self.state_at(binary, self.cutoff)
        remote = self.force_prefix(binary)
        actual = quadrupole_tail_sources(binary, np.concatenate((position, velocity)))
        expected = remote.at_cutoff()
        residual = tuple(np.linalg.norm(actual - expected, axis=(-2, -1)) / np.linalg.norm(expected, axis=(-2, -1)))
        history = PanelHistory(self.cutoff, (3, 3, 3), incoming_identity(binary, self), incoming=remote)

        def combined_flux(parameter):
            cut = self.time_origin+tau*(parameter+parameter**3/3)
            xy = state(parameter) if self.asymptotic_order == 1 else self.state_at(binary, min(cut, self.cutoff))
            leading, next_order, center = recoil_fluxes(binary, *xy)
            lag = self.cutoff - self.time_origin - tau * (parameter + parameter**3 / 3)
            return np.stack((leading, next_order, center + lag * (leading + next_order)))

        if binary.m1 == binary.m2:
            integrated = np.zeros((3, 3))
        else:
            row_scale = np.linalg.norm(combined_flux(barker), axis=1)
            row_scale = np.where(row_scale > 0, row_scale, 1.)

            def integrand(q):
                parameter = barker / q
                weight = -barker * tau * (1 + parameter**2) / q**2
                return combined_flux(parameter) / row_scale[:, None] * weight

            scaled, _, info = quad_vec(integrand, 0, 1, epsabs=0, epsrel=1e-10, full_output=True)
            if not info.success:
                raise RuntimeError(f"Parabolic recoil matching quadrature failed: {info.message}")
            integrated = row_scale[:, None] * scaled
        return InitialState(tuple(position), tuple(velocity)), TailInitialData(
            history, *(tuple(row) for row in integrated), asymptotic_join_residual=residual, incoming_orbit=self,
        )


@dataclass(frozen=True)
class IncomingMatch:
    incoming: IncomingOrbit
    prefix: object
    scaled_target_residual: float
    maximum_compactness: float
    maximum_speed_ratio: float
    forward_evaluations: int


def _dense_weak_field(binary, solution):
    from numpy.polynomial import polynomial

    min_distance = np.inf
    max_speed = 0.0
    for coefficients in solution.dense:
        for components, name in ((coefficients[:, :3], "position"), (coefficients[:, 3:6], "velocity")):
            norm_squared = sum(np.convolve(axis, axis) for axis in components.T)
            turning_points = polynomial.polyroots(polynomial.polyder(norm_squared))
            internal = [real for real, imag in zip(turning_points.real, turning_points.imag)
                        if abs(imag) <= 1e-11 and 0 < real < 1]
            values = np.linalg.norm(polynomial.polyval(np.array([0.0, 1.0, *internal]), components), axis=0)
            if name == "position":
                min_distance = min(min_distance, values.min())
            else:
                max_speed = max(max_speed, values.max())
    return binary.G * binary.mass / (binary.c**2 * min_distance), max_speed / binary.c


def match_incoming_state(binary, guess, target, reference_time, accuracy, *, target_tolerance=1e-8, max_evaluations=35):
    """Shoot from a hyperbolic asymptote to a target using direct retarded evolution.

    Changes to the incoming velocity and impact use the supplied cutoff and
    origin convention. The achieved residual and weak-field diagnostics refer
    to this retained prefix and the chosen incoming asymptotic order.
    """
    from scipy.optimize import least_squares

    from .conventions import SH_MODELS
    from .dynamics import evolve

    if accuracy.dynamics not in SH_MODELS or accuracy.canonical_grade is not None:
        raise ValueError("Shooting requires a retarded SH development accuracy with the incoming force history.")
    duration = reference_time - guess.cutoff
    if not np.isfinite([duration, target_tolerance]).all() or min(duration, target_tolerance) <= 0:
        raise ValueError("Supply a target cut after the incoming cutoff and a positive target tolerance.")
    target_state = np.concatenate((target.position, target.velocity))
    if target_state.shape != (6,) or not np.isfinite(target_state).all():
        raise ValueError("The shooting target must have finite Cartesian position and velocity.")
    speed = np.linalg.norm(guess.velocity)
    target_radius = np.linalg.norm(target_state[:3])
    if speed <= 0 or target_radius <= 0:
        raise ValueError("Use a nonzero incoming speed and a nonzero target separation.")
    parameters = np.concatenate((guess.velocity, guess.impact))
    parameter_scale = np.repeat([speed, max(target_radius, duration * speed)], 3)
    residual_scale = np.repeat([target_radius, max(speed, np.linalg.norm(target_state[3:]))], 3)
    evaluations = 0

    def forward(scaled_change):
        nonlocal evaluations
        evaluations += 1
        change = parameters + scaled_change * parameter_scale
        orbit = IncomingOrbit(
            tuple(change[:3]), tuple(change[3:]), guess.cutoff,
            guess.log_reference, guess.time_origin, guess.asymptotic_order,
        )
        initial, past = orbit.prepare(binary)
        prefix = evolve(binary, initial, past, (orbit.cutoff, reference_time), accuracy)
        return orbit, prefix

    def residual(scaled_change):
        _, prefix = forward(scaled_change)
        return (np.concatenate(prefix.state(reference_time)) - target_state) / residual_scale

    result = least_squares(
        residual, np.zeros(6), max_nfev=max_evaluations,
        ftol=target_tolerance * 0.1, xtol=target_tolerance * 0.1, gtol=target_tolerance * 0.1,
    )
    achieved = np.max(np.abs(result.fun))
    if achieved > target_tolerance:
        raise RuntimeError(f"The forward incoming solution missed the scaled target: {achieved:.6g}.")
    orbit, prefix = forward(result.x)
    compactness, speed_ratio = _dense_weak_field(binary, prefix._solution)
    return IncomingMatch(orbit, prefix, achieved, compactness, speed_ratio, evaluations)
