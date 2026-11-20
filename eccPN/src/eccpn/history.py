"""Causal logarithmic integrals of polynomial source panels."""

from bisect import bisect_left
from dataclasses import dataclass
from functools import cache, cached_property
from math import comb, fsum, log

import numpy as np

from .incoming import IncomingOrbit, IncomingQuadrupole, ParabolicOrbit, ParabolicQuadrupole, SampledIncomingTail


@dataclass(frozen=True)
class SourcePanel:
    """One committed polynomial with coefficients in the panel coordinate."""

    start: float
    end: float
    coefficients: np.ndarray

    def __post_init__(self):
        coefficients = np.array(self.coefficients, dtype=float, copy=True)
        if not np.isfinite([self.start, self.end]).all() or self.end <= self.start:
            raise ValueError("Source panels require a finite, positive time interval.")
        if coefficients.ndim == 0 or not coefficients.shape[0] or not np.isfinite(coefficients).all():
            raise ValueError("Source panels require a nonempty array of finite coefficients.")
        coefficients.flags.writeable = False
        object.__setattr__(self, "coefficients", coefficients)

    def integral(self, time, scale, power):
        return integrate_log_panel(self.start, self.end, self.coefficients, time, scale, power)


@dataclass(frozen=True)
class PanelHistory:
    """An immutable finite history for a named source and initialization."""

    start: float
    source_shape: tuple
    identity: str
    panels: tuple[SourcePanel, ...] = ()
    incoming: IncomingQuadrupole | ParabolicQuadrupole | SampledIncomingTail | None = None

    def __post_init__(self):
        panels = tuple(self.panels)
        if not np.isfinite(self.start) or not self.identity:
            raise ValueError("A history requires a finite start and explicit identity.")
        if self.incoming is not None and (
            self.start != self.incoming.cutoff or self.source_shape != getattr(self.incoming,"source_shape",(3,3,3))
        ):
            raise ValueError("The incoming quadrupole must join the time and shape of the panel history.")
        previous = self.start
        for panel in panels:
            if panel.start != previous or panel.coefficients.shape[1:] != self.source_shape:
                raise ValueError("The source panels must cover adjacent intervals with one source shape.")
            previous = panel.end
        object.__setattr__(self, "panels", panels)

    @property
    def end(self):
        return self.panels[-1].end if self.panels else self.start

    def append(self, *panels):
        return PanelHistory(
            self.start, self.source_shape, self.identity, self.panels + panels, self.incoming,
        )

    @cached_property
    def _panel_arrays(self):
        if not self.panels or max(len(panel.coefficients) for panel in self.panels) > 3:
            return None
        starts = np.array([panel.start for panel in self.panels])
        ends = np.array([panel.end for panel in self.panels])
        coefficients = np.zeros((len(self.panels), 3) + self.source_shape)
        for index, panel in enumerate(self.panels):
            coefficients[index, :len(panel.coefficients)] = panel.coefficients
        return starts, ends, coefficients

    @cached_property
    def _panel_ends(self):
        return tuple(panel.end for panel in self.panels)

    @cached_property
    def _parts_cache(self):
        return {}

    def integral(self, time, scale, power=1, trial=None):
        """Integrate up to time, optionally including the clipped trial interval."""
        if not np.isfinite([time, scale]).all() or scale <= 0 or power not in (0, 1, 2):
            raise ValueError("History integrals require a finite time, positive scale and power 0–2.")
        if trial is not None and (
            trial.start != self.end or trial.coefficients.shape[1:] != self.source_shape
        ):
            raise ValueError("The provisional panel must join the committed source history.")
        end = trial.end if trial is not None else self.end
        if not self.start <= time <= end:
            raise ValueError("The represented history must cover the requested time.")
        initial = self._committed_integral(time, scale, power)
        if trial is not None and trial.start < time:
            initial = initial + trial.integral(time, scale, power)
        return initial

    def _committed_integral(self, time, scale, power):
        """Integrate the committed prefix up to min(time, end), including for a trial cut."""
        initial = np.zeros(self.source_shape)
        if self.incoming is not None:
            initial = self.incoming.integral(time, scale, power)
        if self._panel_arrays is not None:
            starts, ends, coefficients = self._panel_arrays
            count = np.searchsorted(starts, time, side="left")
            if count:
                weights = batched_log_weights(starts[:count], ends[:count], time, scale, power)
                initial = initial + np.tensordot(weights, coefficients[:count], axes=((0, 1), (0, 1)))
        else:
            initial = sum(
                (panel.integral(time, scale, power) for panel in self.panels if panel.start < time), initial,
            )
        return initial

    def value(self, time):
        """Evaluate the represented source on a covered finite panel."""
        if not np.isfinite(time) or not self.start <= time <= self.end or not self.panels:
            raise ValueError("A finite source panel must cover the evaluation time.")
        panel = self.panels[bisect_left(self._panel_ends, time)]
        fraction = (time - panel.start) / (panel.end - panel.start)
        return np.polynomial.polynomial.polyval(fraction, panel.coefficients)

    def integrate_jet(self, time, scale, power, orders):
        """Integrate selected derivatives of a canonical jet."""
        return self.integral(time, scale, power)[orders]

    def start_jump(self, order):
        """Initial source jump from the explicitly finite prefix convention."""
        return self.value(self.start)[order]


@dataclass(frozen=True)
class TailInitialData:
    """Finite quadrupole histories and independent inherited recoil integrals."""

    quadrupole: PanelHistory
    leading_recoil: tuple[float, float, float]
    next_recoil: tuple[float, float, float]
    center_integral: tuple[float, float, float]
    asymptotic_join_residual: tuple[float, float, float] | None = None
    orbital_prefix: "FiniteOrbitalPast | None" = None
    incoming_orbit: IncomingOrbit | ParabolicOrbit | None = None

    def recoil_state(self):
        state = np.asarray((self.leading_recoil, self.next_recoil, self.center_integral), dtype=float)
        if state.shape != (3, 3) or not np.isfinite(state).all():
            raise ValueError("Inherit three finite Cartesian vectors for recoil and CM history.")
        return state.reshape(-1)


@dataclass(frozen=True)
class FiniteOrbitalPast:
    """A caller-declared finite SH orbital prefix and its conserved inputs."""

    boundaries: tuple[float, ...]
    state: object
    mass_grades: dict
    dipole_grades: dict
    infrared_scale: float
    identity: str

    def __post_init__(self):
        cuts = tuple(self.boundaries)
        if len(cuts) < 2 or not np.isfinite(cuts).all() or not (np.diff(cuts) > 0).all():
            raise ValueError("The finite prefix must supply increasing source-time boundaries.")
        if not callable(self.state) or not self.identity:
            raise ValueError("The prefix must supply a state sampler and an explicit identity.")
        if not np.isfinite(self.infrared_scale) or self.infrared_scale <= 0:
            raise ValueError("The infrared regulator is a finite positive length.")
        object.__setattr__(self, "boundaries", cuts)

    def initial_data(self, binary, leading_recoil, next_recoil, center_integral):
        """Build dynamical tail panels from the same worldline as the wave source."""
        from .hereditary import _interpolate
        from .tails import quadrupole_tail_sources

        def sample(time):
            return quadrupole_tail_sources(binary, np.concatenate(self.state(time)))

        panels = tuple(
            SourcePanel(start, stop, _interpolate([sample(time) for time in (start, (start + stop) / 2, stop)]))
            for start, stop in zip(self.boundaries[:-1], self.boundaries[1:])
        )
        history = PanelHistory(self.boundaries[0], (3, 3, 3), self.identity, panels)
        return TailInitialData(history, leading_recoil, next_recoil, center_integral, orbital_prefix=self)


def _endpoint_moments(degree):
    harmonic = fsum(1 / index for index in range(1, degree + 2))
    harmonic_squared = fsum(1 / index**2 for index in range(1, degree + 2))
    divisor = degree + 1
    return -harmonic / divisor, (harmonic**2 + harmonic_squared) / divisor


def _primitive(distance, degree, power):
    if distance == 0:
        return 0.0
    divisor = degree + 1
    logarithm = log(distance)
    if power == 1:
        factor = logarithm - 1 / divisor
    else:
        factor = logarithm**2 - 2 * logarithm / divisor + 2 / divisor**2
    return distance**divisor / divisor * factor


def _log_moments(degree, ratio):
    """Integrate xi**degree times log(1-ratio*xi) and its square."""
    if ratio == 1:
        return _endpoint_moments(degree)
    if ratio > 0.5:
        return tuple(
            fsum(
                (-1)**index * comb(degree, index) * (
                    _primitive(1, index, power) - _primitive(1 - ratio, index, power)
                )
                for index in range(degree + 1)
            ) / ratio**(degree + 1)
            for power in (1, 2)
        )

    # For a distant panel, a convergent power series avoids subtracting
    # nearly equal primitives. Summing 64 terms is adequate for ratio <= 1/2.
    linear_terms = []
    squared_terms = []
    harmonic = 0.0
    ratio_power = ratio
    for index in range(1, 65):
        term = ratio_power / (index * (degree + index + 1))
        linear_terms.append(-term)
        squared_terms.append(2 * harmonic * term)
        harmonic += 1 / index
        ratio_power *= ratio
    return fsum(linear_terms), fsum(squared_terms)


@cache
def _far_series(degree, power):
    orders = np.arange(1, 65)
    denominator = orders * (degree + orders + 1)
    if power == 1:
        coefficients = -1 / denominator
    else:
        harmonic = np.concatenate(([0], np.cumsum(1 / orders[:-1])))
        coefficients = 2 * harmonic / denominator
    return np.concatenate(([0], coefficients))


def _batch_log_moment(degree, ratio, power):
    result = np.empty_like(ratio)
    far = ratio <= 0.5
    end = ratio == 1
    result[end] = _endpoint_moments(degree)[power - 1]
    result[far] = np.polynomial.polynomial.polyval(ratio[far], _far_series(degree, power))
    partial = ~(far | end)
    fraction = ratio[partial]
    if fraction.size:
        remaining = 1 - fraction
        logarithm = np.log(remaining)
        numerator = 0
        for index in range(degree + 1):
            divisor = index + 1
            if power == 1:
                lower = remaining**divisor / divisor * (logarithm - 1 / divisor)
                upper = -1 / divisor**2
            else:
                lower = remaining**divisor / divisor * (
                    logarithm**2 - 2 * logarithm / divisor + 2 / divisor**2
                )
                upper = 2 / divisor**3
            numerator += (-1)**index * comb(degree, index) * (upper - lower)
        result[partial] = numerator / fraction**(degree + 1)
    return result


def batched_log_weights(starts,ends,time,scale,power,degree_count=3):
    """Weights for a batch of covered, possibly clipped polynomial panels."""
    width = np.minimum(ends, time) - starts
    delay = time - starts
    fraction = width / (ends - starts)
    ratio = width / delay
    logarithm = np.log(delay) - log(scale)
    weights = []
    for degree in range(degree_count):
        constant = 1 / (degree + 1)
        if power == 0:
            weight = constant
        else:
            first_moment = _batch_log_moment(degree, ratio, 1)
            if power == 1:
                weight = logarithm * constant + first_moment
            else:
                second_moment = _batch_log_moment(degree, ratio, 2)
                weight = logarithm**2 * constant + 2 * logarithm * first_moment + second_moment
        weights.append(width * fraction**degree * weight)
    return np.stack(weights, axis=1)


def integrate_log_panel(start, end, coefficients, time, scale, power=1):
    """Integrate a panel against log((time - t)/scale)**power.

    Coefficients are ascending powers of (t - start)/(end - start).
    The leading array axis holds the polynomial degree, with optional tensor
    axes for the source. Only the part of the panel earlier than time enters.
    """
    if not np.isfinite([start, end, time, scale]).all() or end <= start or scale <= 0:
        raise ValueError("The panel must have a finite positive width and logarithmic scale.")
    if power not in (0, 1, 2):
        raise ValueError("The logarithm power must be 0, 1 or 2.")
    coefficients = np.asarray(coefficients, dtype=float)
    if coefficients.ndim == 0 or not coefficients.shape[0] or not np.isfinite(coefficients).all():
        raise ValueError("Supply a finite, nonempty coefficient array.")
    if time <= start:
        return np.zeros(coefficients.shape[1:], dtype=float)

    delay = time - start
    width = min(end, time) - start
    ratio = width / delay
    fraction = width / (end - start)
    logarithm = log(delay) - log(scale)
    weights = []
    for degree in range(coefficients.shape[0]):
        constant = 1 / (degree + 1)
        if power == 0:
            weight = constant
        else:
            linear, squared = _log_moments(degree, ratio)
            if power == 1:
                weight = logarithm * constant + linear
            else:
                weight = logarithm**2 * constant + 2 * logarithm * linear + squared
        weights.append(width * fraction**degree * weight)
    return np.tensordot(weights, coefficients, axes=1)
