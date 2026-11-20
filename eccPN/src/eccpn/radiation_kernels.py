"""Rank-wise MPM tail and tail-of-tail kernels on finite canonical histories."""

from math import comb, fsum

import numpy as np


def _integral_and_derivative(history, time, order, terms, include_endpoint):
    remote_cut = getattr(getattr(history, "past", None), "orbit", None)
    first = history.start if remote_cut is None else remote_cut.cutoff
    if not first <= time <= history.end:
        raise ValueError("Supply a covered canonical history for the tail kernel.")
    if include_endpoint and time == history.start:
        raise ValueError("A finite-start tail derivative requires a cut strictly after its start.")
    include_endpoint = include_endpoint and time > history.start
    lag = time - history.start
    result = 0
    endpoint = 0
    for coefficient, scale, power in terms:
        integrated = history.integrate_jet(time, scale, power, slice(order, order+2))
        if len(integrated) != 2:
            raise ValueError("The canonical history omits a required tail derivative.")
        result = result + coefficient * integrated
        if include_endpoint:
            endpoint += coefficient * np.log(lag / scale)**power
    if include_endpoint:
        result[1] += endpoint * history.start_jump(order)
    return result


def _ordinary_terms(rank, kind, radiative_scale):
    if rank < 2 or kind not in ("mass", "current"):
        raise ValueError("Specify the mass or current tail of a radiating rank >= 2.")
    if kind == "mass":
        offset = (
            (2 * rank**2 + 5 * rank + 4) / (rank * (rank + 1) * (rank + 2))
            + fsum(1 / order for order in range(1, rank - 1))
        )
    else:
        offset = (rank - 1) / (rank * (rank + 1)) + fsum(1 / order for order in range(1, rank))
    scale = 2 * radiative_scale
    return (1, scale, 1), (offset, scale, 0)


def _tail_of_tail_terms(binary, rank, kind, radiative_scale, infrared_scale):
    offsets = {
        ("mass", 2): (11 / 6, -107 / 105, 124627 / 44100),
        ("mass", 3): (97 / 30, -13 / 21, 13283 / 8820),
        ("current", 2): (7 / 3, -107 / 105, 13127 / 11025),
    }
    if (kind, rank) not in offsets:
        raise ValueError("This tail of tail covers the quadrupoles and mass octupole.")
    coefficients = offsets[kind, rank]
    scale = 2 * radiative_scale
    return (
        (1, scale, 2), (coefficients[0], scale, 1),
        (coefficients[1], 2 * infrared_scale / binary.c, 1), (coefficients[2], scale, 0),
    )


def _news_endpoint(history, order, terms, scale):
    """Expand each finite-start term in log((u-u_start)/scale)."""
    initial = history.start_jump(order)
    coefficients = {}
    for k, entry_scale, degree in terms:
        offset = np.log(scale / entry_scale)
        for power in range(degree + 1):
            coefficients[power] = coefficients.get(power, 0) + k * comb(degree, power) * offset**(degree - power)
    return {power: k * initial for power, k in coefficients.items()}


def ordinary_tail(binary, history, time, rank, kind, canonical_mass, radiative_scale, *, endpoint=None, split_jet=None):
    """Return the ordinary MPM tail and its news for one canonical PN coefficient."""
    terms = _ordinary_terms(rank, kind, radiative_scale)
    grade_factor = 2 * binary.G * canonical_mass / binary.c**3
    if split_jet is None:
        result = _integral_and_derivative(history,time,rank+2,terms,endpoint is None)
    else:
        result = split_jet(time,terms,rank+2,2)
        if endpoint is False:
            lag = time-history.start
            result[1] -= history.start_jump(rank+2)*sum(
                coefficient*np.log(lag/scale)**power for coefficient,scale,power in terms)
    result *= grade_factor
    if isinstance(endpoint, dict):
        endpoint.update({power: grade_factor * tensor
                         for power, tensor in _news_endpoint(history, rank + 2, terms, 2 * radiative_scale).items()})
    return result


def tail_of_tail(
    binary, history, time, rank, kind, canonical_mass, radiative_scale, infrared_scale, *, endpoint=None, split_jet=None,
):
    """Return the 3PN mass octupole or mass/current quadrupole tail of tail."""
    terms = _tail_of_tail_terms(binary, rank, kind, radiative_scale, infrared_scale)
    grade_factor = 2 * (binary.G * canonical_mass)**2 / binary.c**6
    if split_jet is None:
        result = _integral_and_derivative(history,time,rank+3,terms,endpoint is None)
    else:
        result = split_jet(time,terms,rank+3,2)
        if endpoint is False:
            lag = time-history.start
            result[1] -= history.start_jump(rank+3)*sum(
                coefficient*np.log(lag/scale)**power for coefficient,scale,power in terms)
    result *= grade_factor
    if isinstance(endpoint, dict):
        endpoint.update({power: grade_factor * tensor
                         for power, tensor in _news_endpoint(history, rank + 3, terms, 2 * radiative_scale).items()})
    return result
