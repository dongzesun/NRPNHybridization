"""Analytic Newtonian Cartesian mass and current multipole jets."""

from math import factorial

import numpy as np

from .leading_modes import newtonian_position_series
from .projection import stf_tensor


def _outer_series(left, right, order):
    shape = left.shape[1:] + (3,)
    result = np.zeros((order + 1,) + shape)
    for n in range(order + 1):
        for k in range(n + 1):
            result[n] += left[k][..., None] * right[n - k]
    return result


def newtonian_source_jets(binary, position, velocity, highest_derivative=8, highest_rank=6):
    """Return fixed-history source dictionaries with grade zero for each STF tensor."""
    if highest_rank not in range(2, 10) or highest_derivative < 1:
        raise ValueError("Newtonian jets require a positive order and a mass/current rank from 2 to 9.")
    series = newtonian_position_series(binary, position, velocity, highest_derivative)
    angular = np.cross(series[0], series[1])
    x1, x2 = binary.m1 / binary.mass, binary.m2 / binary.mass
    factors = [factorial(n) for n in range(highest_derivative + 1)]
    power = series
    result = {"S1": {0: np.array([binary.mu * angular, *([np.zeros(3)] * highest_derivative)])}}
    for rank in range(2, highest_rank + 1):
        mass = binary.mu * (x2**(rank - 1) + (-1)**rank * x1**(rank - 1))
        current = binary.mu * (x2**rank - (-1)**rank * x1**rank)
        physical_factors = np.asarray(factors).reshape((-1,) + (1,) * rank)
        current_outer = power[..., None] * angular
        power = _outer_series(power, series, highest_derivative)
        result[f"M{rank}"] = {0: mass * physical_factors * stf_tensor(power, rank)}
        result[f"S{rank}"] = {0: current * physical_factors * stf_tensor(current_outer, rank)}
    return result
