"""All leading mass and current modes on an instantaneous Newtonian jet."""

from functools import cache
from math import factorial, pi, prod, sqrt

import numpy as np


def _times(left, right):
    return np.convolve(left, right)[:len(left)]


def _power(base, exponent):
    result = np.zeros_like(base)
    result[0] = 1
    for _ in range(exponent):
        result = _times(result, base)
    return result


def newtonian_position_series(binary, position, velocity, order):
    """Taylor coefficients of x(t + dt) based on direct Cartesian Newtonian dynamics."""
    if order < 1:
        raise ValueError("The Newtonian jet includes at least the position and velocity.")
    position, velocity = np.asarray(position, dtype=float), np.asarray(velocity, dtype=float)
    radius = np.linalg.norm(position)
    if (position.shape != (3,) or velocity.shape != (3,) or not np.isfinite(velocity).all()
            or not np.isfinite(radius) or radius <= 0):
        raise ValueError("Finite Cartesian position and velocity with positive separation are required.")
    coefficients = np.zeros((order + 1, 3))
    coefficients[0], coefficients[1] = position, velocity
    radii = np.zeros(order + 1)
    inverse_radius_cubed = np.zeros(order + 1)
    radii[0], inverse_radius_cubed[0] = radius**2, radius**-3
    for n in range(order - 1):
        if n:
            radii[n] = sum(coefficients[k] @ coefficients[n - k] for k in range(n + 1))
            inverse_radius_cubed[n] = sum(
                (-k / 2 - n) * radii[k] * inverse_radius_cubed[n - k]
                for k in range(1, n + 1)
            ) / (n * radii[0])
        acceleration = -binary.G * binary.mass * sum(
            inverse_radius_cubed[k] * coefficients[n - k] for k in range(n + 1)
        )
        coefficients[n + 2] = acceleration / ((n + 1) * (n + 2))
    return coefficients


@cache
def _solid_coefficients(ell, m):
    angular = (-1)**m * sqrt(
        (2 * ell + 1) / (4 * pi) * factorial(ell - m) / factorial(ell + m)
    )
    return tuple(
        angular * (-1)**k * factorial(2 * ell - 2 * k)
        / (2**ell * factorial(k) * factorial(ell - k) * factorial(ell - 2 * k - m))
        for k in range((ell - m) // 2 + 1)
    )


def _solid_harmonic(ell, m, position, direction):
    """Conjugated solid harmonic and its Cartesian derivative along direction."""
    _, _, z = position.T
    horizontal = position[:, 0] - 1j * position[:, 1]
    radii = sum(_times(axis, axis) for axis in position.T)
    polynomial = np.zeros_like(z)
    dz = np.zeros_like(z)
    # The direction is the conserved Newtonian orbital vector x cross v,
    # perpendicular to the position at every time of this jet.
    for k, amplitude in enumerate(_solid_coefficients(ell, m)):
        power = ell - m - 2 * k
        radius_factor = _power(radii, k)
        polynomial += amplitude * _times(_power(z, power), radius_factor)
        if power:
            dz += amplitude * power * _times(_power(z, power - 1), radius_factor)
    transverse = _power(horizontal, m)
    solid = _times(transverse, polynomial)
    derivative = direction[2] * _times(transverse, dz)
    if m:
        derivative += m * (direction[0] - 1j * direction[1]) * _times(_power(horizontal, m - 1), polynomial)
    return solid, derivative


def leading_mode_coefficients(binary, position, velocity, ell_max=9):
    """Return leading H and news contributions keyed by (grade, ell, m).

    The grade is the absolute 1/c power relative to the Newtonian mass
    quadrupole. Individual leading terms do not provide completeness at
    all lower grades for the hierarchy.
    """
    if ell_max not in range(2, 10):
        raise ValueError("The leading inventory covers degrees 2–9.")
    x = newtonian_position_series(binary, position, velocity, ell_max + 1)
    orbital = np.cross(x[0], x[1])
    X1, X2 = binary.m1 / binary.mass, binary.m2 / binary.mass
    result = ({}, {})
    for ell in range(2, ell_max + 1):
        mass = binary.mu * (X2**(ell - 1) + (-1)**ell * X1**(ell - 1))
        current = binary.mu * (X2**ell - (-1)**ell * X1**ell)
        projection = 16 * pi / prod(range(1, 2 * ell + 2, 2))
        common = binary.G / (sqrt(2) * binary.c**(ell + 2))
        mass_scale = common * mass * projection * sqrt((ell + 1) * (ell + 2) / (2 * ell * (ell - 1)))
        current_scale = (
            -2 * common * current * projection / (ell * binary.c)
            * sqrt(ell * (ell + 2) / (2 * (ell + 1) * (ell - 1)))
        )
        for m in range(ell + 1):
            solid, directed = _solid_harmonic(ell, m, x, orbital)
            for index, output in enumerate(result):
                derivative = ell + index
                factor = factorial(derivative)
                mass_mode, current_mode = mass_scale * factor * solid[derivative], current_scale * factor * directed[derivative]
                output[ell - 2, ell, m] = mass_mode
                if ell < 9:
                    output[ell - 1, ell, m] = -1j * current_mode
                if m:
                    output[ell - 2, ell, -m] = (-1)**m * np.conj(mass_mode)
                    if ell < 9:
                        output[ell - 1, ell, -m] = -1j * (-1)**m * np.conj(current_mode)
    return result
