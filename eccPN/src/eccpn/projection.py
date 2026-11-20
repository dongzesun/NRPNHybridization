"""Analytic STF-to-spherical projections from solid harmonics."""

from functools import cache
from itertools import product
from math import comb, factorial, pi, prod, sqrt

import numpy as np

from .leading_modes import _solid_coefficients


@cache
def _projector(ell, m):
    if m < 0:
        return (-1)**m * np.conj(_projector(ell, -m))
    polynomial = {}
    for k, coefficient in enumerate(_solid_coefficients(ell, m)):
        for y in range(m + 1):
            transverse = comb(m, y) * (-1j)**y
            for kx in range(k + 1):
                for ky in range(k - kx + 1):
                    kz = k - kx - ky
                    powers = (m - y + 2 * kx, y + 2 * ky, ell - m - 2 * k + 2 * kz)
                    factor = comb(k, kx) * comb(k - kx, ky)
                    polynomial[powers] = polynomial.get(powers, 0) + coefficient * transverse * factor
    tensor = np.zeros((3,) * ell, dtype=complex)
    for indices in product(range(3), repeat=ell):
        powers = tuple(indices.count(axis) for axis in range(3))
        permutation_count = factorial(ell) / prod(factorial(power) for power in powers)
        tensor[indices] = polynomial.get(powers, 0) / permutation_count
    tensor.flags.writeable = False
    return tensor


def project_stf(binary, rank, tensor, kind="mass"):
    """Return asymptotic strain contributions of a single radiative tensor."""
    normalization = 16 * pi * binary.G / (sqrt(2) * binary.c**(rank + 2) * prod(range(1, 2 * rank + 2, 2)))
    if kind == "mass":
        scale = normalization * sqrt((rank + 1) * (rank + 2) / (2 * rank * (rank - 1)))
    elif kind == "current":
        scale = 2j * normalization / binary.c * sqrt(rank * (rank + 2) / (2 * (rank + 1) * (rank - 1)))
    else:
        raise ValueError("Multipole type must be mass or current.")
    return {
        (rank, m): scale * np.tensordot(tensor, _projector(rank, m), axes=rank)
        for m in range(-rank, rank + 1)
    }


def stf_of_vectors(vectors):
    """Exactly project an outer product into its Cartesian STF subspace."""
    rank = len(vectors)
    axes = "ijklmnopq"[:rank]
    terms = ",".join("..." + axis for axis in axes)
    raw = np.einsum(terms + "->..." + axes, *vectors)
    return stf_tensor(raw, rank)


def stf_tensor(raw, rank):
    """Project the last rank Cartesian axes of a real tensor to STF form."""
    normalization = 4 * pi * factorial(rank) / prod(range(1, 2 * rank + 2, 2))
    result = np.zeros_like(raw, dtype=float)
    for m in range(rank + 1):
        harmonic = _projector(rank, m)
        contraction = np.tensordot(raw, harmonic, axes=rank)
        field = np.multiply.outer(contraction, harmonic.conj()).real
        result += normalization * (2 if m else 1) * field
    return result
