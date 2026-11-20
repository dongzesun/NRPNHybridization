"""Exact product-rule assembly of the canonical mass quadrupole."""

from math import comb

import numpy as np

from .radiative_quadrupole import _spin_contract, _stf


_TERMS = (
    # QC2 Eq. (2.39). Fully dimensionally regularized source convention.
    (5, 1, "W", 2, "I2", 0, "scale"),
    (5, -1, "W", 1, "I2", 1, "scale"),
    (7, 4 / 7, "W2", 1, "I2", 3, "matrix"),
    (7, 6 / 7, "W2", 0, "I2", 4, "matrix"),
    (7, -1 / 7, "Y2", 3, "I2", 0, "matrix"),
    (7, -1, "Y2", 0, "I2", 3, "matrix"),
    (7, -2, "X", 0, "I2", 3, "scale"),
    (7, -5 / 21, "W1", 4, "I3", 0, "rank3"),
    (7, 1 / 63, "W1", 3, "I3", 1, "rank3"),
    (7, -25 / 21, "Y1", 3, "I3", 0, "rank3"),
    (7, -22 / 63, "Y1", 2, "I3", 1, "rank3"),
    (7, 5 / 63, "Y1", 1, "I3", 2, "rank3"),
    (7, 2, "W", 3, "W2", 0, "scale"),
    (7, 2, "W", 2, "W2", 1, "scale"),
    (7, -4 / 3, "W1", 0, "W1", 3, "outer"),
    (7, 2, "W", 2, "Y2", 0, "scale"),
    (7, -4, "W1", 0, "Y1", 2, "outer"),
    (7, 1 / 3, "I2", 0, "Z1", 3, "epsilon"),
    (7, -1, "I2", 3, "Z1", 0, "epsilon"),
    (7, 4 / 9, "J2", 0, "W1", 3, "epsilon"),
    (7, -4 / 9, "J2", 0, "Y1", 2, "epsilon"),
    (7, 8 / 9, "J2", 1, "Y1", 1, "epsilon"),
)

_RANKS = {
    "I2": 2, "I3": 3, "J2": 2, "W": 0, "W1": 1, "W2": 2,
    "X": 0, "Y1": 1, "Y2": 2, "Z1": 1,
}


def _product(left, right, kind):
    if kind == "scale":
        return left[..., None, None] * right
    if kind == "matrix":
        return _stf(left @ right)
    if kind == "rank3":
        return _stf(np.einsum("...a,...ija->...ij", left, right))
    if kind == "outer":
        return _stf(left[..., :, None] * right[..., None, :])
    return -_spin_contract(right, left)


def canonical_mass_quadrupole(binary, moments, derivative=0, max_grade=8):
    """Assemble M_ij^(derivative) by explicit 1/c grade through Eq. (2.39).

    Each moment maps its relative 1/c grade to an array whose first axis holds
    actual time derivatives starting at zero. Coefficients already include
    their own powers of 1/c. The mass input I2 includes all available
    renormalized source terms in the same chart.
    """
    if derivative not in range(6) or max_grade not in range(9) or (derivative >= 4 and max_grade > 5):
        raise ValueError("Canonical quadrupole covers derivatives 0–3 through grade 8, and 4–5 through grade 5.")
    required = {"I2"}
    for grade, _, left, _, right, _, _ in _TERMS:
        if grade <= max_grade:
            required.update((left, right))
    jets = {}
    for name in required:
        if name not in moments:
            raise ValueError(f"Explicit source jet {name} is required.")
        rank = _RANKS[name]
        jets[name] = {grade: np.asarray(jet) for grade, jet in moments[name].items()}
        if 0 not in jets[name] or any(
            grade not in range(9) or jet.ndim < 1 + rank or (rank and jet.shape[-rank:] != (3,) * rank)
            for grade, jet in jets[name].items()
        ):
            raise ValueError(f"Invalid grade inventory or Cartesian shape for source jet {name}.")
    result = {
        grade: _derivative("I2", jet, derivative).copy()
        for grade, jet in jets["I2"].items() if grade <= max_grade
    }
    for grade, coefficient, left, l_order, right, r_order, kind in _TERMS:
        if grade > max_grade:
            continue
        scale = 4 * binary.G / binary.c**grade * coefficient
        for left_grade, left_jet in jets[left].items():
            for right_grade, right_jet in jets[right].items():
                total = grade + left_grade + right_grade
                if total > max_grade:
                    continue
                value = sum(
                    comb(derivative, index) * _product(
                        _derivative(left, left_jet, l_order + index),
                        _derivative(right, right_jet, r_order + derivative - index), kind,
                    )
                    for index in range(derivative + 1)
                )
                result[total] = result.get(total, 0) + scale * value
    return result


def _derivative(name, jet, order):
    if order >= len(jet):
        raise ValueError(f"Source jet {name} needs time derivative {order}.")
    return jet[order]
