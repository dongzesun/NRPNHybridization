"""PN-graded canonical, instantaneous and formal memory STF maps."""

from functools import cache, lru_cache
from math import comb

import numpy as np

from .generated import mpm_tables
from .projection import stf_tensor


_LEVI_CIVITA = np.zeros((3, 3, 3))
_LEVI_CIVITA[0, 1, 2] = _LEVI_CIVITA[1, 2, 0] = _LEVI_CIVITA[2, 0, 1] = 1
_LEVI_CIVITA[0, 2, 1] = _LEVI_CIVITA[2, 1, 0] = _LEVI_CIVITA[1, 0, 2] = -1
_LEVI_CIVITA.flags.writeable = False


@cache
def _contraction_equation(tensors, rank):
    return ",".join("..." + indices for _, _, indices in tensors) + "->..." + "ijklmn"[:rank]


@lru_cache(maxsize=256)
def _contraction_path(equation, shapes):
    # NumPy's planner uses operand dimensions. Reuse its chosen order so
    # identical contractions keep the same floating-point reduction.
    operands = [np.broadcast_to(0., shape) for shape in shapes]
    return np.einsum_path(equation, *operands, optimize=True)[0]


def _input(name, jet, index):
    if index >= len(jet):
        raise ValueError(f"{name} is missing its time derivative {index}.")
    return jet[index]


def _apply(binary, target, rows, moments, derivative_order, max_grade):
    rank = int(target[1])
    result = {}
    for grade, numerator, denominator, tensors in rows:
        if grade > max_grade:
            continue
        first, second = [(symbol.rstrip("0"), order, indices)
                         for symbol, order, indices in tensors if symbol != "epsilon"]
        if first[0] not in moments or second[0] not in moments:
            raise ValueError(f"{target} requires moment jets {first[0]} and {second[0]}.")
        prefactor = binary.G * numerator / (binary.c**grade * denominator)
        equation = _contraction_equation(tensors, rank)
        for grade_a, a in moments[first[0]].items():
            for grade_b, b in moments[second[0]].items():
                total = grade + grade_a + grade_b
                if total > max_grade:
                    continue
                row = 0
                for index in range(derivative_order + 1):
                    factors = iter((
                        _input(first[0], a, first[1] + index),
                        _input(second[0], b, second[1] + derivative_order - index),
                    ))
                    values = [_LEVI_CIVITA if symbol == "epsilon" else next(factors)
                              for symbol, _, _ in tensors]
                    path = _contraction_path(equation, tuple(value.shape for value in values))
                    row += comb(derivative_order, index) * np.einsum(equation, *values, optimize=path)
                result[total] = result.get(total, 0) + prefactor * row
    return {grade: stf_tensor(value, rank) for grade, value in result.items()}


def higher_canonical_corrections(binary, target, source_jets, derivative=0, max_grade=7):
    """Nonlinear map for M3, M4, S2 and S3 from their SH source moments."""
    target = target.lower()
    return _apply(binary, target, mpm_tables.CANONICAL[target], source_jets, derivative, max_grade)


def instantaneous_corrections(binary, target, canonical_jets, derivative=0, maximum_absolute_grade=7):
    """Return corrections to the leading radiative derivative by relative 1/c grade."""
    target = target.lower()
    offset = int(target[1]) - 2 + (target[0] == "v")
    rows = mpm_tables.INSTANTANEOUS.get(target, ())
    return _apply(binary, target, rows, canonical_jets, derivative, maximum_absolute_grade - offset)


def formal_memory_rates(binary, target, canonical_jets, maximum_absolute_grade=7):
    """Return the formal MPM integrand in native coordinates, before history integration."""
    target = target.lower()
    offset = int(target[1]) - 2
    return _apply(
        binary, target, mpm_tables.MEMORY.get(target, ()), canonical_jets, 0,
        maximum_absolute_grade - offset,
    )
