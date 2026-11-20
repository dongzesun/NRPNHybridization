"""Standard harmonic source and gauge derivatives for the MPM quadrupole."""

import numpy as np

from .generated import gauge_moments as generated
from .multipoles import _state_quantities, _stf, source_mass_quadrupole


def _rank3(first, second, third):
    polynomial = np.einsum("...i,...j,...k->...ijk", first, second, third)
    dimensions = tuple(range(polynomial.ndim - 3))
    axes = tuple(range(polynomial.ndim - 3, polynomial.ndim))
    symmetric = sum(
        polynomial.transpose(dimensions + tuple(axes[i] for i in permutation))
        for permutation in ((0, 1, 2), (0, 2, 1), (1, 0, 2), (1, 2, 0), (2, 0, 1), (2, 1, 0))
    ) / 6
    trace = np.trace(symmetric, axis1=-2, axis2=-1)
    tail = trace[..., :, None, None] * np.eye(3)
    return symmetric - (tail + tail.swapaxes(-2, -3) + tail.swapaxes(-1, -3)) / 5


def canonical_source_jets(binary, position, velocity, harmonic_scale, infrared_scale):
    """Source/gauge inventory for the canonical mass quadrupole through 3.5PN."""
    if not np.isfinite([harmonic_scale, infrared_scale]).all() or min(harmonic_scale, infrared_scale) <= 0:
        raise ValueError("The harmonic and infrared scales must be finite positive lengths.")
    radius, normal, radial, speed_squared, potential, velocity = _state_quantities(binary, position, velocity)
    axial = np.cross(normal, velocity)
    delta = (binary.m1 - binary.m2) / binary.mass
    scalar = (np.ones_like(radius),)
    basis = {
        (0, 0): scalar,
        (0, 1): (axial,),
        (1, 0): (normal, velocity),
        (1, 1): (_stf(normal, axial), _stf(velocity, axial)),
        (2, 0): (_stf(normal, normal), _stf(normal, velocity), _stf(velocity, velocity)),
        (3, 0): tuple(_rank3(*([normal] * (3 - n_v) + [velocity] * n_v)) for n_v in range(4)),
    }
    inputs = (
        potential, radial, speed_squared, binary.nu,
        np.log(radius / harmonic_scale), np.log(radius / infrared_scale),
    )
    result = {}
    for name, (power, rank, axial_order, last, max_grade) in generated.INVENTORY.items():
        tensors = basis[rank, axial_order]
        entries = getattr(generated, name)(*inputs)
        grades = range(0, max_grade + 1, 2)
        count = len(grades)
        slots = {grade: [] for grade in grades}
        extra_axes = (None,) * (rank + axial_order)
        mass = binary.mu * (delta if name in ("w1", "y1", "i3", "j2") else 1)
        for order in range(last + 1):
            factor = (mass * radius**(power - order))[(...,) + extra_axes]
            for index, grade in enumerate(grades):
                coefficients = entries[order * count + index]
                field = sum(
                    np.asarray(coefficient)[(...,) + extra_axes] * tensor
                    for coefficient, tensor in zip(coefficients, tensors)
                )
                slots[grade].append(factor * field / binary.c**grade)
        result[name.upper()] = {grade: np.stack(field) for grade, field in slots.items()}

    high_grades = (4, 5, 6, 7)
    jets = [
        source_mass_quadrupole(binary, position, velocity, order, harmonic_scale, infrared_scale)
        for order in range(6)
    ]
    for grade in high_grades:
        result["I2"][grade] = np.stack([jet[grade] for jet in jets if grade in jet])
    return result
