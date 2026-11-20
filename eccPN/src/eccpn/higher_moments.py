"""Published general-orbit higher source moments and local derivatives."""

import numpy as np

from .generated import higher_multipoles as higher, octupole_current_quadrupole as rank_two_three
from .multipoles import _state_quantities
from .projection import stf_of_vectors


def higher_source_moment(binary, name, position, velocity, derivative, harmonic_scale, infrared_scale):
    """Return the available SH coefficients of an individual source moment.

    A key is the explicit correction grade above the moment's leading
    tensor. The returned derivatives are analytically order-reduced. The
    inventory is deliberately limited to source coefficients supplied by
    Sections 6 and 8 and their pinned published references.
    """
    if derivative < 0:
        raise ValueError("The time derivative must be nonnegative.")
    jets = higher_source_jet(binary, name, position, velocity, derivative, harmonic_scale, infrared_scale)
    return {grade: value[derivative] for grade, value in jets.items()}


def higher_source_jet(binary, name, position, velocity, highest, harmonic_scale, infrared_scale):
    """Assemble all available coefficients up to the requested time derivative."""
    if not np.isfinite([harmonic_scale, infrared_scale]).all() or min(harmonic_scale, infrared_scale) <= 0:
        raise ValueError("The harmonic and infrared scales must be finite positive lengths.")
    if name == "i3":
        power, rank, axial_order, last, grades = 3, 3, 0, 4, rank_two_three.GRADES
        evaluator = rank_two_three.i3
    elif name == "j2":
        power, rank, axial_order, last, grades = 2, 1, 1, 3, rank_two_three.GRADES
        evaluator = rank_two_three.j2
    elif name in higher.INVENTORY:
        power, rank, axial_order, last, grades = higher.INVENTORY[name]
        evaluator = getattr(higher, name)
    else:
        raise ValueError(f"No higher source jet is available for {name}.")
    extended_order = getattr(higher, "TAIL_ORDERS", {}).get(name, last)
    if highest not in range(extended_order + 1):
        raise ValueError(f"{name} supports derivatives up to order {extended_order}.")
    radius, normal, radial, speed_squared, potential, velocity = _state_quantities(binary, position, velocity)
    axial = [np.cross(normal, velocity)] if axial_order else []
    basis = [
        stf_of_vectors([normal] * (rank - replacements) + [velocity] * replacements + axial)
        for replacements in range(rank + 1)
    ]
    arguments = (
        potential, radial, speed_squared, binary.nu,
        np.log(radius / harmonic_scale), np.log(radius / infrared_scale),
    )
    values = evaluator(*arguments)
    derivatives = [(grade, min(highest, last), values, len(grades), index)
                   for index, grade in enumerate(grades)]
    if highest > last:
        tail = getattr(higher, name + "_tail")(*arguments)
        derivatives[:2] = [(grade, highest, tail, 2, index) for index, grade in enumerate((0, 2))]
    delta = (binary.m1 - binary.m2) / binary.mass
    body_factor = -delta if name in ("i5", "j4") else delta if name in ("i3", "i7", "j2", "j6") else 1
    axes = (...,) + (None,) * (rank + axial_order)
    scales = [(binary.mu * body_factor * radius**(power - derivative))[axes] for derivative in range(highest + 1)]
    return {
        grade: np.stack([
            scale / binary.c**grade * sum(
                np.asarray(coefficient)[axes] * tensor
                for coefficient, tensor in zip(values[derivative * count + index], basis)
            )
            for derivative, scale in enumerate(scales[:limit + 1])
        ])
        for grade, limit, values, count, index in derivatives
    }
