"""Assemble the canonical waveform moment inventory from SH phase space."""

import numpy as np

from .canonical import canonical_mass_quadrupole
from .gauge_moments import canonical_source_jets
from .higher_moments import higher_source_jet
from .mpm import higher_canonical_corrections
from .newtonian_moments import newtonian_source_jets


def _add_grade(moment, grade, correction, order):
    if grade not in moment:
        moment[grade] = np.zeros((len(moment[0]),) + correction.shape)
    moment[grade][order] += correction


def canonical_waveform_jets(
    binary, position, velocity, harmonic_scale, infrared_scale, dipole_grades, max_absolute_grade=6,
):
    """Return canonical Cartesian jets through absolute 3PN for ranks 2–6.

    The caller supplies the conserved orbital current dipole in its incoming
    frame. A key in dipole_grades is the correction grade of that charge. The
    optional grade-seven inventory carries the available 3.5PN terms.
    """
    if max_absolute_grade not in (6, 7) or 0 not in dipole_grades or not dipole_grades.keys() <= {0, 2}:
        raise ValueError("Choose the supported grade 6 or 7 inventory and supply the conserved current dipole.")
    extended = max_absolute_grade == 7
    leading = newtonian_source_jets(binary, position, velocity, highest_derivative=9, highest_rank=6 + extended)
    # The instantaneous source-to-canonical equations take the locally
    # reduced Newtonian dipole, as does the published SH canonical flux.
    # The independently supplied conserved dipole belongs in propagation
    # and may differ from this fixed-grade representative at higher order.
    local_current_dipole = dict(leading["S1"])
    for grade, dipole in dipole_grades.items():
        dipole = np.asarray(dipole, dtype=float)
        if dipole.shape != (3,) or not np.isfinite(dipole).all():
            raise ValueError("Each supplied current dipole coefficient must be a finite Cartesian vector.")
        array = np.zeros_like(leading["S1"][0])
        array[0] = dipole
        leading["S1"][grade] = array
    source = canonical_source_jets(binary, position, velocity, harmonic_scale, infrared_scale)
    source["J1"] = local_current_dipole
    orders = [("i3", 4), ("j2", 3), ("i4", 5 + 2 * extended), ("i5", 6), ("i6", 7),
              ("j3", 4 + 2 * extended), ("j4", 5), ("j5", 6), ("j6", 7)]
    if extended:
        orders.append(("i7", 8))
    for name, highest in orders:
        higher = higher_source_jet(binary, name, position, velocity, highest, harmonic_scale, infrared_scale)
        available = source.setdefault(name.upper(), {})
        for grade, jet in higher.items():
            if grade and (grade not in available or len(jet) > len(available[grade])):
                available[grade] = jet
    for name, values in tuple(source.items()):
        if name[0] not in "IJ" or name == "J1":
            continue
        rank = int(name[1])
        target = ("M" if name[0] == "I" else "S") + str(rank)
        for grade, jet in values.items():
            if grade:
                leading[target][grade] = jet.copy()
        source[name][0] = leading[target][0]

    # Eq. (2.39) includes the full source I2. Higher canonical maps return
    # only their nonlinear terms, so each source coefficient enters once.
    for order in range(6 if extended else 4):
        budget = max_absolute_grade if order < 4 else 5
        values = canonical_mass_quadrupole(binary, source, order, max_grade=budget)
        for grade, tensor in values.items():
            if grade:
                leading["M2"][grade][order] = tensor

    for target in ("m3", "m4", "s2", "s3"):
        rank = int(target[1])
        correction_grade = max_absolute_grade - (rank - 2 + (target[0] == "s"))
        for order in range(rank + 2):
            corrections = higher_canonical_corrections(binary, target, source, order, correction_grade)
            for grade, tensor in corrections.items():
                _add_grade(leading[target.upper()], grade, tensor, order)
    return leading
