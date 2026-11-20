"""Three-piece general-orbit luminosity ingredients of QC2, Section VI A."""

from math import comb, factorial, prod

import numpy as np

from .generated import local_luminosity as coefficient


def canonical_local_luminosity(binary, position, velocity, harmonic_scale, infrared_scale):
    """Return the canonical local flux by quadrupole-relative 1/c grade."""
    position = np.asarray(position, dtype=float)
    velocity = np.asarray(velocity, dtype=float)
    if position.shape != velocity.shape or position.shape[-1:] != (3,):
        raise ValueError("Matching Cartesian position and velocity arrays are required.")
    radius = np.linalg.norm(position, axis=-1)
    if (not np.isfinite([harmonic_scale, infrared_scale]).all()
            or min(harmonic_scale, infrared_scale) <= 0 or not (radius > 0).all()
            or not np.isfinite(radius).all() or not np.isfinite(velocity).all()):
        raise ValueError("The scales, state and separation must be finite with positive lengths.")
    potential = binary.G * binary.mass / radius
    radial_speed = np.sum(position * velocity, axis=-1) / radius
    speed_squared = np.sum(velocity**2, axis=-1)
    args = (potential, radial_speed, speed_squared, binary.nu)
    with_logs = (*args, np.log(radius / harmonic_scale), np.log(radius / infrared_scale))
    values = {
        0: coefficient.flux_0(*args), 2: coefficient.flux_2(*args),
        4: coefficient.flux_4(*args), 5: coefficient.flux_5(*args),
        6: coefficient.flux_6(*with_logs), 7: coefficient.flux_7(*args),
        8: coefficient.flux_8(*with_logs),
    }
    return {
        grade: binary.mass / (radius * binary.c**(5 + grade)) * raw
        for grade, raw in values.items()
    }


def _log_kernel_derivative(lag, scale, derivative):
    if derivative == 0:
        return np.log(lag / scale)
    return (-1)**(derivative - 1) * factorial(derivative - 1) / lag**derivative


def canonical_nonlocal_luminosity(binary, history, time, position, velocity, canonical_mass, infrared_scale, current=None):
    """Evaluate QC2 Eq. (6.4) as printed, with finite-start source boundary terms."""
    _, source, printed, _ = _canonical_nonlocal_parts(
        binary, history, time, position, velocity, canonical_mass, infrared_scale, current,
    )
    return source+printed


def _canonical_nonlocal_parts(
    binary, history, time, position, velocity, canonical_mass, infrared_scale, current, *,
    dynamical_tail_jet=None, split_jet=None,
):
    if (history.source_shape not in ((9, 3, 3), (10, 3, 3))
            or not history.start < time <= history.end):
        raise ValueError("The canonical nonlocal luminosity requires at least nine finite covered quadrupole jets.")
    position, velocity = np.asarray(position, dtype=float), np.asarray(velocity, dtype=float)
    radius = np.linalg.norm(position)
    if (position.shape != (3,) or velocity.shape != (3,) or not np.isfinite(velocity).all()
            or not np.isfinite([radius, canonical_mass, infrared_scale]).all()
            or min(radius, canonical_mass, infrared_scale) <= 0):
        raise ValueError("Provide a finite relative state, positive canonical mass and positive IR scale.")
    moments = history.value(time) if current is None else np.asarray(current)
    if moments.shape != history.source_shape or not np.isfinite(moments).all():
        raise ValueError("All Newtonian quadrupole derivatives must be supplied at the current cut.")
    scale = 2 * infrared_scale / binary.c
    lag = time - history.start
    tails = history.integrate_jet(time,scale,1,slice(5,9)).copy()
    for k in range(1,4):
        for j in range(k):
            tails[k] += history.start_jump(5+k-1-j)*_log_kernel_derivative(lag,scale,j)

    def source_contraction(entries):
        return sum(comb(3,k)*moments[3-k]@entries[k].T for k in range(4))/7

    printed_source = source_contraction(tails)
    source = (printed_source if split_jet is None
              else source_contraction(split_jet(time,((1,scale,1),),5,4)))
    scale = 2*radius/binary.c
    tail7, tail8 = history.integrate_jet(time, scale, 1, slice(7, 9))
    radial_fraction = np.dot(position, velocity)/radius**2
    lag_rate = (history.start_jump(7)*np.log(lag/scale)
                - radial_fraction*history.integrate_jet(time, scale, 0, 7))
    dynamics = -(binary.mu / 15) * (
        np.outer(position, position) @ tail8.T
        + (3 * np.outer(velocity, position) + np.outer(position, velocity)) @ tail7.T
    )
    extra = -(2*binary.mu/15)*np.outer(position,position) @ lag_rate.T
    pre = 48*binary.G**3*canonical_mass/(5*binary.c**13)
    retarded = 2*dynamics+extra
    if dynamical_tail_jet is not None:
        tail, rate = dynamical_tail_jet
        retarded = -(2*binary.mu/15) * (
            np.outer(position, position) @ rate.T
            + (3*np.outer(velocity, position)+np.outer(position, velocity)) @ tail.T
        )
    return tuple(pre*np.sum(moments[3]*value) for value in (source,printed_source,dynamics,retarded))


def causal_canonical_nonlocal_luminosity(binary, position, velocity, history, time, infrared_scale, current=None):
    """Contract the actual retarded source and acceleration jets with Newtonian news."""
    source, _, _, retarded = _canonical_nonlocal_parts(
        binary, history, time, position, velocity, binary.mass, infrared_scale, current,
    )
    return source+retarded


def _stf_prefactor(binary, kind, rank):
    double_factorial = prod(range(1, 2 * rank + 2, 2))
    if kind == "mass":
        coefficient = (rank + 1) * (rank + 2) / ((rank - 1) * rank * factorial(rank) * double_factorial)
        grade = 2 * rank + 1
    else:
        coefficient = 4 * rank * (rank + 2) / ((rank - 1) * factorial(rank + 1) * double_factorial)
        grade = 2 * rank + 3
    return binary.G * coefficient / binary.c**grade, grade - 5


def canonical_luminosity_from_jets(binary, jets, max_grade=8):
    """Independently contract the canonical news tensors by explicit flux grade."""
    result = {}
    for name, dictionary in jets.items():
        rank = int(name[1])
        if rank < 2:
            continue
        kind = "mass" if name.startswith("M") else "current"
        scale, offset = _stf_prefactor(binary, kind, rank)
        axes = tuple(range(-rank, 0))
        for grade_a, first in dictionary.items():
            for grade_b, second in dictionary.items():
                grade = offset + grade_a + grade_b
                if grade <= max_grade:
                    contribution = scale * np.sum(first[rank + 1] * second[rank + 1], axis=axes)
                    result[grade] = result.get(grade, 0) + contribution
    return result


def mode_luminosity(binary, news, max_grade):
    """Sky-integrated mode power, with products reexpanded before summing."""
    grouped = {}
    for (grade, ell, m), coefficient in news.items():
        grouped.setdefault((ell, m), []).append((grade, coefficient))
    result = {}
    scale = binary.c**3 / (16 * np.pi * binary.G)
    for terms in grouped.values():
        for grade_a, first in terms:
            for grade_b, second in terms:
                grade = grade_a + grade_b
                if grade <= max_grade:
                    product = scale * np.real(first * np.conj(second))
                    result[grade] = result.get(grade, 0) + product
    return result


def angular_momentum_flux(binary, radiative_moments, max_grade):
    """Native angular-momentum flux from STF strain/news by explicit grade."""
    result = {}
    for (kind, rank), (strain, news) in radiative_moments.items():
        scale, offset = _stf_prefactor(binary, kind, rank)
        for first_grade, first in strain.items():
            for second_grade, second in news.items():
                grade = offset + first_grade + second_grade
                if grade <= max_grade:
                    torque = np.cross(first.reshape(3, -1).T, second.reshape(3, -1).T).sum(axis=0)
                    result[grade] = result.get(grade, 0) + rank * scale * torque
    return result


def nonlinear_luminosity(binary, canonical_news, corrections, max_grade=8):
    """Compute the graded radiative correction to the canonical STF power.

    The dictionaries map (kind, rank) to a mapping of additional 1/c grade
    above that moment's leading term and the associated Cartesian tensor.
    Corrections hold U-I or V-J time derivatives. Their grading excludes the
    common projection factors, which are supplied by the rank and parity.
    """
    pieces = {}
    for (kind, rank), radiative in corrections.items():
        if kind not in ("mass", "current") or rank not in range(2, 8):
            raise ValueError("Expected mass or current STF tensors of rank 2–7.")
        canonical = canonical_news[kind, rank]
        scale, offset = _stf_prefactor(binary, kind, rank)
        axes = tuple(range(-rank, 0))
        for first_grade, first in radiative.items():
            for multiplier, other in ((2, canonical), (1, radiative)):
                for second_grade, second in other.items():
                    grade = offset + first_grade + second_grade
                    if grade <= max_grade:
                        contribution = multiplier * scale * np.sum(first * second, axis=axes)
                        pieces[grade] = pieces.get(grade, 0) + contribution
    return pieces


def finite_luminosity_grade_eight(binary, history, time, position, velocity, harmonic_scales, radiative_scale):
    """Compute the three 4PN scalar ingredients on the finite canonical history.

    Return both the published nonlocal term and the contraction obtained
    from the retarded source and acceleration used by the evolution.
    """
    from .waveform_history import finite_quadrupole_grade_eight, finite_radiation_available

    _, _, radiative, current = finite_radiation_available(binary, history, time, position, velocity, radiative_scale)
    parts = finite_quadrupole_grade_eight(
        binary, history, time, position, velocity, harmonic_scales, radiative_scale, leading_jet=current["M2"][0],
    )
    return _flux_grade_eight(binary, history, time, position, velocity, radiative, current, parts)


def _flux_grade_eight(
    binary, history, time, position, velocity, radiative, current, parts, *, dynamical_tail_jet=None, split_jet=None,
):
    canonical_news, corrections = {}, {}
    for (kind, rank), (_, rates) in radiative.items():
        name = ("M" if kind == "mass" else "S") + str(rank)
        canonical_news[kind, rank] = {grade: jet[rank + 1] for grade, jet in current[name].items()}
        corrections[kind, rank] = {
            grade: tensor - canonical_news[kind, rank].get(grade, 0)
            for grade, tensor in rates.items()
        }
    corrections["mass", 2][8] = sum(parts[name][1] for name in ("ordinary_tail", "tail_of_tail", "cubic"))
    reference = canonical_local_luminosity(binary, position, velocity, history.harmonic_scale, history.infrared_scale)
    prefix = history.quadrupole_for_cubic_terms
    leading = current["M2"][0]
    source, printed_source, printed, retarded = _canonical_nonlocal_parts(
        binary, prefix, time, position, velocity, history.mass_grades[0], history.infrared_scale, leading,
        dynamical_tail_jet=dynamical_tail_jet,split_jet=split_jet,
    )
    return {
        "canonical_local": reference,
        "canonical_nonlocal_published": {8: printed_source+printed},
        "canonical_nonlocal_retarded": {8: source+retarded},
        "propagation": nonlinear_luminosity(binary, canonical_news, corrections),
    }
