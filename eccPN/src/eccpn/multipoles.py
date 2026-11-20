"""Fixed-history source moments on the local standard harmonic trajectory."""

import numpy as np
from math import comb, factorial

from .generated import source_mass_quadrupole as generated


def _stf(left, right):
    product = left[..., :, None] * right[..., None, :]
    symmetric = (product + np.swapaxes(product, -1, -2)) / 2
    return symmetric - np.trace(symmetric, axis1=-1, axis2=-2)[..., None, None] * np.eye(3) / 3


def _state_quantities(binary, position, velocity):
    position = np.asarray(position, dtype=float)
    velocity = np.asarray(velocity, dtype=float)
    if position.shape != velocity.shape or position.shape[-1:] != (3,):
        raise ValueError("Matching Cartesian position and velocity arrays are required.")
    radius = np.linalg.norm(position, axis=-1)
    if (not np.isfinite(radius).all() or not (radius > 0).all() or not np.isfinite(velocity).all()):
        raise ValueError("Position and velocity must be finite and separation must be positive.")
    normal = position / radius[..., None]
    radial_speed = np.sum(normal * velocity, axis=-1)
    speed_squared = np.sum(velocity**2, axis=-1)
    return radius, normal, radial_speed, speed_squared, binary.G * binary.mass / radius, velocity


def source_mass_quadrupole(binary, position, velocity, derivative, harmonic_scale, infrared_scale):
    """Return the source STF tensor by fixed-history grade through 3.5PN.

    A dictionary key is the explicit power of 1/c. Values include that power
    and the dimensional scale. This source moment still requires the canonical
    and radiative transformations to enter strain.
    """
    evaluators = (
        generated.derivative_0, generated.derivative_1,
        generated.derivative_2, generated.derivative_3,
        generated.derivative_4, generated.derivative_5,
    )
    if derivative not in range(len(evaluators)):
        raise ValueError("Only derivatives from zero to five are available.")
    scales = np.asarray([harmonic_scale, infrared_scale], dtype=float)
    if not (np.isfinite(scales) & (scales > 0)).all():
        raise ValueError("The source and infrared length scales must be positive and finite.")
    radius, normal, radial, speed_squared, potential, velocity = _state_quantities(binary, position, velocity)
    coefficients = evaluators[derivative](
        potential, radial, speed_squared, binary.nu,
        np.log(radius / harmonic_scale), np.log(radius / infrared_scale),
    )
    basis = (_stf(normal, normal), _stf(normal, velocity), _stf(velocity, velocity))
    factor = binary.mu * radius[..., None, None]**(2 - derivative)
    grades = generated.GRADES if derivative <= 3 else (0, 2, 4, 5)
    return {
        grade: factor / binary.c**grade * sum(
            np.asarray(value)[..., None, None] * tensor for value, tensor in zip(values, basis)
        )
        for grade, values in zip(grades, coefficients)
    }


def quadrupole_infrared_commutator(binary, position, velocity):
    """Return the standalone 4PN infrared commutator of QC2, Sec. III C."""
    radius, normal, radial, speed_squared, potential, _ = _state_quantities(binary, position, velocity)
    coefficient = (
        -4 * binary.mu * binary.nu * radius**2 * potential**3 / binary.c**8
        * (speed_squared - 3 * radial**2 - potential)
    )
    return coefficient[..., None, None] * _stf(normal, normal)


def quadrupole_renormalization(binary, position, velocity):
    """Finite local difference relative to the Hadamard source at 4PN."""
    _, normal, radial, speed_squared, potential, velocity = _state_quantities(binary, position, velocity)
    nu = binary.nu
    radial_coefficient = (
        potential**2 * (1072 / 105 - 571826 * nu / 14175)
        + potential * radial**2 * (3004 / 147 + 1800322 * nu / 11025)
        + potential * speed_squared * (-30004 / 2205 - 9985783 * nu / 99225)
    )
    mixed_coefficient = potential * radial * (-7516 / 245 - 2437018 * nu / 19845)
    velocity_coefficient = (
        speed_squared * (7516 / 735 - 123944 * nu / 6615)
        + potential * (-15056 / 2205 + 161096 * nu / 6615)
    )
    factor = binary.G**2 * binary.mass**3 * nu / binary.c**8
    return factor * (
        radial_coefficient[..., None, None] * _stf(normal, normal)
        + mixed_coefficient[..., None, None] * _stf(normal, velocity)
        + velocity_coefficient[..., None, None] * _stf(velocity, velocity)
    )


def quadrupole_nonlocal_source(trajectory, time, infrared_scale):
    """Retarded 4PN nonlocal source moment before its radiative derivatives."""
    if trajectory.history is None:
        raise ValueError("The source tail requires a covered quadrupole history.")
    if np.ndim(time) != 0:
        raise ValueError("Request the source tail on one cut.")
    binary = trajectory.binary
    position, _ = trajectory.state(time)
    quadrupole = binary.mu * _stf(position, position)
    lagged = trajectory.history.integral(time, 2 * infrared_scale / binary.c)[1]
    product = quadrupole @ lagged
    tensor = (product + product.T) / 2 - np.eye(3) * np.trace(product) / 3
    return (24 * binary.G**2 * binary.mass / (7 * binary.c**8)) * tensor


def known_local_quadrupole_grade_eight(binary, position, velocity, harmonic_scales, infrared_scale):
    """Known grade-eight local source jets, including lower-source differentiation."""
    from .generated import partial_quadrupole_4pn as local

    scales = np.asarray((*harmonic_scales, infrared_scale))
    if scales.shape != (3,) or not (np.isfinite(scales) & (scales > 0)).all():
        raise ValueError("Supply two positive SH scales and one positive infrared scale.")
    radius, normal, radial, speed_squared, potential, velocity = _state_quantities(binary, position, velocity)
    arguments = (potential, radial, speed_squared, binary.nu,
                 np.log(radius / scales[0]), np.log(radius / scales[2]), np.log(radius / scales[1]))
    basis = (_stf(normal, normal), _stf(normal, velocity), _stf(velocity, velocity))
    result = []
    for order, evaluate in enumerate((local.derivative_0, local.derivative_1, local.derivative_2, local.derivative_3)):
        factor = binary.mu * radius[..., None, None]**(2 - order) / binary.c**8
        result.append(factor * sum(
            np.asarray(value)[..., None, None] * tensor for value, tensor in zip(evaluate(*arguments)[0], basis)
        ))
    return np.stack(result)


def additional_local_quadrupole_grade_eight(binary, position, velocity, harmonic_scale):
    """The remaining local Hadamard term after the explicit 4PN corrections."""
    from .generated import local_quadrupole_4pn as source

    radius, normal, s, w, u, velocity = _state_quantities(binary, position, velocity)
    if not np.isfinite(harmonic_scale) or harmonic_scale <= 0:
        raise ValueError("The common SH ultraviolet scale must be positive.")
    basis = (_stf(normal, normal), _stf(normal, velocity), _stf(velocity, velocity))
    args = (u, s, w, binary.nu, np.log(radius/harmonic_scale), 0.)
    return np.stack([
        (binary.mu*radius**(2-order)/binary.c**8)[...,None,None]
        * sum(np.asarray(coefficient)[...,None,None]*tensor for coefficient,tensor in zip(evaluate(*args)[0],basis))
        for order,evaluate in enumerate((source.derivative_0, source.derivative_1,
                                         source.derivative_2, source.derivative_3))
    ])


def _constant_scale_log_jet(history, time, scale, initial_order, last_derivative):
    """Derivatives of a finite-prefix log integral at fixed scale."""
    result = history.integrate_jet(time, scale, 1, slice(initial_order, initial_order+last_derivative+1)).copy()
    lag = time - history.start
    kernel = [np.log(lag / scale)] + [
        (-1)**(order - 1) * factorial(order - 1) / lag**order for order in range(1, last_derivative)
    ]
    for order in range(1, last_derivative + 1):
        result[order] += sum(history.start_jump(initial_order+order-1-j)*kernel[j] for j in range(order))
    return result


def nonlocal_quadrupole_grade_eight(
    binary, newtonian_history, time, position, velocity, infrared_scale, current=None, *,
    dynamical_tail_jet=None, split_jet=None,
):
    """Retarded source and order reduction of its local Newtonian term, through news."""
    if (newtonian_history.source_shape != (10, 3, 3)
            or not newtonian_history.start < time <= newtonian_history.end):
        raise ValueError("Supply a finite ten-derivative Newtonian quadrupole history covering the cut.")
    radius, _, radial, _, _, velocity = _state_quantities(binary, position, velocity)
    if radius.ndim or not np.isfinite(infrared_scale) or infrared_scale <= 0:
        raise ValueError("Supply one Cartesian relative state and a positive infrared scale.")
    current = newtonian_history.value(time) if current is None else np.asarray(current, dtype=float)
    if current.shape != (10, 3, 3) or not np.isfinite(current).all():
        raise ValueError("Supply ten Cartesian Newtonian quadrupole derivatives at this cut.")
    prefactor = binary.G**2 * binary.mass / binary.c**8
    source_scale = 2*infrared_scale/binary.c
    if split_jet is None:
        lagged = _constant_scale_log_jet(newtonian_history,time,source_scale,5,3)
    else:
        lagged = split_jet(time,((1,source_scale,1),),5,4)
    # Q4b/QC2 nonlocal source moment. The product uses the instantaneous jet.
    result = np.array([
        24 / 7 * prefactor * sum(
            comb(order, j) * sum(_stf(a, b) for a, b in zip(current[order - j], lagged[j].T))
            for j in range(order + 1)
        )
        for order in range(4)
    ])
    # The scale in the dynamical tail changes with the separation.
    if dynamical_tail_jet is None:
        lag, lag_rate = _constant_scale_log_jet(newtonian_history, time, 2 * radius / binary.c, 7, 1)
        lag_rate -= radial / radius * newtonian_history.integrate_jet(time, 1.0, 0, 7)
    else:
        lag, lag_rate = dynamical_tail_jet
    acceleration = -8 / 5 * prefactor * (lag @ position)
    jerk = -8 / 5 * prefactor * (lag_rate @ position + lag @ velocity)
    result[2] += 2 * binary.mu * _stf(position, acceleration)
    result[3] += 2 * binary.mu * (_stf(position, jerk) + 3 * _stf(velocity, acceleration))
    return result
