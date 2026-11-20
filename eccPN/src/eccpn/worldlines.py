"""Standard harmonic body worldlines on the relative orbit and inherited recoil."""

import numpy as np

from .generated import center_of_mass
from .multipoles import _state_quantities
from .recoil import center_flux_rate, recoil_fluxes


def local_common_shift(binary, position, velocity, harmonic_scales, maximum_grade=9):
    """Return the grade-tagged common displacement and velocity of D1 and D2."""
    scales = np.asarray(harmonic_scales)
    if scales.shape != (2,) or not (np.isfinite(scales) & (scales > 0)).all():
        raise ValueError("Supply the two positive SH scale lengths.")
    if maximum_grade not in range(10):
        raise ValueError("The local worldline inventory ends at explicit order nine.")
    radius, normal, s, w, u, velocity = _state_quantities(binary, position, velocity)
    terms = center_of_mass.local_common_shift(
        u, s, w, binary.nu, np.log(radius / scales[0]), np.log(radius / scales[1]),
    )
    delta = (binary.m1 - binary.m2) / binary.mass
    result = {}
    for grade, (p, q, v_n, v_v) in zip(center_of_mass.GRADES, terms):
        if grade <= maximum_grade:
            scale = binary.nu * delta / binary.c**grade
            displacement = scale * radius[..., None] * (np.asarray(p)[..., None] * normal + np.asarray(q)[..., None] * velocity)
            rate = scale * (np.asarray(v_n)[..., None] * normal + np.asarray(v_v)[..., None] * velocity)
            result[grade] = (displacement, rate)
    return result


def body_states(trajectory, times):
    """Return both SH positions and velocities in the trajectory's inherited frame."""
    position, velocity = trajectory.state(times)
    return harmonic_body_states(
        trajectory.binary, position, velocity, *trajectory.recoil(times), trajectory.accuracy.harmonic_scales,
    )


def harmonic_body_states(binary, position, velocity, pi_leading, pi_next, gamma, harmonic_scales):
    """Recover both worldlines from saved SH phase space and inherited recoil."""
    leading_flux, next_flux, center_flux = recoil_fluxes(binary, position, velocity)
    gamma_rate = pi_leading + pi_next + center_flux
    gamma_acceleration = leading_flux + next_flux + center_flux_rate(binary, position, velocity)
    radius, normal, _, speed_squared, potential, _ = _state_quantities(binary, position, velocity)
    acceleration = -(potential / radius)[..., None] * normal
    energy = (speed_squared / 2 - potential)[..., None]
    speed_gamma = np.sum(velocity * gamma_rate, axis=-1, keepdims=True)
    derivative = np.sum(acceleration * gamma_rate + velocity * gamma_acceleration, axis=-1, keepdims=True)
    relativistic = binary.nu / (binary.mass * binary.c**2)
    common_position = -gamma / binary.mass + relativistic * (energy * gamma + speed_gamma * position)
    common_velocity = -gamma_rate / binary.mass + relativistic * (
        energy * gamma_rate + derivative * position + speed_gamma * velocity
    )
    for displacement, rate in local_common_shift(binary, position, velocity, harmonic_scales).values():
        common_position = common_position + displacement
        common_velocity = common_velocity + rate
    weights = (binary.m2 / binary.mass, -binary.m1 / binary.mass)
    return tuple((factor * position + common_position, factor * velocity + common_velocity) for factor in weights)
