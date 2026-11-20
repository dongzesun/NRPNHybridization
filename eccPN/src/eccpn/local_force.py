"""Higher-order local force contributions in standard harmonic coordinates."""

import numpy as np

from .generated import sh_local as coefficients


def higher_local_acceleration(binary, position, velocity, harmonic_scales):
    """Return separate 1–4PN conservative and 3.5/4.5PN matter accelerations."""
    radius = np.linalg.norm(position, axis=-1, keepdims=True)
    normal = position / radius
    potential = binary.G * binary.mass / radius
    radial_speed = np.sum(normal * velocity, axis=-1, keepdims=True)
    speed_squared = np.sum(velocity**2, axis=-1, keepdims=True)
    nu, c = binary.nu, binary.c
    log_one, log_two = (np.log(radius / scale) for scale in harmonic_scales)
    arguments = (potential, radial_speed, speed_squared, nu)
    with_logs = (*arguments, log_one, log_two)
    radial4 = (
        coefficients.a4_0(radial_speed, speed_squared, nu)
        + coefficients.a4_1(*arguments)
        + coefficients.a4_2(*with_logs)
        + coefficients.a4_3(*with_logs)
        + coefficients.a4_4(potential, nu, log_one, log_two)
    )
    velocity4 = (
        coefficients.b4_0(radial_speed, speed_squared, nu)
        + coefficients.b4_1(*arguments)
        + coefficients.b4_2(*with_logs)
        + coefficients.b4_3(potential, radial_speed, nu, log_one, log_two)
    )
    radial = (
        coefficients.a1(*arguments) / c**2 + coefficients.a2(*arguments) / c**4
        + coefficients.a3(*arguments, log_one) / c**6 + radial4 / c**8
    )
    tangential = (
        coefficients.b1(radial_speed, nu) / c**2 + coefficients.b2(*arguments) / c**4
        + coefficients.b3(*arguments, log_one) / c**6 + velocity4 / c**8
    )
    conservative = -potential / radius * (radial * normal + tangential * velocity)
    matter = nu * potential**2 / radius * (
        (coefficients.rr35n(*arguments) * normal + coefficients.rr35v(*arguments) * velocity) / c**7
        + (coefficients.rr45mattern(*arguments) * normal
           + coefficients.rr45matterv(*arguments) * velocity) / c**9
    )
    return conservative, matter
