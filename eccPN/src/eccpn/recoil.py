"""Instantaneous momentum fluxes and SH recoil forces of [D2,D3]."""

import numpy as np


def _orbital_quantities(binary, position, velocity):
    radius = np.linalg.norm(position, axis=-1, keepdims=True)
    normal = position / radius
    radial = np.sum(normal * velocity, axis=-1, keepdims=True)
    speed_squared = np.sum(velocity**2, axis=-1, keepdims=True)
    potential = binary.G * binary.mass / radius
    delta = (binary.m1 - binary.m2) / binary.mass
    return radius, normal, radial, speed_squared, potential, delta


def recoil_fluxes(binary, position, velocity):
    """Return the leading and 1PN momentum fluxes and the leading CM flux."""
    radius, normal, s, w, u, delta = _orbital_quantities(binary, position, velocity)
    nu = binary.nu
    tangential = velocity - s * normal
    transverse_squared = np.sum(tangential**2, axis=-1, keepdims=True)
    q = transverse_squared
    c = binary.c
    prefactor = binary.mass * delta * nu**2 / c**7
    momentum_scale = prefactor * u**3 / radius
    # Bini et al., arXiv:2611.12093v1, Eqs. (2.10a,b).
    leading = momentum_scale * (
        s * normal * (32 * u / 35 - 24 * s**2 / 7 + 88 * w / 21)
        + velocity * (-64 * u / 105 + 304 * s**2 / 105 - 80 * w / 21)
    )
    center = prefactor * u**2 * (
        normal * (272 * u * s**2 / 105 - 64 * u * w / 15 + 48 * s**2 * w / 35 - 32 * w**2 / 35)
        + s * velocity * (16 * u / 21 - 24 * s**2 / 7 + 24 * w / 7)
    )
    # Eq. (5.17), projected onto n and the transverse velocity. Keeping
    # the transverse vector avoids a division at radial or near-radial motion.
    radial_coefficient = s * (
        2 * (74 * nu - 129) * u**2
        + u * (4 * (53 * nu - 197) * s**2 + 6 * (235 * nu - 832) * q)
        + 12 * (13 - 4 * nu) * s**4 + 6 * (70 * nu - 97) * s**2 * q
        + 9 * (259 - 174 * nu) * q**2
    )
    transverse_coefficient = (
        -8 * (189 + 17 * nu) * u**2
        + u * (6 * (248 * nu - 11) * s**2 + 9 * (162 * nu - 907) * q)
        + 12 * (35 * nu - 83) * s**4 + 18 * (83 * nu - 345) * s**2 * q
        + 111 * (25 - 28 * nu) * q**2
    )
    next_order = (8 * momentum_scale / (1890 * c**2)) * (
        radial_coefficient * normal - transverse_coefficient * tangential
    )
    return leading, next_order, center


def recoil_accelerations(binary, position, velocity, leading_recoil, next_recoil, center_flux):
    """Return the formal 4.5PN and precession-enhanced 5.5PN forces."""
    radius, normal, s, w, u, delta = _orbital_quantities(binary, position, velocity)
    c, nu = binary.c, binary.nu
    complete_recoil = leading_recoil + next_recoil + center_flux
    n_complete = np.sum(normal * complete_recoil, axis=-1, keepdims=True)
    v_complete = np.sum(velocity * complete_recoil, axis=-1, keepdims=True)
    prefactor = binary.G * delta / radius**2
    leading = prefactor / c**2 * (2 * v_complete * normal + n_complete * velocity)
    n_leading = np.sum(normal * leading_recoil, axis=-1, keepdims=True)
    v_leading = np.sum(velocity * leading_recoil, axis=-1, keepdims=True)
    # Eq. (5.32b), independently reproduced by A$RR$5p5PN$CM$rad$nonlocal.
    radial = (
        (3 * nu / 2) * s * n_leading * (3 * w - 5 * s**2)
        + v_leading * (3 * nu * s**2 - 4 * nu * w - u * (8 - nu))
    )
    along_velocity = (
        n_leading * ((9 * nu / 2) * s**2 + (2 - nu) * (w / 2 - 2 * u))
        - (4 + nu) * s * v_leading
    )
    promotion = prefactor / c**4 * (normal * radial + velocity * along_velocity)
    return leading, promotion


def center_flux_rate(binary, position, velocity):
    """Newtonian order reduction of the derivative of the leading CM flux."""
    radius, normal, s, w, u, delta = _orbital_quantities(binary, position, velocity)
    ds, dw, du = w - s**2 - u, -2 * u * s, -u * s
    prefactor = binary.mass * delta * binary.nu**2 * u**2 / (binary.c**7 * radius)
    a = 272 * u * s**2 / 105 - 64 * u * w / 15 + 48 * s**2 * w / 35 - 32 * w**2 / 35
    b = 16 * u / 21 - 24 * s**2 / 7 + 24 * w / 7
    a_rate = (
        (272 * s**2 / 105 - 64 * w / 15) * du
        + (544 * u / 105 + 96 * w / 35) * s * ds
        + (-64 * u / 15 + 48 * s**2 / 35 - 64 * w / 35) * dw
    )
    b_rate = 16 * du / 21 - 48 * s * ds / 7 + 24 * dw / 7
    return prefactor * (
        normal * (a_rate - 3 * s * a - u * s * b)
        + velocity * (a + ds * b + s * b_rate - 2 * s**2 * b)
    )
