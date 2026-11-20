"""Bondi cut retiming, polarization transport and harmonic reprojection."""

from dataclasses import dataclass, field
from functools import cache
from math import factorial, pi, sqrt

import numpy as np
from scipy.optimize import brentq


@cache
def angular_grid(polar_points, azimuthal_points):
    z, weights = np.polynomial.legendre.leggauss(polar_points)
    azimuth = np.arange(azimuthal_points) * (2 * pi / azimuthal_points)
    z, azimuth = np.meshgrid(z, azimuth, indexing="ij")
    sin_theta = np.sqrt(1 - z**2)
    vectors = np.stack((sin_theta * np.cos(azimuth), sin_theta * np.sin(azimuth), z), axis=-1)
    return vectors.reshape(-1, 3), np.repeat(weights * (2 * pi / azimuthal_points), azimuthal_points)


def _angles(vectors):
    x, y, z = np.moveaxis(vectors, -1, 0)
    return np.arctan2(np.hypot(x, y), z), np.arctan2(y, x)


def spin_harmonic(ell, m, spin, theta, phi):
    """Spin-weighted Y_lm in the spherical Condon–Shortley convention."""
    if abs(spin) > ell or abs(m) > ell:
        return np.zeros(np.broadcast_shapes(np.shape(theta), np.shape(phi)), dtype=complex)
    other = -spin
    normalization = sqrt(
        (2 * ell + 1) / (4 * pi) * factorial(ell + m) * factorial(ell - m)
        * factorial(ell + other) * factorial(ell - other)
    )
    cosine, sine = np.cos(theta / 2), np.sin(theta / 2)
    lower, upper = max(0, other - m), min(ell + other, ell - m)
    radial = sum(
        (-1)**(index - m + other) * normalization
        * cosine**(2 * ell + other - m - 2 * index) * sine**(2 * index - other + m)
        / (factorial(ell + other - index) * factorial(index) * factorial(ell - m - index)
           * factorial(index - other + m))
        for index in range(lower, upper + 1)
    )
    return (-1)**spin * radial * np.exp(1j * m * phi)


def _basis(theta, phi):
    polar = np.stack((np.cos(theta) * np.cos(phi), np.cos(theta) * np.sin(phi), -np.sin(theta)), axis=-1)
    azimuthal = np.stack((-np.sin(phi), np.cos(phi), np.zeros_like(phi)), axis=-1)
    return polar, azimuthal


def _boost(beta):
    squared = beta @ beta
    gamma = (1 - squared)**-0.5
    result = np.empty((4, 4))
    result[0, 0] = gamma
    result[0, 1:] = result[1:, 0] = -gamma * beta
    result[1:, 1:] = np.eye(3) + gamma**2 / (1 + gamma) * np.outer(beta, beta)
    return result


@dataclass(frozen=True)
class BondiFrame:
    """Passive boost followed by a spatial rotation and a source-frame alpha."""

    boost: tuple[float, float, float] = (0, 0, 0)
    rotation: np.ndarray = field(default_factory=lambda: np.eye(3))
    supertranslation: dict = field(default_factory=dict)
    ell_max: int = 4
    polar_points: int = 14
    azimuthal_points: int = 28

    def __post_init__(self):
        beta, rotation = np.asarray(self.boost, dtype=float), np.array(self.rotation, dtype=float, copy=True)
        if beta.shape != (3,) or not np.isfinite(beta).all() or beta @ beta >= 1:
            raise ValueError("The frame boost must be a finite subluminal Cartesian vector.")
        if (rotation.shape != (3, 3) or not np.isfinite(rotation).all()
                or not np.allclose(rotation @ rotation.T, np.eye(3), rtol=0, atol=2e-13)
                or not np.isclose(np.linalg.det(rotation), 1, rtol=0, atol=2e-13)):
            raise ValueError("The frame rotation must be a proper orthogonal matrix.")
        modes = dict(self.supertranslation)
        for (ell, m), value in modes.items():
            if ell < 0 or abs(m) > ell or not np.isfinite(value):
                raise ValueError("Invalid supertranslation mode.")
            if not np.isclose(modes.get((ell, -m), 0), (-1)**m * np.conj(value), atol=0, rtol=2e-13):
                raise ValueError("The supertranslation coefficients must describe a real scalar field.")
        if self.ell_max < 2:
            raise ValueError("The output harmonic degree must be at least two.")
        self._check_grid(self.ell_max)
        rotation.flags.writeable = False
        object.__setattr__(self, "rotation", rotation)
        object.__setattr__(self, "boost", tuple(beta))
        object.__setattr__(self, "supertranslation", modes)

    def _check_grid(self, degree):
        if self.polar_points < degree + 1 or self.azimuthal_points < 2*degree + 1:
            raise ValueError(f"Angular grid {self.polar_points} x {self.azimuthal_points} requires at least "
                             f"{degree + 1} x {2*degree + 1} points for harmonic degree {degree}.")

    def _grid(self):
        out_vectors, weights = angular_grid(self.polar_points, self.azimuthal_points)
        beta = np.array(self.boost)
        rotation = np.eye(4)
        rotation[1:, 1:] = self.rotation
        inverse = _boost(-beta) @ rotation.T
        reciprocal = inverse[0, 0] + out_vectors @ inverse[0, 1:]
        vectors = (inverse[1:, 0] + out_vectors @ inverse[1:, 1:].T) / reciprocal[:, None]
        theta, phi = _angles(vectors)
        theta_out, phi_out = _angles(out_vectors)
        source_polar, _ = _basis(theta, phi)
        output_polar, output_azimuthal = _basis(theta_out, phi_out)
        forward = rotation @ _boost(beta)
        transported = source_polar @ forward[1:, 1:].T
        transported -= (source_polar @ forward[0, 1:])[:, None] * out_vectors
        phase = np.sum(transported * output_polar, axis=-1) + 1j * np.sum(
            transported * output_azimuthal, axis=-1,
        )
        return theta, phi, theta_out, phi_out, reciprocal, phase.conj()**2, weights


def retime_and_project(frame, times, speed_of_light, sample_native, native_modes):
    """Return transformed H and news modes from a retimed native sampler.

    The sampler receives a flat array of native times in time units and
    returns dictionaries of asymptotic strain and news modes.
    """
    native_modes = tuple(native_modes)
    degree = max([frame.ell_max, *(ell for ell, _ in native_modes),
                  *(ell for (ell, _), value in frame.supertranslation.items() if value != 0)])
    frame._check_grid(degree)
    theta, phi, theta_out, phi_out, conformal, phase, weights = frame._grid()
    alpha = np.zeros_like(theta, dtype=complex)
    shear_offset = np.zeros_like(theta, dtype=complex)
    for (ell, m), coefficient in frame.supertranslation.items():
        alpha += coefficient * spin_harmonic(ell, m, 0, theta, phi)
        if ell >= 2:
            spin_factor = sqrt((ell - 1) * ell * (ell + 1) * (ell + 2))
            shear_offset += coefficient * spin_factor * spin_harmonic(ell, m, 2, theta, phi)
    source_times = times[..., None] / conformal + alpha.real / speed_of_light
    distinct_times,inverse = np.unique(source_times,return_inverse=True)
    channels = sample_native(distinct_times)
    inverse = inverse.reshape(source_times.shape)
    native_basis = {
        (ell, m): spin_harmonic(ell, m, -2, theta, phi) for ell, m in native_modes
    }
    fields = [
        sum(series[key][inverse]*native_basis[key] for key in native_modes if key in series)
        for series in channels
    ]
    fields[0] = phase * (fields[0] - shear_offset.conj()) / conformal
    fields[1] = phase * fields[1] / conformal**2
    transformed = ({}, {})
    for ell in range(2, frame.ell_max + 1):
        for m in range(-ell, ell + 1):
            projector = weights * spin_harmonic(ell, m, -2, theta_out, phi_out).conj()
            for output, field in zip(transformed, fields):
                output[ell, m] = field @ projector
    return transformed


@dataclass(frozen=True)
class IncomingBondiData:
    """Inertial past charges, electric shear and the initial PN Bondi frame."""

    mass_lengths: tuple
    betas: np.ndarray
    shear_strain: dict
    frame: BondiFrame

    def mass_aspect(self, directions):
        """Evaluate the anisotropic incoming mass aspect in native directions."""
        rays = np.asarray(directions, dtype=float)
        if rays.shape[-1:] != (3,) or not np.isfinite(rays).all() or not np.allclose(
            np.linalg.norm(rays, axis=-1), 1.0, atol=2e-13, rtol=0.0,
        ):
            raise ValueError("Use Cartesian unit directions for the native mass aspect.")
        inverse_gamma = np.sqrt(1 - np.sum(self.betas**2, axis=-1))
        return np.sum(np.asarray(self.mass_lengths) * (inverse_gamma / (1 - rays @ self.betas.T))**3, axis=-1)


def incoming_bondi_data(binary, incoming, *, ell_max=9, polar_points=24, azimuthal_points=48):
    """Fix the electric incoming shear of a hyperbolic or parabolic rest frame."""
    from .incoming import IncomingOrbit, ParabolicOrbit

    if isinstance(incoming, IncomingOrbit):
        incoming.conserved_charges(binary)
        relative = np.asarray(incoming.velocity)
        speed_ratio = np.linalg.norm(relative) / binary.c
        tangent = relative / np.linalg.norm(relative)
        transverse = np.asarray(incoming.impact) - tangent * np.dot(incoming.impact, tangent)
        transverse_length = np.linalg.norm(transverse)
        # The inertial momentum is equal and opposite, whereas the specified
        # speed is the difference of the velocities, in the rest frame.
        masses = np.asarray((binary.m1, binary.m2))
        limit = min(masses) * speed_ratio / np.sqrt(1 - speed_ratio**2)
        p = brentq(lambda p: np.sum(p / np.hypot(masses, p)) - speed_ratio, 0., limit,
                   xtol=np.finfo(float).tiny, rtol=1e-14)
        betas = (np.array([1., -1.]) * p / np.hypot(masses, p))[:, None] * tangent
        # Head-on incoming data leave one rotation angle arbitrary.
        if transverse_length <= 1e-13 * np.linalg.norm(incoming.impact):
            transverse = np.eye(3)[np.argmin(np.abs(tangent))]
            transverse -= tangent * np.dot(transverse, tangent)
            transverse_length = np.linalg.norm(transverse)
        radial = transverse / transverse_length
    elif isinstance(incoming, ParabolicOrbit):
        radial, tangent = np.asarray((incoming.pericenter_direction, incoming.transverse_direction), dtype=float)
        betas = np.zeros((2, 3))
    else:
        raise ValueError("A Bondi rest frame requires physical incoming orbit data.")

    rotation = np.stack((radial, tangent, np.cross(radial, tangent)))
    # Construct the frame first to validate coverage and the chosen axes.
    frame = BondiFrame(rotation=rotation, ell_max=ell_max, polar_points=polar_points,
                       azimuthal_points=azimuthal_points)
    rays, weights = angular_grid(polar_points, azimuthal_points)
    theta, phi = _angles(rays)
    polar, azimuthal = _basis(theta, phi)
    lengths = tuple(binary.G / binary.c**2 * np.asarray((binary.m1, binary.m2)))
    gamma = 1 / np.sqrt(1 - np.sum(betas**2, axis=-1))
    strain = 2 * np.sum(
        np.asarray(lengths) * gamma * ((polar - 1j * azimuthal) @ betas.T)**2
        / (1 - rays @ betas.T), axis=-1,
    )
    modes, alpha = {}, {}
    for ell in range(2, ell_max + 1):
        normalization = sqrt((ell - 1) * ell * (ell + 1) * (ell + 2))
        for m in range(ell + 1):
            positive = np.dot(weights * spin_harmonic(ell, m, -2, theta, phi).conj(), strain)
            negative = np.dot(weights * spin_harmonic(ell, -m, -2, theta, phi).conj(), strain)
            electric = (positive + (-1)**m * negative.conj()) / 2
            if abs(positive - electric) > 2e-12 * max(np.max(np.abs(strain)), np.finfo(float).tiny):
                raise ValueError("The projected incoming shear is not electric at the requested angular resolution.")
            modes[ell, m], modes[ell, -m] = electric, (-1)**m * electric.conj()
            alpha[ell, m], alpha[ell, -m] = modes[ell, m] / normalization, modes[ell, -m] / normalization
    frame = BondiFrame(rotation=rotation, supertranslation=alpha, ell_max=ell_max, polar_points=polar_points,
                       azimuthal_points=azimuthal_points)
    betas.flags.writeable = False
    return IncomingBondiData(lengths, betas, modes, frame)
