"""Named nonlinear 5PN acceleration sectors in source harmonic coordinates."""

import numpy as np

from .conventions import NONLINEAR_5PN_SECTORS


def sh_nonlinear_5pn_terms(binary,position,velocity,sectors=NONLINEAR_5PN_SECTORS):
    """Evaluate the three matched radiative sectors on one or more SH states."""
    from .generated import radiative_reaction_5pn as coefficients

    position = np.asarray(position,dtype=float)
    velocity = np.asarray(velocity,dtype=float)
    if position.shape!=velocity.shape or position.shape[-1:]!=(3,) or not np.isfinite(velocity).all():
        raise ValueError("Supply matching finite relative Cartesian states.")
    radius = np.linalg.norm(position,axis=-1,keepdims=True)
    if not (np.isfinite(radius)&(radius>0)).all():
        raise ValueError("The relative separation must be positive and finite.")
    normal = position/radius
    arguments = (binary.G*binary.mass/radius,
                 np.sum(normal*velocity,axis=-1,keepdims=True),
                 np.sum(velocity**2,axis=-1,keepdims=True),binary.nu)
    functions = {
        "failed_tail":(coefficients.failed_tail_n,coefficients.failed_tail_v),
        "memory":(coefficients.memory_n,coefficients.memory_v),
        "reaction_square":(coefficients.reaction_square_n,coefficients.reaction_square_v),
    }
    return {name:(functions[name][0](*arguments)*normal+functions[name][1](*arguments)*velocity)
                 /(radius*binary.c**10) for name in sectors}
