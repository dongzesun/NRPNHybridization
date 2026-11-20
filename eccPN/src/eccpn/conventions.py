"""Units and accuracy for the local development models."""

from dataclasses import dataclass

import numpy as np


TAIL_MODEL = "newtonian+2.5PN+4PN-tail"
SH_LOCAL_TAIL_MODEL = "sh-through-4.5PN+recoil"
SH_RADIATIVE_5PN_MODEL = "sh-through-5PN-radiative-sectors"
COMPLETE_SECTOR_MODELS = ("full",SH_RADIATIVE_5PN_MODEL)
SH_MODELS = (SH_LOCAL_TAIL_MODEL,*COMPLETE_SECTOR_MODELS)
HEREDITARY_MODELS = (TAIL_MODEL,*SH_MODELS)
NONLINEAR_5PN_SECTORS = ("failed_tail","memory","reaction_square")


@dataclass(frozen=True)
class Binary:
    """Rest masses and constants in a consistent unit system."""

    m1: float
    m2: float
    G: float = 1.0
    c: float = 1.0

    def __post_init__(self):
        values = np.array([self.m1, self.m2, self.G, self.c, self.mass], dtype=float)
        if not (np.isfinite(values) & (values > 0)).all():
            raise ValueError("Rest masses, total mass, G and c must be finite and positive.")

    @property
    def mass(self):
        return self.m1 + self.m2

    @property
    def nu(self):
        return (self.m1 / self.mass) * (self.m2 / self.mass)

    @property
    def mu(self):
        return self.mass * self.nu


@dataclass(frozen=True)
class InitialState:
    """Relative position and velocity in standard harmonic coordinates."""

    position: tuple[float, float, float]
    velocity: tuple[float, float, float]


@dataclass(frozen=True)
class Accuracy:
    """Explicit development model selection and integration settings."""

    dynamics: str = "full"
    rtol: float = 1e-10
    atol: float = 1e-12
    max_step: float = np.inf
    radiative_scale: float = 1.0
    harmonic_scales: tuple[float, float] = (1.0, 1.0)
    waveform_pn: float | None = None
    retain_memory_excess: bool | None = None
    incoming_infrared_scale: float = 1.0
    incoming_waveform_nodes: int = 32
    compress_history: bool | None = None
    causal_tail_interpolation: str | None = None
    causal_tail_split_time: float | None = None
    radiative_tail_split_time: float | None = None
    nonlinear_5pn_sectors: tuple[str, ...] | None = None

    @property
    def selected_nonlinear_5pn_sectors(self):
        if self.dynamics not in COMPLETE_SECTOR_MODELS:
            return ()
        return NONLINEAR_5PN_SECTORS if self.nonlinear_5pn_sectors is None else self.nonlinear_5pn_sectors

    @property
    def canonical_grade(self):
        """Explicit 1/c grade of the selected canonical-history waveform."""
        return {3.0: 6, 3.5: 7, 4.0: 7}.get(self.waveform_pn)

    def __post_init__(self):
        if self.dynamics not in ("newtonian","newtonian+2.5PN",*HEREDITARY_MODELS):
            raise ValueError("Unknown dynamics model.")
        full = self.dynamics=="full"
        if full and self.waveform_pn is None:
            object.__setattr__(self,"waveform_pn",4.)
        if self.retain_memory_excess is None:
            object.__setattr__(self,"retain_memory_excess",full and self.canonical_grade==7)
        if self.compress_history is None:
            object.__setattr__(self,"compress_history",full)
        if not full and self.causal_tail_interpolation is None:
            object.__setattr__(self,"causal_tail_interpolation","pointwise")
        tolerances = np.array([self.rtol, self.atol], dtype=float)
        if not (np.isfinite(tolerances) & (tolerances > 0)).all():
            raise ValueError("Integration tolerances must be finite and positive.")
        if not self.max_step > 0:
            raise ValueError("max_step must be positive.")
        if not np.isfinite(self.radiative_scale) or self.radiative_scale <= 0:
            raise ValueError("radiative_scale must be a positive time.")
        scales = np.asarray(self.harmonic_scales, dtype=float)
        if scales.shape != (2,) or not (np.isfinite(scales) & (scales > 0)).all():
            raise ValueError("harmonic_scales holds the two positive effective CM scale lengths.")
        object.__setattr__(self, "harmonic_scales", tuple(scales))
        if self.waveform_pn is not None and (
            self.waveform_pn not in (1.5, 3.0, 3.5, 4.0) or self.dynamics not in SH_MODELS
        ):
            raise ValueError("Select a supported staged waveform order and the higher-order SH trajectory.")
        if self.waveform_pn == 4.0 and scales[0] != scales[1]:
            raise ValueError("The 4PN quadrupole uses the common harmonic ultraviolet scale.")
        if self.retain_memory_excess and self.canonical_grade != 7:
            raise ValueError("Retained higher formal memory requires the complete 3.5PN news inventory.")
        if self.compress_history and self.dynamics not in HEREDITARY_MODELS:
            raise ValueError("Causal history compression requires a hereditary evolution.")
        if (self.causal_tail_interpolation not in (None,"pointwise","primitive","integration_by_parts")
                or (self.causal_tail_interpolation != "pointwise" and self.dynamics not in SH_MODELS)):
            raise ValueError("Choose pointwise, primitive or integration_by_parts for the SH tail.")
        if self.causal_tail_interpolation == "integration_by_parts":
            if (self.causal_tail_split_time is None or not np.isfinite(self.causal_tail_split_time)
                    or self.causal_tail_split_time<=0):
                raise ValueError("Integration by parts requires a positive causal split lag.")
        elif self.causal_tail_split_time is not None and self.causal_tail_interpolation is not None:
            raise ValueError("A split lag is used only for the integration-by-parts tail.")
        if self.causal_tail_interpolation is None and self.causal_tail_split_time is not None:
            if not np.isfinite(self.causal_tail_split_time) or self.causal_tail_split_time<=0:
                raise ValueError("The causal split lag must be finite and positive.")
        if self.nonlinear_5pn_sectors is not None:
            sectors = tuple(self.nonlinear_5pn_sectors)
            if (self.dynamics!=SH_RADIATIVE_5PN_MODEL or len(set(sectors))!=len(sectors)
                    or not set(sectors)<=set(NONLINEAR_5PN_SECTORS)):
                raise ValueError("Select distinct radiative 5PN sectors with the SH radiative model.")
            object.__setattr__(self,"nonlinear_5pn_sectors",sectors)
        if self.radiative_tail_split_time is not None and (
            not np.isfinite(self.radiative_tail_split_time) or self.radiative_tail_split_time<=0
            or self.canonical_grade is None
        ):
            raise ValueError("Select a positive radiative split lag for a native canonical waveform.")
        if (not np.isfinite(self.incoming_infrared_scale) or self.incoming_infrared_scale <= 0
                or not isinstance(self.incoming_waveform_nodes, int) or self.incoming_waveform_nodes < 8):
            raise ValueError("The incoming waveform requires a positive infrared scale and at least eight nodes.")
