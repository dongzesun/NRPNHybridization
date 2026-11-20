"""Canonical source quadrature on the incoming SH phase expansion."""

from dataclasses import asdict, dataclass
from functools import cached_property

import numpy as np

from .conventions import Binary
from .incoming import IncomingOrbit, ParabolicOrbit, incoming_grid
from .waveform_jets import canonical_waveform_jets


@dataclass(frozen=True)
class IncomingWaveform:
    """One compactified quadrature grid shared by the remote mode inventory."""

    binary: Binary
    orbit: IncomingOrbit | ParabolicOrbit
    harmonic_scale: float
    infrared_scale: float
    dipoles: dict
    maximum_grade: int
    order: int

    @cached_property
    def grid(self):
        # q=0 is past infinity, q=1 is the numerical join. The trigonometric
        # map clusters nodes at both ends to resolve powers and the lag log.
        return incoming_grid(self.order,self.orbit.cutoff,self.orbit.time_origin)

    def jet(self, time):
        return canonical_waveform_jets(
            self.binary, *self.orbit.state_at(self.binary, time), self.harmonic_scale,
            self.infrared_scale, self.dipoles, self.maximum_grade,
        )

    @cached_property
    def join(self):
        return self.jet(self.orbit.cutoff)

    @cached_property
    def sampled(self):
        data = [self.jet(time) for time in self.grid[0]]
        return {name: {grade: np.stack([cut[name][grade] for cut in data]) for grade in grades}
                for name, grades in data[0].items()}

    def to_record(self):
        return {
            "binary": asdict(self.binary), "orbit": asdict(self.orbit), "harmonic": self.harmonic_scale,
            "kind": "hyperbolic" if isinstance(self.orbit, IncomingOrbit) else "parabolic",
            "infrared": self.infrared_scale, "dipoles": {g: np.asarray(v).tolist() for g,v in self.dipoles.items()},
            "grade": self.maximum_grade, "nodes": self.order,
        }

    @classmethod
    def from_record(cls, record):
        source = record["orbit"]
        kind = record.get("kind", "hyperbolic")
        if kind == "hyperbolic":
            orbit = IncomingOrbit(**{**source, "velocity": tuple(source["velocity"]), "impact": tuple(source["impact"])})
        elif kind == "parabolic":
            orbit = ParabolicOrbit(**{**source, "pericenter_direction": tuple(source["pericenter_direction"]),
                                     "transverse_direction": tuple(source["transverse_direction"])})
        else:
            raise ValueError("The checkpoint has an unknown incoming waveform model.")
        return cls(
            Binary(**record["binary"]), orbit,
            record["harmonic"], record["infrared"], {int(g): np.asarray(v) for g,v in record["dipoles"].items()},
            record["grade"], record["nodes"],
        )


@dataclass(frozen=True)
class IncomingMomentHistory:
    """One canonical mode's committed panels and source-time prefix."""

    finite: object
    past: IncomingWaveform
    name: str
    grade: int

    @property
    def start(self):
        return self.finite.start

    @property
    def end(self):
        return self.finite.end

    @property
    def source_shape(self):
        return self.finite.source_shape

    @property
    def panels(self):
        return self.finite.panels

    def value(self, time):
        return self.finite.value(time)

    def integrate_jet(self, time, scale, power, orders):
        times, weights = self.past.grid
        if time >= self.start:
            finite = self.finite.integrate_jet(time, scale, power, orders)
        elif time >= self.past.orbit.cutoff:
            finite = 0.
        else:
            raise ValueError("The remote quadrature must reach the requested cut.")
        kernel = weights * np.log((time-times)/scale)**power
        return finite + np.tensordot(kernel, self.past.sampled[self.name][self.grade][:, orders], axes=1)

    def start_jump(self, order):
        return self.finite.start_jump(order)-self.past.join[self.name][self.grade][order]
