"""Sample source moments on accepted orbital panels, in chronological order."""

from collections import deque

import numpy as np

from .compact_moments import CompactMomentHistory, pack_stf
from .hereditary import _interpolate
from .history import PanelHistory, SourcePanel
from .waveform_jets import canonical_waveform_jets


def _stored_rank(name):
    if name.startswith(("M", "S")) and int(name[1:]) >= 3:
        return int(name[1:])


def _sample_batch(binary, harmonic_scale, infrared_scale, dipoles, grade, samples):
    output = []
    for time, position, velocity in samples:
        jets = canonical_waveform_jets(binary, position, velocity, harmonic_scale,
                                       infrared_scale, dipoles, grade)
        for name, inventory in jets.items():
            rank = _stored_rank(name)
            if rank is not None:
                jets[name] = {g: pack_stf(jet, rank) for g, jet in inventory.items()}
        jets["orbital"] = {0: np.stack((position, velocity))}
        output.append((time, jets))
    return output


class CanonicalSampler:
    """Stream alternating boundary and midpoint samples into quadratic panels."""

    def __init__(self, binary, mass_grades, dipole_grades, harmonic_scale, infrared_scale,
                 maximum_absolute_grade, incoming, executor):
        masses = {0, 2} if maximum_absolute_grade == 6 else {0, 2, 4}
        dipoles = {0} if maximum_absolute_grade == 6 else {0, 2}
        if (maximum_absolute_grade not in (6, 7) or not masses <= mass_grades.keys()
                or not dipoles <= dipole_grades.keys() or mass_grades[0] != binary.mass
                or any(g not in (0, 2, 4) or not np.isfinite(v) for g, v in mass_grades.items())):
            raise ValueError("Supply all finite conserved mass and dipole grades required by the waveform budget.")
        self.settings = (binary, harmonic_scale, infrared_scale, dipole_grades, maximum_absolute_grade)
        self.history_args = (dict(mass_grades), dict(dipole_grades), harmonic_scale, infrared_scale,
                             maximum_absolute_grade, incoming)
        self.executor = executor
        self.batch = []
        self.pending = deque()
        self.max_pending = max(2, 2 * (getattr(executor, "_max_workers", None) or 16))
        self.panels = {}
        self.left = None
        self.middle = None
        self.midpoint_at_next_boundary = False

    def __enter__(self):
        return self

    def __exit__(self, *exc):
        for future in self.pending:
            future.cancel()

    def _collect(self, samples):
        for time, jets in samples:
            if self.left is None:
                self.left = (time, jets)
            elif self.middle is None:
                self.middle = jets
            else:
                start, first = self.left
                for name, inventory in first.items():
                    panels = self.panels.setdefault(name, {})
                    for grade, value in inventory.items():
                        coefficients = _interpolate((value, self.middle[name][grade], jets[name][grade]))
                        panels.setdefault(grade, []).append(SourcePanel(start, time, coefficients))
                self.left = (time, jets)
                self.middle = None

    def _submit(self):
        future = self.executor.submit(_sample_batch, *self.settings, self.batch)
        self.batch = []
        self.pending.append(future)
        if len(self.pending) >= self.max_pending:
            self._collect(self.pending.popleft().result())

    def _sample(self, time, state):
        position, velocity = state
        if self.executor is None:
            self._collect(_sample_batch(*self.settings, [(time, position, velocity)]))
        else:
            # A caller's state function may reuse its own output buffer.
            self.batch.append((time, np.copy(position), np.copy(velocity)))
            if len(self.batch) == 8:
                self._submit()

    def boundary(self, time, state):
        if self.midpoint_at_next_boundary:
            self._sample(time, state)
            self.midpoint_at_next_boundary = False
        self._sample(time, state)

    def interval(self, start, stop, state):
        self.boundary(start, state(start))
        middle = (start + stop) / 2
        if middle == stop:
            # The dense solution selects the following panel at an internal boundary.
            self.midpoint_at_next_boundary = True
        else:
            self._sample(middle, state(middle))

    def accepted(self, left, right):
        for _, panel, dense in (left, right):
            def state(time):
                u = (time - panel.start) / (panel.end - panel.start)
                value = dense[0] + u * (dense[1] + u * (dense[2] + u * dense[3]))
                return value[:3], value[3:6]

            self.interval(panel.start, panel.end, state)

    def finish(self):
        from .waveform_history import CanonicalHistory

        if self.batch:
            self._submit()
        while self.pending:
            self._collect(self.pending.popleft().result())

        def history(name, grade, panels):
            finite = PanelHistory(panels[0].start, panels[0].coefficients.shape[1:],
                                  f"{name}-grade-{grade}", tuple(panels))
            rank = _stored_rank(name)
            return finite if rank is None else CompactMomentHistory(finite, rank)

        moments = {name: {g: history(name, g, panels) for g, panels in inventory.items()}
                   for name, inventory in self.panels.items()}
        return CanonicalHistory(moments, *self.history_args)
