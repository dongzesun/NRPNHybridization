"""Adaptive causal integration for a state coupled to a retarded source."""

from contextlib import contextmanager, nullcontext
from dataclasses import dataclass
from multiprocessing import get_context
from traceback import format_exc
from hashlib import sha256
import json
import os
from tempfile import NamedTemporaryFile
from pathlib import Path

import numpy as np

from .compressed_history import CompressedPanelHistory
from .history import PanelHistory, SourcePanel
from .conventions import Binary
from .incoming import (IncomingOrbit, IncomingQuadrupole, ParabolicOrbit, ParabolicQuadrupole,
                       SampledIncomingTail, incoming_identity, incoming_parameters)
from .tail_parts import PartsIncomingTail
from .tails import PrimitiveIncomingTail


# Three-node Lobatto IIIA integration matrix. The stage derivative is quadratic.
_LOBATTO = np.array([[0, 0, 0], [5 / 24, 1 / 3, -1 / 24], [1 / 6, 2 / 3, 1 / 6]])
_NODES = np.array([0, 0.5, 1])


def _interpolate(values):
    first, middle, last = values
    return np.stack((first, 4 * middle - 3 * first - last, 2 * (first - 2 * middle + last)))


def _collocate(rhs, source, start, stop, initial, history, rtol, atol, interpolate_source, first_rate):
    width = stop - start
    times = start + width * _NODES
    times[-1] = stop
    stages = initial + (times - start)[:, None] * first_rate
    for _ in range(16):
        panel = SourcePanel(start, stop, interpolate_source(source, stages, width))
        rates = np.stack((first_rate, *(
            rhs(time, state, history, panel) for time, state in zip(times[1:], stages[1:])
        )))
        updated = initial + width * (_LOBATTO @ rates)
        scale = atol + rtol * np.maximum(np.abs(updated), np.abs(stages))
        if np.max(np.abs(updated - stages) / scale) <= 0.03:
            constant, linear, quadratic = _interpolate(rates)
            dense = np.stack((initial, width * constant, width * linear / 2, width * quadratic / 3))
            source_panel = SourcePanel(start, stop, interpolate_source(source, updated, width))
            return updated[-1], source_panel, dense
        stages = updated
    return None


def _incoming_record_digest(record):
    """Bind the saved force and waveform prefix parameters, including their quadrature grids."""
    keys = ("physical_incoming", "incoming", "parabolic", "sampled_incoming",
            "parts_incoming", "primitive_incoming")
    metadata = {key: record[key] for key in keys if key in record}
    metadata["waveform"] = record.get("canonical", {}).get("incoming")
    data = json.dumps(metadata, sort_keys=True, separators=(",", ":"), allow_nan=False)
    return sha256(data.encode("utf-8")).hexdigest()


@dataclass(frozen=True)
class AcceptedHistory:
    """Last accepted state and its history, sufficient to continue integration."""

    state: np.ndarray
    history: PanelHistory | CompressedPanelHistory
    next_step: float
    configuration: str | None = None
    canonical_history: object = None

    def save(self, path):
        """Store numeric state and committed source panels in a versioned file."""
        record = {
            "format": "eccpn.accepted-history.v2",
            "state": self.state.tolist(),
            "next_step": self.next_step,
            "start": self.history.start,
            "source_shape": self.history.source_shape,
            "identity": self.history.identity,
            "configuration": self.configuration,
            "compressed": isinstance(self.history,CompressedPanelHistory),
            "panels": [
                [panel.start, panel.end, panel.coefficients.tolist()]
                for panel in self.history.panels
            ],
        }
        incoming = self.history.incoming
        if getattr(incoming, "orbit", None) is not None:
            record["physical_incoming"] = incoming_parameters(incoming.binary, incoming.orbit)
        if isinstance(incoming,PrimitiveIncomingTail):
            record["primitive_incoming"] = dict(record["physical_incoming"], nodes=incoming.nodes)
            incoming = incoming.base
        elif isinstance(incoming,PartsIncomingTail):
            record["parts_incoming"] = record["physical_incoming"]
            incoming = incoming.base
        if isinstance(incoming, ParabolicQuadrupole):
            record["parabolic"] = [incoming.cutoff, incoming.time_origin, incoming.tensor.tolist()]
        elif isinstance(incoming,SampledIncomingTail):
            record["sampled_incoming"] = {
                **incoming_parameters(incoming.binary, incoming.orbit), "nodes": incoming.nodes,
            }
        elif incoming is not None:
            correction = None if incoming.subleading is None else incoming.subleading.tolist()
            record["incoming"] = [
                incoming.cutoff, incoming.time_origin, incoming.tensor.tolist(), correction, incoming.log_reference,
            ]
        if self.canonical_history is not None:
            record["canonical"] = self.canonical_history.to_record()
        if "physical_incoming" in record:
            record["incoming_parameters_sha256"] = _incoming_record_digest(record)
        # Publish only a complete serialization. link() retains the exclusive
        # create contract if another writer claims this destination first.
        with NamedTemporaryFile(mode="w", dir=Path(path).parent, encoding="utf-8") as stream:
            json.dump(record, stream, allow_nan=False)
            stream.flush()
            os.link(stream.name, path)

    @classmethod
    def load(cls, path, *, identity):
        """Resume only a checkpoint with the caller's expected problem identity."""
        with Path(path).open(encoding="utf-8") as stream:
            record = json.load(stream)
        physical = record.get("physical_incoming")
        has_remote = any(key in record for key in ("incoming", "parabolic", "sampled_incoming"))
        if physical is None and (record["identity"].startswith("eccpn.incoming.v1:") or
                                 (has_remote and record["format"] == "eccpn.accepted-history.v1")):
            raise ValueError("Checkpoint lacks verified incoming parameters. Replay with its original version "
                             "or recompute from the physical inputs.")
        if record["format"] not in ("eccpn.accepted-history.v1", "eccpn.accepted-history.v2") or record["identity"] != identity:
            raise ValueError("The checkpoint format or problem identity does not match.")
        if (physical is not None and record["format"] == "eccpn.accepted-history.v2"
                and record.get("incoming_parameters_sha256") != _incoming_record_digest(record)):
            raise ValueError("The saved incoming force or waveform parameters failed their checksum.")
        incoming = record.get("incoming")
        parabolic = record.get("parabolic")
        sampled = record.get("sampled_incoming")
        if sum(item is not None for item in (incoming,parabolic,sampled)) > 1:
            raise ValueError("A checkpoint can contain only one incoming model.")
        if incoming is not None:
            remote = IncomingQuadrupole(*incoming)
        elif parabolic is not None:
            remote = ParabolicQuadrupole(*parabolic)
        elif sampled is not None:
            remote = SampledIncomingTail(IncomingOrbit(**sampled["orbit"]), Binary(**sampled["binary"]),
                                         sampled["nodes"])
        else:
            remote = None
        if physical is not None:
            types = {"IncomingOrbit": IncomingOrbit, "ParabolicOrbit": ParabolicOrbit}
            binary = Binary(**physical["binary"])
            if physical["type"] not in types:
                raise ValueError("The checkpoint has an unknown incoming orbit type.")
            orbit = types[physical["type"]](**physical["orbit"])
            if incoming_identity(binary, orbit) != identity:
                raise ValueError("The checkpoint incoming parameters do not match its identity.")
            expected = orbit.force_prefix(binary)
            if type(expected) is not type(remote):
                raise ValueError("The saved force prefix does not match the incoming orbit.")
            if isinstance(remote, SampledIncomingTail):
                if incoming_identity(remote.binary, remote.orbit) != identity:
                    raise ValueError("The sampled force orbit does not match the checkpoint identity.")
            else:
                fields = ("cutoff", "time_origin", "tensor")
                if isinstance(remote, IncomingQuadrupole):
                    fields += ("subleading", "log_reference")
                if any(not np.array_equal(getattr(remote, key), getattr(expected, key)) for key in fields):
                    raise ValueError("The saved force coefficients do not match the incoming parameters.")
                remote = expected
        if "primitive_incoming" in record or "parts_incoming" in record:
            primitive = record.get("primitive_incoming",record.get("parts_incoming"))
            types = {"IncomingOrbit":IncomingOrbit,"ParabolicOrbit":ParabolicOrbit}
            orbit = types[primitive["type"]](**primitive["orbit"])
            binary = Binary(**primitive["binary"])
            if physical is None or incoming_identity(binary, orbit) != identity:
                raise ValueError("The saved force orbit does not match the checkpoint identity.")
            remote = (PartsIncomingTail(remote,orbit,binary) if "parts_incoming" in record
                      else PrimitiveIncomingTail(remote,orbit,binary,primitive["nodes"]))
        history = PanelHistory(
            record["start"], tuple(record["source_shape"]), identity,
            tuple(SourcePanel(*entry) for entry in record["panels"]),
            remote,
        )
        if record.get("compressed",False):
            history = CompressedPanelHistory.from_history(history)
        canonical = None
        if "canonical" in record:
            from .waveform_history import CanonicalHistory

            canonical = CanonicalHistory.from_record(record["canonical"])
            wave = canonical.incoming
            if (wave is None) != (physical is None) or (
                wave is not None and incoming_identity(wave.binary, wave.orbit) != identity
            ):
                raise ValueError("The saved waveform orbit does not match the force identity.")
            if wave is not None:
                masses, dipoles = wave.orbit.conserved_charges(wave.binary)
                if (wave.harmonic_scale != canonical.harmonic_scale or
                        wave.infrared_scale != canonical.infrared_scale or
                        wave.maximum_grade != canonical.maximum_absolute_grade or
                        masses != canonical.mass_grades or
                        any(saved.keys() != dipoles.keys() or
                            any(not np.array_equal(saved[g], dipoles[g]) for g in dipoles)
                            for saved in (wave.dipoles, canonical.dipole_grades))):
                    raise ValueError("The incoming waveform scales, grades or charges differ from the saved history.")
            quadrupole = canonical.moments["M2"][0]
            if quadrupole.start != history.start or quadrupole.end != history.end:
                raise ValueError("Saved wave and force histories must share their first and last cuts.")
        return cls(
            np.array(record["state"], dtype=float), history, float(record["next_step"]),
            record.get("configuration"), canonical,
        )


@dataclass(frozen=True)
class HereditarySolution:
    checkpoint: AcceptedHistory
    boundaries: np.ndarray
    dense: np.ndarray
    rejected_steps: int

    def __call__(self, times):
        times = np.atleast_1d(times)
        if not (np.isfinite(times) & (times >= self.boundaries[0]) & (times <= self.boundaries[-1])).all():
            raise ValueError("Evaluation times must be covered by these accepted steps.")
        indices = np.clip(np.searchsorted(self.boundaries, times, side="right") - 1, 0, len(self.dense) - 1)
        coordinate = ((times - self.boundaries[indices])
                      / (self.boundaries[indices + 1] - self.boundaries[indices]))
        values = self.dense[indices]
        coordinate = coordinate[:, None]
        result = values[:, 0] + coordinate * (
            values[:, 1] + coordinate * (values[:, 2] + coordinate * values[:, 3])
        )
        result[times == self.checkpoint.history.end] = self.checkpoint.state
        return result.T


def _pointwise_source(source, stages, width):
    return _interpolate([source(state) for state in stages])


def _coarse_worker(channel, history, rhs, source, rtol, atol, interpolate_source):
    with channel:
        try:
            while True:
                try:
                    request = channel.recv()
                except EOFError:
                    return
                if request is None:
                    return
                start, end, state, first_rate, committed = request
                if committed:
                    history = history.append(*committed)
                step = _collocate(rhs, source, start, end, state, history, rtol, atol,
                                  interpolate_source, first_rate)
                channel.send((True, None if step is None else step[0]))
        except Exception:
            channel.send((False, format_exc()))


@contextmanager
def _coarse_process(*settings):
    context = get_context("spawn")
    local, remote = context.Pipe()
    process = context.Process(target=_coarse_worker, args=(remote, *settings))
    try:
        process.start()
        remote.close()
        try:
            yield local
            local.send(None)
        except BaseException:
            process.terminate()
            raise
        finally:
            process.join()
            process.close()
    finally:
        local.close()
        remote.close()


def integrate_hereditary(rhs, source, checkpoint, stop, rtol, atol, max_step, *, max_steps=100000,
                         interpolate_source=_pointwise_source, on_accept=None, parallel_orbit=False):
    """Advance with step doubling, committing the two fine source panels on acceptance.

    on_accept(left, right) receives each accepted half as (end_state,
    source_panel, dense_coefficients). The halves are in chronological order.
    parallel_orbit overlaps the coarse error estimate with both fine steps
    in one dedicated process and requires picklable RHS and source functions.
    """
    pool = (_coarse_process(checkpoint.history, rhs, source, rtol, atol, interpolate_source)
            if parallel_orbit else nullcontext(None))
    with pool as worker:
        return _integrate_steps(rhs, source, checkpoint, stop, rtol, atol, max_step, max_steps,
                                interpolate_source, on_accept, worker)


def _integrate_steps(rhs, source, checkpoint, stop, rtol, atol, max_step, max_steps,
                     interpolate_source, on_accept, worker):
    history = checkpoint.history
    state = np.asarray(checkpoint.state, dtype=float)
    current = history.end
    step = min(checkpoint.next_step, max_step)
    if (state.ndim != 1 or not state.size or not np.isfinite(state).all()
            or not np.isfinite(stop) or stop <= current
            or not np.isfinite([rtol, atol, step]).all() or min(rtol, atol, step) <= 0):
        raise ValueError("An advance requires a finite state, later stop and positive finite initial step and tolerances.")
    boundaries = [current]
    dense_steps = []
    rejected = 0
    first_rate = rhs(current, state, history, None)
    committed = ()
    for _ in range(max_steps):
        end = min(stop, current + step)
        middle = current + (end - current) / 2
        if not current < middle < end:
            raise RuntimeError("Hereditary integration step underflow.")
        if worker is None:
            coarse = _collocate(rhs, source, current, end, state, history, rtol, atol, interpolate_source, first_rate)
            if coarse is None:
                step /= 2
                rejected += 1
                continue
            coarse_endpoint = coarse[0]
        else:
            worker.send((current, end, state, first_rate, committed))
            committed = ()
        left = _collocate(rhs, source, current, middle, state, history, rtol, atol, interpolate_source, first_rate)
        right = None
        if left is not None:
            midpoint_history = history.append(left[1])
            midpoint_rate = rhs(middle, left[0], midpoint_history, None)
            right = _collocate(rhs, source, middle, end, left[0], midpoint_history,
                               rtol, atol, interpolate_source, midpoint_rate)
        if worker is not None:
            try:
                ok, coarse_endpoint = worker.recv()
            except EOFError as exc:
                raise RuntimeError("The coarse orbital worker stopped before returning its estimate.") from exc
            if not ok:
                raise RuntimeError(f"Coarse orbital step failed:\n{coarse_endpoint}")
        if right is None or coarse_endpoint is None:
            step /= 2
            rejected += 1
            continue
        error = np.max(np.abs(right[0] - coarse_endpoint) / (atol + rtol * np.maximum(np.abs(state), np.abs(right[0]))))
        if error > 1:
            step /= 2
            rejected += 1
            continue
        state = right[0]
        history = midpoint_history.append(right[1])
        if worker is not None:
            committed = (left[1], right[1])
        if on_accept is not None:
            on_accept(left, right)
        dense_steps.extend((left[2], right[2]))
        boundaries.extend((middle, end))
        current = end
        step = min(max_step, (end - boundaries[-3]) * (2 if error < 0.04 else 1))
        if current == stop:
            return HereditarySolution(
                AcceptedHistory(state, history, step, checkpoint.configuration),
                np.array(boundaries), np.array(dense_steps), rejected,
            )
        first_rate = rhs(current, state, history, None)
    raise RuntimeError("Hereditary integration exceeded its step limit.")
