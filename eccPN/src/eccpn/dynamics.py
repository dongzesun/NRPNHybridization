"""Direct integration in SH coordinates for the local development models."""

from dataclasses import asdict, dataclass, field, replace
from functools import partial
from hashlib import sha256
from importlib.resources import files
import json
from typing import Callable

import numpy as np
from scipy.integrate import solve_ivp

from .conventions import Accuracy, Binary, COMPLETE_SECTOR_MODELS, HEREDITARY_MODELS, SH_MODELS
from .hereditary import AcceptedHistory, integrate_hereditary
from .history import PanelHistory, TailInitialData
from .local_force import higher_local_acceleration
from .nonlinear_force import sh_nonlinear_5pn_terms
from .recoil import recoil_accelerations, recoil_fluxes
from .compressed_history import CompressedPanelHistory
from .tail_parts import PartsIncomingTail, parts_jet
from .tails import PrimitiveIncomingTail, dynamical_tail, primitive_tail_panel, quadrupole_tail_sources


def acceleration(binary, position, velocity, model):
    """Return Newtonian and selected 2.5PN reaction accelerations separately."""
    radius = np.linalg.norm(position, axis=-1, keepdims=True)
    if not (np.isfinite(radius) & (radius > 0)).all():
        raise ValueError("The relative separation must be finite and nonzero.")
    normal = position / radius
    potential = binary.G * binary.mass / radius
    newtonian = -potential * normal / radius
    if model == "newtonian":
        return newtonian, np.zeros_like(newtonian)
    radial_speed = np.sum(normal * velocity, axis=-1, keepdims=True)
    speed_squared = np.sum(velocity**2, axis=-1, keepdims=True)
    # Blanchet et al., arXiv:2601.06743v3, Eq. (6.6).
    coefficient = (8 / 5) * binary.nu * potential**2 / (radius * binary.c**5)
    reaction = coefficient * (
        radial_speed * (3 * speed_squared + (17 / 3) * potential) * normal
        - (speed_squared + 3 * potential) * velocity
    )
    return newtonian, reaction


@dataclass(frozen=True)
class Trajectory:
    binary: Binary
    time_span: tuple[float, float]
    accuracy: Accuracy
    _solution: Callable = field(repr=False)
    history: PanelHistory | None = None
    asymptotic_join_residual: tuple[float, float, float] | None = None
    _canonical_history: object = field(default=None, repr=False)

    def _values(self, times):
        times = np.asarray(times, dtype=float)
        if times.size == 0:
            raise ValueError("At least one evaluation time is required.")
        if not (
            np.isfinite(times) & (times >= self.time_span[0]) & (times <= self.time_span[1])
        ).all():
            raise ValueError("Evaluation times must lie within the trajectory interval.")
        values = np.moveaxis(self._solution(times.reshape(-1)), 0, -1)
        return values.reshape(times.shape + (values.shape[-1],))

    def state(self, times):
        """Evaluate accepted dense steps in source time without extrapolation."""
        values = self._values(times)
        return values[..., :3], values[..., 3:6]

    def recoil(self, times):
        """Return the two recoil contributions and CM integral at the given times."""
        if self.accuracy.dynamics not in SH_MODELS:
            raise ValueError("This trajectory does not carry recoil or the center-of-mass integral.")
        values = self._values(times)
        return values[..., 6:9], values[..., 9:12], values[..., 12:15]

    def nonlinear_5pn_terms(self,times):
        """Return the selected source-harmonic acceleration contributions."""
        position,velocity = self.state(times)
        sectors = self.accuracy.selected_nonlinear_5pn_sectors
        return sh_nonlinear_5pn_terms(self.binary,position,velocity,sectors) if sectors else {}

    def body_states(self, times):
        """Evaluate the two standard harmonic worldlines on covered source times."""
        from .worldlines import body_states

        return body_states(self, times)

    @property
    def checkpoint(self):
        """Return the last committed hereditary state, including its source history."""
        if self.history is None:
            raise ValueError("The selected local trajectory does not have a hereditary checkpoint.")
        return replace(self._solution.checkpoint, canonical_history=self._canonical_history)


def _resolve_full_accuracy(binary,accuracy,prehistory):
    if accuracy.dynamics!="full":
        return accuracy
    if prehistory is None:
        raise NotImplementedError("The full SH trajectory requires a matched causal past and recoil state.")
    geometric_time = binary.G*binary.mass/binary.c**3
    settings = {}
    if accuracy.causal_tail_interpolation is None:
        settings.update(causal_tail_interpolation="integration_by_parts",
                        causal_tail_split_time=accuracy.causal_tail_split_time or 480*geometric_time)
    physical_incoming = (getattr(prehistory,"incoming_orbit",None) is not None
                         or getattr(getattr(prehistory,"history",None),"incoming",None) is not None)
    if accuracy.radiative_tail_split_time is None and accuracy.canonical_grade is not None and physical_incoming:
        settings["radiative_tail_split_time"] = 7680*geometric_time
    return replace(accuracy,**settings)


def evolve(binary, initial_state, prehistory, time_span, accuracy=None):
    """Integrate the selected local model from the explicit SH initial state."""
    accuracy = Accuracy() if accuracy is None else accuracy
    accuracy = _resolve_full_accuracy(binary,accuracy,prehistory)
    if accuracy.dynamics not in HEREDITARY_MODELS and prehistory is not None:
        raise NotImplementedError("The local development models accept prehistory=None.")
    start, stop = map(float, time_span)
    if not np.isfinite([start, stop]).all() or stop <= start:
        raise ValueError("time_span must contain finite, increasing times.")
    position = np.asarray(initial_state.position, dtype=float)
    velocity = np.asarray(initial_state.velocity, dtype=float)
    if position.shape != (3,) or velocity.shape != (3,):
        raise ValueError("Position and velocity must each have three Cartesian components.")
    initial = np.concatenate((position, velocity))
    if not np.isfinite(initial).all():
        raise ValueError("The initial state must be finite.")
    acceleration(binary, position, velocity, accuracy.dynamics)

    if accuracy.dynamics in HEREDITARY_MODELS:
        return _evolve_tail(binary, initial, prehistory, start, stop, accuracy)

    def rhs(_, state):
        velocity = state[3:]
        newtonian, reaction = acceleration(binary, state[:3], velocity, accuracy.dynamics)
        return np.concatenate((velocity, newtonian + reaction))

    solution = solve_ivp(
        rhs, (start, stop), initial, method="DOP853", dense_output=True,
        rtol=accuracy.rtol, atol=accuracy.atol, max_step=accuracy.max_step,
    )
    if not solution.success:
        raise RuntimeError(f"Orbital integration failed: {solution.message}")
    return Trajectory(binary, (start, stop), accuracy, solution.sol)


def _configuration(binary, accuracy):
    payload = sha256(json.dumps((asdict(binary), asdict(accuracy)), sort_keys=True).encode())
    package = files("eccpn")
    names = ("dynamics.py", "tails.py", "history.py", "incoming.py", "recoil.py",
             "hereditary.py", "local_force.py", "generated/sh_local.py")
    if accuracy.dynamics in COMPLETE_SECTOR_MODELS:
        names += ("nonlinear_force.py","generated/radiative_reaction_5pn.py",
                  "conventions.py","radiation.py","consistent_tail.py","radiation_kernels.py",
                  "radiative_quadrupole.py","luminosity.py","memory.py")
    if accuracy.compress_history:
        names += ("compressed_history.py",)
    if accuracy.causal_tail_interpolation == "integration_by_parts":
        names += ("tail_parts.py",)
    if accuracy.canonical_grade is not None:
        names += ("waveform_history.py", "compact_moments.py", "waveform_jets.py", "multipoles.py", "higher_moments.py",
                  "gauge_moments.py", "canonical.py", "mpm.py", "newtonian_moments.py", "leading_modes.py",
                  "projection.py", "asymptotic_waveform.py", "incoming_waveform.py",
                  "generated/mpm_tables.py", "generated/source_mass_quadrupole.py",
                  "generated/gauge_moments.py", "generated/octupole_current_quadrupole.py",
                  "generated/higher_multipoles.py")
        if accuracy.waveform_pn == 4.0:
            names += ("generated/local_quadrupole_4pn.py", "generated/partial_quadrupole_4pn.py")
    for name in names:
        payload.update(package.joinpath(*name.split("/")).read_bytes())
    return payload.hexdigest()


def _evolve_tail(binary, initial, prehistory, start, stop, accuracy):
    has_recoil = accuracy.dynamics in SH_MODELS
    nonlinear_sectors = accuracy.selected_nonlinear_5pn_sectors
    uses_parts = accuracy.causal_tail_interpolation == "integration_by_parts"
    asymptotic_residual = None
    orbital_prefix = None
    incoming_orbit = None
    configuration = _configuration(binary, accuracy)
    restart = prehistory if isinstance(prehistory, AcceptedHistory) else None
    if restart is not None:
        if restart.configuration != configuration:
            raise ValueError("Checkpoint evolution configuration differs. Replay with its original version "
                             "and settings or recompute from the original inputs.")
        if (restart.history.end != start
                or not np.array_equal(initial, restart.state[:6])
                or (accuracy.canonical_grade is not None) != (restart.canonical_history is not None)):
            raise ValueError("The checkpoint, start state and supported waveform configuration must match.")
        initial = restart.state
        prehistory = restart.history
    elif has_recoil:
        if not isinstance(prehistory, TailInitialData):
            raise ValueError("This model requires a matched quadrupole and recoil prehistory.")
        initial = np.concatenate((initial, prehistory.recoil_state()))
        asymptotic_residual = prehistory.asymptotic_join_residual
        orbital_prefix = prehistory.orbital_prefix
        incoming_orbit = prehistory.incoming_orbit
        prehistory = prehistory.quadrupole
        if incoming_orbit is not None:
            from .incoming import incoming_identity

            force = prehistory.incoming
            if (getattr(force, "orbit", None) is None or
                    incoming_identity(binary, incoming_orbit) != prehistory.identity or
                    incoming_identity(force.binary, force.orbit) != prehistory.identity):
                raise ValueError("The physical force and orbital prefixes must share the incoming identity.")
    source_shape = (5,3,3) if uses_parts and restart is not None else (3,3,3)
    if (not isinstance(prehistory, (PanelHistory,CompressedPanelHistory)) or not (prehistory.panels or prehistory.incoming)
            or prehistory.end != start or prehistory.source_shape != source_shape):
        raise ValueError("The tail development model requires finite source panels ending at the initial time.")
    if accuracy.canonical_grade is not None and restart is None:
        finite_match = (orbital_prefix is not None and orbital_prefix.boundaries[-1] == start
                        and orbital_prefix.identity == prehistory.identity and prehistory.incoming is None)
        remote_match = (incoming_orbit is not None and incoming_orbit.cutoff == start
                        and prehistory.incoming is not None and orbital_prefix is None)
        if not (finite_match or remote_match):
            raise ValueError("The canonical waveform requires a shared finite or incoming orbital history.")

    def source(state):
        if uses_parts and prehistory.source_shape == (5,3,3):
            return parts_jet(binary,state)
        return quadrupole_tail_sources(binary, state[:6])

    expected = source(initial)
    if prehistory.panels:
        joined_source = np.sum(prehistory.panels[-1].coefficients, axis=0)
        if np.max(np.abs(joined_source - expected)) > 1e-10 * np.max(np.abs(expected)):
            raise ValueError("The tail prehistory must match all three quadrupole derivatives at the join.")
    else:
        remote = prehistory.incoming.at_cutoff()
        relative = np.linalg.norm(remote - expected, axis=(-2, -1)) / np.linalg.norm(remote, axis=(-2, -1))
        if (asymptotic_residual is None or not np.allclose(relative, asymptotic_residual, rtol=1e-10, atol=1e-14)):
            raise ValueError("An analytic remote join requires the measured asymptotic residuals.")

    solver_options = {}
    if accuracy.causal_tail_interpolation in ("primitive","integration_by_parts"):
        solver_options["interpolate_source"] = partial(primitive_tail_panel,binary)
        if accuracy.causal_tail_interpolation == "primitive" and restart is None and incoming_orbit is not None:
            prehistory = replace(prehistory,incoming=PrimitiveIncomingTail(prehistory.incoming,incoming_orbit,binary))
    if uses_parts and restart is None:
        if incoming_orbit is not None:
            remote = PartsIncomingTail(prehistory.incoming,incoming_orbit,binary)
            prehistory = PanelHistory(prehistory.start,(5,3,3),prehistory.identity,incoming=remote)
        elif orbital_prefix is not None:
            from .history import SourcePanel
            def sample(t):
                return np.concatenate(orbital_prefix.state(t))
            panels = []
            for a,b in zip(orbital_prefix.boundaries[:-1],orbital_prefix.boundaries[1:]):
                panel = primitive_tail_panel(binary,partial(parts_jet,binary),(sample(a),sample((a+b)/2),sample(b)),b-a)
                panels.append(SourcePanel(a,b,panel))
            prehistory = PanelHistory(prehistory.start,(5,3,3),prehistory.identity,tuple(panels))
        else:
            raise ValueError("The integration-by-parts tail needs a supplied phase-space prehistory.")

    def rhs(time, state, history, trial):
        position, velocity = state[:3], state[3:6]
        newtonian, reaction = acceleration(binary, position, velocity, accuracy.dynamics)
        tail = dynamical_tail(binary,time,state[:6],history,trial,split_time=accuracy.causal_tail_split_time)
        if has_recoil:
            conservative, matter = higher_local_acceleration(
                binary, position, velocity, accuracy.harmonic_scales,
            )
            leading, next_order, center = recoil_fluxes(binary, position, velocity)
            leading_recoil, next_recoil = state[6:9], state[9:12]
            recoil, promoted = recoil_accelerations(
                binary, position, velocity, leading_recoil, next_recoil, center,
            )
            relative = newtonian + reaction + tail + conservative + matter + recoil + promoted
            if nonlinear_sectors:
                relative += sum(sh_nonlinear_5pn_terms(binary,position,velocity,nonlinear_sectors).values())
            return np.concatenate((velocity, relative, leading, next_order, leading_recoil + next_recoil + center))
        return np.concatenate((velocity, newtonian + reaction + tail))

    if accuracy.compress_history and isinstance(prehistory,PanelHistory):
        prehistory = CompressedPanelHistory.from_history(prehistory)
    step = min(accuracy.max_step, (stop - start) / 10)
    solution = integrate_hereditary(
        rhs, source, restart if restart is not None else AcceptedHistory(initial, prehistory, step, configuration), stop,
        accuracy.rtol, accuracy.atol, accuracy.max_step, **solver_options,
    )
    trajectory = Trajectory(
        binary, (start, stop), accuracy, solution, solution.checkpoint.history, asymptotic_residual,
    )
    if accuracy.canonical_grade is not None:
        from .waveform_history import build_canonical_history

        if restart is None and incoming_orbit is not None:
            from .incoming_waveform import IncomingWaveform

            if not np.array_equal(np.concatenate(incoming_orbit.state_at(binary, start)), initial[:6]):
                raise ValueError("The incoming waveform must share the force's phase-space join.")
            mass_grades, dipole_grades = incoming_orbit.conserved_charges(binary)
            remote = IncomingWaveform(
                binary, incoming_orbit, accuracy.harmonic_scales[0], accuracy.incoming_infrared_scale,
                dipole_grades, accuracy.canonical_grade, accuracy.incoming_waveform_nodes,
            )
            history = build_canonical_history(
                binary, solution.boundaries, trajectory.state, mass_grades, dipole_grades,
                remote.harmonic_scale, remote.infrared_scale, accuracy.canonical_grade, incoming=remote,
            )
        elif restart is None:
            def sampled_state(time):
                return orbital_prefix.state(time) if time < start else trajectory.state(time)

            if not np.allclose(np.concatenate(orbital_prefix.state(start)), initial[:6], rtol=1e-12, atol=1e-14):
                raise ValueError("The supplied orbital prefix and trajectory must share the phase-space join.")
            history = build_canonical_history(
                binary, (*orbital_prefix.boundaries[:-1], *solution.boundaries), sampled_state,
                orbital_prefix.mass_grades, orbital_prefix.dipole_grades,
                accuracy.harmonic_scales[0], orbital_prefix.infrared_scale, accuracy.canonical_grade,
            )
        else:
            past = restart.canonical_history
            history = past.extend(build_canonical_history(
                binary, solution.boundaries, trajectory.state, past.mass_grades, past.dipole_grades,
                past.harmonic_scale, past.infrared_scale, accuracy.canonical_grade,
            ))
        trajectory = replace(trajectory, _canonical_history=history)
    return trajectory
