"""Leading mass quadrupole in the native frame of the trajectory."""

from dataclasses import dataclass, replace

import numpy as np

from .conventions import Accuracy, SH_MODELS
from .dynamics import acceleration
from .frame import BondiFrame, retime_and_project
from .local_force import higher_local_acceleration
from .lower_waveform import waveform_15
from .luminosity import canonical_local_luminosity, mode_luminosity, nonlinear_luminosity
from .recoil import recoil_accelerations, recoil_fluxes
from .tails import dynamical_tail, radiative_tail


def _outer(left, right):
    return left[..., :, None] * right[..., None, :]


def _stf(tensor):
    symmetric = (tensor + np.swapaxes(tensor, -1, -2)) / 2
    trace = np.trace(tensor, axis1=-2, axis2=-1)
    return symmetric - trace[..., None, None] * np.eye(3) / 3


def _quadrupole_jets(trajectory, times, position, velocity):
    """Return the reduced quadrupole and its Newtonian and model derivatives."""
    binary = trajectory.binary
    newtonian, reaction = acceleration(
        binary, position, velocity, trajectory.accuracy.dynamics,
    )
    if trajectory.accuracy.dynamics in SH_MODELS:
        conservative, matter = higher_local_acceleration(
            binary, position, velocity, trajectory.accuracy.harmonic_scales,
        )
        leading_recoil, next_recoil, _ = trajectory.recoil(times)
        _, _, center = recoil_fluxes(binary, position, velocity)
        recoil, promoted = recoil_accelerations(
            binary, position, velocity, leading_recoil, next_recoil, center,
        )
        reaction = reaction + conservative + matter + recoil + promoted
        nonlinear = trajectory.nonlinear_5pn_terms(times)
        if nonlinear:
            reaction += sum(nonlinear.values())
    if trajectory.history is not None:
        states = np.concatenate((position, velocity), axis=-1).reshape(-1, 6)
        tail_acceleration = np.stack([
            dynamical_tail(binary,time,state,trajectory.history,split_time=trajectory.accuracy.causal_tail_split_time)
            for time, state in zip(times.flat, states)
        ])
        reaction = reaction + tail_acceleration.reshape(velocity.shape)
    radius = np.linalg.norm(position, axis=-1, keepdims=True)
    radial_rate = np.sum(position * velocity, axis=-1, keepdims=True) / radius**2
    newtonian_jerk = -binary.G * binary.mass * (velocity - 3 * radial_rate * position) / radius**3
    quadrupole = 2 * binary.mu * _stf(
        _outer(velocity, velocity) + _outer(position, newtonian),
    )
    leading_news = 2 * binary.mu * _stf(
        3 * _outer(velocity, newtonian) + _outer(position, newtonian_jerk),
    )
    derivative_correction = 4 * binary.mu * _stf(_outer(velocity, reaction))
    return quadrupole, leading_news, leading_news + derivative_correction


def _mode_projection(tensor, binary):
    """Project the STF tensor using the Section 2 polarization with H = R h."""
    normalization = binary.G / binary.c**4 * np.sqrt(np.pi / 5)
    mode22 = normalization * (
        tensor[..., 0, 0] - tensor[..., 1, 1] - 2j * tensor[..., 0, 1]
    ) * 2
    mode21 = -4 * normalization * (tensor[..., 0, 2] - 1j * tensor[..., 1, 2])
    mode20 = np.sqrt(24) * normalization * tensor[..., 2, 2]
    return {
        (2, -2): np.conj(mode22),
        (2, -1): -np.conj(mode21),
        (2, 0): mode20,
        (2, 1): mode21,
        (2, 2): mode22,
    }


@dataclass(frozen=True)
class Radiation:
    retarded_times: np.ndarray
    frame: str | BondiFrame
    radiative_mass_quadrupole: np.ndarray
    radiative_mass_quadrupole_derivative: np.ndarray
    strain: dict
    news: dict
    energy_flux_infinity: np.ndarray
    angular_momentum_flux_infinity: np.ndarray
    strain_normalization: str = "H = R h"
    fixed_history_amplitude_pn: float = 0.0
    fixed_history_luminosity_relative_pn: float = 0.0
    fixed_history_quadrupole_pn: float = 0.0
    leading_quadrupole_tail: bool = False
    native_strain: dict | None = None
    native_news: dict | None = None
    energy_flux_from_modes: np.ndarray | None = None
    energy_flux_parts: dict | None = None
    retained_memory_strain: dict | None = None
    retained_memory_news: dict | None = None
    incoming_asymptotic_order: int | None = None
    incoming_quadrature_nodes: int | None = None
    native_recoil_quadrupole: dict | None = None
    radiative_tail_split_time: float | None = None

    @property
    def native_retarded_times(self):
        """Times of the native STF and luminosity sampled beside the output modes."""
        return self.retarded_times

    @property
    def h22(self):
        """Native (2,2) mode in the strain normalization recorded above."""
        series = self.native_strain if self.native_strain is not None else self.strain
        return series[2, 2]

    @property
    def news22(self):
        series = self.native_news if self.native_news is not None else self.news
        return series[2, 2]


def compute_radiation(trajectory, output_times, frame="native", accuracy=None):
    """Return H = R h and news on the requested cuts, plus native total power."""
    if not isinstance(frame, BondiFrame) and frame != "native":
        raise NotImplementedError("Select frame='native' or supply a BondiFrame.")
    if accuracy is not None:
        if not isinstance(accuracy, Accuracy) or accuracy.dynamics != trajectory.accuracy.dynamics:
            raise ValueError("Radiation accuracy must match the local trajectory model.")
    times = np.array(output_times,dtype=float,copy=True)
    if times.size==0 or not (np.isfinite(times)&(times>=trajectory.time_span[0])
                             &(times<=trajectory.time_span[1])).all():
        raise ValueError("Radiation cuts must lie within the covered trajectory.")
    split = None
    if trajectory.accuracy.radiative_tail_split_time is not None:
        from .consistent_tail import RadiativeLogHistory

        split = RadiativeLogHistory(trajectory._canonical_history,trajectory.accuracy.radiative_tail_split_time)
    memory = None
    if (accuracy or trajectory.accuracy).retain_memory_excess:
        from .waveform_history import build_retained_memory

        if trajectory.accuracy.canonical_grade != 7:
            raise ValueError("Retained higher memory requires a trajectory with the 3.5PN canonical history.")
        memory = build_retained_memory(
            trajectory.binary,trajectory._canonical_history,trajectory.accuracy.radiative_scale,
            log_history=split,max_time=float(times.max()) if frame=="native" else None,
        )
    # Both native cuts and BMS-resampled cuts share the same ordered past.
    prefix_cache = {}
    native = _native_radiation(trajectory,times,memory,prefix_cache,split)
    if isinstance(frame, BondiFrame):
        def sample_native(sample_times):
            result = _native_radiation(trajectory,sample_times,memory,prefix_cache,split)
            return result.strain, result.news

        strain, news = retime_and_project(frame, times, trajectory.binary.c, sample_native, native.strain)
        return replace(
            native, frame=frame, strain=strain, news=news,
            native_strain=native.strain, native_news=native.news,
        )
    return native


def _native_radiation(trajectory,times,memory=None,prefix_cache=None,split=None):
    position, velocity = trajectory.state(times)
    if trajectory.accuracy.canonical_grade is not None:
        return _native_canonical(trajectory,times,position,velocity,memory,prefix_cache,split)
    quadrupole, leading_news, derivative = _quadrupole_jets(trajectory, times, position, velocity)
    binary = trajectory.binary
    # Use the formally leading news in the leading flux contractions.
    power = binary.G / (5 * binary.c**5) * np.sum(leading_news**2, axis=(-2, -1))
    torque = 2 * binary.G / (5 * binary.c**5) * np.sum(
        np.cross(quadrupole, leading_news), axis=-2,
    )
    has_tail = trajectory.history is not None
    if has_tail:
        tail = np.stack([
            radiative_tail(binary, time, trajectory.history, trajectory.accuracy.radiative_scale)
            for time in times.flat
        ]).reshape(times.shape + (2, 3, 3))
        if trajectory.accuracy.waveform_pn == 1.5:
            modes, news, quad, quad_news = waveform_15(trajectory, position, velocity, tail)
            return Radiation(
                times, "native", quad, quad_news, modes, news, power, torque,
                fixed_history_amplitude_pn=1.5, fixed_history_quadrupole_pn=1.5, leading_quadrupole_tail=True,
            )
        quadrupole = quadrupole + tail[..., 0, :, :]
        derivative = derivative + tail[..., 1, :, :]
    return Radiation(
        times, "native", quadrupole, derivative, _mode_projection(quadrupole, binary),
        _mode_projection(derivative, binary), power, torque, leading_quadrupole_tail=has_tail,
    )


def _native_canonical(trajectory,times,position,velocity,memory,prefix_cache,split):
    from .luminosity import _flux_grade_eight, angular_momentum_flux
    from .multipoles import additional_local_quadrupole_grade_eight
    from .tail_derivative import fourfold_tail_jet
    from .waveform_history import finite_radiation_available, finite_radiation_grade_eight

    binary = trajectory.binary
    history = trajectory._canonical_history
    maximum_grade = trajectory.accuracy.canonical_grade
    upgraded = trajectory.accuracy.waveform_pn == 4.0
    flux_grade = 8 if upgraded else maximum_grade
    radiative_scale = trajectory.accuracy.radiative_scale
    samples = []
    channels_by_cut = []
    tail_lag = trajectory.accuracy.causal_tail_split_time
    for time, x, v in zip(times.flat, position.reshape(-1, 3), velocity.reshape(-1, 3)):
        logs = split if split is not None and split.active(time) else None
        if upgraded:
            tail_jet = None
            if tail_lag is not None:
                radius = np.linalg.norm(x)
                tail_jet = fourfold_tail_jet(
                    trajectory.history, time, 2*radius/binary.c, x@v/radius**2, tail_lag,
                )
            local = additional_local_quadrupole_grade_eight(binary, x, v, history.harmonic_scale)[2:4]
            strain, news, tensors, canonical, quadrupole_parts = finite_radiation_grade_eight(
                binary, history, time, x, v, trajectory.accuracy.harmonic_scales, radiative_scale, local,
                prefix_cache=prefix_cache, dynamical_tail_jet=tail_jet, log_history=logs,
            )
            parts = _flux_grade_eight(
                binary, history, time, x, v, tensors, canonical, quadrupole_parts, dynamical_tail_jet=tail_jet,
                split_jet=None if logs is None else logs.leading,
            )
            power = sum(sum(parts[key].values())
                        for key in ("canonical_local", "canonical_nonlocal_retarded", "propagation"))
            channels_by_cut.append(parts)
        else:
            strain, news, tensors, canonical = finite_radiation_available(
                binary,history,time,x,v,radiative_scale,log_history=logs,prefix_cache=prefix_cache,
            )
            power = _canonical_power_through_seven(binary, history, x, v, tensors, canonical, maximum_grade)
        sky_power = sum(mode_luminosity(binary, news, flux_grade).values())
        quad, quad_news = tensors["mass", 2]
        torque = sum(angular_momentum_flux(binary, tensors, flux_grade).values())
        samples.append((strain, news, sum(quad.values()), sum(quad_news.values()), power, sky_power, torque))

    channels = ({}, {})
    for index, channel in enumerate(channels):
        keys = {key[1:] for key in samples[0][index]}
        for ell, m in keys:
            values = [sum(value for key, value in sample[index].items() if key[1:] == (ell, m))
                      for sample in samples]
            channel[ell, m] = np.array(values).reshape(times.shape)
    values = [np.stack([sample[index] for sample in samples]) for index in range(2, 7)]
    mass_quadrupole, quadrupole_news, power, sky, torque = (
        value.reshape(times.shape + value.shape[1:]) for value in values
    )
    retained = (None, None)
    if memory is not None:
        increments = [memory.value(t) for t in times.flat]
        retained = tuple(
            {key: np.array([value[index][key] for value in increments]).reshape(times.shape)
             for key in increments[0][index]}
            for index in range(2)
        )
        for channel, graded in zip(channels, retained):
            for (_, ell, m), value in graded.items():
                channel[ell, m] = channel.get((ell, m), 0) + value
    parts = None
    if upgraded:
        parts = {
            name: {grade: np.array([cut[name][grade] for cut in channels_by_cut]).reshape(times.shape)
                   for grade in values}
            for name, values in channels_by_cut[0].items()
        }
    from .recoil_radiation import recoil_quadrupole_terms

    leading_recoil,_,center = trajectory.recoil(times)
    recoil = [
        recoil_quadrupole_terms(binary,x,v,momentum,displacement)
        for x,v,momentum,displacement in zip(*(value.reshape(-1,3)
                                              for value in (position,velocity,leading_recoil,center)))
    ]
    native_recoil = {
        key:np.stack([part[key] for part in recoil]).reshape(times.shape+(2,3,3))
        for key in recoil[0]
    }
    return Radiation(
        times, "native", mass_quadrupole, quadrupole_news, *channels, power, torque,
        fixed_history_amplitude_pn=maximum_grade/2, fixed_history_luminosity_relative_pn=flux_grade/2,
        fixed_history_quadrupole_pn=flux_grade/2,
        leading_quadrupole_tail=True, energy_flux_from_modes=sky, energy_flux_parts=parts,
        retained_memory_strain=retained[0], retained_memory_news=retained[1],
        incoming_asymptotic_order=None if history.incoming is None else history.incoming.orbit.asymptotic_order,
        incoming_quadrature_nodes=None if history.incoming is None else history.incoming.order,
        native_recoil_quadrupole=native_recoil,
        radiative_tail_split_time=None if split is None else split.lag,
    )


def _canonical_power_through_seven(binary, history, x, v, tensors, canonical, maximum_grade):
    canonical_news, corrections = {}, {}
    for (kind, rank), (_, grades) in tensors.items():
        name = ("M" if kind == "mass" else "S")+str(rank)
        baseline = {grade: field[rank+1] for grade,field in canonical[name].items() if grade in grades}
        canonical_news[kind, rank] = baseline
        corrections[kind, rank] = {grade: value-baseline.get(grade,0) for grade,value in grades.items()}
    local = canonical_local_luminosity(binary,x,v,history.harmonic_scale,history.infrared_scale)
    nonlinear = nonlinear_luminosity(binary,canonical_news,corrections,max_grade=maximum_grade)
    return sum(value for grade,value in local.items() if grade <= maximum_grade)+sum(nonlinear.values())
