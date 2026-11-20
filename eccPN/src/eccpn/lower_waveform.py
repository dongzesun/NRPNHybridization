"""Complete mode assembly through absolute 1.5PN on a supplied SH history."""

import numpy as np

from .gauge_moments import canonical_source_jets
from .leading_modes import leading_mode_coefficients
from .projection import project_stf


def _add(output, contributions):
    for mode, field in contributions.items():
        output[mode] = output.get(mode, 0) + field


def waveform_15(trajectory, position, velocity, quadrupole_tail):
    """Return all native 1.5PN modes, analytic graded news and the STF pair."""
    binary = trajectory.binary
    sample_shape = position.shape[:-1]
    shape = sample_shape + (3, 3)
    harmonic_scale = trajectory.accuracy.harmonic_scales[0]
    # No infrared regulator enters through the first PN correction.
    inputs = canonical_source_jets(binary, position, velocity, harmonic_scale, 1.0)
    complete_quadrupole = inputs["I2"][0][2:4] + inputs["I2"][2][2:4] + np.moveaxis(quadrupole_tail, -3, 0)

    waveforms = [
        leading_mode_coefficients(binary, x, v, ell_max=5)
        for x, v in zip(position.reshape(-1, 3), velocity.reshape(-1, 3))
    ]
    result = ({}, {})
    for index, output in enumerate(result):
        for key in waveforms[0][index]:
            grade, ell, m = key
            if (ell == 2 and grade == 0) or grade > 3:
                continue
            field = np.array([waveform[index][key] for waveform in waveforms]).reshape(sample_shape)
            _add(output, {(ell, m): field})
        _add(output, project_stf(binary, 2, complete_quadrupole[index].reshape(shape)))
        _add(output, project_stf(binary, 3, inputs["I3"][2][3 + index]))
        _add(output, project_stf(binary, 2, inputs["J2"][2][2 + index], "current"))
    return *result, complete_quadrupole[0], complete_quadrupole[1]
