"""Time derivative of the quadrature used in the dynamical 4PN tail."""

import numpy as np

from .newtonian_moments import newtonian_source_jets
from .tail_parts import _grid, _inverse_five_panel, fourfold_tail_integral, near_log_panel, parts_jet


def _slope(panel,time):
    width = panel.end-panel.start
    polynomial = np.polynomial.polynomial.polyder(panel.coefficients,axis=0)
    return np.polynomial.polynomial.polyval((time-panel.start)/width,polynomial)/width


def _incoming_slope(incoming,time):
    binary = incoming.binary
    position,velocity,acceleration = incoming.orbit.state_at(binary,time,include_acceleration=True)
    newtonian = -binary.G*binary.mass*position/np.linalg.norm(position)**3
    nominal = newtonian_source_jets(binary,position,velocity,highest_derivative=8,highest_rank=2)["M2"][0]
    # Only the velocity derivative of the nominal Newtonian jet changes.
    # The source is polynomial in velocity, so its directional derivative
    # is evaluated without a finite difference in retarded time.
    step = 1e-24
    displaced = velocity+1j*step*(acceleration-newtonian)
    correction = np.imag(parts_jet(binary,np.r_[position,displaced]))/step
    return nominal[[5,6,8,4,7]]+correction


def _recent_derivative(panels,time,split,scale):
    result = np.zeros((3,3))
    previous = None
    for panel in panels:
        start = max(split,panel.start)
        result += near_log_panel(panel,time,split,scale,differentiate=True)
        if start>split:
            jump = panel.coefficients[0,2].copy()
            if previous is not None:
                jump -= previous.coefficients[:,2].sum(axis=0)
            result += jump*np.log((time-start)/scale)
        previous = panel
    return result


def fourfold_tail_jet(history,time,scale,log_scale_rate,lag):
    """Covered acceleration-tail logarithmic integral and its cut derivative."""
    if not history.start<time<=history.end:
        raise ValueError("The tail jet needs a covered cut after the start of the finite panels.")
    split = time-lag
    incoming = history.incoming
    value = fourfold_tail_integral(history,time,scale,lag,None)
    if incoming is None:
        panels = (history.intersecting(history.start,time) if hasattr(history,"intersecting") else
                  (panel for panel in history.panels if panel.start < time))
        # Translation of a piecewise polynomial gives its interior slopes,
        # the onset and each right-minus-left jump at a panel boundary.
        derivative = history.value(history.start)[2] * np.log((time-history.start)/scale)
        derivative += _recent_derivative(panels,time,history.start,scale)
        return value,derivative-log_scale_rate*history._parts_cache[time][1]
    recent = (history.intersecting(split,time) if hasattr(history,"intersecting") else
              (panel for panel in history.panels if panel.end>split and panel.start<time))
    derivative = _recent_derivative(recent,time,split,scale)

    if split>=history.start:
        boundary = history.value(split)
        panel = (history.covering_panel(split) if hasattr(history,"covering_panel") else
                 next(panel for panel in history.panels if panel.start<=split<=panel.end))
        slope = _slope(panel,split)
    else:
        boundary = incoming.value(split)
        slope = _incoming_slope(incoming,split)
        width = incoming.cutoff-split
        nodes,weights = _grid(24)
        squared = nodes**2
        sample_times = incoming.cutoff-width*squared
        samples = np.stack([incoming.value(t)[2] for t in sample_times])
        slopes = np.stack([_incoming_slope(incoming,t)[2] for t in sample_times])
        delays = time-incoming.cutoff+width*squared
        logarithms = np.log(delays/scale)
        derivative += np.tensordot(
            2*weights*nodes*(width*(1-squared)/delays-logarithms),samples,axes=1,
        )
        derivative += np.tensordot(2*width*weights*nodes*squared*logarithms,slopes,axes=1)
        # At a moving remote split the integration nodes have fixed
        # lags. Differentiate that quadrature with the same source grid.
        source_times = time-lag/nodes
        far_slopes = np.stack([_incoming_slope(incoming,t)[3] for t in source_times])
        derivative -= 6*np.tensordot(weights*nodes**2/lag**3,far_slopes,axes=1)
    if split>=history.start:
        if hasattr(history,"inverse_power"):
            far = history.inverse_power(time,split,3,_inverse_five_panel,5)
        else:
            far = sum((_inverse_five_panel(panel,time,split,3)
                       for panel in history.panels if panel.start<split),np.zeros((3,3)))
        derivative += 24*far-6*boundary[3]/lag**4
        if incoming is not None:
            delay = time-incoming.cutoff
            if delay>=(incoming.orbit.time_origin-incoming.cutoff)/4:
                derivative += 24*incoming.inverse_power(time,split,5)
            else:
                nodes,weights = _grid(24)
                source_times = time-delay/nodes
                data = np.stack([incoming.value(t)[3] for t in source_times])
                slopes = np.stack([_incoming_slope(incoming,t)[3] for t in source_times])
                terms = (1-1/nodes)[:,None,None]*slopes-3*data/delay
                derivative -= 6*np.tensordot(weights*nodes**2/delay**3,terms,axes=1)
    derivative += (slope[4]*np.log(lag/scale)+slope[1]/lag-slope[0]/lag**2+2*slope[3]/lag**3)
    moments = history._parts_cache[(time,lag)]
    return value,derivative-log_scale_rate*moments[1]
