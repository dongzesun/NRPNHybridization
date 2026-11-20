"""Causal far-domain integration by parts of the standard harmonic tail."""

from dataclasses import dataclass
from functools import cache, cached_property

import numpy as np

from .history import integrate_log_panel
from .tails import _stf, primitive_tail_panel, quadrupole_tail_sources, sixth_quadrupole


@cache
def _grid(count):
    nodes, weights = np.polynomial.legendre.leggauss(count)
    return (1+nodes)/2, weights/2


def parts_jet(binary,state):
    """Return Q4,Q5,Q7,Q3,Q6 in the same source gauge."""
    position,velocity = state[:3],state[3:6]
    radius = np.linalg.norm(position)
    rate = position@velocity/radius**2
    g = binary.G*binary.mass
    acceleration = -g*position/radius**3
    jerk = -g*(velocity-3*rate*position)/radius**3
    third = 2*binary.mu*(3*_stf(velocity,acceleration)+_stf(position,jerk))
    return np.stack((*quadrupole_tail_sources(binary,state[:6]),third,sixth_quadrupole(binary,state)))


def inverse_four_panel(panel,time,stop,source_index):
    return _inverse_panel(panel,time,stop,source_index,4)


def _inverse_five_panel(panel,time,stop,source_index):
    return _inverse_panel(panel,time,stop,source_index,5)


def _inverse_panel(panel,time,stop,source_index,power):
    """Exact inverse-power moment of a quadratic panel using y=1/(u-t)."""
    a,b = panel.start,min(stop,panel.end)
    y_start,y_stop = 1/(time-a),1/(time-b)
    nodes,weights = _grid(2)
    width = y_stop-y_start
    angular = y_start+width*nodes
    # Evaluate y² times the source without subtracting nearly equal times.
    period = panel.end-panel.start
    unit = (time-a)*(angular-y_start)/period
    basis = np.stack((angular**2,unit*angular,unit**2))
    coefficients = panel.coefficients[:,source_index]
    return np.tensordot(width*(basis[:len(coefficients)]@(weights*angular**(power-4))),coefficients,axes=1)


def near_log_panel(panel,time,start,scale,*,differentiate=False):
    a,b = max(start,panel.start),min(time,panel.end)
    if b<=a:
        return np.zeros((3,3))
    original_width = panel.end-panel.start
    coefficients = np.zeros((3,3,3))
    values = panel.coefficients[:,2]
    if differentiate:
        values = np.polynomial.polynomial.polyder(values,axis=0)/original_width
    coefficients[:len(values)] = values
    offset = (a-panel.start)/original_width
    fraction = (b-a)/original_width
    c0,c1,c2 = coefficients
    local = np.stack((c0+offset*(c1+offset*c2),fraction*(c1+2*offset*c2),fraction**2*c2))
    return integrate_log_panel(a,b,local,time,scale,1)


def near_moments(panel,time,start):
    """Return the unscaled log and its exact constant-kernel coefficient."""
    a,b = max(start,panel.start),min(time,panel.end)
    if b<=a:
        return np.zeros((2,3,3))
    original_width = panel.end-panel.start
    x,y = (a-panel.start)/original_width,(b-panel.start)/original_width
    polynomial = panel.coefficients[:,2]
    orders = np.arange(1,len(polynomial)+1)
    coefficient = np.tensordot(original_width*(y**orders-x**orders)/orders,polynomial,axes=1)
    return np.stack((near_log_panel(panel,time,start,1.),coefficient))


@dataclass(frozen=True)
class PartsIncomingTail:
    """Incoming moments supplying the same SH phase-space join as the force."""
    base: object
    orbit: object
    binary: object

    source_shape = (5,3,3)

    @property
    def cutoff(self):
        return self.orbit.cutoff

    def at_cutoff(self):
        return self.value(self.cutoff)

    def value(self,time):
        return parts_jet(self.binary,np.concatenate(self.orbit.state_at(self.binary,time)))

    def integral(self,time,scale,power=1):
        ordinary = self.base.integral(time,scale,power)
        return np.concatenate((ordinary,np.zeros((2,3,3))))

    @cached_property
    def _remote_source(self):
        nodes,weights = _grid(48)
        distance = self.orbit.time_origin-self.cutoff
        elapsed = distance*nodes/(1-nodes)
        data = np.stack([self.value(self.cutoff-delay)[3] for delay in elapsed])
        return elapsed, weights*distance/(1-nodes)**2, data

    def inverse_four(self,time,stop):
        return self.inverse_power(time,stop,4)

    def inverse_power(self,time,stop,power):
        nearest = min(stop,self.cutoff)
        delay = time-nearest
        natural_distance = self.orbit.time_origin-self.cutoff
        if nearest==self.cutoff and delay>=natural_distance/4:
            offsets,weights,source = self._remote_source
            return np.tensordot(weights/(delay+offsets)**power,source,axes=1)
        nodes,weights = _grid(24)
        moments = np.stack([self.value(time-delay/node)[3] for node in nodes])
        return np.tensordot(weights*nodes**(power-2)/delay**(power-1),moments,axes=1)

    def near_moments(self,time,start):
        # Integrate after the substitution t=time-(time-cutoff)-d z².
        # The weight z renders the possible right-endpoint log integrable.
        width = self.cutoff-start
        nodes,weights = _grid(24)
        shifts = width*nodes**2
        source = np.stack([self.value(self.cutoff-lag)[2] for lag in shifts])
        mass = 2*width*weights*nodes
        factor = mass*np.log((time-self.cutoff)+shifts)
        return np.tensordot(np.stack((factor,mass)),source,axes=1)


def _parts_moments(history,time,lag,trial):
    split = time-lag
    remote = history.incoming
    if split<=history.start and remote is None:
        return np.stack([history.integral(time,1.,power,trial=trial)[2] for power in (1,0)])
    if split>=history.start:
        if trial is not None and split>trial.start:
            fraction = (split-trial.start)/(trial.end-trial.start)
            boundary = np.polynomial.polynomial.polyval(fraction,trial.coefficients)
        elif split==history.start and remote is not None and not history.end>history.start:
            boundary = remote.at_cutoff()
        else:
            boundary = history.value(split)
    else:
        boundary = remote.value(split)
    result = np.stack((boundary[4]*np.log(lag)+boundary[1]/lag-boundary[0]/lag**2
                       +2*boundary[3]/lag**3,boundary[4]))
    if hasattr(history,"inverse_four"):
        far = history.inverse_four(time,split,3,inverse_four_panel)
        recent = history.intersecting(split,time)
    else:
        far = sum((inverse_four_panel(panel,time,split,3)
                   for panel in history.panels if panel.start<split),np.zeros((3,3)))
        recent = (panel for panel in history.panels if panel.start<time and panel.end>split)
    if remote is not None:
        far += remote.inverse_four(time,split)
        if split<remote.cutoff:
            result += remote.near_moments(time,split)
    if trial is not None and trial.start<split:
        far += inverse_four_panel(trial,time,split,3)
    result[0] -= 6*far
    result += sum((near_moments(panel,time,split) for panel in recent),np.zeros((2,3,3)))
    if trial is not None:
        result += near_moments(trial,time,split)
    return result


def fourfold_tail_integral(history,time,scale,lag,trial):
    """Use the original log at recent times and the inverse quartic on the earlier past."""
    split = time-lag
    if split>history.end or (time>history.end and split<=history.start and history.incoming is None):
        # The provisional interval contains the split, or its finite
        # start is the boundary. Its contributions must remain provisional.
        moments = _parts_moments(history,time,lag,trial)
        return moments[0]-np.log(scale)*moments[1]
    key = (time,lag)
    cache = history._parts_cache
    if key not in cache:
        cache[key] = _parts_moments(history,time,lag,None)
    moments = cache[key]
    result = moments[0]-np.log(scale)*moments[1]
    if trial is not None:
        result += near_log_panel(trial,time,split,scale)
    return result

