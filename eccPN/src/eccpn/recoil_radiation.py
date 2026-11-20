"""Native quadrupole terms from motion of the inherited center of mass."""

from math import factorial

import numpy as np

from .leading_modes import newtonian_position_series
from .multipoles import _stf
from .recoil import recoil_fluxes


def _product(left,right,operation):
    return np.stack([
        sum(operation(left[k],right[power-k]) for k in range(power+1))
        for power in range(len(left))
    ])


def _scale(scalar,tensor):
    return _product(scalar,tensor,lambda left,right:left*right)


def _compact_terms(binary,position,velocity,recoil,center):
    """Differentiate the compact 1PN source with a linear common worldline."""
    position_jet = newtonian_position_series(binary,position,velocity,6)
    speed_jet = position_jet[1:]*np.arange(1,7)[:,None]
    position_jet = position_jet[:-1]
    center_jet = np.zeros_like(position_jet)
    center_jet[0],center_jet[1] = -center/binary.mass,-recoil/binary.mass
    drift_jet = np.zeros_like(center_jet)
    drift_jet[0] = center_jet[1]
    radius_squared = _product(position_jet,position_jet,np.dot)
    inverse_radius = np.zeros_like(radius_squared)
    inverse_radius[0] = radius_squared[0]**-0.5
    for order in range(1,len(inverse_radius)):
        inverse_radius[order] = sum(
            (k/2-order)*radius_squared[k]*inverse_radius[order-k] for k in range(1,order+1)
        )/(order*radius_squared[0])
    potential = binary.G*binary.mass*inverse_radius
    vv = _product(speed_jet,speed_jet,np.dot)
    xx = _product(position_jet,position_jet,_stf)
    xcenter = _product(position_jet,center_jet,_stf)
    xv = _product(position_jet,speed_jet,_stf)
    xdot = _product(position_jet,drift_jet,_stf)
    center_v = _product(center_jet,speed_jet,_stf)
    vh = _product(speed_jet,drift_jet,np.dot)
    xc = _product(position_jet,center_jet,np.dot)
    xvh = _product(position_jet,drift_jet,np.dot)+_product(center_jet,speed_jet,np.dot)
    x_v = _product(position_jet,speed_jet,np.dot)
    density = 3*_scale(vh,xx)+_scale(2*vv-potential,xcenter)
    retardation = 2*(_scale(xc,xx)+_scale(radius_squared,xcenter))
    current = (_scale(xvh,xx)+2*_scale(x_v,xcenter)
               -0.4*(2*_scale(xc,xv)+_scale(radius_squared,center_v+xdot)))
    prefactor = -binary.mu*(binary.m1-binary.m2)/(binary.mass*binary.c**2)
    return prefactor*np.stack([
        factorial(order)*density[order]+factorial(order+2)*retardation[order+2]/14
        -20*factorial(order+1)*current[order+1]/21 for order in (2,3)
    ])


def recoil_quadrupole_terms(binary,position,velocity,recoil,center):
    """Return strain and news tensors linear in the leading inherited recoil.

    The acceleration term holds the inherited momentum fixed in its main
    contribution. The separate flux derivative restores its ordinary
    time dependence in the news. The compact source term keeps a constant
    recoil velocity through the local 1PN differentiation. These tensors
    carry their explicit prefactors and remain separately recoverable.
    """
    radius = np.linalg.norm(position)
    normal = position/radius
    radial = normal@velocity
    normal_rate = (velocity-radial*normal)/radius
    newtonian = -binary.G*binary.mass/radius**2*normal
    strength = binary.G*(binary.m1-binary.m2)/(binary.mass*radius**2*binary.c**2)

    def acceleration(momentum):
        return strength*(2*(velocity@momentum)*normal+(normal@momentum)*velocity)

    local = acceleration(recoil)
    jerk = strength*(2*(newtonian@recoil)*normal+2*(velocity@recoil)*normal_rate
                     +(normal_rate@recoil)*velocity+(normal@recoil)*newtonian)
    jerk -= 2*radial/radius*local
    held = 2*binary.mu*np.stack((_stf(position,local),3*_stf(velocity,local)+_stf(position,jerk)))
    flux_jerk = acceleration(recoil_fluxes(binary,position,velocity)[0])
    flux = 2*binary.mu*np.stack((np.zeros((3,3)),_stf(position,flux_jerk)))
    return {
        "held_momentum_force":held,
        "momentum_flux_derivative":flux,
        "linear_compact_source":_compact_terms(binary,position,velocity,recoil,center),
    }
