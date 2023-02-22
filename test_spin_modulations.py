#############################################################################
##
##      Filename: test_spin_modulations.py
##
##      Author: Vijay Varma
##
##      Created: 22-02-2023
##
##      Description: Check if there is a relation between spin modulations in
##                   NR and orbital frequency.
##
#############################################################################

import numpy as np
import matplotlib.pyplot as P
from scipy.signal import find_peaks

from test_vs_lal import load_NR, generate_lal_pn

#-------------------------------------------------------------------------
#-------------------------------------------------------------------------
if __name__ == "__main__":


    base_dir = 'data_HybTest'
    case = '0011'
    Res = 'HiRes'
    NRDir =  f'{base_dir}/{Res}/NR'


    t_nr, h_nr, q_nr, chiA_nr, chiB_nr, omega_nr, orbphase_nr \
            = load_NR(case, NRDir)

    omega_pn, orbphase_pn, chiA_pn, chiB_pn, LNhat_pn, t_pn, h_pn \
            = generate_lal_pn(q_nr, chiA_nr[0], chiB_nr[0], omega_nr[0],
                          omega_nr[0], 0, t_nr[0])

    # Title tag for plots
    chiA_tag = '$\chi_A$=['
    chiB_tag = '$\chi_B$=['
    for idx in range(3):
        chiA_tag += f"{chiA_nr[0][idx]:.2f}"
        chiB_tag += f"{chiB_nr[0][idx]:.2f}"
        if idx == 2:
            chiA_tag += ']'
            chiB_tag += ']'
        else:
            chiA_tag += ', '
            chiB_tag += ', '

    title_tag = f"{base_dir}/{case} $q={q_nr:.2f}$ {chiA_tag} {chiB_tag} at " \
                f"t={t_nr[0]:.2f}"


    # Truncate nr and pn to common times. Do this for other quantities if you
    # use them.
    tmax = -50
    keep_nr = t_nr <= tmax
    t_nr = t_nr[keep_nr]
    chiA_nr = chiA_nr[keep_nr]
    orbphase_nr = orbphase_nr[keep_nr]
    keep_pn = t_pn <= tmax
    t_pn = t_pn[keep_pn]
    chiA_pn = chiA_pn[keep_pn]

    P.figure(figsize=(20, 12))
    P.subplots_adjust(hspace=0)

    # Let's consider x-components of the spins
    ax = P.subplot(3,1,1, aspect='auto')
    ax.plot(t_nr, chiA_nr.T[0], label='nr', lw=2)
    ax.plot(t_pn, chiA_pn.T[0], label='lal pn', ls='--')
    ax.legend(loc='best')
    ax.set_ylabel('chiAx')

    # Use the difference in NR and PN spins to find the period of the NR spin
    # oscillations. This works because PN seems to be missing the oscillations.
    ax = P.subplot(3,1,2, aspect='auto')
    residual = chiA_nr.T[0] - chiA_pn.T[0]
    peak_idx = find_peaks(residual)[0]
    ax.plot(t_nr, residual, color='C4', label='chiAx_nr - chiAx_pn')
    ax.plot(t_nr[peak_idx], residual[peak_idx], marker='o', lw=0, color='C3',
            label='peak locations')
    ax.set_ylabel('Spin residual')
    ax.legend(loc='best')

    # Find the periods between orbital phase changes of 2pi.
    # I'm sure there's a more efficient way, but I'm just finding times closest
    # to multiples of 2pi with a dumb for loop.
    peak_times_diff = np.diff(t_nr[peak_idx])
    orbit_steps = np.arange(orbphase_nr[0], orbphase_nr[-1], 2*np.pi)
    orbit_times = []
    for i in range(len(orbit_steps)):
        orbit_times.append(t_nr[np.argmin(np.abs(orbphase_nr-orbit_steps[i]))])
    # Time difference between orbphase steps of 2pi
    orbit_times_diff = np.diff(orbit_times)

    ax = P.subplot(3,1,3, aspect='auto')
    ax.plot(t_nr[peak_idx][1:], peak_times_diff, marker='o', lw=0,
            label='Betwen peaks of spin residual (markers above)', color='C2')
    ax.plot(orbit_times[1:], orbit_times_diff,
            label='Between orbphase change of 2pi', color='royalblue')
    # Dividing by two seems to give a great match.
    ax.plot(orbit_times[1:], orbit_times_diff/2,
            label='Between orbphase change of pi', color='tomato', lw=2)
    ax.set_ylabel('Time difference [M]')
    ax.legend(loc='best')

    P.suptitle(title_tag, y=0.91)
    P.savefig('test_spin_period.png', bbox_inches='tight')
    P.close()
