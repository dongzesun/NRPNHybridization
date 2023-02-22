import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as P
import lal
import lalsimulation as lalsim
import quaternion

from PYPostNewtonian.Code import PostNewtonian

#-------------------------------------------------------------------------
def lMax_from_num_modes(num_modes):
    """Assuming h is provided as an array of modes, with num_modes = len(h),
    and the modes are sorted as
    [(2, -2), (2, -1), (2, 0), (2, 1), (2, 2), (3, -3), (3, -2), ...],
    we take num_modes as input and return the lMax of the mode array h.
    """
    # num_modes given lMax, given by (lMax -1) * (lMax + 3)
    nModesVsLMax = [-3, 0, 5, 12, 21, 32, 45, 60, 77]
    # Fancy way to invert the relation. Jon comes up with these stuff, don't
    # ask me.
    lMax = nModesVsLMax.index(num_modes)
    return lMax

# ----------------------------------------------------------------------------
def mode_dict_from_list(h_modes):
    """ Given a list of modes, h_modes return a dict of modes, h_dict.
        h_modes should have order: [(2,-2),(2,-1),(2,0),(2,1),(2,2),(3,-3),...]
        Usage: h_31 = h_dict['h_l3m1']
    """

    lmax = lMax_from_num_modes(len(h_modes))

    h_dict = {}
    idx = 0
    for ell in range(2, lmax+1):
        for emm in range(-ell, ell+1):
            h_dict[f'h_l{ell}m{emm}'] = h_modes[idx]
            idx += 1

    if idx != len(h_modes):
        raise Exception('Mismatch between length of h_modes and lmax.')

    return h_dict


#-----------------------------------------------------------------------------
def generate_lal_pn(q, chiA0, chiB0, omega_ref, omega_start, omega_end, t_ref,
                **kwargs):
    """
    Wrapper for PN waveform, spin and dynamics evolution in LAL.
    Evolves both forward and backwards in time.
    Inputs:
        - q:              Mass ratio (q>=1)
        - chiA0:          Dimensionless spin of BhA at initial freq.
        - chiB0:          Dimensionless spin of BhB at initial freq.
        - omega_ref:      Reference orbital frequency in dimensionless units.
        - omega_start:    Initial orbital frequency in dimensionless units at
                          which to terminate the backward evolution.
                          Set omega_start=omega_ref for only forward evolution.
        - omega_end:      Final orbital frequency in dimensionless units at
                          which to terminate the forward evolution.
                          Set omega_end=omega_ref for only backward evolution.
                          Set omega_end=0 to evolve until ISCO (I think).
        - t_ref:          Time to associate with omega_ref.
        - kwargs:         See kwargs of generate_lal_pn_one_way()

    Outputs (all are time series):
        - Omega:          Dimensionless orbital frequency.
        - Phi:            Orbital phase (radians).
        - ChiA:           Dimensionless spin of BhA.
        - ChiB:           Dimensionless spin of BhB.
        - LNhat:          Orbital angular momentum direction.
        - t:              Time array.
        - hlm_dict:       Waveform modes dictionary.
                          To get a specific mode, do h22 = hlm_dict['h_l2m2']

    The frame is defined at the reference frequency omega_ref, as follows: \n
        - z-axis is set by the orbital angular momentum direction.
        - x-axis is the separation vector from the lighter BH to the heavier BH.
        - y-axis completes the triad by right-hand rule. \n
        All quantities are defined in this fixed frame, including the waveform,
        initial spins, returned spins, other vectors like LNhat, etc.
    """

    if omega_start > omega_ref:
        raise Exception('omega_start cannot be larger than omega_ref')
    if omega_end != 0 and omega_end < omega_ref:
        raise Exception('omega_end cannot be smaller than omega_ref, unless '
                        'omega_end=0')


    # Do the backwards evolution if necessary
    data_backwards = None
    if omega_start < omega_ref:
        data_backwards = generate_lal_pn_one_way(q, chiA0, chiB0, omega_ref,
                                                 omega_start, t_ref, **kwargs)

    # Do the forwards evolution if necessary
    data_forwards = None
    if omega_end == 0 or omega_end > omega_ref:
        data_forwards = generate_lal_pn_one_way(q, chiA0, chiB0, omega_ref,
                                                omega_end, t_ref, **kwargs)

    if data_backwards == None and data_forwards == None:
        raise Exception("Both data_forwards and data_backwards are None!")

    if data_backwards == None:
        return data_forwards
    elif data_forwards == None:
        return data_backwards
    else:
        omega_backwards, phi_backwards, chiA_backwards, chiB_backwards, \
            LNhat_backwards, t_backwards, hlm_dict_backwards = data_backwards
        omega_forwards, phi_forwards, chiA_forwards, chiB_forwards, \
            LNhat_forwards, t_forwards, hlm_dict_forwards = data_forwards

        ## Combine backwards and forward evolutions
        if not np.allclose([omega_backwards[-1],
                            phi_backwards[-1],
                            t_backwards[-1]],
                            [omega_forwards[0],
                            phi_forwards[0],
                            t_forwards[0]]) \
               or not np.allclose(chiA_backwards[-1], chiA_forwards[0]) \
               or not np.allclose(chiB_backwards[-1], chiB_forwards[0]) \
               or not np.allclose(LNhat_backwards[-1], LNhat_forwards[0]):
            raise Exception("Expected params at attachment point to match.")
        else:
            omega_full = np.concatenate((omega_backwards, omega_forwards[1:]))
            phi_full = np.concatenate((phi_backwards, phi_forwards[1:]))
            chiA_full = np.concatenate((chiA_backwards, chiA_forwards[1:]))
            chiB_full = np.concatenate((chiB_backwards, chiB_forwards[1:]))
            LNhat_full = np.concatenate((LNhat_backwards, LNhat_forwards[1:]))
            t_full = np.concatenate((t_backwards, t_forwards[1:]))

        if not np.all(np.diff(t_full) > 0):
            raise Exception('Non monotonic time array.')

        hlm_dict_full = {}
        for key in hlm_dict_backwards:
            hlm_back = hlm_dict_backwards[key]
            hlm_forw = hlm_dict_forwards[key]
            if not np.allclose(hlm_back[-1], hlm_forw[0]):
                raise Exception("Expected WF at attachment point to match.")
            hlm_dict_full[key] = np.concatenate((hlm_back, hlm_forw[1:]))

        return omega_full, phi_full, chiA_full, chiB_full, LNhat_full, \
               t_full, hlm_dict_full


#-----------------------------------------------------------------------------
def generate_lal_pn_one_way(q, chiA0, chiB0, omega_ref, omega_final, t_ref,
        approximant='SpinTaylorT4', dt=0.1, spinO=6, phaseO=7, lmax=4):
    """
    Wrapper for PN waveform, spin and dynamics evolution in LAL.
    Evolves only forward or backwards in time.
    Inputs:
        - q:              Mass ratio (q>=1)
        - chiA0:          Dimensionless spin of BhA at initial freq.
        - chiB0:          Dimensionless spin of BhB at initial freq.
        - omega_ref:      Reference orbital frequency in dimensionless units.
        - omega_final:    Final orbital frequency in dimensionless units.
                          If omega_final > omega_ref or omega_final = 0, the
                          evolution proceeds forwards in time, else backwards.
        - t_ref:          Time to associate with omega_ref.
        - approximant:    'SpinTaylorT1/T4/T5'.
        - dt:             Dimensionless step time for evolution.
        - spinO:          Twice PN order of spin effects.
        - phaseO:         Twice PN order in phase.
        - lmax:           Max ell.

    Outputs (all are time series):
        - Omega:          Dimensionless orbital frequency.
        - Phi:            Orbital phase (radians).
        - ChiA:           Dimensionless spin of BhA.
        - ChiB:           Dimensionless spin of BhB.
        - LNhat:          Orbital angular momentum direction.
        - t:              Time array.
        - hlm_dict:       Waveform modes dictionary.
                          To get a specific mode, do h22 = hlm_dict['h_l2m2']

    The frame is defined at the reference frequency omega_ref, as follows: \n
        - z-axis is set by the orbital angular momentum direction.
        - x-axis is the separation vector from the lighter BH to the heavier BH.
        - y-axis completes the triad by right-hand rule. \n
        All quantities are defined in this fixed frame, including the waveform,
        initial spins, returned spins, other vectors like LNhat, etc.
    """

    approxTag = lalsim.SimInspiralGetApproximantFromString(approximant)

    # These get scaled out in the end
    M = 100                         # Total mass in solar masses
    dist_SI = 1.0e6 * lal.PC_SI     # 1 Mpc in SI

    # time step and initial GW freq in SI units
    MT = M*lal.MTSUN_SI
    deltaT = dt*MT
    fStart = omega_ref/np.pi/MT

    # component masses of the binary
    m1_SI =  M*lal.MSUN_SI*q/(1.+q)
    m2_SI =  M*lal.MSUN_SI/(1.+q)

    # spins at fStart
    s1x, s1y, s1z = chiA0
    s2x, s2y, s2z = chiB0

    # Final frequency (can go backwards as well)
    fEnd = omega_final/np.pi/MT

    # initial value of orbital angular momentum unit vector, i.e at fStart
    lnhatx, lnhaty, lnhatz = 0,0,1

    # initial value of orbital plane basis vector, i.e at fStart
    e1x, e1y, e1z = 1, 0, 0

    # tidal deformability parameters
    lambda1, lambda2 = 0, 0

    # spin-induced quadrupole moments
    quadparam1, quadparam2 = 1, 1

    # twice PN order of tidal effects
    tideO = 0

    # include some known L-S terms
    lscorr = 1

    # evolve spins and collect data into a nice format. The input values start
    # with lower-case letters, while output values start with upper-case
    # letters.
    V, Phi, S1x, S1y, S1z, S2x, S2y, S2z, LNhatx, LNhaty, LNhatz, \
        E1x, E1y, E1z = lalsim.SimInspiralSpinTaylorPNEvolveOrbit(deltaT, \
        m1_SI, m2_SI, fStart, fEnd, s1x, s1y, s1z, s2x, s2y, s2z, \
        lnhatx, lnhaty, lnhatz, e1x, e1y, e1z, lambda1, lambda2, \
        quadparam1, quadparam2, spinO, tideO, phaseO, lscorr, approxTag)

    modearray=lalsim.SimInspiralCreateModeArray();
    for ell in range(2, lmax+1):
        modearray = lalsim.SimInspiralModeArrayActivateAllModesAtL(modearray, ell)

    ampO = 3    # 3 means 1.5 PN
    _, Hlms = lalsim.SimInspiralSpinTaylorHlmModesFromOrbit(V, Phi,
            LNhatx, LNhaty, LNhatz, E1x, E1y, E1z,
            S1x, S1y, S1z, S2x, S2y, S2z,
            m1_SI, m2_SI, dist_SI, ampO, modearray)

    # Convert to human format
    V = np.array(V.data.data)
    Phi = np.array(Phi.data.data)
    ChiA = np.array([S1x.data.data, S1y.data.data, S1z.data.data]).T
    ChiB = np.array([S2x.data.data, S2y.data.data, S2z.data.data]).T
    LNhat = np.array([LNhatx.data.data, LNhaty.data.data, LNhatz.data.data]).T
    E1 = np.array([E1x.data.data, E1y.data.data, E1z.data.data]).T
    Omega = V**3  # orbital frequency

    hlm_dict = {}
    for ell in range(2, lmax+1):
        for emm in range(-ell, ell+1):
            # Convert to human format, and make dimless
            # Flip sign to match convention of Dongze's code
            hlm_dict[f'h_l{ell}m{emm}'] = -lalsim.SphHarmTimeSeriesGetMode(Hlms, ell, emm).data.data * dist_SI/MT/lal.C_SI

    # Time array
    t = np.arange(len(Omega)) * dt
    if omega_final == 0 or omega_final > omega_ref:
        # In this case, the evolution goes forward, and we set t=t_ref at the
        # first index
        t +=  t_ref - t[0]
    else:
        # Else, the evolving goes backward, and we set t=t_ref at the last index
        t +=  t_ref - t[-1]

    return Omega, Phi, ChiA, ChiB, LNhat, t, hlm_dict


#-----------------------------------------------------------------------------
def generate_dongze_pn(q, chiA0, chiB0, omega0, t_ref, omega_start, omega_end):

    if omega_end == 0:
        omega_end = None

    # At this point w is in corotating frame
    w, chiA, chiB = PostNewtonian.PNWaveform(q, omega0,
            chiA0, chiB0, t_0=t_ref, omega_start=omega_start,
            omega_end=omega_end, return_chi=True)

    # Now in inertial frame
    w.to_inertial_frame()

    t = w.t
    h = w.data.T
    chiA = quaternion.as_float_array(chiA)[:,1:]
    chiB = quaternion.as_float_array(chiB)[:,1:]

    # Convert to dict format
    h_dict = mode_dict_from_list(h)

    return t, h_dict, chiA, chiB


# Params corresponding to q8_7d_0701
q = 2.8270343151377024
chiA0 = [-0.31006488,  0.0156336,  -0.51537769]
chiB0 = [0.15821173, -0.15421887, -0.25157237]
omega_ref = 0.013182673547263367
omega_start = omega_ref/1.2    # something smaller
omega_end = 0       # Go forward until PN dies
t_ref = 0


omega_pn, orbphase_pn, chiA_pn, chiB_pn, LNhat_pn, t_pn, h_pn \
        = generate_lal_pn(q, chiA0, chiB0, omega_ref, omega_start, omega_end,
                t_ref, approximant='SpinTaylorT4')

t_pn_d, h_pn_d, chiA_pn_d, chiB_pn_d  \
    = generate_dongze_pn(q, chiA0, chiB0, omega_ref, t_ref,
                         omega_start, omega_end)

# Title tag for plots
chiA_tag = '$\chi_A$=['
chiB_tag = '$\chi_B$=['
for idx in range(3):
    chiA_tag += f"{chiA0[idx]:.2f}"
    chiB_tag += f"{chiB0[idx]:.2f}"
    if idx == 2:
        chiA_tag += ']'
        chiB_tag += ']'
    else:
        chiA_tag += ', '
        chiB_tag += ', '

title_tag = f"$q={q:.2f}$ {chiA_tag} {chiB_tag} at t={t_ref:.2f}"

P.figure(figsize=(12, 12))
P.subplots_adjust(hspace=0)

# Plot the waveform
mode_keys_list = ['h_l2m2', 'h_l2m1']
for idx, mode_key in enumerate(mode_keys_list):
    ax = P.subplot(4,1,idx+1, aspect='auto')

    ax.plot(t_pn, np.real(h_pn[mode_key]), label='lal_pn')
    ax.plot(t_pn_d, np.real(h_pn_d[mode_key]), label='dongze_pn', ls='--')

    ax.set_ylabel(mode_key)
    ax.set_xlabel('t')
    if idx == 0:
        ax.legend(loc='best')

# Plot spins
ax = P.subplot(4,1,3, aspect='auto')
spin_idx = 0
ax.plot(t_pn, chiA_pn.T[spin_idx], label='lal_pn')
ax.plot(t_pn_d, chiA_pn_d.T[spin_idx], label='dongze_pn', ls='--')
ax.set_ylabel('chiA'+['x','y','z'][spin_idx])
ax.set_xlabel('t')
ax = P.subplot(4,1,4, aspect='auto')
ax.plot(t_pn, chiB_pn.T[spin_idx], label='lal_pn')
ax.plot(t_pn_d, chiB_pn_d.T[spin_idx], label='dongze_pn', ls='--')
ax.set_ylabel('chiB'+['x','y','z'][spin_idx])
ax.set_xlabel('t')

P.suptitle(title_tag, y=0.93)
P.savefig('test.png', bbox_inches='tight')
P.close()
