import scri
from PYPostNewtonian.Code import PNEvolution
from PYPostNewtonian.Code import PNWaveformModes
import numpy as np
import quaternionic
import quaternion

def PNWaveform(delta,omega_0,chi1_0,chi2_0,frame_0,t_0=0.0, t_PNStart=False, t_PNEnd=False, PNEvolutionOrder=3.5, PNWaveformModeOrder=3.5, TaylorTn=1, StepsPerOrbit=32, ForwardInTime=True, tol=1e-12, MinStep=1e-7):
    """
    delta = (m1-m2)/(m1+m2), float number,
    omega_0: magnititude of angular velocity at t_0, float number,
    chi1_0 and chi2_0: spin vectors at t_0, 3-d vectors,
    frame_0: the frame quaternion at t_0, quaternionic_array object,
    t_0: the corresponding time of the above given initial values, float number,
    t_PNStart: the start time of PN relative to t_0: t_PNStart=t_real_start-t_0, float number. If false, default is t_0-t_real_start=3(t_merger-t_0),
    t_PNEnd: the end time of PN relative to t_0: t_PNEnd=t_real_end-t_0, float number. If false, default is merger time,
    PNEvolutionOrder: float number in [0,0.5,1,1.5,2,2.5,3,3.5,4], default is 3.5,
    PNWaveformModeOrder: float number in [0,0.5,1,1.5,2,2.5,3,3.5,4], default is 3.5,
    TaylorTn: int number in [1,4,5], default is 1,
    StepsPerOrbit: float number,
    ForwardInTime: whether to evlove PN forward in time, bool number, default if True,
    tol: tolerance of the integrator, float number,
    MinStep: minimal time interval for the PN waveform, float number.
    """
    
    if not PNEvolutionOrder in [0,0.5,1,1.5,2,2.5,3,3.5,4]:
        message=("PNEvolutionOrder must be a float number in [0,0.5,1,1.5,2,2.5,3,3.5,4].")
        raise ValueError(message)
    if not PNWaveformModeOrder in [0,0.5,1,1.5,2,2.5,3,3.5,4]:
        message=("PNWaveformModeOrder must be a float number in [0,0.5,1,1.5,2,2.5,3,3.5,4].")
        raise ValueError(message)
    if not TaylorTn in [1,4,5]:
        message=("TaylorTn must be an int number in [1,4,5].")
        raise ValueError(message)          
        
    xHat=quaternionic.array([0.0,1.0,0.0,0.0])
    yHat=quaternionic.array([0.0,0.0,1.0,0.0])
    zHat=quaternionic.array([0.0,0.0,0.0,1.0])
    m1=(1+delta)/2.0
    m2=(1-delta)/2.0
    v_0=omega_0**(1/3)
    chi1Mag=np.sqrt(quaternionic.array([0,chi1_0[0],chi1_0[1],chi1_0[2]]).norm)
    chi2Mag=np.sqrt(quaternionic.array([0,chi2_0[0],chi2_0[1],chi2_0[2]]).norm)
    
    # Quaternions that rotate z-axis to spin vectors
    S_chi1_0=quaternionic.array([0.0,0.0,0.0,0.0])
    S_chi2_0=quaternionic.array([0.0,0.0,0.0,0.0])
    if chi1Mag>1e-12:
        S_chi1_0=np.sqrt(chi1Mag)*np.sqrt(\
            -quaternionic.array([0,chi1_0[0],chi1_0[1],chi1_0[2]]).normalized*zHat).normalized
    if chi2Mag>1e-12:
        S_chi2_0=np.sqrt(chi2Mag)*np.sqrt(\
            -quaternionic.array([0,chi2_0[0],chi2_0[1],chi2_0[2]]).normalized*zHat).normalized
        
    rfrak_frame_0=np.log(frame_0).vec # logarithm of frame quaternion
    PN=PNEvolution.PNEv.Evolution(xHat, yHat, zHat, m1, m2, v_0,S_chi1_0, S_chi2_0, rfrak_frame_0, t_PNStart, t_PNEnd,\
        PNEvolutionOrder, TaylorTn, StepsPerOrbit, ForwardInTime, tol, MinStep)# Evolve PN parameters, PN.t is PN time, PN.y=[v, chi1_x, chi1_y
        # chi2_x, chi2_y, rfrak_frame_x, rfrak_frame_y, rfrak_frame_z]

    W_PN_corot=scri.WaveformModes()
    W_PN_corot.t=PN.t+t_0
    W_PN_corot.frame=np.empty(len(PN.t), dtype=quaternion.quaternion)
    for i in range (len(PN.t)):
        W_PN_corot.frame[i]=np.exp(quaternion.quaternion(0.0,PN.y[5,i],PN.y[6,i],PN.y[7,i]))
    W_PN_corot.data, W_PN_corot.ells = PNWaveformModes.Modes(xHat, yHat, zHat, m1, m2, v_0,S_chi1_0, S_chi2_0, rfrak_frame_0, PN.y, PNWaveformModeOrder)

    return W_PN_corot
