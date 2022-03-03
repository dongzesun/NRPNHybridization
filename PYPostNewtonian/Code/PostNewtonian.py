import scri
from PYPostNewtonian.Code import PNEvolution
from PYPostNewtonian.Code import PNWaveformModes
import numpy as np
import quaternion
import quaternionic

def PNWaveform(m1,m2,omega_0,chi1_0,chi2_0,frame_0,t_0=0.0, t_PNStart=False, t_PNEnd=False, PNEvolutionOrder=3.5, PNWaveformModeOrder=3.5, TaylorTn=1, StepsPerOrbit=32, ForwardInTime=True, tol=1e-12, MinStep=1e-7):
    xHat=quaternionic.array([0.0,1.0,0.0,0.0])
    yHat=quaternionic.array([0.0,0.0,1.0,0.0])
    zHat=quaternionic.array([0.0,0.0,0.0,1.0])
    delta=(m1-m2)/(m1+m2)
    m1=(1+delta)/2.0
    m2=(1-delta)/2.0
    v_0=omega_0**(1/3)
    chi1Mag=quaternionic.array([0,chi1_0[0],chi1_0[1],chi1_0[2]]).norm
    chi2Mag=quaternionic.array([0,chi2_0[0],chi2_0[1],chi2_0[2]]).norm
    S_chi1_0=quaternionic.array([0.0,0.0,0.0,0.0])
    S_chi2_0=S_chi1_0
    if chi1Mag>1e-12:
        S_chi1_0=np.sqrt(chi1Mag)*np.sqrt(\
            -quaternionic.array([0,chi1_0[0],chi1_0[1],chi1_0[2]]).normalized*zHat).normalized
    if chi2Mag>1e-12:
        S_chi2_0=np.sqrt(chi2Mag)*np.sqrt(\
            -quaternionic.array([0,chi2_0[0],chi2_0[1],chi2_0[2]]).normalized*zHat).normalized
    rfrak_frame_0=np.log(frame_0).vec
    PN=PNEvolution.PNEv.Evolution(xHat, yHat, zHat, m1, m2, v_0,S_chi1_0, S_chi2_0,\
        0.0, 0.0, 0.0, 0.0,rfrak_frame_0[0], rfrak_frame_0[1], rfrak_frame_0[2], t_PNStart, t_PNEnd,\
        PNEvolutionOrder, TaylorTn, StepsPerOrbit, ForwardInTime, tol, MinStep)

    W_PN_corot=scri.WaveformModes()
    W_PN_corot.t=PN.t+t_0
    frame=np.empty(len(PN.t), dtype=quaternion.quaternion)
    for i in range (len(PN.t)):
        frame[i]=np.exp(quaternion.quaternion(0.0,PN.y[5,i],PN.y[6,i],PN.y[7,i]))
    W_PN_corot.frame=frame
    W_PN_corot.data=PNWaveformModes.Modes(xHat, yHat, zHat, m1, m2, v_0,S_chi1_0, S_chi2_0, 0.0, 0.0, 0.0, 0.0,\
        rfrak_frame_0[0], rfrak_frame_0[1], rfrak_frame_0[2], PN.y, PNWaveformModeOrder)
    ell_min, ell_max = 2, 8
    W_PN_corot.ells = ell_min, ell_max

    return W_PN_corot
