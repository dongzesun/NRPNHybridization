# File produced automatically by PNCodeGen.ipynb
import numpy as np
from numpy import conjugate, dot, exp, log, sqrt, pi
from numpy import euler_gamma as EulerGamma
import quaternion

class WaveformModes_0PN :
    def WaveformModes_0PN(xHat_i, yHat_i, zHat_i, M1_i, M2_i, v_i, S_chi1_i, S_chi2_i, rfrak_chi1_x_i, rfrak_chi1_y_i, rfrak_chi2_x_i, rfrak_chi2_y_i, rfrak_frame_x_i, rfrak_frame_y_i, rfrak_frame_z_i,y) :            
        ModeData=np.empty((len(y[0]),77), dtype=complex)
        def WaveformModes(M1_i, M2_i, v_i):
            I=1j
            global rhOverM_coeff
            M1=M1_i
            M2=M2_i
            v=v_i
            M=M1 + M2
            nu=M1*M2/M**2
            rhOverM_coeff=8*sqrt(5)*nu*v**2*sqrt(pi)/5
            hHat_2_0_0=-5*sqrt(6)/84
            hHat_2_2_0=1
            hHat_4_0_0=-sqrt(2)/1008
            


            rhOverM_coeff = 8*sqrt(5)*nu*v**2*sqrt(pi)/5

            Modes=np.empty(77, dtype=complex)
            # (ell, m) = (2, +/- 0)
            Symm = hHat_2_0_0*rhOverM_coeff
            Asymm = 0
            Modes[2] = Symm + Asymm
            # (ell, m) = (2, +/- 1)
            Symm = 0
            Asymm = 0
            Modes[3] = Symm + Asymm
            Modes[1] = conjugate(Symm - Asymm)
            # (ell, m) = (2, +/- 2)
            Symm = hHat_2_2_0*rhOverM_coeff
            Asymm = 0
            Modes[4] = Symm + Asymm
            Modes[0] = conjugate(Symm - Asymm)
            # (ell, m) = (3, +/- 0)
            Symm = 0
            Asymm = 0
            Modes[8] = Symm + Asymm
            # (ell, m) = (3, +/- 1)
            Symm = 0
            Asymm = 0
            Modes[9] = Symm + Asymm
            Modes[7] = -conjugate(Symm - Asymm)
            # (ell, m) = (3, +/- 2)
            Symm = 0
            Asymm = 0
            Modes[10] = Symm + Asymm
            Modes[6] = -conjugate(Symm - Asymm)
            # (ell, m) = (3, +/- 3)
            Symm = 0
            Asymm = 0
            Modes[11] = Symm + Asymm
            Modes[5] = -conjugate(Symm - Asymm)
            # (ell, m) = (4, +/- 0)
            Symm = hHat_4_0_0*rhOverM_coeff
            Asymm = 0
            Modes[16] = Symm + Asymm
            # (ell, m) = (4, +/- 1)
            Symm = 0
            Asymm = 0
            Modes[17] = Symm + Asymm
            Modes[15] = conjugate(Symm - Asymm)
            # (ell, m) = (4, +/- 2)
            Symm = 0
            Asymm = 0
            Modes[18] = Symm + Asymm
            Modes[14] = conjugate(Symm - Asymm)
            # (ell, m) = (4, +/- 3)
            Symm = 0
            Asymm = 0
            Modes[19] = Symm + Asymm
            Modes[13] = conjugate(Symm - Asymm)
            # (ell, m) = (4, +/- 4)
            Symm = 0
            Asymm = 0
            Modes[20] = Symm + Asymm
            Modes[12] = conjugate(Symm - Asymm)
            # (ell, m) = (5, +/- 0)
            Symm = 0
            Asymm = 0
            Modes[26] = Symm + Asymm
            # (ell, m) = (5, +/- 1)
            Symm = 0
            Asymm = 0
            Modes[27] = Symm + Asymm
            Modes[25] = -conjugate(Symm - Asymm)
            # (ell, m) = (5, +/- 2)
            Symm = 0
            Asymm = 0
            Modes[28] = Symm + Asymm
            Modes[24] = -conjugate(Symm - Asymm)
            # (ell, m) = (5, +/- 3)
            Symm = 0
            Asymm = 0
            Modes[29] = Symm + Asymm
            Modes[23] = -conjugate(Symm - Asymm)
            # (ell, m) = (5, +/- 4)
            Symm = 0
            Asymm = 0
            Modes[30] = Symm + Asymm
            Modes[22] = -conjugate(Symm - Asymm)
            # (ell, m) = (5, +/- 5)
            Symm = 0
            Asymm = 0
            Modes[31] = Symm + Asymm
            Modes[21] = -conjugate(Symm - Asymm)
            # (ell, m) = (6, +/- 0)
            Symm = 0
            Asymm = 0
            Modes[38] = Symm + Asymm
            # (ell, m) = (6, +/- 1)
            Symm = 0
            Asymm = 0
            Modes[39] = Symm + Asymm
            Modes[37] = conjugate(Symm - Asymm)
            # (ell, m) = (6, +/- 2)
            Symm = 0
            Asymm = 0
            Modes[40] = Symm + Asymm
            Modes[36] = conjugate(Symm - Asymm)
            # (ell, m) = (6, +/- 3)
            Symm = 0
            Asymm = 0
            Modes[41] = Symm + Asymm
            Modes[35] = conjugate(Symm - Asymm)
            # (ell, m) = (6, +/- 4)
            Symm = 0
            Asymm = 0
            Modes[42] = Symm + Asymm
            Modes[34] = conjugate(Symm - Asymm)
            # (ell, m) = (6, +/- 5)
            Symm = 0
            Asymm = 0
            Modes[43] = Symm + Asymm
            Modes[33] = conjugate(Symm - Asymm)
            # (ell, m) = (6, +/- 6)
            Symm = 0
            Asymm = 0
            Modes[44] = Symm + Asymm
            Modes[32] = conjugate(Symm - Asymm)
            # (ell, m) = (7, +/- 0)
            Symm = 0
            Asymm = 0
            Modes[52] = Symm + Asymm
            # (ell, m) = (7, +/- 1)
            Symm = 0
            Asymm = 0
            Modes[53] = Symm + Asymm
            Modes[51] = -conjugate(Symm - Asymm)
            # (ell, m) = (7, +/- 2)
            Symm = 0
            Asymm = 0
            Modes[54] = Symm + Asymm
            Modes[50] = -conjugate(Symm - Asymm)
            # (ell, m) = (7, +/- 3)
            Symm = 0
            Asymm = 0
            Modes[55] = Symm + Asymm
            Modes[49] = -conjugate(Symm - Asymm)
            # (ell, m) = (7, +/- 4)
            Symm = 0
            Asymm = 0
            Modes[56] = Symm + Asymm
            Modes[48] = -conjugate(Symm - Asymm)
            # (ell, m) = (7, +/- 5)
            Symm = 0
            Asymm = 0
            Modes[57] = Symm + Asymm
            Modes[47] = -conjugate(Symm - Asymm)
            # (ell, m) = (7, +/- 6)
            Symm = 0
            Asymm = 0
            Modes[58] = Symm + Asymm
            Modes[46] = -conjugate(Symm - Asymm)
            # (ell, m) = (7, +/- 7)
            Symm = 0
            Asymm = 0
            Modes[59] = Symm + Asymm
            Modes[45] = -conjugate(Symm - Asymm)
            # (ell, m) = (8, +/- 0)
            Symm = 0
            Asymm = 0
            Modes[68] = Symm + Asymm
            # (ell, m) = (8, +/- 1)
            Symm = 0
            Asymm = 0
            Modes[69] = Symm + Asymm
            Modes[67] = conjugate(Symm - Asymm)
            # (ell, m) = (8, +/- 2)
            Symm = 0
            Asymm = 0
            Modes[70] = Symm + Asymm
            Modes[66] = conjugate(Symm - Asymm)
            # (ell, m) = (8, +/- 3)
            Symm = 0
            Asymm = 0
            Modes[71] = Symm + Asymm
            Modes[65] = conjugate(Symm - Asymm)
            # (ell, m) = (8, +/- 4)
            Symm = 0
            Asymm = 0
            Modes[72] = Symm + Asymm
            Modes[64] = conjugate(Symm - Asymm)
            # (ell, m) = (8, +/- 5)
            Symm = 0
            Asymm = 0
            Modes[73] = Symm + Asymm
            Modes[63] = conjugate(Symm - Asymm)
            # (ell, m) = (8, +/- 6)
            Symm = 0
            Asymm = 0
            Modes[74] = Symm + Asymm
            Modes[62] = conjugate(Symm - Asymm)
            # (ell, m) = (8, +/- 7)
            Symm = 0
            Asymm = 0
            Modes[75] = Symm + Asymm
            Modes[61] = conjugate(Symm - Asymm)
            # (ell, m) = (8, +/- 8)
            Symm = 0
            Asymm = 0
            Modes[76] = Symm + Asymm
            Modes[60] = conjugate(Symm - Asymm)

            return Modes
            
        for i in range(len(y[0])):
            v = y[0,i]
            rfrak_chi1_x = y[1,i]
            rfrak_chi1_y = y[2,i]
            rfrak_chi2_x = y[3,i]
            rfrak_chi2_y = y[4,i]
            rfrak_frame_x = y[5,i]
            rfrak_frame_y = y[6,i]
            rfrak_frame_z = y[7,i]
            Phi = y[8,i]
            xHat=xHat_i
            yHat=yHat_i
            zHat=zHat_i
            M1=M1_i
            M2=M2_i 
            S_chi1=S_chi1_i
            S_chi2=S_chi2_i
            M=M1 + M2
            delta=(M1 - M2)/M
            nu=M1*M2/M**2
            R=exp(rfrak_frame_x*xHat + rfrak_frame_y*yHat + rfrak_frame_z*zHat)
            nHat=R*xHat*conjugate(R)
            lambdaHat=R*yHat*conjugate(R)
            ellHat=R*zHat*conjugate(R)
            R_S1=exp(rfrak_chi1_x*xHat + rfrak_chi1_y*yHat)
            R_S2=exp(rfrak_chi2_x*xHat + rfrak_chi2_y*yHat)
            chiVec1=S_chi1*R_S1*zHat*conjugate(R_S1)*conjugate(S_chi1)
            chiVec2=S_chi2*R_S2*zHat*conjugate(R_S2)*conjugate(S_chi2)
            v_i=v
            rfrak_frame_x_i=rfrak_frame_x
            rfrak_frame_y_i=rfrak_frame_y
            rfrak_frame_z_i=rfrak_frame_z
            chiVec1_i=quaternion.as_float_array(chiVec1)
            chiVec2_i=quaternion.as_float_array(chiVec2)
            ModeData[i,:]=WaveformModes(M1_i, M2_i, v_i)
        return ModeData

class WaveformModes_0p50PN :
    def WaveformModes_0p50PN(xHat_i, yHat_i, zHat_i, M1_i, M2_i, v_i, S_chi1_i, S_chi2_i, rfrak_chi1_x_i, rfrak_chi1_y_i, rfrak_chi2_x_i, rfrak_chi2_y_i, rfrak_frame_x_i, rfrak_frame_y_i, rfrak_frame_z_i,y) :            
        ModeData=np.empty((len(y[0]),77), dtype=complex)
        def WaveformModes(M1_i, M2_i, v_i):
            I=1j
            global rhOverM_coeff
            M1=M1_i
            M2=M2_i
            v=v_i
            M=M1 + M2
            delta=(M1 - M2)/M
            nu=M1*M2/M**2
            rhOverM_coeff=8*sqrt(5)*nu*v**2*sqrt(pi)/5
            hHat_2_0_0=-5*sqrt(6)/84
            hHat_2_1_1=I*delta/3
            hHat_2_2_0=1
            hHat_3_1_1=sqrt(14)*I*delta/168
            hHat_3_3_1=-3*sqrt(210)*I*delta/56
            hHat_4_0_0=-sqrt(2)/1008
            


            rhOverM_coeff = 8*sqrt(5)*nu*v**2*sqrt(pi)/5

            Modes=np.empty(77, dtype=complex)
            # (ell, m) = (2, +/- 0)
            Symm = hHat_2_0_0*rhOverM_coeff
            Asymm = 0
            Modes[2] = Symm + Asymm
            # (ell, m) = (2, +/- 1)
            Symm = hHat_2_1_1*rhOverM_coeff*v
            Asymm = 0
            Modes[3] = Symm + Asymm
            Modes[1] = conjugate(Symm - Asymm)
            # (ell, m) = (2, +/- 2)
            Symm = hHat_2_2_0*rhOverM_coeff
            Asymm = 0
            Modes[4] = Symm + Asymm
            Modes[0] = conjugate(Symm - Asymm)
            # (ell, m) = (3, +/- 0)
            Symm = 0
            Asymm = 0
            Modes[8] = Symm + Asymm
            # (ell, m) = (3, +/- 1)
            Symm = hHat_3_1_1*rhOverM_coeff*v
            Asymm = 0
            Modes[9] = Symm + Asymm
            Modes[7] = -conjugate(Symm - Asymm)
            # (ell, m) = (3, +/- 2)
            Symm = 0
            Asymm = 0
            Modes[10] = Symm + Asymm
            Modes[6] = -conjugate(Symm - Asymm)
            # (ell, m) = (3, +/- 3)
            Symm = hHat_3_3_1*rhOverM_coeff*v
            Asymm = 0
            Modes[11] = Symm + Asymm
            Modes[5] = -conjugate(Symm - Asymm)
            # (ell, m) = (4, +/- 0)
            Symm = hHat_4_0_0*rhOverM_coeff
            Asymm = 0
            Modes[16] = Symm + Asymm
            # (ell, m) = (4, +/- 1)
            Symm = 0
            Asymm = 0
            Modes[17] = Symm + Asymm
            Modes[15] = conjugate(Symm - Asymm)
            # (ell, m) = (4, +/- 2)
            Symm = 0
            Asymm = 0
            Modes[18] = Symm + Asymm
            Modes[14] = conjugate(Symm - Asymm)
            # (ell, m) = (4, +/- 3)
            Symm = 0
            Asymm = 0
            Modes[19] = Symm + Asymm
            Modes[13] = conjugate(Symm - Asymm)
            # (ell, m) = (4, +/- 4)
            Symm = 0
            Asymm = 0
            Modes[20] = Symm + Asymm
            Modes[12] = conjugate(Symm - Asymm)
            # (ell, m) = (5, +/- 0)
            Symm = 0
            Asymm = 0
            Modes[26] = Symm + Asymm
            # (ell, m) = (5, +/- 1)
            Symm = 0
            Asymm = 0
            Modes[27] = Symm + Asymm
            Modes[25] = -conjugate(Symm - Asymm)
            # (ell, m) = (5, +/- 2)
            Symm = 0
            Asymm = 0
            Modes[28] = Symm + Asymm
            Modes[24] = -conjugate(Symm - Asymm)
            # (ell, m) = (5, +/- 3)
            Symm = 0
            Asymm = 0
            Modes[29] = Symm + Asymm
            Modes[23] = -conjugate(Symm - Asymm)
            # (ell, m) = (5, +/- 4)
            Symm = 0
            Asymm = 0
            Modes[30] = Symm + Asymm
            Modes[22] = -conjugate(Symm - Asymm)
            # (ell, m) = (5, +/- 5)
            Symm = 0
            Asymm = 0
            Modes[31] = Symm + Asymm
            Modes[21] = -conjugate(Symm - Asymm)
            # (ell, m) = (6, +/- 0)
            Symm = 0
            Asymm = 0
            Modes[38] = Symm + Asymm
            # (ell, m) = (6, +/- 1)
            Symm = 0
            Asymm = 0
            Modes[39] = Symm + Asymm
            Modes[37] = conjugate(Symm - Asymm)
            # (ell, m) = (6, +/- 2)
            Symm = 0
            Asymm = 0
            Modes[40] = Symm + Asymm
            Modes[36] = conjugate(Symm - Asymm)
            # (ell, m) = (6, +/- 3)
            Symm = 0
            Asymm = 0
            Modes[41] = Symm + Asymm
            Modes[35] = conjugate(Symm - Asymm)
            # (ell, m) = (6, +/- 4)
            Symm = 0
            Asymm = 0
            Modes[42] = Symm + Asymm
            Modes[34] = conjugate(Symm - Asymm)
            # (ell, m) = (6, +/- 5)
            Symm = 0
            Asymm = 0
            Modes[43] = Symm + Asymm
            Modes[33] = conjugate(Symm - Asymm)
            # (ell, m) = (6, +/- 6)
            Symm = 0
            Asymm = 0
            Modes[44] = Symm + Asymm
            Modes[32] = conjugate(Symm - Asymm)
            # (ell, m) = (7, +/- 0)
            Symm = 0
            Asymm = 0
            Modes[52] = Symm + Asymm
            # (ell, m) = (7, +/- 1)
            Symm = 0
            Asymm = 0
            Modes[53] = Symm + Asymm
            Modes[51] = -conjugate(Symm - Asymm)
            # (ell, m) = (7, +/- 2)
            Symm = 0
            Asymm = 0
            Modes[54] = Symm + Asymm
            Modes[50] = -conjugate(Symm - Asymm)
            # (ell, m) = (7, +/- 3)
            Symm = 0
            Asymm = 0
            Modes[55] = Symm + Asymm
            Modes[49] = -conjugate(Symm - Asymm)
            # (ell, m) = (7, +/- 4)
            Symm = 0
            Asymm = 0
            Modes[56] = Symm + Asymm
            Modes[48] = -conjugate(Symm - Asymm)
            # (ell, m) = (7, +/- 5)
            Symm = 0
            Asymm = 0
            Modes[57] = Symm + Asymm
            Modes[47] = -conjugate(Symm - Asymm)
            # (ell, m) = (7, +/- 6)
            Symm = 0
            Asymm = 0
            Modes[58] = Symm + Asymm
            Modes[46] = -conjugate(Symm - Asymm)
            # (ell, m) = (7, +/- 7)
            Symm = 0
            Asymm = 0
            Modes[59] = Symm + Asymm
            Modes[45] = -conjugate(Symm - Asymm)
            # (ell, m) = (8, +/- 0)
            Symm = 0
            Asymm = 0
            Modes[68] = Symm + Asymm
            # (ell, m) = (8, +/- 1)
            Symm = 0
            Asymm = 0
            Modes[69] = Symm + Asymm
            Modes[67] = conjugate(Symm - Asymm)
            # (ell, m) = (8, +/- 2)
            Symm = 0
            Asymm = 0
            Modes[70] = Symm + Asymm
            Modes[66] = conjugate(Symm - Asymm)
            # (ell, m) = (8, +/- 3)
            Symm = 0
            Asymm = 0
            Modes[71] = Symm + Asymm
            Modes[65] = conjugate(Symm - Asymm)
            # (ell, m) = (8, +/- 4)
            Symm = 0
            Asymm = 0
            Modes[72] = Symm + Asymm
            Modes[64] = conjugate(Symm - Asymm)
            # (ell, m) = (8, +/- 5)
            Symm = 0
            Asymm = 0
            Modes[73] = Symm + Asymm
            Modes[63] = conjugate(Symm - Asymm)
            # (ell, m) = (8, +/- 6)
            Symm = 0
            Asymm = 0
            Modes[74] = Symm + Asymm
            Modes[62] = conjugate(Symm - Asymm)
            # (ell, m) = (8, +/- 7)
            Symm = 0
            Asymm = 0
            Modes[75] = Symm + Asymm
            Modes[61] = conjugate(Symm - Asymm)
            # (ell, m) = (8, +/- 8)
            Symm = 0
            Asymm = 0
            Modes[76] = Symm + Asymm
            Modes[60] = conjugate(Symm - Asymm)

            return Modes
            
        for i in range(len(y[0])):
            v = y[0,i]
            rfrak_chi1_x = y[1,i]
            rfrak_chi1_y = y[2,i]
            rfrak_chi2_x = y[3,i]
            rfrak_chi2_y = y[4,i]
            rfrak_frame_x = y[5,i]
            rfrak_frame_y = y[6,i]
            rfrak_frame_z = y[7,i]
            Phi = y[8,i]
            xHat=xHat_i
            yHat=yHat_i
            zHat=zHat_i
            M1=M1_i
            M2=M2_i 
            S_chi1=S_chi1_i
            S_chi2=S_chi2_i
            M=M1 + M2
            delta=(M1 - M2)/M
            nu=M1*M2/M**2
            R=exp(rfrak_frame_x*xHat + rfrak_frame_y*yHat + rfrak_frame_z*zHat)
            nHat=R*xHat*conjugate(R)
            lambdaHat=R*yHat*conjugate(R)
            ellHat=R*zHat*conjugate(R)
            R_S1=exp(rfrak_chi1_x*xHat + rfrak_chi1_y*yHat)
            R_S2=exp(rfrak_chi2_x*xHat + rfrak_chi2_y*yHat)
            chiVec1=S_chi1*R_S1*zHat*conjugate(R_S1)*conjugate(S_chi1)
            chiVec2=S_chi2*R_S2*zHat*conjugate(R_S2)*conjugate(S_chi2)
            v_i=v
            rfrak_frame_x_i=rfrak_frame_x
            rfrak_frame_y_i=rfrak_frame_y
            rfrak_frame_z_i=rfrak_frame_z
            chiVec1_i=quaternion.as_float_array(chiVec1)
            chiVec2_i=quaternion.as_float_array(chiVec2)
            ModeData[i,:]=WaveformModes(M1_i, M2_i, v_i)
        return ModeData

class WaveformModes_1p0PN :
    def WaveformModes_1p0PN(xHat_i, yHat_i, zHat_i, M1_i, M2_i, v_i, S_chi1_i, S_chi2_i, rfrak_chi1_x_i, rfrak_chi1_y_i, rfrak_chi2_x_i, rfrak_chi2_y_i, rfrak_frame_x_i, rfrak_frame_y_i, rfrak_frame_z_i,y) :            
        ModeData=np.empty((len(y[0]),77), dtype=complex)
        def WaveformModes(xHat_i, yHat_i, zHat_i, M1_i, M2_i, v_i, S_chi1_i, rfrak_chi1_x_i, rfrak_chi1_y_i, rfrak_frame_x_i, rfrak_frame_y_i, rfrak_frame_z_i, chiVec1_i, chiVec2_i):
            I=1j
            global R
            global R_S1
            global chi1_n
            global chi1_lambda
            global chi1_ell
            global chi2_n
            global chi2_lambda
            global chi2_ell
            global Sigma_ell
            global Sigma_n
            global Sigma_lambda
            global rhOverM_coeff
            global hHat_spin_Symm_2_1_2
            global hHat_spin_Asymm_2_2_2
            global hHat_spin_Asymm_2_0_2
            xHat=xHat_i
            yHat=yHat_i
            zHat=zHat_i
            M1=M1_i
            M2=M2_i
            v=v_i
            S_chi1=S_chi1_i
            rfrak_chi1_x=rfrak_chi1_x_i
            rfrak_chi1_y=rfrak_chi1_y_i
            rfrak_frame_x=rfrak_frame_x_i
            rfrak_frame_y=rfrak_frame_y_i
            rfrak_frame_z=rfrak_frame_z_i
            M=M1 + M2
            delta=(M1 - M2)/M
            nu=M1*M2/M**2
            R=exp(rfrak_frame_x*xHat + rfrak_frame_y*yHat + rfrak_frame_z*zHat)
            R_S1=exp(rfrak_chi1_x*xHat + rfrak_chi1_y*yHat)
            chiVec1=chiVec1_i
            chiVec2=chiVec2_i
            chi1_n=chiVec1[1]
            chi1_lambda=chiVec1[2]
            chi1_ell=chiVec1[3]
            chi2_n=chiVec2[1]
            chi2_lambda=chiVec2[2]
            chi2_ell=chiVec2[3]
            Sigma_ell=M*(-M1*chi1_ell + M2*chi2_ell)
            Sigma_n=M*(-M1*chi1_n + M2*chi2_n)
            Sigma_lambda=M*(-M1*chi1_lambda + M2*chi2_lambda)
            rhOverM_coeff=8*sqrt(5)*nu*v**2*sqrt(pi)/5
            hHat_2_0_0=-5*sqrt(6)/84
            hHat_2_1_1=I*delta/3
            hHat_2_2_0=1
            hHat_2_2_2=55*nu/42 - 107/42
            hHat_3_1_1=sqrt(14)*I*delta/168
            hHat_3_2_2=sqrt(35)*(1 - 3*nu)/21
            hHat_3_3_1=-3*sqrt(210)*I*delta/56
            hHat_4_0_0=-sqrt(2)/1008
            hHat_4_2_2=sqrt(5)*(1 - 3*nu)/63
            hHat_4_4_2=8*sqrt(35)*(3*nu - 1)/63
            hHat_spin_Symm_2_1_2=I*Sigma_ell/(2*M**2)
            hHat_spin_Asymm_2_2_2=(-Sigma_lambda - I*Sigma_n)/(2*M**2)
            hHat_spin_Asymm_2_0_2=sqrt(6)*I*Sigma_n/(6*M**2)
            
            chiVec1 = [0., chi1_n, chi1_lambda, chi1_ell]
            chiVec2 = [0., chi2_n, chi2_lambda, chi2_ell]


            R = exp(rfrak_frame_x*xHat + rfrak_frame_y*yHat + rfrak_frame_z*zHat)
            R_S1 = exp(rfrak_chi1_x*xHat + rfrak_chi1_y*yHat)
            chi1_n = chiVec1[1]
            chi1_lambda = chiVec1[2]
            chi1_ell = chiVec1[3]
            chi2_n = chiVec2[1]
            chi2_lambda = chiVec2[2]
            chi2_ell = chiVec2[3]
            Sigma_ell = M*(-M1*chi1_ell + M2*chi2_ell)
            Sigma_n = M*(-M1*chi1_n + M2*chi2_n)
            Sigma_lambda = M*(-M1*chi1_lambda + M2*chi2_lambda)
            rhOverM_coeff = 8*sqrt(5)*nu*v**2*sqrt(pi)/5
            hHat_spin_Symm_2_1_2 = I*Sigma_ell/(2*M**2)
            hHat_spin_Asymm_2_2_2 = (-Sigma_lambda - I*Sigma_n)/(2*M**2)
            hHat_spin_Asymm_2_0_2 = sqrt(6)*I*Sigma_n/(6*M**2)

            Modes=np.empty(77, dtype=complex)
            # (ell, m) = (2, +/- 0)
            Symm = hHat_2_0_0*rhOverM_coeff
            Asymm = hHat_spin_Asymm_2_0_2*rhOverM_coeff*v**2
            Modes[2] = Symm + Asymm
            # (ell, m) = (2, +/- 1)
            Symm = rhOverM_coeff*v*(hHat_2_1_1 + hHat_spin_Symm_2_1_2*v)
            Asymm = 0
            Modes[3] = Symm + Asymm
            Modes[1] = conjugate(Symm - Asymm)
            # (ell, m) = (2, +/- 2)
            Symm = rhOverM_coeff*(hHat_2_2_0 + hHat_2_2_2*v**2)
            Asymm = hHat_spin_Asymm_2_2_2*rhOverM_coeff*v**2
            Modes[4] = Symm + Asymm
            Modes[0] = conjugate(Symm - Asymm)
            # (ell, m) = (3, +/- 0)
            Symm = 0
            Asymm = 0
            Modes[8] = Symm + Asymm
            # (ell, m) = (3, +/- 1)
            Symm = hHat_3_1_1*rhOverM_coeff*v
            Asymm = 0
            Modes[9] = Symm + Asymm
            Modes[7] = -conjugate(Symm - Asymm)
            # (ell, m) = (3, +/- 2)
            Symm = hHat_3_2_2*rhOverM_coeff*v**2
            Asymm = 0
            Modes[10] = Symm + Asymm
            Modes[6] = -conjugate(Symm - Asymm)
            # (ell, m) = (3, +/- 3)
            Symm = hHat_3_3_1*rhOverM_coeff*v
            Asymm = 0
            Modes[11] = Symm + Asymm
            Modes[5] = -conjugate(Symm - Asymm)
            # (ell, m) = (4, +/- 0)
            Symm = hHat_4_0_0*rhOverM_coeff
            Asymm = 0
            Modes[16] = Symm + Asymm
            # (ell, m) = (4, +/- 1)
            Symm = 0
            Asymm = 0
            Modes[17] = Symm + Asymm
            Modes[15] = conjugate(Symm - Asymm)
            # (ell, m) = (4, +/- 2)
            Symm = hHat_4_2_2*rhOverM_coeff*v**2
            Asymm = 0
            Modes[18] = Symm + Asymm
            Modes[14] = conjugate(Symm - Asymm)
            # (ell, m) = (4, +/- 3)
            Symm = 0
            Asymm = 0
            Modes[19] = Symm + Asymm
            Modes[13] = conjugate(Symm - Asymm)
            # (ell, m) = (4, +/- 4)
            Symm = hHat_4_4_2*rhOverM_coeff*v**2
            Asymm = 0
            Modes[20] = Symm + Asymm
            Modes[12] = conjugate(Symm - Asymm)
            # (ell, m) = (5, +/- 0)
            Symm = 0
            Asymm = 0
            Modes[26] = Symm + Asymm
            # (ell, m) = (5, +/- 1)
            Symm = 0
            Asymm = 0
            Modes[27] = Symm + Asymm
            Modes[25] = -conjugate(Symm - Asymm)
            # (ell, m) = (5, +/- 2)
            Symm = 0
            Asymm = 0
            Modes[28] = Symm + Asymm
            Modes[24] = -conjugate(Symm - Asymm)
            # (ell, m) = (5, +/- 3)
            Symm = 0
            Asymm = 0
            Modes[29] = Symm + Asymm
            Modes[23] = -conjugate(Symm - Asymm)
            # (ell, m) = (5, +/- 4)
            Symm = 0
            Asymm = 0
            Modes[30] = Symm + Asymm
            Modes[22] = -conjugate(Symm - Asymm)
            # (ell, m) = (5, +/- 5)
            Symm = 0
            Asymm = 0
            Modes[31] = Symm + Asymm
            Modes[21] = -conjugate(Symm - Asymm)
            # (ell, m) = (6, +/- 0)
            Symm = 0
            Asymm = 0
            Modes[38] = Symm + Asymm
            # (ell, m) = (6, +/- 1)
            Symm = 0
            Asymm = 0
            Modes[39] = Symm + Asymm
            Modes[37] = conjugate(Symm - Asymm)
            # (ell, m) = (6, +/- 2)
            Symm = 0
            Asymm = 0
            Modes[40] = Symm + Asymm
            Modes[36] = conjugate(Symm - Asymm)
            # (ell, m) = (6, +/- 3)
            Symm = 0
            Asymm = 0
            Modes[41] = Symm + Asymm
            Modes[35] = conjugate(Symm - Asymm)
            # (ell, m) = (6, +/- 4)
            Symm = 0
            Asymm = 0
            Modes[42] = Symm + Asymm
            Modes[34] = conjugate(Symm - Asymm)
            # (ell, m) = (6, +/- 5)
            Symm = 0
            Asymm = 0
            Modes[43] = Symm + Asymm
            Modes[33] = conjugate(Symm - Asymm)
            # (ell, m) = (6, +/- 6)
            Symm = 0
            Asymm = 0
            Modes[44] = Symm + Asymm
            Modes[32] = conjugate(Symm - Asymm)
            # (ell, m) = (7, +/- 0)
            Symm = 0
            Asymm = 0
            Modes[52] = Symm + Asymm
            # (ell, m) = (7, +/- 1)
            Symm = 0
            Asymm = 0
            Modes[53] = Symm + Asymm
            Modes[51] = -conjugate(Symm - Asymm)
            # (ell, m) = (7, +/- 2)
            Symm = 0
            Asymm = 0
            Modes[54] = Symm + Asymm
            Modes[50] = -conjugate(Symm - Asymm)
            # (ell, m) = (7, +/- 3)
            Symm = 0
            Asymm = 0
            Modes[55] = Symm + Asymm
            Modes[49] = -conjugate(Symm - Asymm)
            # (ell, m) = (7, +/- 4)
            Symm = 0
            Asymm = 0
            Modes[56] = Symm + Asymm
            Modes[48] = -conjugate(Symm - Asymm)
            # (ell, m) = (7, +/- 5)
            Symm = 0
            Asymm = 0
            Modes[57] = Symm + Asymm
            Modes[47] = -conjugate(Symm - Asymm)
            # (ell, m) = (7, +/- 6)
            Symm = 0
            Asymm = 0
            Modes[58] = Symm + Asymm
            Modes[46] = -conjugate(Symm - Asymm)
            # (ell, m) = (7, +/- 7)
            Symm = 0
            Asymm = 0
            Modes[59] = Symm + Asymm
            Modes[45] = -conjugate(Symm - Asymm)
            # (ell, m) = (8, +/- 0)
            Symm = 0
            Asymm = 0
            Modes[68] = Symm + Asymm
            # (ell, m) = (8, +/- 1)
            Symm = 0
            Asymm = 0
            Modes[69] = Symm + Asymm
            Modes[67] = conjugate(Symm - Asymm)
            # (ell, m) = (8, +/- 2)
            Symm = 0
            Asymm = 0
            Modes[70] = Symm + Asymm
            Modes[66] = conjugate(Symm - Asymm)
            # (ell, m) = (8, +/- 3)
            Symm = 0
            Asymm = 0
            Modes[71] = Symm + Asymm
            Modes[65] = conjugate(Symm - Asymm)
            # (ell, m) = (8, +/- 4)
            Symm = 0
            Asymm = 0
            Modes[72] = Symm + Asymm
            Modes[64] = conjugate(Symm - Asymm)
            # (ell, m) = (8, +/- 5)
            Symm = 0
            Asymm = 0
            Modes[73] = Symm + Asymm
            Modes[63] = conjugate(Symm - Asymm)
            # (ell, m) = (8, +/- 6)
            Symm = 0
            Asymm = 0
            Modes[74] = Symm + Asymm
            Modes[62] = conjugate(Symm - Asymm)
            # (ell, m) = (8, +/- 7)
            Symm = 0
            Asymm = 0
            Modes[75] = Symm + Asymm
            Modes[61] = conjugate(Symm - Asymm)
            # (ell, m) = (8, +/- 8)
            Symm = 0
            Asymm = 0
            Modes[76] = Symm + Asymm
            Modes[60] = conjugate(Symm - Asymm)

            return Modes
            
        for i in range(len(y[0])):
            v = y[0,i]
            rfrak_chi1_x = y[1,i]
            rfrak_chi1_y = y[2,i]
            rfrak_chi2_x = y[3,i]
            rfrak_chi2_y = y[4,i]
            rfrak_frame_x = y[5,i]
            rfrak_frame_y = y[6,i]
            rfrak_frame_z = y[7,i]
            Phi = y[8,i]
            xHat=xHat_i
            yHat=yHat_i
            zHat=zHat_i
            M1=M1_i
            M2=M2_i 
            S_chi1=S_chi1_i
            S_chi2=S_chi2_i
            M=M1 + M2
            delta=(M1 - M2)/M
            nu=M1*M2/M**2
            R=exp(rfrak_frame_x*xHat + rfrak_frame_y*yHat + rfrak_frame_z*zHat)
            nHat=R*xHat*conjugate(R)
            lambdaHat=R*yHat*conjugate(R)
            ellHat=R*zHat*conjugate(R)
            R_S1=exp(rfrak_chi1_x*xHat + rfrak_chi1_y*yHat)
            R_S2=exp(rfrak_chi2_x*xHat + rfrak_chi2_y*yHat)
            chiVec1=S_chi1*R_S1*zHat*conjugate(R_S1)*conjugate(S_chi1)
            chiVec2=S_chi2*R_S2*zHat*conjugate(R_S2)*conjugate(S_chi2)
            v_i=v
            rfrak_frame_x_i=rfrak_frame_x
            rfrak_frame_y_i=rfrak_frame_y
            rfrak_frame_z_i=rfrak_frame_z
            chiVec1_i=quaternion.as_float_array(chiVec1)
            chiVec2_i=quaternion.as_float_array(chiVec2)
            ModeData[i,:]=WaveformModes(xHat_i, yHat_i, zHat_i, M1_i, M2_i, v_i, S_chi1_i, rfrak_chi1_x_i, rfrak_chi1_y_i, rfrak_frame_x_i, rfrak_frame_y_i, rfrak_frame_z_i, chiVec1_i, chiVec2_i)
        return ModeData

class WaveformModes_1p5PN :
    def WaveformModes_1p5PN(xHat_i, yHat_i, zHat_i, M1_i, M2_i, v_i, S_chi1_i, S_chi2_i, rfrak_chi1_x_i, rfrak_chi1_y_i, rfrak_chi2_x_i, rfrak_chi2_y_i, rfrak_frame_x_i, rfrak_frame_y_i, rfrak_frame_z_i,y) :            
        ModeData=np.empty((len(y[0]),77), dtype=complex)
        def WaveformModes(M1_i, M2_i, v_i, chiVec1_i, chiVec2_i):
            I=1j
            global chi1_n
            global chi1_lambda
            global chi1_ell
            global chi2_n
            global chi2_lambda
            global chi2_ell
            global S_ell
            global S_n
            global S_lambda
            global Sigma_ell
            global Sigma_n
            global Sigma_lambda
            global rhOverM_coeff
            global hHat_spin_Symm_2_2_3
            global hHat_spin_Symm_2_1_2
            global hHat_spin_Symm_3_2_3
            global hHat_spin_Asymm_2_2_2
            global hHat_spin_Asymm_2_1_3
            global hHat_spin_Asymm_2_0_2
            global hHat_spin_Asymm_3_3_3
            global hHat_spin_Asymm_3_1_3
            M1=M1_i
            M2=M2_i
            v=v_i
            M=M1 + M2
            delta=(M1 - M2)/M
            nu=M1*M2/M**2
            chiVec1=chiVec1_i
            chiVec2=chiVec2_i
            chi1_n=chiVec1[1]
            chi1_lambda=chiVec1[2]
            chi1_ell=chiVec1[3]
            chi2_n=chiVec2[1]
            chi2_lambda=chiVec2[2]
            chi2_ell=chiVec2[3]
            S_ell=M1**2*chi1_ell + M2**2*chi2_ell
            S_n=M1**2*chi1_n + M2**2*chi2_n
            S_lambda=M1**2*chi1_lambda + M2**2*chi2_lambda
            Sigma_ell=M*(-M1*chi1_ell + M2*chi2_ell)
            Sigma_n=M*(-M1*chi1_n + M2*chi2_n)
            Sigma_lambda=M*(-M1*chi1_lambda + M2*chi2_lambda)
            rhOverM_coeff=8*sqrt(5)*nu*v**2*sqrt(pi)/5
            hHat_2_0_0=-5*sqrt(6)/84
            hHat_2_1_1=I*delta/3
            hHat_2_1_3=I*delta*(20*nu - 17)/84
            hHat_2_2_0=1
            hHat_2_2_2=55*nu/42 - 107/42
            hHat_2_2_3=2*pi
            hHat_3_1_1=sqrt(14)*I*delta/168
            hHat_3_1_3=-sqrt(14)*I*delta*(nu + 4)/252
            hHat_3_2_2=sqrt(35)*(1 - 3*nu)/21
            hHat_3_3_1=-3*sqrt(210)*I*delta/56
            hHat_3_3_3=-3*sqrt(210)*I*delta*(nu - 2)/28
            hHat_4_0_0=-sqrt(2)/1008
            hHat_4_1_3=sqrt(10)*I*delta*(1 - 2*nu)/840
            hHat_4_2_2=sqrt(5)*(1 - 3*nu)/63
            hHat_4_3_3=9*sqrt(70)*I*delta*(2*nu - 1)/280
            hHat_4_4_2=8*sqrt(35)*(3*nu - 1)/63
            hHat_5_1_3=sqrt(385)*I*delta*(1 - 2*nu)/110880
            hHat_5_3_3=9*sqrt(330)*I*delta*(2*nu - 1)/3520
            hHat_5_5_3=-625*sqrt(66)*I*delta*(2*nu - 1)/6336
            hHat_spin_Symm_2_2_3=(-6*S_ell - 2*Sigma_ell*delta)/(3*M**2)
            hHat_spin_Symm_2_1_2=I*Sigma_ell/(2*M**2)
            hHat_spin_Symm_3_2_3=2*sqrt(35)*(S_ell + Sigma_ell*delta)/(21*M**2)
            hHat_spin_Asymm_2_2_2=(-Sigma_lambda - I*Sigma_n)/(2*M**2)
            hHat_spin_Asymm_2_1_3=(4*I*S_lambda + 25*S_n + 4*I*Sigma_lambda*delta + 13*Sigma_n*delta)/(6*M**2)
            hHat_spin_Asymm_2_0_2=sqrt(6)*I*Sigma_n/(6*M**2)
            hHat_spin_Asymm_3_3_3=sqrt(210)*I*(S_lambda + I*S_n + delta*(Sigma_lambda + I*Sigma_n))/(21*M**2)
            hHat_spin_Asymm_3_1_3=sqrt(14)*(I*S_lambda + S_n + delta*(I*Sigma_lambda + Sigma_n))/(21*M**2)
            
            chiVec1 = [0., chi1_n, chi1_lambda, chi1_ell]
            chiVec2 = [0., chi2_n, chi2_lambda, chi2_ell]


            chi1_n = chiVec1[1]
            chi1_lambda = chiVec1[2]
            chi1_ell = chiVec1[3]
            chi2_n = chiVec2[1]
            chi2_lambda = chiVec2[2]
            chi2_ell = chiVec2[3]
            S_ell = M1**2*chi1_ell + M2**2*chi2_ell
            S_n = M1**2*chi1_n + M2**2*chi2_n
            S_lambda = M1**2*chi1_lambda + M2**2*chi2_lambda
            Sigma_ell = M*(-M1*chi1_ell + M2*chi2_ell)
            Sigma_n = M*(-M1*chi1_n + M2*chi2_n)
            Sigma_lambda = M*(-M1*chi1_lambda + M2*chi2_lambda)
            rhOverM_coeff = 8*sqrt(5)*nu*v**2*sqrt(pi)/5
            hHat_spin_Symm_2_2_3 = (-6*S_ell - 2*Sigma_ell*delta)/(3*M**2)
            hHat_spin_Symm_2_1_2 = I*Sigma_ell/(2*M**2)
            hHat_spin_Symm_3_2_3 = 2*sqrt(35)*(S_ell + Sigma_ell*delta)/(21*M**2)
            hHat_spin_Asymm_2_2_2 = (-Sigma_lambda - I*Sigma_n)/(2*M**2)
            hHat_spin_Asymm_2_1_3 = (4*I*S_lambda + 25*S_n + 4*I*Sigma_lambda*delta + 13*Sigma_n*delta)/(6*M**2)
            hHat_spin_Asymm_2_0_2 = sqrt(6)*I*Sigma_n/(6*M**2)
            hHat_spin_Asymm_3_3_3 = sqrt(210)*I*(S_lambda + I*S_n + delta*(Sigma_lambda + I*Sigma_n))/(21*M**2)
            hHat_spin_Asymm_3_1_3 = sqrt(14)*(I*S_lambda + S_n + delta*(I*Sigma_lambda + Sigma_n))/(21*M**2)

            Modes=np.empty(77, dtype=complex)
            # (ell, m) = (2, +/- 0)
            Symm = hHat_2_0_0*rhOverM_coeff
            Asymm = hHat_spin_Asymm_2_0_2*rhOverM_coeff*v**2
            Modes[2] = Symm + Asymm
            # (ell, m) = (2, +/- 1)
            Symm = rhOverM_coeff*v*(hHat_2_1_1 + v*(hHat_2_1_3*v + hHat_spin_Symm_2_1_2))
            Asymm = hHat_spin_Asymm_2_1_3*rhOverM_coeff*v**3
            Modes[3] = Symm + Asymm
            Modes[1] = conjugate(Symm - Asymm)
            # (ell, m) = (2, +/- 2)
            Symm = rhOverM_coeff*(hHat_2_2_0 + v**2*(hHat_2_2_2 + v*(hHat_2_2_3 + hHat_spin_Symm_2_2_3)))
            Asymm = hHat_spin_Asymm_2_2_2*rhOverM_coeff*v**2
            Modes[4] = Symm + Asymm
            Modes[0] = conjugate(Symm - Asymm)
            # (ell, m) = (3, +/- 0)
            Symm = 0
            Asymm = 0
            Modes[8] = Symm + Asymm
            # (ell, m) = (3, +/- 1)
            Symm = rhOverM_coeff*v*(hHat_3_1_1 + hHat_3_1_3*v**2)
            Asymm = hHat_spin_Asymm_3_1_3*rhOverM_coeff*v**3
            Modes[9] = Symm + Asymm
            Modes[7] = -conjugate(Symm - Asymm)
            # (ell, m) = (3, +/- 2)
            Symm = rhOverM_coeff*v**2*(hHat_3_2_2 + hHat_spin_Symm_3_2_3*v)
            Asymm = 0
            Modes[10] = Symm + Asymm
            Modes[6] = -conjugate(Symm - Asymm)
            # (ell, m) = (3, +/- 3)
            Symm = rhOverM_coeff*v*(hHat_3_3_1 + hHat_3_3_3*v**2)
            Asymm = hHat_spin_Asymm_3_3_3*rhOverM_coeff*v**3
            Modes[11] = Symm + Asymm
            Modes[5] = -conjugate(Symm - Asymm)
            # (ell, m) = (4, +/- 0)
            Symm = hHat_4_0_0*rhOverM_coeff
            Asymm = 0
            Modes[16] = Symm + Asymm
            # (ell, m) = (4, +/- 1)
            Symm = hHat_4_1_3*rhOverM_coeff*v**3
            Asymm = 0
            Modes[17] = Symm + Asymm
            Modes[15] = conjugate(Symm - Asymm)
            # (ell, m) = (4, +/- 2)
            Symm = hHat_4_2_2*rhOverM_coeff*v**2
            Asymm = 0
            Modes[18] = Symm + Asymm
            Modes[14] = conjugate(Symm - Asymm)
            # (ell, m) = (4, +/- 3)
            Symm = hHat_4_3_3*rhOverM_coeff*v**3
            Asymm = 0
            Modes[19] = Symm + Asymm
            Modes[13] = conjugate(Symm - Asymm)
            # (ell, m) = (4, +/- 4)
            Symm = hHat_4_4_2*rhOverM_coeff*v**2
            Asymm = 0
            Modes[20] = Symm + Asymm
            Modes[12] = conjugate(Symm - Asymm)
            # (ell, m) = (5, +/- 0)
            Symm = 0
            Asymm = 0
            Modes[26] = Symm + Asymm
            # (ell, m) = (5, +/- 1)
            Symm = hHat_5_1_3*rhOverM_coeff*v**3
            Asymm = 0
            Modes[27] = Symm + Asymm
            Modes[25] = -conjugate(Symm - Asymm)
            # (ell, m) = (5, +/- 2)
            Symm = 0
            Asymm = 0
            Modes[28] = Symm + Asymm
            Modes[24] = -conjugate(Symm - Asymm)
            # (ell, m) = (5, +/- 3)
            Symm = hHat_5_3_3*rhOverM_coeff*v**3
            Asymm = 0
            Modes[29] = Symm + Asymm
            Modes[23] = -conjugate(Symm - Asymm)
            # (ell, m) = (5, +/- 4)
            Symm = 0
            Asymm = 0
            Modes[30] = Symm + Asymm
            Modes[22] = -conjugate(Symm - Asymm)
            # (ell, m) = (5, +/- 5)
            Symm = hHat_5_5_3*rhOverM_coeff*v**3
            Asymm = 0
            Modes[31] = Symm + Asymm
            Modes[21] = -conjugate(Symm - Asymm)
            # (ell, m) = (6, +/- 0)
            Symm = 0
            Asymm = 0
            Modes[38] = Symm + Asymm
            # (ell, m) = (6, +/- 1)
            Symm = 0
            Asymm = 0
            Modes[39] = Symm + Asymm
            Modes[37] = conjugate(Symm - Asymm)
            # (ell, m) = (6, +/- 2)
            Symm = 0
            Asymm = 0
            Modes[40] = Symm + Asymm
            Modes[36] = conjugate(Symm - Asymm)
            # (ell, m) = (6, +/- 3)
            Symm = 0
            Asymm = 0
            Modes[41] = Symm + Asymm
            Modes[35] = conjugate(Symm - Asymm)
            # (ell, m) = (6, +/- 4)
            Symm = 0
            Asymm = 0
            Modes[42] = Symm + Asymm
            Modes[34] = conjugate(Symm - Asymm)
            # (ell, m) = (6, +/- 5)
            Symm = 0
            Asymm = 0
            Modes[43] = Symm + Asymm
            Modes[33] = conjugate(Symm - Asymm)
            # (ell, m) = (6, +/- 6)
            Symm = 0
            Asymm = 0
            Modes[44] = Symm + Asymm
            Modes[32] = conjugate(Symm - Asymm)
            # (ell, m) = (7, +/- 0)
            Symm = 0
            Asymm = 0
            Modes[52] = Symm + Asymm
            # (ell, m) = (7, +/- 1)
            Symm = 0
            Asymm = 0
            Modes[53] = Symm + Asymm
            Modes[51] = -conjugate(Symm - Asymm)
            # (ell, m) = (7, +/- 2)
            Symm = 0
            Asymm = 0
            Modes[54] = Symm + Asymm
            Modes[50] = -conjugate(Symm - Asymm)
            # (ell, m) = (7, +/- 3)
            Symm = 0
            Asymm = 0
            Modes[55] = Symm + Asymm
            Modes[49] = -conjugate(Symm - Asymm)
            # (ell, m) = (7, +/- 4)
            Symm = 0
            Asymm = 0
            Modes[56] = Symm + Asymm
            Modes[48] = -conjugate(Symm - Asymm)
            # (ell, m) = (7, +/- 5)
            Symm = 0
            Asymm = 0
            Modes[57] = Symm + Asymm
            Modes[47] = -conjugate(Symm - Asymm)
            # (ell, m) = (7, +/- 6)
            Symm = 0
            Asymm = 0
            Modes[58] = Symm + Asymm
            Modes[46] = -conjugate(Symm - Asymm)
            # (ell, m) = (7, +/- 7)
            Symm = 0
            Asymm = 0
            Modes[59] = Symm + Asymm
            Modes[45] = -conjugate(Symm - Asymm)
            # (ell, m) = (8, +/- 0)
            Symm = 0
            Asymm = 0
            Modes[68] = Symm + Asymm
            # (ell, m) = (8, +/- 1)
            Symm = 0
            Asymm = 0
            Modes[69] = Symm + Asymm
            Modes[67] = conjugate(Symm - Asymm)
            # (ell, m) = (8, +/- 2)
            Symm = 0
            Asymm = 0
            Modes[70] = Symm + Asymm
            Modes[66] = conjugate(Symm - Asymm)
            # (ell, m) = (8, +/- 3)
            Symm = 0
            Asymm = 0
            Modes[71] = Symm + Asymm
            Modes[65] = conjugate(Symm - Asymm)
            # (ell, m) = (8, +/- 4)
            Symm = 0
            Asymm = 0
            Modes[72] = Symm + Asymm
            Modes[64] = conjugate(Symm - Asymm)
            # (ell, m) = (8, +/- 5)
            Symm = 0
            Asymm = 0
            Modes[73] = Symm + Asymm
            Modes[63] = conjugate(Symm - Asymm)
            # (ell, m) = (8, +/- 6)
            Symm = 0
            Asymm = 0
            Modes[74] = Symm + Asymm
            Modes[62] = conjugate(Symm - Asymm)
            # (ell, m) = (8, +/- 7)
            Symm = 0
            Asymm = 0
            Modes[75] = Symm + Asymm
            Modes[61] = conjugate(Symm - Asymm)
            # (ell, m) = (8, +/- 8)
            Symm = 0
            Asymm = 0
            Modes[76] = Symm + Asymm
            Modes[60] = conjugate(Symm - Asymm)

            return Modes
            
        for i in range(len(y[0])):
            v = y[0,i]
            rfrak_chi1_x = y[1,i]
            rfrak_chi1_y = y[2,i]
            rfrak_chi2_x = y[3,i]
            rfrak_chi2_y = y[4,i]
            rfrak_frame_x = y[5,i]
            rfrak_frame_y = y[6,i]
            rfrak_frame_z = y[7,i]
            Phi = y[8,i]
            xHat=xHat_i
            yHat=yHat_i
            zHat=zHat_i
            M1=M1_i
            M2=M2_i 
            S_chi1=S_chi1_i
            S_chi2=S_chi2_i
            M=M1 + M2
            delta=(M1 - M2)/M
            nu=M1*M2/M**2
            R=exp(rfrak_frame_x*xHat + rfrak_frame_y*yHat + rfrak_frame_z*zHat)
            nHat=R*xHat*conjugate(R)
            lambdaHat=R*yHat*conjugate(R)
            ellHat=R*zHat*conjugate(R)
            R_S1=exp(rfrak_chi1_x*xHat + rfrak_chi1_y*yHat)
            R_S2=exp(rfrak_chi2_x*xHat + rfrak_chi2_y*yHat)
            chiVec1=S_chi1*R_S1*zHat*conjugate(R_S1)*conjugate(S_chi1)
            chiVec2=S_chi2*R_S2*zHat*conjugate(R_S2)*conjugate(S_chi2)
            v_i=v
            rfrak_frame_x_i=rfrak_frame_x
            rfrak_frame_y_i=rfrak_frame_y
            rfrak_frame_z_i=rfrak_frame_z
            chiVec1_i=quaternion.as_float_array(chiVec1)
            chiVec2_i=quaternion.as_float_array(chiVec2)
            ModeData[i,:]=WaveformModes(M1_i, M2_i, v_i, chiVec1_i, chiVec2_i)
        return ModeData

class WaveformModes_2p0PN :
    def WaveformModes_2p0PN(xHat_i, yHat_i, zHat_i, M1_i, M2_i, v_i, S_chi1_i, S_chi2_i, rfrak_chi1_x_i, rfrak_chi1_y_i, rfrak_chi2_x_i, rfrak_chi2_y_i, rfrak_frame_x_i, rfrak_frame_y_i, rfrak_frame_z_i,y) :            
        ModeData=np.empty((len(y[0]),77), dtype=complex)
        def WaveformModes(M1_i, M2_i, v_i, chiVec1_i, chiVec2_i):
            I=1j
            global chi1_n
            global chi1_lambda
            global chi1_ell
            global chi2_n
            global chi2_lambda
            global chi2_ell
            global S_ell
            global S_n
            global S_lambda
            global Sigma_ell
            global Sigma_n
            global Sigma_lambda
            global S1_ell
            global S1_n
            global S1_lambda
            global S2_ell
            global S2_n
            global S2_lambda
            global rhOverM_coeff
            global hHat_spin_Symm_2_2_3
            global hHat_spin_Symm_2_2_4
            global hHat_spin_Symm_2_1_2
            global hHat_spin_Symm_2_1_4
            global hHat_spin_Symm_2_0_4
            global hHat_spin_Symm_3_3_4
            global hHat_spin_Symm_3_2_3
            global hHat_spin_Symm_3_1_4
            global hHat_spin_Symm_4_3_4
            global hHat_spin_Symm_4_1_4
            global hHat_spin_Asymm_2_2_2
            global hHat_spin_Asymm_2_2_4
            global hHat_spin_Asymm_2_1_3
            global hHat_spin_Asymm_2_1_4
            global hHat_spin_Asymm_2_0_2
            global hHat_spin_Asymm_2_0_4
            global hHat_spin_Asymm_3_3_3
            global hHat_spin_Asymm_3_2_4
            global hHat_spin_Asymm_3_1_3
            global hHat_spin_Asymm_3_0_4
            global hHat_spin_Asymm_4_4_4
            global hHat_spin_Asymm_4_2_4
            global hHat_spin_Asymm_4_0_4
            M1=M1_i
            M2=M2_i
            v=v_i
            M=M1 + M2
            delta=(M1 - M2)/M
            nu=M1*M2/M**2
            chiVec1=chiVec1_i
            chiVec2=chiVec2_i
            chi1_n=chiVec1[1]
            chi1_lambda=chiVec1[2]
            chi1_ell=chiVec1[3]
            chi2_n=chiVec2[1]
            chi2_lambda=chiVec2[2]
            chi2_ell=chiVec2[3]
            S_ell=M1**2*chi1_ell + M2**2*chi2_ell
            S_n=M1**2*chi1_n + M2**2*chi2_n
            S_lambda=M1**2*chi1_lambda + M2**2*chi2_lambda
            Sigma_ell=M*(-M1*chi1_ell + M2*chi2_ell)
            Sigma_n=M*(-M1*chi1_n + M2*chi2_n)
            Sigma_lambda=M*(-M1*chi1_lambda + M2*chi2_lambda)
            S1_ell=M1**2*chi1_ell
            S1_n=M1**2*chi1_n
            S1_lambda=M1**2*chi1_lambda
            S2_ell=M2**2*chi2_ell
            S2_n=M2**2*chi2_n
            S2_lambda=M2**2*chi2_lambda
            rhOverM_coeff=8*sqrt(5)*nu*v**2*sqrt(pi)/5
            hHat_2_0_0=-5*sqrt(6)/84
            hHat_2_1_1=I*delta/3
            hHat_2_1_3=I*delta*(20*nu - 17)/84
            hHat_2_1_4=delta*(1 + log(16) + 2*I*pi)/6
            hHat_2_2_0=1
            hHat_2_2_2=55*nu/42 - 107/42
            hHat_2_2_3=2*pi
            hHat_2_2_4=nu*(2047*nu - 7483)/1512 - 2173/1512
            hHat_3_1_1=sqrt(14)*I*delta/168
            hHat_3_1_3=-sqrt(14)*I*delta*(nu + 4)/252
            hHat_3_1_4=sqrt(14)*delta*(log(1024) + 7 + 5*I*pi)/840
            hHat_3_2_2=sqrt(35)*(1 - 3*nu)/21
            hHat_3_2_4=sqrt(35)*(nu*(725 - 365*nu) - 193)/1890
            hHat_3_3_1=-3*sqrt(210)*I*delta/56
            hHat_3_3_3=-3*sqrt(210)*I*delta*(nu - 2)/28
            hHat_3_3_4=9*sqrt(210)*delta*(-7 + 10*log(3/2) - 5*I*pi)/280
            hHat_4_0_0=-sqrt(2)/1008
            hHat_4_1_3=sqrt(10)*I*delta*(1 - 2*nu)/840
            hHat_4_2_2=sqrt(5)*(1 - 3*nu)/63
            hHat_4_2_4=sqrt(5)*(nu*(4025 - 285*nu) - 1311)/20790
            hHat_4_3_3=9*sqrt(70)*I*delta*(2*nu - 1)/280
            hHat_4_4_2=8*sqrt(35)*(3*nu - 1)/63
            hHat_4_4_4=sqrt(35)*(20*nu*(525*nu - 1273) + 7116)/10395
            hHat_5_1_3=sqrt(385)*I*delta*(1 - 2*nu)/110880
            hHat_5_2_4=sqrt(55)*(2*nu*(5*nu - 5) + 2)/1485
            hHat_5_3_3=9*sqrt(330)*I*delta*(2*nu - 1)/3520
            hHat_5_4_4=sqrt(165)*(-32*nu*(5*nu - 5) - 32)/1485
            hHat_5_5_3=-625*sqrt(66)*I*delta*(2*nu - 1)/6336
            hHat_6_2_4=sqrt(65)*(2*nu*(5*nu - 5) + 2)/19305
            hHat_6_4_4=-128*sqrt(78)*(nu*(5*nu - 5) + 1)/19305
            hHat_6_6_4=sqrt(143)*(54*nu*(5*nu - 5) + 54)/715
            hHat_spin_Symm_2_2_3=(-6*S_ell - 2*Sigma_ell*delta)/(3*M**2)
            hHat_spin_Symm_2_2_4=(12*S1_ell*S2_ell + 10*S1_lambda*S2_lambda - 15*I*S1_lambda*S2_n - 15*I*S1_n*S2_lambda - 22*S1_n*S2_n)/(6*M**4*nu)
            hHat_spin_Symm_2_1_2=I*Sigma_ell/(2*M**2)
            hHat_spin_Symm_2_1_4=I*(-86*S_ell*delta + Sigma_ell*(139*nu - 79))/(42*M**2)
            hHat_spin_Symm_2_0_4=sqrt(6)*(-S1_lambda*S2_lambda + S1_n*S2_n)/(3*M**4*nu)
            hHat_spin_Symm_3_3_4=3*sqrt(210)*I*(7*S_ell*delta - 3*Sigma_ell*(3*nu - 1))/(112*M**2)
            hHat_spin_Symm_3_2_3=2*sqrt(35)*(S_ell + Sigma_ell*delta)/(21*M**2)
            hHat_spin_Symm_3_1_4=sqrt(14)*I*(S_ell*delta - 5*Sigma_ell*(3*nu - 1))/(336*M**2)
            hHat_spin_Symm_4_3_4=9*sqrt(70)*I*(-S_ell*delta + 3*Sigma_ell*nu - Sigma_ell)/(112*M**2)
            hHat_spin_Symm_4_1_4=sqrt(10)*I*(S_ell*delta - 3*Sigma_ell*nu + Sigma_ell)/(336*M**2)
            hHat_spin_Asymm_2_2_2=(-Sigma_lambda - I*Sigma_n)/(2*M**2)
            hHat_spin_Asymm_2_2_4=(19*S_lambda*delta + 182*I*S_n*delta - 43*Sigma_lambda*nu + 5*Sigma_lambda - 280*I*Sigma_n*nu + 98*I*Sigma_n)/(84*M**2)
            hHat_spin_Asymm_2_1_3=(4*I*S_lambda + 25*S_n + 4*I*Sigma_lambda*delta + 13*Sigma_n*delta)/(6*M**2)
            hHat_spin_Asymm_2_1_4=(-3*S1_ell*S2_n - 3*S1_n*S2_ell)/(2*M**4*nu)
            hHat_spin_Asymm_2_0_2=sqrt(6)*I*Sigma_n/(6*M**2)
            hHat_spin_Asymm_2_0_4=sqrt(6)*I*(255*S_n*delta - Sigma_n*(506*nu - 45))/(126*M**2)
            hHat_spin_Asymm_3_3_3=sqrt(210)*I*(S_lambda + I*S_n + delta*(Sigma_lambda + I*Sigma_n))/(21*M**2)
            hHat_spin_Asymm_3_2_4=sqrt(35)*(-Sigma_lambda*(83*nu - 17) + 4*I*Sigma_n*(55*nu - 13) + 25*delta*(S_lambda - 4*I*S_n))/(168*M**2)
            hHat_spin_Asymm_3_1_3=sqrt(14)*(I*S_lambda + S_n + delta*(I*Sigma_lambda + Sigma_n))/(21*M**2)
            hHat_spin_Asymm_3_0_4=sqrt(42)*(-17*S_lambda*delta + Sigma_lambda*(35*nu - 9))/(168*M**2)
            hHat_spin_Asymm_4_4_4=9*sqrt(35)*(-3*Sigma_lambda*nu + Sigma_lambda - I*Sigma_n*(3*nu - 1) + delta*(S_lambda + I*S_n))/(56*M**2)
            hHat_spin_Asymm_4_2_4=sqrt(5)*(-13*Sigma_lambda*(3*nu - 1) + 14*I*Sigma_n*(3*nu - 1) + delta*(13*S_lambda - 14*I*S_n))/(168*M**2)
            hHat_spin_Asymm_4_0_4=sqrt(2)*I*(S_n*delta - 3*Sigma_n*nu + Sigma_n)/(168*M**2)
            
            chiVec1 = [0., chi1_n, chi1_lambda, chi1_ell]
            chiVec2 = [0., chi2_n, chi2_lambda, chi2_ell]


            chi1_n = chiVec1[1]
            chi1_lambda = chiVec1[2]
            chi1_ell = chiVec1[3]
            chi2_n = chiVec2[1]
            chi2_lambda = chiVec2[2]
            chi2_ell = chiVec2[3]
            S_ell = M1**2*chi1_ell + M2**2*chi2_ell
            S_n = M1**2*chi1_n + M2**2*chi2_n
            S_lambda = M1**2*chi1_lambda + M2**2*chi2_lambda
            Sigma_ell = M*(-M1*chi1_ell + M2*chi2_ell)
            Sigma_n = M*(-M1*chi1_n + M2*chi2_n)
            Sigma_lambda = M*(-M1*chi1_lambda + M2*chi2_lambda)
            S1_ell = M1**2*chi1_ell
            S1_n = M1**2*chi1_n
            S1_lambda = M1**2*chi1_lambda
            S2_ell = M2**2*chi2_ell
            S2_n = M2**2*chi2_n
            S2_lambda = M2**2*chi2_lambda
            rhOverM_coeff = 8*sqrt(5)*nu*v**2*sqrt(pi)/5
            hHat_spin_Symm_2_2_3 = (-6*S_ell - 2*Sigma_ell*delta)/(3*M**2)
            hHat_spin_Symm_2_2_4 = (12*S1_ell*S2_ell + 10*S1_lambda*S2_lambda - 15*I*S1_lambda*S2_n - 15*I*S1_n*S2_lambda - 22*S1_n*S2_n)/(6*M**4*nu)
            hHat_spin_Symm_2_1_2 = I*Sigma_ell/(2*M**2)
            hHat_spin_Symm_2_1_4 = I*(-86*S_ell*delta + Sigma_ell*(139*nu - 79))/(42*M**2)
            hHat_spin_Symm_2_0_4 = sqrt(6)*(-S1_lambda*S2_lambda + S1_n*S2_n)/(3*M**4*nu)
            hHat_spin_Symm_3_3_4 = 3*sqrt(210)*I*(7*S_ell*delta - 3*Sigma_ell*(3*nu - 1))/(112*M**2)
            hHat_spin_Symm_3_2_3 = 2*sqrt(35)*(S_ell + Sigma_ell*delta)/(21*M**2)
            hHat_spin_Symm_3_1_4 = sqrt(14)*I*(S_ell*delta - 5*Sigma_ell*(3*nu - 1))/(336*M**2)
            hHat_spin_Symm_4_3_4 = 9*sqrt(70)*I*(-S_ell*delta + 3*Sigma_ell*nu - Sigma_ell)/(112*M**2)
            hHat_spin_Symm_4_1_4 = sqrt(10)*I*(S_ell*delta - 3*Sigma_ell*nu + Sigma_ell)/(336*M**2)
            hHat_spin_Asymm_2_2_2 = (-Sigma_lambda - I*Sigma_n)/(2*M**2)
            hHat_spin_Asymm_2_2_4 = (19*S_lambda*delta + 182*I*S_n*delta - 43*Sigma_lambda*nu + 5*Sigma_lambda - 280*I*Sigma_n*nu + 98*I*Sigma_n)/(84*M**2)
            hHat_spin_Asymm_2_1_3 = (4*I*S_lambda + 25*S_n + 4*I*Sigma_lambda*delta + 13*Sigma_n*delta)/(6*M**2)
            hHat_spin_Asymm_2_1_4 = (-3*S1_ell*S2_n - 3*S1_n*S2_ell)/(2*M**4*nu)
            hHat_spin_Asymm_2_0_2 = sqrt(6)*I*Sigma_n/(6*M**2)
            hHat_spin_Asymm_2_0_4 = sqrt(6)*I*(255*S_n*delta - Sigma_n*(506*nu - 45))/(126*M**2)
            hHat_spin_Asymm_3_3_3 = sqrt(210)*I*(S_lambda + I*S_n + delta*(Sigma_lambda + I*Sigma_n))/(21*M**2)
            hHat_spin_Asymm_3_2_4 = sqrt(35)*(-Sigma_lambda*(83*nu - 17) + 4*I*Sigma_n*(55*nu - 13) + 25*delta*(S_lambda - 4*I*S_n))/(168*M**2)
            hHat_spin_Asymm_3_1_3 = sqrt(14)*(I*S_lambda + S_n + delta*(I*Sigma_lambda + Sigma_n))/(21*M**2)
            hHat_spin_Asymm_3_0_4 = sqrt(42)*(-17*S_lambda*delta + Sigma_lambda*(35*nu - 9))/(168*M**2)
            hHat_spin_Asymm_4_4_4 = 9*sqrt(35)*(-3*Sigma_lambda*nu + Sigma_lambda - I*Sigma_n*(3*nu - 1) + delta*(S_lambda + I*S_n))/(56*M**2)
            hHat_spin_Asymm_4_2_4 = sqrt(5)*(-13*Sigma_lambda*(3*nu - 1) + 14*I*Sigma_n*(3*nu - 1) + delta*(13*S_lambda - 14*I*S_n))/(168*M**2)
            hHat_spin_Asymm_4_0_4 = sqrt(2)*I*(S_n*delta - 3*Sigma_n*nu + Sigma_n)/(168*M**2)

            Modes=np.empty(77, dtype=complex)
            # (ell, m) = (2, +/- 0)
            Symm = rhOverM_coeff*(hHat_2_0_0 + hHat_spin_Symm_2_0_4*v**4)
            Asymm = rhOverM_coeff*v**2*(hHat_spin_Asymm_2_0_2 + hHat_spin_Asymm_2_0_4*v**2)
            Modes[2] = Symm + Asymm
            # (ell, m) = (2, +/- 1)
            Symm = rhOverM_coeff*v*(hHat_2_1_1 + v*(hHat_spin_Symm_2_1_2 + v*(hHat_2_1_3 + v*(hHat_2_1_4 + hHat_spin_Symm_2_1_4))))
            Asymm = rhOverM_coeff*v**3*(hHat_spin_Asymm_2_1_3 + hHat_spin_Asymm_2_1_4*v)
            Modes[3] = Symm + Asymm
            Modes[1] = conjugate(Symm - Asymm)
            # (ell, m) = (2, +/- 2)
            Symm = rhOverM_coeff*(hHat_2_2_0 + v**2*(hHat_2_2_2 + v*(hHat_2_2_3 + hHat_spin_Symm_2_2_3 + v*(hHat_2_2_4 + hHat_spin_Symm_2_2_4))))
            Asymm = rhOverM_coeff*v**2*(hHat_spin_Asymm_2_2_2 + hHat_spin_Asymm_2_2_4*v**2)
            Modes[4] = Symm + Asymm
            Modes[0] = conjugate(Symm - Asymm)
            # (ell, m) = (3, +/- 0)
            Symm = 0
            Asymm = hHat_spin_Asymm_3_0_4*rhOverM_coeff*v**4
            Modes[8] = Symm + Asymm
            # (ell, m) = (3, +/- 1)
            Symm = rhOverM_coeff*v*(hHat_3_1_1 + v**2*(hHat_3_1_3 + v*(hHat_3_1_4 + hHat_spin_Symm_3_1_4)))
            Asymm = hHat_spin_Asymm_3_1_3*rhOverM_coeff*v**3
            Modes[9] = Symm + Asymm
            Modes[7] = -conjugate(Symm - Asymm)
            # (ell, m) = (3, +/- 2)
            Symm = rhOverM_coeff*v**2*(hHat_3_2_2 + v*(hHat_3_2_4*v + hHat_spin_Symm_3_2_3))
            Asymm = hHat_spin_Asymm_3_2_4*rhOverM_coeff*v**4
            Modes[10] = Symm + Asymm
            Modes[6] = -conjugate(Symm - Asymm)
            # (ell, m) = (3, +/- 3)
            Symm = rhOverM_coeff*v*(hHat_3_3_1 + v**2*(hHat_3_3_3 + v*(hHat_3_3_4 + hHat_spin_Symm_3_3_4)))
            Asymm = hHat_spin_Asymm_3_3_3*rhOverM_coeff*v**3
            Modes[11] = Symm + Asymm
            Modes[5] = -conjugate(Symm - Asymm)
            # (ell, m) = (4, +/- 0)
            Symm = hHat_4_0_0*rhOverM_coeff
            Asymm = hHat_spin_Asymm_4_0_4*rhOverM_coeff*v**4
            Modes[16] = Symm + Asymm
            # (ell, m) = (4, +/- 1)
            Symm = rhOverM_coeff*v**3*(hHat_4_1_3 + hHat_spin_Symm_4_1_4*v)
            Asymm = 0
            Modes[17] = Symm + Asymm
            Modes[15] = conjugate(Symm - Asymm)
            # (ell, m) = (4, +/- 2)
            Symm = rhOverM_coeff*v**2*(hHat_4_2_2 + hHat_4_2_4*v**2)
            Asymm = hHat_spin_Asymm_4_2_4*rhOverM_coeff*v**4
            Modes[18] = Symm + Asymm
            Modes[14] = conjugate(Symm - Asymm)
            # (ell, m) = (4, +/- 3)
            Symm = rhOverM_coeff*v**3*(hHat_4_3_3 + hHat_spin_Symm_4_3_4*v)
            Asymm = 0
            Modes[19] = Symm + Asymm
            Modes[13] = conjugate(Symm - Asymm)
            # (ell, m) = (4, +/- 4)
            Symm = rhOverM_coeff*v**2*(hHat_4_4_2 + hHat_4_4_4*v**2)
            Asymm = hHat_spin_Asymm_4_4_4*rhOverM_coeff*v**4
            Modes[20] = Symm + Asymm
            Modes[12] = conjugate(Symm - Asymm)
            # (ell, m) = (5, +/- 0)
            Symm = 0
            Asymm = 0
            Modes[26] = Symm + Asymm
            # (ell, m) = (5, +/- 1)
            Symm = hHat_5_1_3*rhOverM_coeff*v**3
            Asymm = 0
            Modes[27] = Symm + Asymm
            Modes[25] = -conjugate(Symm - Asymm)
            # (ell, m) = (5, +/- 2)
            Symm = hHat_5_2_4*rhOverM_coeff*v**4
            Asymm = 0
            Modes[28] = Symm + Asymm
            Modes[24] = -conjugate(Symm - Asymm)
            # (ell, m) = (5, +/- 3)
            Symm = hHat_5_3_3*rhOverM_coeff*v**3
            Asymm = 0
            Modes[29] = Symm + Asymm
            Modes[23] = -conjugate(Symm - Asymm)
            # (ell, m) = (5, +/- 4)
            Symm = hHat_5_4_4*rhOverM_coeff*v**4
            Asymm = 0
            Modes[30] = Symm + Asymm
            Modes[22] = -conjugate(Symm - Asymm)
            # (ell, m) = (5, +/- 5)
            Symm = hHat_5_5_3*rhOverM_coeff*v**3
            Asymm = 0
            Modes[31] = Symm + Asymm
            Modes[21] = -conjugate(Symm - Asymm)
            # (ell, m) = (6, +/- 0)
            Symm = 0
            Asymm = 0
            Modes[38] = Symm + Asymm
            # (ell, m) = (6, +/- 1)
            Symm = 0
            Asymm = 0
            Modes[39] = Symm + Asymm
            Modes[37] = conjugate(Symm - Asymm)
            # (ell, m) = (6, +/- 2)
            Symm = hHat_6_2_4*rhOverM_coeff*v**4
            Asymm = 0
            Modes[40] = Symm + Asymm
            Modes[36] = conjugate(Symm - Asymm)
            # (ell, m) = (6, +/- 3)
            Symm = 0
            Asymm = 0
            Modes[41] = Symm + Asymm
            Modes[35] = conjugate(Symm - Asymm)
            # (ell, m) = (6, +/- 4)
            Symm = hHat_6_4_4*rhOverM_coeff*v**4
            Asymm = 0
            Modes[42] = Symm + Asymm
            Modes[34] = conjugate(Symm - Asymm)
            # (ell, m) = (6, +/- 5)
            Symm = 0
            Asymm = 0
            Modes[43] = Symm + Asymm
            Modes[33] = conjugate(Symm - Asymm)
            # (ell, m) = (6, +/- 6)
            Symm = hHat_6_6_4*rhOverM_coeff*v**4
            Asymm = 0
            Modes[44] = Symm + Asymm
            Modes[32] = conjugate(Symm - Asymm)
            # (ell, m) = (7, +/- 0)
            Symm = 0
            Asymm = 0
            Modes[52] = Symm + Asymm
            # (ell, m) = (7, +/- 1)
            Symm = 0
            Asymm = 0
            Modes[53] = Symm + Asymm
            Modes[51] = -conjugate(Symm - Asymm)
            # (ell, m) = (7, +/- 2)
            Symm = 0
            Asymm = 0
            Modes[54] = Symm + Asymm
            Modes[50] = -conjugate(Symm - Asymm)
            # (ell, m) = (7, +/- 3)
            Symm = 0
            Asymm = 0
            Modes[55] = Symm + Asymm
            Modes[49] = -conjugate(Symm - Asymm)
            # (ell, m) = (7, +/- 4)
            Symm = 0
            Asymm = 0
            Modes[56] = Symm + Asymm
            Modes[48] = -conjugate(Symm - Asymm)
            # (ell, m) = (7, +/- 5)
            Symm = 0
            Asymm = 0
            Modes[57] = Symm + Asymm
            Modes[47] = -conjugate(Symm - Asymm)
            # (ell, m) = (7, +/- 6)
            Symm = 0
            Asymm = 0
            Modes[58] = Symm + Asymm
            Modes[46] = -conjugate(Symm - Asymm)
            # (ell, m) = (7, +/- 7)
            Symm = 0
            Asymm = 0
            Modes[59] = Symm + Asymm
            Modes[45] = -conjugate(Symm - Asymm)
            # (ell, m) = (8, +/- 0)
            Symm = 0
            Asymm = 0
            Modes[68] = Symm + Asymm
            # (ell, m) = (8, +/- 1)
            Symm = 0
            Asymm = 0
            Modes[69] = Symm + Asymm
            Modes[67] = conjugate(Symm - Asymm)
            # (ell, m) = (8, +/- 2)
            Symm = 0
            Asymm = 0
            Modes[70] = Symm + Asymm
            Modes[66] = conjugate(Symm - Asymm)
            # (ell, m) = (8, +/- 3)
            Symm = 0
            Asymm = 0
            Modes[71] = Symm + Asymm
            Modes[65] = conjugate(Symm - Asymm)
            # (ell, m) = (8, +/- 4)
            Symm = 0
            Asymm = 0
            Modes[72] = Symm + Asymm
            Modes[64] = conjugate(Symm - Asymm)
            # (ell, m) = (8, +/- 5)
            Symm = 0
            Asymm = 0
            Modes[73] = Symm + Asymm
            Modes[63] = conjugate(Symm - Asymm)
            # (ell, m) = (8, +/- 6)
            Symm = 0
            Asymm = 0
            Modes[74] = Symm + Asymm
            Modes[62] = conjugate(Symm - Asymm)
            # (ell, m) = (8, +/- 7)
            Symm = 0
            Asymm = 0
            Modes[75] = Symm + Asymm
            Modes[61] = conjugate(Symm - Asymm)
            # (ell, m) = (8, +/- 8)
            Symm = 0
            Asymm = 0
            Modes[76] = Symm + Asymm
            Modes[60] = conjugate(Symm - Asymm)

            return Modes
            
        for i in range(len(y[0])):
            v = y[0,i]
            rfrak_chi1_x = y[1,i]
            rfrak_chi1_y = y[2,i]
            rfrak_chi2_x = y[3,i]
            rfrak_chi2_y = y[4,i]
            rfrak_frame_x = y[5,i]
            rfrak_frame_y = y[6,i]
            rfrak_frame_z = y[7,i]
            Phi = y[8,i]
            xHat=xHat_i
            yHat=yHat_i
            zHat=zHat_i
            M1=M1_i
            M2=M2_i 
            S_chi1=S_chi1_i
            S_chi2=S_chi2_i
            M=M1 + M2
            delta=(M1 - M2)/M
            nu=M1*M2/M**2
            R=exp(rfrak_frame_x*xHat + rfrak_frame_y*yHat + rfrak_frame_z*zHat)
            nHat=R*xHat*conjugate(R)
            lambdaHat=R*yHat*conjugate(R)
            ellHat=R*zHat*conjugate(R)
            R_S1=exp(rfrak_chi1_x*xHat + rfrak_chi1_y*yHat)
            R_S2=exp(rfrak_chi2_x*xHat + rfrak_chi2_y*yHat)
            chiVec1=S_chi1*R_S1*zHat*conjugate(R_S1)*conjugate(S_chi1)
            chiVec2=S_chi2*R_S2*zHat*conjugate(R_S2)*conjugate(S_chi2)
            v_i=v
            rfrak_frame_x_i=rfrak_frame_x
            rfrak_frame_y_i=rfrak_frame_y
            rfrak_frame_z_i=rfrak_frame_z
            chiVec1_i=quaternion.as_float_array(chiVec1)
            chiVec2_i=quaternion.as_float_array(chiVec2)
            ModeData[i,:]=WaveformModes(M1_i, M2_i, v_i, chiVec1_i, chiVec2_i)
        return ModeData

class WaveformModes_2p5PN :
    def WaveformModes_2p5PN(xHat_i, yHat_i, zHat_i, M1_i, M2_i, v_i, S_chi1_i, S_chi2_i, rfrak_chi1_x_i, rfrak_chi1_y_i, rfrak_chi2_x_i, rfrak_chi2_y_i, rfrak_frame_x_i, rfrak_frame_y_i, rfrak_frame_z_i,y) :            
        ModeData=np.empty((len(y[0]),77), dtype=complex)
        def WaveformModes(M1_i, M2_i, v_i, chiVec1_i, chiVec2_i):
            I=1j
            global chi1_n
            global chi1_lambda
            global chi1_ell
            global chi2_n
            global chi2_lambda
            global chi2_ell
            global S_ell
            global S_n
            global S_lambda
            global Sigma_ell
            global Sigma_n
            global Sigma_lambda
            global S1_ell
            global S1_n
            global S1_lambda
            global S2_ell
            global S2_n
            global S2_lambda
            global rhOverM_coeff
            global hHat_spin_Symm_2_2_3
            global hHat_spin_Symm_2_2_4
            global hHat_spin_Symm_2_1_2
            global hHat_spin_Symm_2_1_4
            global hHat_spin_Symm_2_0_4
            global hHat_spin_Symm_3_3_4
            global hHat_spin_Symm_3_2_3
            global hHat_spin_Symm_3_1_4
            global hHat_spin_Symm_4_3_4
            global hHat_spin_Symm_4_1_4
            global hHat_spin_Asymm_2_2_2
            global hHat_spin_Asymm_2_2_4
            global hHat_spin_Asymm_2_1_3
            global hHat_spin_Asymm_2_1_4
            global hHat_spin_Asymm_2_0_2
            global hHat_spin_Asymm_2_0_4
            global hHat_spin_Asymm_3_3_3
            global hHat_spin_Asymm_3_2_4
            global hHat_spin_Asymm_3_1_3
            global hHat_spin_Asymm_3_0_4
            global hHat_spin_Asymm_4_4_4
            global hHat_spin_Asymm_4_2_4
            global hHat_spin_Asymm_4_0_4
            M1=M1_i
            M2=M2_i
            v=v_i
            M=M1 + M2
            delta=(M1 - M2)/M
            nu=M1*M2/M**2
            chiVec1=chiVec1_i
            chiVec2=chiVec2_i
            chi1_n=chiVec1[1]
            chi1_lambda=chiVec1[2]
            chi1_ell=chiVec1[3]
            chi2_n=chiVec2[1]
            chi2_lambda=chiVec2[2]
            chi2_ell=chiVec2[3]
            S_ell=M1**2*chi1_ell + M2**2*chi2_ell
            S_n=M1**2*chi1_n + M2**2*chi2_n
            S_lambda=M1**2*chi1_lambda + M2**2*chi2_lambda
            Sigma_ell=M*(-M1*chi1_ell + M2*chi2_ell)
            Sigma_n=M*(-M1*chi1_n + M2*chi2_n)
            Sigma_lambda=M*(-M1*chi1_lambda + M2*chi2_lambda)
            S1_ell=M1**2*chi1_ell
            S1_n=M1**2*chi1_n
            S1_lambda=M1**2*chi1_lambda
            S2_ell=M2**2*chi2_ell
            S2_n=M2**2*chi2_n
            S2_lambda=M2**2*chi2_lambda
            rhOverM_coeff=8*sqrt(5)*nu*v**2*sqrt(pi)/5
            hHat_2_0_0=-5*sqrt(6)/84
            hHat_2_1_1=I*delta/3
            hHat_2_1_3=I*delta*(20*nu - 17)/84
            hHat_2_1_4=delta*(1 + log(16) + 2*I*pi)/6
            hHat_2_1_5=I*delta*(nu*(237*nu - 2036) - 172)/1512
            hHat_2_2_0=1
            hHat_2_2_2=55*nu/42 - 107/42
            hHat_2_2_3=2*pi
            hHat_2_2_4=nu*(2047*nu - 7483)/1512 - 2173/1512
            hHat_2_2_5=-24*I*nu + pi*(34*nu - 107)/21
            hHat_3_0_5=-2*sqrt(42)*I*nu/35
            hHat_3_1_1=sqrt(14)*I*delta/168
            hHat_3_1_3=-sqrt(14)*I*delta*(nu + 4)/252
            hHat_3_1_4=sqrt(14)*delta*(log(1024) + 7 + 5*I*pi)/840
            hHat_3_1_5=-sqrt(14)*I*delta*(nu*(247*nu + 272) - 607)/33264
            hHat_3_2_2=sqrt(35)*(1 - 3*nu)/21
            hHat_3_2_4=sqrt(35)*(nu*(725 - 365*nu) - 193)/1890
            hHat_3_2_5=sqrt(35)*(-30*nu*pi + 66*I*nu + 10*pi - 15*I)/105
            hHat_3_3_1=-3*sqrt(210)*I*delta/56
            hHat_3_3_3=-3*sqrt(210)*I*delta*(nu - 2)/28
            hHat_3_3_4=9*sqrt(210)*delta*(-7 + 10*log(3/2) - 5*I*pi)/280
            hHat_3_3_5=-sqrt(210)*I*delta*(nu*(887*nu - 3676) + 369)/6160
            hHat_4_0_0=-sqrt(2)/1008
            hHat_4_1_3=sqrt(10)*I*delta*(1 - 2*nu)/840
            hHat_4_1_5=-sqrt(10)*I*delta*(nu*(332*nu - 1011) + 404)/110880
            hHat_4_2_2=sqrt(5)*(1 - 3*nu)/63
            hHat_4_2_4=sqrt(5)*(nu*(4025 - 285*nu) - 1311)/20790
            hHat_4_2_5=sqrt(5)*(nu*(-30*pi + 84*I) + 10*pi - 21*I)/315
            hHat_4_3_3=9*sqrt(70)*I*delta*(2*nu - 1)/280
            hHat_4_3_5=3*sqrt(70)*I*delta*(nu*(524*nu - 1267) + 468)/12320
            hHat_4_4_2=8*sqrt(35)*(3*nu - 1)/63
            hHat_4_4_4=sqrt(35)*(20*nu*(525*nu - 1273) + 7116)/10395
            hHat_4_4_5=sqrt(35)*(pi*(480*nu - 160) + I*(-1193*nu + (960*nu - 320)*log(2) + 336))/315
            hHat_5_1_3=sqrt(385)*I*delta*(1 - 2*nu)/110880
            hHat_5_1_5=-sqrt(385)*I*delta*(nu*(4*nu - 352) + 179)/4324320
            hHat_5_2_4=sqrt(55)*(2*nu*(5*nu - 5) + 2)/1485
            hHat_5_3_3=9*sqrt(330)*I*delta*(2*nu - 1)/3520
            hHat_5_3_5=3*sqrt(330)*I*delta*(8*nu*(11*nu - 58) + 207)/45760
            hHat_5_4_4=sqrt(165)*(-32*nu*(5*nu - 5) - 32)/1485
            hHat_5_5_3=-625*sqrt(66)*I*delta*(2*nu - 1)/6336
            hHat_5_5_5=-625*sqrt(66)*I*delta*(16*nu*(16*nu - 43) + 263)/247104
            hHat_6_1_5=sqrt(26)*I*delta*(nu - 1)*(3*nu - 1)/216216
            hHat_6_2_4=sqrt(65)*(2*nu*(5*nu - 5) + 2)/19305
            hHat_6_3_5=-81*sqrt(65)*I*delta*(nu - 1)*(3*nu - 1)/40040
            hHat_6_4_4=-128*sqrt(78)*(nu*(5*nu - 5) + 1)/19305
            hHat_6_5_5=3125*sqrt(429)*I*delta*(nu - 1)*(3*nu - 1)/216216
            hHat_6_6_4=sqrt(143)*(54*nu*(5*nu - 5) + 54)/715
            hHat_7_1_5=sqrt(2)*I*delta*(nu - 1)*(3*nu - 1)/1729728
            hHat_7_3_5=-243*sqrt(6)*I*delta*(nu - 1)*(3*nu - 1)/320320
            hHat_7_5_5=15625*sqrt(66)*I*delta*(nu - 1)*(3*nu - 1)/1729728
            hHat_7_7_5=-16807*sqrt(6006)*I*delta*(nu - 1)*(3*nu - 1)/1235520
            hHat_spin_Symm_2_2_3=(-6*S_ell - 2*Sigma_ell*delta)/(3*M**2)
            hHat_spin_Symm_2_2_4=(12*S1_ell*S2_ell + 10*S1_lambda*S2_lambda - 15*I*S1_lambda*S2_n - 15*I*S1_n*S2_lambda - 22*S1_n*S2_n)/(6*M**4*nu)
            hHat_spin_Symm_2_1_2=I*Sigma_ell/(2*M**2)
            hHat_spin_Symm_2_1_4=I*(-86*S_ell*delta + Sigma_ell*(139*nu - 79))/(42*M**2)
            hHat_spin_Symm_2_0_4=sqrt(6)*(-S1_lambda*S2_lambda + S1_n*S2_n)/(3*M**4*nu)
            hHat_spin_Symm_3_3_4=3*sqrt(210)*I*(7*S_ell*delta - 3*Sigma_ell*(3*nu - 1))/(112*M**2)
            hHat_spin_Symm_3_2_3=2*sqrt(35)*(S_ell + Sigma_ell*delta)/(21*M**2)
            hHat_spin_Symm_3_1_4=sqrt(14)*I*(S_ell*delta - 5*Sigma_ell*(3*nu - 1))/(336*M**2)
            hHat_spin_Symm_4_3_4=9*sqrt(70)*I*(-S_ell*delta + 3*Sigma_ell*nu - Sigma_ell)/(112*M**2)
            hHat_spin_Symm_4_1_4=sqrt(10)*I*(S_ell*delta - 3*Sigma_ell*nu + Sigma_ell)/(336*M**2)
            hHat_spin_Asymm_2_2_2=(-Sigma_lambda - I*Sigma_n)/(2*M**2)
            hHat_spin_Asymm_2_2_4=(19*S_lambda*delta + 182*I*S_n*delta - 43*Sigma_lambda*nu + 5*Sigma_lambda - 280*I*Sigma_n*nu + 98*I*Sigma_n)/(84*M**2)
            hHat_spin_Asymm_2_1_3=(4*I*S_lambda + 25*S_n + 4*I*Sigma_lambda*delta + 13*Sigma_n*delta)/(6*M**2)
            hHat_spin_Asymm_2_1_4=(-3*S1_ell*S2_n - 3*S1_n*S2_ell)/(2*M**4*nu)
            hHat_spin_Asymm_2_0_2=sqrt(6)*I*Sigma_n/(6*M**2)
            hHat_spin_Asymm_2_0_4=sqrt(6)*I*(255*S_n*delta - Sigma_n*(506*nu - 45))/(126*M**2)
            hHat_spin_Asymm_3_3_3=sqrt(210)*I*(S_lambda + I*S_n + delta*(Sigma_lambda + I*Sigma_n))/(21*M**2)
            hHat_spin_Asymm_3_2_4=sqrt(35)*(-Sigma_lambda*(83*nu - 17) + 4*I*Sigma_n*(55*nu - 13) + 25*delta*(S_lambda - 4*I*S_n))/(168*M**2)
            hHat_spin_Asymm_3_1_3=sqrt(14)*(I*S_lambda + S_n + delta*(I*Sigma_lambda + Sigma_n))/(21*M**2)
            hHat_spin_Asymm_3_0_4=sqrt(42)*(-17*S_lambda*delta + Sigma_lambda*(35*nu - 9))/(168*M**2)
            hHat_spin_Asymm_4_4_4=9*sqrt(35)*(-3*Sigma_lambda*nu + Sigma_lambda - I*Sigma_n*(3*nu - 1) + delta*(S_lambda + I*S_n))/(56*M**2)
            hHat_spin_Asymm_4_2_4=sqrt(5)*(-13*Sigma_lambda*(3*nu - 1) + 14*I*Sigma_n*(3*nu - 1) + delta*(13*S_lambda - 14*I*S_n))/(168*M**2)
            hHat_spin_Asymm_4_0_4=sqrt(2)*I*(S_n*delta - 3*Sigma_n*nu + Sigma_n)/(168*M**2)
            
            chiVec1 = [0., chi1_n, chi1_lambda, chi1_ell]
            chiVec2 = [0., chi2_n, chi2_lambda, chi2_ell]


            chi1_n = chiVec1[1]
            chi1_lambda = chiVec1[2]
            chi1_ell = chiVec1[3]
            chi2_n = chiVec2[1]
            chi2_lambda = chiVec2[2]
            chi2_ell = chiVec2[3]
            S_ell = M1**2*chi1_ell + M2**2*chi2_ell
            S_n = M1**2*chi1_n + M2**2*chi2_n
            S_lambda = M1**2*chi1_lambda + M2**2*chi2_lambda
            Sigma_ell = M*(-M1*chi1_ell + M2*chi2_ell)
            Sigma_n = M*(-M1*chi1_n + M2*chi2_n)
            Sigma_lambda = M*(-M1*chi1_lambda + M2*chi2_lambda)
            S1_ell = M1**2*chi1_ell
            S1_n = M1**2*chi1_n
            S1_lambda = M1**2*chi1_lambda
            S2_ell = M2**2*chi2_ell
            S2_n = M2**2*chi2_n
            S2_lambda = M2**2*chi2_lambda
            rhOverM_coeff = 8*sqrt(5)*nu*v**2*sqrt(pi)/5
            hHat_spin_Symm_2_2_3 = (-6*S_ell - 2*Sigma_ell*delta)/(3*M**2)
            hHat_spin_Symm_2_2_4 = (12*S1_ell*S2_ell + 10*S1_lambda*S2_lambda - 15*I*S1_lambda*S2_n - 15*I*S1_n*S2_lambda - 22*S1_n*S2_n)/(6*M**4*nu)
            hHat_spin_Symm_2_1_2 = I*Sigma_ell/(2*M**2)
            hHat_spin_Symm_2_1_4 = I*(-86*S_ell*delta + Sigma_ell*(139*nu - 79))/(42*M**2)
            hHat_spin_Symm_2_0_4 = sqrt(6)*(-S1_lambda*S2_lambda + S1_n*S2_n)/(3*M**4*nu)
            hHat_spin_Symm_3_3_4 = 3*sqrt(210)*I*(7*S_ell*delta - 3*Sigma_ell*(3*nu - 1))/(112*M**2)
            hHat_spin_Symm_3_2_3 = 2*sqrt(35)*(S_ell + Sigma_ell*delta)/(21*M**2)
            hHat_spin_Symm_3_1_4 = sqrt(14)*I*(S_ell*delta - 5*Sigma_ell*(3*nu - 1))/(336*M**2)
            hHat_spin_Symm_4_3_4 = 9*sqrt(70)*I*(-S_ell*delta + 3*Sigma_ell*nu - Sigma_ell)/(112*M**2)
            hHat_spin_Symm_4_1_4 = sqrt(10)*I*(S_ell*delta - 3*Sigma_ell*nu + Sigma_ell)/(336*M**2)
            hHat_spin_Asymm_2_2_2 = (-Sigma_lambda - I*Sigma_n)/(2*M**2)
            hHat_spin_Asymm_2_2_4 = (19*S_lambda*delta + 182*I*S_n*delta - 43*Sigma_lambda*nu + 5*Sigma_lambda - 280*I*Sigma_n*nu + 98*I*Sigma_n)/(84*M**2)
            hHat_spin_Asymm_2_1_3 = (4*I*S_lambda + 25*S_n + 4*I*Sigma_lambda*delta + 13*Sigma_n*delta)/(6*M**2)
            hHat_spin_Asymm_2_1_4 = (-3*S1_ell*S2_n - 3*S1_n*S2_ell)/(2*M**4*nu)
            hHat_spin_Asymm_2_0_2 = sqrt(6)*I*Sigma_n/(6*M**2)
            hHat_spin_Asymm_2_0_4 = sqrt(6)*I*(255*S_n*delta - Sigma_n*(506*nu - 45))/(126*M**2)
            hHat_spin_Asymm_3_3_3 = sqrt(210)*I*(S_lambda + I*S_n + delta*(Sigma_lambda + I*Sigma_n))/(21*M**2)
            hHat_spin_Asymm_3_2_4 = sqrt(35)*(-Sigma_lambda*(83*nu - 17) + 4*I*Sigma_n*(55*nu - 13) + 25*delta*(S_lambda - 4*I*S_n))/(168*M**2)
            hHat_spin_Asymm_3_1_3 = sqrt(14)*(I*S_lambda + S_n + delta*(I*Sigma_lambda + Sigma_n))/(21*M**2)
            hHat_spin_Asymm_3_0_4 = sqrt(42)*(-17*S_lambda*delta + Sigma_lambda*(35*nu - 9))/(168*M**2)
            hHat_spin_Asymm_4_4_4 = 9*sqrt(35)*(-3*Sigma_lambda*nu + Sigma_lambda - I*Sigma_n*(3*nu - 1) + delta*(S_lambda + I*S_n))/(56*M**2)
            hHat_spin_Asymm_4_2_4 = sqrt(5)*(-13*Sigma_lambda*(3*nu - 1) + 14*I*Sigma_n*(3*nu - 1) + delta*(13*S_lambda - 14*I*S_n))/(168*M**2)
            hHat_spin_Asymm_4_0_4 = sqrt(2)*I*(S_n*delta - 3*Sigma_n*nu + Sigma_n)/(168*M**2)

            Modes=np.empty(77, dtype=complex)
            # (ell, m) = (2, +/- 0)
            Symm = rhOverM_coeff*(hHat_2_0_0 + hHat_spin_Symm_2_0_4*v**4)
            Asymm = rhOverM_coeff*v**2*(hHat_spin_Asymm_2_0_2 + hHat_spin_Asymm_2_0_4*v**2)
            Modes[2] = Symm + Asymm
            # (ell, m) = (2, +/- 1)
            Symm = rhOverM_coeff*v*(hHat_2_1_1 + v*(hHat_spin_Symm_2_1_2 + v*(hHat_2_1_3 + v*(hHat_2_1_4 + hHat_2_1_5*v + hHat_spin_Symm_2_1_4))))
            Asymm = rhOverM_coeff*v**3*(hHat_spin_Asymm_2_1_3 + hHat_spin_Asymm_2_1_4*v)
            Modes[3] = Symm + Asymm
            Modes[1] = conjugate(Symm - Asymm)
            # (ell, m) = (2, +/- 2)
            Symm = rhOverM_coeff*(hHat_2_2_0 + v**2*(hHat_2_2_2 + v*(hHat_2_2_3 + hHat_spin_Symm_2_2_3 + v*(hHat_2_2_4 + hHat_2_2_5*v + hHat_spin_Symm_2_2_4))))
            Asymm = rhOverM_coeff*v**2*(hHat_spin_Asymm_2_2_2 + hHat_spin_Asymm_2_2_4*v**2)
            Modes[4] = Symm + Asymm
            Modes[0] = conjugate(Symm - Asymm)
            # (ell, m) = (3, +/- 0)
            Symm = hHat_3_0_5*rhOverM_coeff*v**5
            Asymm = hHat_spin_Asymm_3_0_4*rhOverM_coeff*v**4
            Modes[8] = Symm + Asymm
            # (ell, m) = (3, +/- 1)
            Symm = rhOverM_coeff*v*(hHat_3_1_1 + v**2*(hHat_3_1_3 + v*(hHat_3_1_4 + hHat_3_1_5*v + hHat_spin_Symm_3_1_4)))
            Asymm = hHat_spin_Asymm_3_1_3*rhOverM_coeff*v**3
            Modes[9] = Symm + Asymm
            Modes[7] = -conjugate(Symm - Asymm)
            # (ell, m) = (3, +/- 2)
            Symm = rhOverM_coeff*v**2*(hHat_3_2_2 + v*(hHat_spin_Symm_3_2_3 + v*(hHat_3_2_4 + hHat_3_2_5*v)))
            Asymm = hHat_spin_Asymm_3_2_4*rhOverM_coeff*v**4
            Modes[10] = Symm + Asymm
            Modes[6] = -conjugate(Symm - Asymm)
            # (ell, m) = (3, +/- 3)
            Symm = rhOverM_coeff*v*(hHat_3_3_1 + v**2*(hHat_3_3_3 + v*(hHat_3_3_4 + hHat_3_3_5*v + hHat_spin_Symm_3_3_4)))
            Asymm = hHat_spin_Asymm_3_3_3*rhOverM_coeff*v**3
            Modes[11] = Symm + Asymm
            Modes[5] = -conjugate(Symm - Asymm)
            # (ell, m) = (4, +/- 0)
            Symm = hHat_4_0_0*rhOverM_coeff
            Asymm = hHat_spin_Asymm_4_0_4*rhOverM_coeff*v**4
            Modes[16] = Symm + Asymm
            # (ell, m) = (4, +/- 1)
            Symm = rhOverM_coeff*v**3*(hHat_4_1_3 + v*(hHat_4_1_5*v + hHat_spin_Symm_4_1_4))
            Asymm = 0
            Modes[17] = Symm + Asymm
            Modes[15] = conjugate(Symm - Asymm)
            # (ell, m) = (4, +/- 2)
            Symm = rhOverM_coeff*v**2*(hHat_4_2_2 + v**2*(hHat_4_2_4 + hHat_4_2_5*v))
            Asymm = hHat_spin_Asymm_4_2_4*rhOverM_coeff*v**4
            Modes[18] = Symm + Asymm
            Modes[14] = conjugate(Symm - Asymm)
            # (ell, m) = (4, +/- 3)
            Symm = rhOverM_coeff*v**3*(hHat_4_3_3 + v*(hHat_4_3_5*v + hHat_spin_Symm_4_3_4))
            Asymm = 0
            Modes[19] = Symm + Asymm
            Modes[13] = conjugate(Symm - Asymm)
            # (ell, m) = (4, +/- 4)
            Symm = rhOverM_coeff*v**2*(hHat_4_4_2 + v**2*(hHat_4_4_4 + hHat_4_4_5*v))
            Asymm = hHat_spin_Asymm_4_4_4*rhOverM_coeff*v**4
            Modes[20] = Symm + Asymm
            Modes[12] = conjugate(Symm - Asymm)
            # (ell, m) = (5, +/- 0)
            Symm = 0
            Asymm = 0
            Modes[26] = Symm + Asymm
            # (ell, m) = (5, +/- 1)
            Symm = rhOverM_coeff*v**3*(hHat_5_1_3 + hHat_5_1_5*v**2)
            Asymm = 0
            Modes[27] = Symm + Asymm
            Modes[25] = -conjugate(Symm - Asymm)
            # (ell, m) = (5, +/- 2)
            Symm = hHat_5_2_4*rhOverM_coeff*v**4
            Asymm = 0
            Modes[28] = Symm + Asymm
            Modes[24] = -conjugate(Symm - Asymm)
            # (ell, m) = (5, +/- 3)
            Symm = rhOverM_coeff*v**3*(hHat_5_3_3 + hHat_5_3_5*v**2)
            Asymm = 0
            Modes[29] = Symm + Asymm
            Modes[23] = -conjugate(Symm - Asymm)
            # (ell, m) = (5, +/- 4)
            Symm = hHat_5_4_4*rhOverM_coeff*v**4
            Asymm = 0
            Modes[30] = Symm + Asymm
            Modes[22] = -conjugate(Symm - Asymm)
            # (ell, m) = (5, +/- 5)
            Symm = rhOverM_coeff*v**3*(hHat_5_5_3 + hHat_5_5_5*v**2)
            Asymm = 0
            Modes[31] = Symm + Asymm
            Modes[21] = -conjugate(Symm - Asymm)
            # (ell, m) = (6, +/- 0)
            Symm = 0
            Asymm = 0
            Modes[38] = Symm + Asymm
            # (ell, m) = (6, +/- 1)
            Symm = hHat_6_1_5*rhOverM_coeff*v**5
            Asymm = 0
            Modes[39] = Symm + Asymm
            Modes[37] = conjugate(Symm - Asymm)
            # (ell, m) = (6, +/- 2)
            Symm = hHat_6_2_4*rhOverM_coeff*v**4
            Asymm = 0
            Modes[40] = Symm + Asymm
            Modes[36] = conjugate(Symm - Asymm)
            # (ell, m) = (6, +/- 3)
            Symm = hHat_6_3_5*rhOverM_coeff*v**5
            Asymm = 0
            Modes[41] = Symm + Asymm
            Modes[35] = conjugate(Symm - Asymm)
            # (ell, m) = (6, +/- 4)
            Symm = hHat_6_4_4*rhOverM_coeff*v**4
            Asymm = 0
            Modes[42] = Symm + Asymm
            Modes[34] = conjugate(Symm - Asymm)
            # (ell, m) = (6, +/- 5)
            Symm = hHat_6_5_5*rhOverM_coeff*v**5
            Asymm = 0
            Modes[43] = Symm + Asymm
            Modes[33] = conjugate(Symm - Asymm)
            # (ell, m) = (6, +/- 6)
            Symm = hHat_6_6_4*rhOverM_coeff*v**4
            Asymm = 0
            Modes[44] = Symm + Asymm
            Modes[32] = conjugate(Symm - Asymm)
            # (ell, m) = (7, +/- 0)
            Symm = 0
            Asymm = 0
            Modes[52] = Symm + Asymm
            # (ell, m) = (7, +/- 1)
            Symm = hHat_7_1_5*rhOverM_coeff*v**5
            Asymm = 0
            Modes[53] = Symm + Asymm
            Modes[51] = -conjugate(Symm - Asymm)
            # (ell, m) = (7, +/- 2)
            Symm = 0
            Asymm = 0
            Modes[54] = Symm + Asymm
            Modes[50] = -conjugate(Symm - Asymm)
            # (ell, m) = (7, +/- 3)
            Symm = hHat_7_3_5*rhOverM_coeff*v**5
            Asymm = 0
            Modes[55] = Symm + Asymm
            Modes[49] = -conjugate(Symm - Asymm)
            # (ell, m) = (7, +/- 4)
            Symm = 0
            Asymm = 0
            Modes[56] = Symm + Asymm
            Modes[48] = -conjugate(Symm - Asymm)
            # (ell, m) = (7, +/- 5)
            Symm = hHat_7_5_5*rhOverM_coeff*v**5
            Asymm = 0
            Modes[57] = Symm + Asymm
            Modes[47] = -conjugate(Symm - Asymm)
            # (ell, m) = (7, +/- 6)
            Symm = 0
            Asymm = 0
            Modes[58] = Symm + Asymm
            Modes[46] = -conjugate(Symm - Asymm)
            # (ell, m) = (7, +/- 7)
            Symm = hHat_7_7_5*rhOverM_coeff*v**5
            Asymm = 0
            Modes[59] = Symm + Asymm
            Modes[45] = -conjugate(Symm - Asymm)
            # (ell, m) = (8, +/- 0)
            Symm = 0
            Asymm = 0
            Modes[68] = Symm + Asymm
            # (ell, m) = (8, +/- 1)
            Symm = 0
            Asymm = 0
            Modes[69] = Symm + Asymm
            Modes[67] = conjugate(Symm - Asymm)
            # (ell, m) = (8, +/- 2)
            Symm = 0
            Asymm = 0
            Modes[70] = Symm + Asymm
            Modes[66] = conjugate(Symm - Asymm)
            # (ell, m) = (8, +/- 3)
            Symm = 0
            Asymm = 0
            Modes[71] = Symm + Asymm
            Modes[65] = conjugate(Symm - Asymm)
            # (ell, m) = (8, +/- 4)
            Symm = 0
            Asymm = 0
            Modes[72] = Symm + Asymm
            Modes[64] = conjugate(Symm - Asymm)
            # (ell, m) = (8, +/- 5)
            Symm = 0
            Asymm = 0
            Modes[73] = Symm + Asymm
            Modes[63] = conjugate(Symm - Asymm)
            # (ell, m) = (8, +/- 6)
            Symm = 0
            Asymm = 0
            Modes[74] = Symm + Asymm
            Modes[62] = conjugate(Symm - Asymm)
            # (ell, m) = (8, +/- 7)
            Symm = 0
            Asymm = 0
            Modes[75] = Symm + Asymm
            Modes[61] = conjugate(Symm - Asymm)
            # (ell, m) = (8, +/- 8)
            Symm = 0
            Asymm = 0
            Modes[76] = Symm + Asymm
            Modes[60] = conjugate(Symm - Asymm)

            return Modes
            
        for i in range(len(y[0])):
            v = y[0,i]
            rfrak_chi1_x = y[1,i]
            rfrak_chi1_y = y[2,i]
            rfrak_chi2_x = y[3,i]
            rfrak_chi2_y = y[4,i]
            rfrak_frame_x = y[5,i]
            rfrak_frame_y = y[6,i]
            rfrak_frame_z = y[7,i]
            Phi = y[8,i]
            xHat=xHat_i
            yHat=yHat_i
            zHat=zHat_i
            M1=M1_i
            M2=M2_i 
            S_chi1=S_chi1_i
            S_chi2=S_chi2_i
            M=M1 + M2
            delta=(M1 - M2)/M
            nu=M1*M2/M**2
            R=exp(rfrak_frame_x*xHat + rfrak_frame_y*yHat + rfrak_frame_z*zHat)
            nHat=R*xHat*conjugate(R)
            lambdaHat=R*yHat*conjugate(R)
            ellHat=R*zHat*conjugate(R)
            R_S1=exp(rfrak_chi1_x*xHat + rfrak_chi1_y*yHat)
            R_S2=exp(rfrak_chi2_x*xHat + rfrak_chi2_y*yHat)
            chiVec1=S_chi1*R_S1*zHat*conjugate(R_S1)*conjugate(S_chi1)
            chiVec2=S_chi2*R_S2*zHat*conjugate(R_S2)*conjugate(S_chi2)
            v_i=v
            rfrak_frame_x_i=rfrak_frame_x
            rfrak_frame_y_i=rfrak_frame_y
            rfrak_frame_z_i=rfrak_frame_z
            chiVec1_i=quaternion.as_float_array(chiVec1)
            chiVec2_i=quaternion.as_float_array(chiVec2)
            ModeData[i,:]=WaveformModes(M1_i, M2_i, v_i, chiVec1_i, chiVec2_i)
        return ModeData

class WaveformModes_3p0PN :
    def WaveformModes_3p0PN(xHat_i, yHat_i, zHat_i, M1_i, M2_i, v_i, S_chi1_i, S_chi2_i, rfrak_chi1_x_i, rfrak_chi1_y_i, rfrak_chi2_x_i, rfrak_chi2_y_i, rfrak_frame_x_i, rfrak_frame_y_i, rfrak_frame_z_i,y) :            
        ModeData=np.empty((len(y[0]),77), dtype=complex)
        def WaveformModes(M1_i, M2_i, v_i, chiVec1_i, chiVec2_i):
            I=1j
            global chi1_n
            global chi1_lambda
            global chi1_ell
            global chi2_n
            global chi2_lambda
            global chi2_ell
            global S_ell
            global S_n
            global S_lambda
            global Sigma_ell
            global Sigma_n
            global Sigma_lambda
            global S1_ell
            global S1_n
            global S1_lambda
            global S2_ell
            global S2_n
            global S2_lambda
            global logv
            global rhOverM_coeff
            global hHat_spin_Symm_2_2_3
            global hHat_spin_Symm_2_2_4
            global hHat_spin_Symm_2_1_2
            global hHat_spin_Symm_2_1_4
            global hHat_spin_Symm_2_0_4
            global hHat_spin_Symm_3_3_4
            global hHat_spin_Symm_3_2_3
            global hHat_spin_Symm_3_1_4
            global hHat_spin_Symm_4_3_4
            global hHat_spin_Symm_4_1_4
            global hHat_spin_Asymm_2_2_2
            global hHat_spin_Asymm_2_2_4
            global hHat_spin_Asymm_2_1_3
            global hHat_spin_Asymm_2_1_4
            global hHat_spin_Asymm_2_0_2
            global hHat_spin_Asymm_2_0_4
            global hHat_spin_Asymm_3_3_3
            global hHat_spin_Asymm_3_2_4
            global hHat_spin_Asymm_3_1_3
            global hHat_spin_Asymm_3_0_4
            global hHat_spin_Asymm_4_4_4
            global hHat_spin_Asymm_4_2_4
            global hHat_spin_Asymm_4_0_4
            M1=M1_i
            M2=M2_i
            v=v_i
            M=M1 + M2
            delta=(M1 - M2)/M
            nu=M1*M2/M**2
            chiVec1=chiVec1_i
            chiVec2=chiVec2_i
            chi1_n=chiVec1[1]
            chi1_lambda=chiVec1[2]
            chi1_ell=chiVec1[3]
            chi2_n=chiVec2[1]
            chi2_lambda=chiVec2[2]
            chi2_ell=chiVec2[3]
            S_ell=M1**2*chi1_ell + M2**2*chi2_ell
            S_n=M1**2*chi1_n + M2**2*chi2_n
            S_lambda=M1**2*chi1_lambda + M2**2*chi2_lambda
            Sigma_ell=M*(-M1*chi1_ell + M2*chi2_ell)
            Sigma_n=M*(-M1*chi1_n + M2*chi2_n)
            Sigma_lambda=M*(-M1*chi1_lambda + M2*chi2_lambda)
            S1_ell=M1**2*chi1_ell
            S1_n=M1**2*chi1_n
            S1_lambda=M1**2*chi1_lambda
            S2_ell=M2**2*chi2_ell
            S2_n=M2**2*chi2_n
            S2_lambda=M2**2*chi2_lambda
            logv=log(v)
            rhOverM_coeff=8*sqrt(5)*nu*v**2*sqrt(pi)/5
            hHat_2_0_0=-5*sqrt(6)/84
            hHat_2_1_1=I*delta/3
            hHat_2_1_3=I*delta*(20*nu - 17)/84
            hHat_2_1_4=delta*(1 + log(16) + 2*I*pi)/6
            hHat_2_1_5=I*delta*(nu*(237*nu - 2036) - 172)/1512
            hHat_2_1_6=delta*(2*nu*(log(4096) + 353 + 6*I*pi) - 17*log(16) - 17 - 34*I*pi)/168
            hHat_2_2_0=1
            hHat_2_2_2=55*nu/42 - 107/42
            hHat_2_2_3=2*pi
            hHat_2_2_4=nu*(2047*nu - 7483)/1512 - 2173/1512
            hHat_2_2_5=-24*I*nu + pi*(34*nu - 107)/21
            hHat_2_2_6=nu*(nu*(114635*nu - 729396) - 834555)/99792 + 41*nu*pi**2/96 - 428*log(16)/105 - 856*EulerGamma/105 + 27027409/646800 + 2*pi*(35*pi + 214*I)/105
            hHat_2_2_lnv_6=-856/105
            hHat_3_0_5=-2*sqrt(42)*I*nu/35
            hHat_3_1_1=sqrt(14)*I*delta/168
            hHat_3_1_3=-sqrt(14)*I*delta*(nu + 4)/252
            hHat_3_1_4=sqrt(14)*delta*(log(1024) + 7 + 5*I*pi)/840
            hHat_3_1_5=-sqrt(14)*I*delta*(nu*(247*nu + 272) - 607)/33264
            hHat_3_1_6=sqrt(14)*delta*(2*nu - 5*I*pi*(7*nu + 16) - 2*(35*nu + 80)*log(2) - 112)/5040
            hHat_3_2_2=sqrt(35)*(1 - 3*nu)/21
            hHat_3_2_4=sqrt(35)*(nu*(725 - 365*nu) - 193)/1890
            hHat_3_2_5=sqrt(35)*(-30*nu*pi + 66*I*nu + 10*pi - 15*I)/105
            hHat_3_2_6=sqrt(35)*(nu*(nu*(100026 - 16023*nu) - 17387) - 1451)/83160
            hHat_3_3_1=-3*sqrt(210)*I*delta/56
            hHat_3_3_3=-3*sqrt(210)*I*delta*(nu - 2)/28
            hHat_3_3_4=9*sqrt(210)*delta*(-7 + 10*log(3/2) - 5*I*pi)/280
            hHat_3_3_5=-sqrt(210)*I*delta*(nu*(887*nu - 3676) + 369)/6160
            hHat_3_3_6=sqrt(210)*delta*(-96206*nu - 3645*I*pi*(3*nu - 8) + 2*(10935*nu - 29160)*log(3/2) + 40824)/45360
            hHat_4_0_0=-sqrt(2)/1008
            hHat_4_1_3=sqrt(10)*I*delta*(1 - 2*nu)/840
            hHat_4_1_5=-sqrt(10)*I*delta*(nu*(332*nu - 1011) + 404)/110880
            hHat_4_1_6=sqrt(10)*delta*(-1661*nu + (60 - 120*nu)*log(2) - 30*I*pi*(2*nu - 1) + 64)/25200
            hHat_4_2_2=sqrt(5)*(1 - 3*nu)/63
            hHat_4_2_4=sqrt(5)*(nu*(4025 - 285*nu) - 1311)/20790
            hHat_4_2_5=sqrt(5)*(nu*(-30*pi + 84*I) + 10*pi - 21*I)/315
            hHat_4_2_6=sqrt(5)*(7*nu*(115*nu*(3363*nu + 34822) - 5460759) + 9342351)/113513400
            hHat_4_3_3=9*sqrt(70)*I*delta*(2*nu - 1)/280
            hHat_4_3_5=3*sqrt(70)*I*delta*(nu*(524*nu - 1267) + 468)/12320
            hHat_4_3_6=sqrt(70)*delta*(16301*nu + (4860 - 9720*nu)*log(3/2) + 2430*I*pi*(2*nu - 1) - 5184)/25200
            hHat_4_4_2=8*sqrt(35)*(3*nu - 1)/63
            hHat_4_4_4=sqrt(35)*(20*nu*(525*nu - 1273) + 7116)/10395
            hHat_4_4_5=sqrt(35)*(pi*(480*nu - 160) + I*(-1193*nu + (960*nu - 320)*log(2) + 336))/315
            hHat_4_4_6=sqrt(35)*(7*nu*(5*nu*(678291*nu - 3231338) + 9793071) - 9618039)/14189175
            hHat_5_1_3=sqrt(385)*I*delta*(1 - 2*nu)/110880
            hHat_5_1_5=-sqrt(385)*I*delta*(nu*(4*nu - 352) + 179)/4324320
            hHat_5_1_6=sqrt(385)*delta*(-28*nu*(log(1024) + 313) - 70*I*pi*(2*nu - 1) + 140*log(2) + 181)/7761600
            hHat_5_2_4=sqrt(55)*(2*nu*(5*nu - 5) + 2)/1485
            hHat_5_2_6=sqrt(55)*(7*nu*(35*nu*(33*nu - 118) + 3079) - 3911)/675675
            hHat_5_3_3=9*sqrt(330)*I*delta*(2*nu - 1)/3520
            hHat_5_3_5=3*sqrt(330)*I*delta*(8*nu*(11*nu - 58) + 207)/45760
            hHat_5_3_6=sqrt(330)*delta*(1171828*nu + 153090*I*pi*(2*nu - 1) - (612360*nu - 306180)*log(3/2) - 395847)/19958400
            hHat_5_4_4=sqrt(165)*(-32*nu*(5*nu - 5) - 32)/1485
            hHat_5_4_6=sqrt(165)*(-112*nu*(5*nu*(339*nu - 1042) + 3619) + 71216)/675675
            hHat_5_5_3=-625*sqrt(66)*I*delta*(2*nu - 1)/6336
            hHat_5_5_5=-625*sqrt(66)*I*delta*(16*nu*(16*nu - 43) + 263)/247104
            hHat_5_5_6=sqrt(66)*delta*(-1481676*nu - 218750*I*pi*(2*nu - 1) + (875000*nu - 437500)*log(5/2) + 565625)/443520
            hHat_6_1_5=sqrt(26)*I*delta*(nu - 1)*(3*nu - 1)/216216
            hHat_6_2_4=sqrt(65)*(2*nu*(5*nu - 5) + 2)/19305
            hHat_6_2_6=sqrt(65)*(7*nu*(nu*(7*nu - 64) + 59) - 81)/135135
            hHat_6_3_5=-81*sqrt(65)*I*delta*(nu - 1)*(3*nu - 1)/40040
            hHat_6_4_4=-128*sqrt(78)*(nu*(5*nu - 5) + 1)/19305
            hHat_6_4_6=-64*sqrt(78)*(7*nu*(nu*(19*nu - 88) + 71) - 93)/135135
            hHat_6_5_5=3125*sqrt(429)*I*delta*(nu - 1)*(3*nu - 1)/216216
            hHat_6_6_4=sqrt(143)*(54*nu*(5*nu - 5) + 54)/715
            hHat_6_6_6=sqrt(143)*(189*nu*(nu*(39*nu - 128) + 91) - 3051)/5005
            hHat_7_1_5=sqrt(2)*I*delta*(nu - 1)*(3*nu - 1)/1729728
            hHat_7_2_6=sqrt(3)*(-7*nu*(nu - 1)**2 + 1)/9009
            hHat_7_3_5=-243*sqrt(6)*I*delta*(nu - 1)*(3*nu - 1)/320320
            hHat_7_4_6=128*sqrt(66)*(7*nu*(nu - 1)**2 - 1)/45045
            hHat_7_5_5=15625*sqrt(66)*I*delta*(nu - 1)*(3*nu - 1)/1729728
            hHat_7_6_6=-81*sqrt(429)*(7*nu*(nu - 1)**2 - 1)/5005
            hHat_7_7_5=-16807*sqrt(6006)*I*delta*(nu - 1)*(3*nu - 1)/1235520
            hHat_8_2_6=sqrt(85)*(-7*nu*(nu - 1)**2 + 1)/765765
            hHat_8_4_6=128*sqrt(374)*(7*nu*(nu - 1)**2 - 1)/765765
            hHat_8_6_6=-243*sqrt(51051)*(7*nu*(nu - 1)**2 - 1)/595595
            hHat_8_8_6=16384*sqrt(170170)*(7*nu*(nu - 1)**2 - 1)/5360355
            hHat_spin_Symm_2_2_3=(-6*S_ell - 2*Sigma_ell*delta)/(3*M**2)
            hHat_spin_Symm_2_2_4=(12*S1_ell*S2_ell + 10*S1_lambda*S2_lambda - 15*I*S1_lambda*S2_n - 15*I*S1_n*S2_lambda - 22*S1_n*S2_n)/(6*M**4*nu)
            hHat_spin_Symm_2_1_2=I*Sigma_ell/(2*M**2)
            hHat_spin_Symm_2_1_4=I*(-86*S_ell*delta + Sigma_ell*(139*nu - 79))/(42*M**2)
            hHat_spin_Symm_2_0_4=sqrt(6)*(-S1_lambda*S2_lambda + S1_n*S2_n)/(3*M**4*nu)
            hHat_spin_Symm_3_3_4=3*sqrt(210)*I*(7*S_ell*delta - 3*Sigma_ell*(3*nu - 1))/(112*M**2)
            hHat_spin_Symm_3_2_3=2*sqrt(35)*(S_ell + Sigma_ell*delta)/(21*M**2)
            hHat_spin_Symm_3_1_4=sqrt(14)*I*(S_ell*delta - 5*Sigma_ell*(3*nu - 1))/(336*M**2)
            hHat_spin_Symm_4_3_4=9*sqrt(70)*I*(-S_ell*delta + 3*Sigma_ell*nu - Sigma_ell)/(112*M**2)
            hHat_spin_Symm_4_1_4=sqrt(10)*I*(S_ell*delta - 3*Sigma_ell*nu + Sigma_ell)/(336*M**2)
            hHat_spin_Asymm_2_2_2=(-Sigma_lambda - I*Sigma_n)/(2*M**2)
            hHat_spin_Asymm_2_2_4=(19*S_lambda*delta + 182*I*S_n*delta - 43*Sigma_lambda*nu + 5*Sigma_lambda - 280*I*Sigma_n*nu + 98*I*Sigma_n)/(84*M**2)
            hHat_spin_Asymm_2_1_3=(4*I*S_lambda + 25*S_n + 4*I*Sigma_lambda*delta + 13*Sigma_n*delta)/(6*M**2)
            hHat_spin_Asymm_2_1_4=(-3*S1_ell*S2_n - 3*S1_n*S2_ell)/(2*M**4*nu)
            hHat_spin_Asymm_2_0_2=sqrt(6)*I*Sigma_n/(6*M**2)
            hHat_spin_Asymm_2_0_4=sqrt(6)*I*(255*S_n*delta - Sigma_n*(506*nu - 45))/(126*M**2)
            hHat_spin_Asymm_3_3_3=sqrt(210)*I*(S_lambda + I*S_n + delta*(Sigma_lambda + I*Sigma_n))/(21*M**2)
            hHat_spin_Asymm_3_2_4=sqrt(35)*(-Sigma_lambda*(83*nu - 17) + 4*I*Sigma_n*(55*nu - 13) + 25*delta*(S_lambda - 4*I*S_n))/(168*M**2)
            hHat_spin_Asymm_3_1_3=sqrt(14)*(I*S_lambda + S_n + delta*(I*Sigma_lambda + Sigma_n))/(21*M**2)
            hHat_spin_Asymm_3_0_4=sqrt(42)*(-17*S_lambda*delta + Sigma_lambda*(35*nu - 9))/(168*M**2)
            hHat_spin_Asymm_4_4_4=9*sqrt(35)*(-3*Sigma_lambda*nu + Sigma_lambda - I*Sigma_n*(3*nu - 1) + delta*(S_lambda + I*S_n))/(56*M**2)
            hHat_spin_Asymm_4_2_4=sqrt(5)*(-13*Sigma_lambda*(3*nu - 1) + 14*I*Sigma_n*(3*nu - 1) + delta*(13*S_lambda - 14*I*S_n))/(168*M**2)
            hHat_spin_Asymm_4_0_4=sqrt(2)*I*(S_n*delta - 3*Sigma_n*nu + Sigma_n)/(168*M**2)
            
            chiVec1 = [0., chi1_n, chi1_lambda, chi1_ell]
            chiVec2 = [0., chi2_n, chi2_lambda, chi2_ell]


            chi1_n = chiVec1[1]
            chi1_lambda = chiVec1[2]
            chi1_ell = chiVec1[3]
            chi2_n = chiVec2[1]
            chi2_lambda = chiVec2[2]
            chi2_ell = chiVec2[3]
            S_ell = M1**2*chi1_ell + M2**2*chi2_ell
            S_n = M1**2*chi1_n + M2**2*chi2_n
            S_lambda = M1**2*chi1_lambda + M2**2*chi2_lambda
            Sigma_ell = M*(-M1*chi1_ell + M2*chi2_ell)
            Sigma_n = M*(-M1*chi1_n + M2*chi2_n)
            Sigma_lambda = M*(-M1*chi1_lambda + M2*chi2_lambda)
            S1_ell = M1**2*chi1_ell
            S1_n = M1**2*chi1_n
            S1_lambda = M1**2*chi1_lambda
            S2_ell = M2**2*chi2_ell
            S2_n = M2**2*chi2_n
            S2_lambda = M2**2*chi2_lambda
            logv = log(v)
            rhOverM_coeff = 8*sqrt(5)*nu*v**2*sqrt(pi)/5
            hHat_spin_Symm_2_2_3 = (-6*S_ell - 2*Sigma_ell*delta)/(3*M**2)
            hHat_spin_Symm_2_2_4 = (12*S1_ell*S2_ell + 10*S1_lambda*S2_lambda - 15*I*S1_lambda*S2_n - 15*I*S1_n*S2_lambda - 22*S1_n*S2_n)/(6*M**4*nu)
            hHat_spin_Symm_2_1_2 = I*Sigma_ell/(2*M**2)
            hHat_spin_Symm_2_1_4 = I*(-86*S_ell*delta + Sigma_ell*(139*nu - 79))/(42*M**2)
            hHat_spin_Symm_2_0_4 = sqrt(6)*(-S1_lambda*S2_lambda + S1_n*S2_n)/(3*M**4*nu)
            hHat_spin_Symm_3_3_4 = 3*sqrt(210)*I*(7*S_ell*delta - 3*Sigma_ell*(3*nu - 1))/(112*M**2)
            hHat_spin_Symm_3_2_3 = 2*sqrt(35)*(S_ell + Sigma_ell*delta)/(21*M**2)
            hHat_spin_Symm_3_1_4 = sqrt(14)*I*(S_ell*delta - 5*Sigma_ell*(3*nu - 1))/(336*M**2)
            hHat_spin_Symm_4_3_4 = 9*sqrt(70)*I*(-S_ell*delta + 3*Sigma_ell*nu - Sigma_ell)/(112*M**2)
            hHat_spin_Symm_4_1_4 = sqrt(10)*I*(S_ell*delta - 3*Sigma_ell*nu + Sigma_ell)/(336*M**2)
            hHat_spin_Asymm_2_2_2 = (-Sigma_lambda - I*Sigma_n)/(2*M**2)
            hHat_spin_Asymm_2_2_4 = (19*S_lambda*delta + 182*I*S_n*delta - 43*Sigma_lambda*nu + 5*Sigma_lambda - 280*I*Sigma_n*nu + 98*I*Sigma_n)/(84*M**2)
            hHat_spin_Asymm_2_1_3 = (4*I*S_lambda + 25*S_n + 4*I*Sigma_lambda*delta + 13*Sigma_n*delta)/(6*M**2)
            hHat_spin_Asymm_2_1_4 = (-3*S1_ell*S2_n - 3*S1_n*S2_ell)/(2*M**4*nu)
            hHat_spin_Asymm_2_0_2 = sqrt(6)*I*Sigma_n/(6*M**2)
            hHat_spin_Asymm_2_0_4 = sqrt(6)*I*(255*S_n*delta - Sigma_n*(506*nu - 45))/(126*M**2)
            hHat_spin_Asymm_3_3_3 = sqrt(210)*I*(S_lambda + I*S_n + delta*(Sigma_lambda + I*Sigma_n))/(21*M**2)
            hHat_spin_Asymm_3_2_4 = sqrt(35)*(-Sigma_lambda*(83*nu - 17) + 4*I*Sigma_n*(55*nu - 13) + 25*delta*(S_lambda - 4*I*S_n))/(168*M**2)
            hHat_spin_Asymm_3_1_3 = sqrt(14)*(I*S_lambda + S_n + delta*(I*Sigma_lambda + Sigma_n))/(21*M**2)
            hHat_spin_Asymm_3_0_4 = sqrt(42)*(-17*S_lambda*delta + Sigma_lambda*(35*nu - 9))/(168*M**2)
            hHat_spin_Asymm_4_4_4 = 9*sqrt(35)*(-3*Sigma_lambda*nu + Sigma_lambda - I*Sigma_n*(3*nu - 1) + delta*(S_lambda + I*S_n))/(56*M**2)
            hHat_spin_Asymm_4_2_4 = sqrt(5)*(-13*Sigma_lambda*(3*nu - 1) + 14*I*Sigma_n*(3*nu - 1) + delta*(13*S_lambda - 14*I*S_n))/(168*M**2)
            hHat_spin_Asymm_4_0_4 = sqrt(2)*I*(S_n*delta - 3*Sigma_n*nu + Sigma_n)/(168*M**2)

            Modes=np.empty(77, dtype=complex)
            # (ell, m) = (2, +/- 0)
            Symm = rhOverM_coeff*(hHat_2_0_0 + hHat_spin_Symm_2_0_4*v**4)
            Asymm = rhOverM_coeff*v**2*(hHat_spin_Asymm_2_0_2 + hHat_spin_Asymm_2_0_4*v**2)
            Modes[2] = Symm + Asymm
            # (ell, m) = (2, +/- 1)
            Symm = rhOverM_coeff*v*(hHat_2_1_1 + v*(hHat_spin_Symm_2_1_2 + v*(hHat_2_1_3 + v*(hHat_2_1_4 + hHat_spin_Symm_2_1_4 + v*(hHat_2_1_5 + hHat_2_1_6*v)))))
            Asymm = rhOverM_coeff*v**3*(hHat_spin_Asymm_2_1_3 + hHat_spin_Asymm_2_1_4*v)
            Modes[3] = Symm + Asymm
            Modes[1] = conjugate(Symm - Asymm)
            # (ell, m) = (2, +/- 2)
            Symm = rhOverM_coeff*(hHat_2_2_0 + v**2*(hHat_2_2_2 + v*(hHat_2_2_3 + hHat_spin_Symm_2_2_3 + v*(hHat_2_2_4 + hHat_spin_Symm_2_2_4 + v*(hHat_2_2_5 + v*(hHat_2_2_6 + hHat_2_2_lnv_6*logv))))))
            Asymm = rhOverM_coeff*v**2*(hHat_spin_Asymm_2_2_2 + hHat_spin_Asymm_2_2_4*v**2)
            Modes[4] = Symm + Asymm
            Modes[0] = conjugate(Symm - Asymm)
            # (ell, m) = (3, +/- 0)
            Symm = hHat_3_0_5*rhOverM_coeff*v**5
            Asymm = hHat_spin_Asymm_3_0_4*rhOverM_coeff*v**4
            Modes[8] = Symm + Asymm
            # (ell, m) = (3, +/- 1)
            Symm = rhOverM_coeff*v*(hHat_3_1_1 + v**2*(hHat_3_1_3 + v*(hHat_3_1_4 + hHat_spin_Symm_3_1_4 + v*(hHat_3_1_5 + hHat_3_1_6*v))))
            Asymm = hHat_spin_Asymm_3_1_3*rhOverM_coeff*v**3
            Modes[9] = Symm + Asymm
            Modes[7] = -conjugate(Symm - Asymm)
            # (ell, m) = (3, +/- 2)
            Symm = rhOverM_coeff*v**2*(hHat_3_2_2 + v*(hHat_spin_Symm_3_2_3 + v*(hHat_3_2_4 + v*(hHat_3_2_5 + hHat_3_2_6*v))))
            Asymm = hHat_spin_Asymm_3_2_4*rhOverM_coeff*v**4
            Modes[10] = Symm + Asymm
            Modes[6] = -conjugate(Symm - Asymm)
            # (ell, m) = (3, +/- 3)
            Symm = rhOverM_coeff*v*(hHat_3_3_1 + v**2*(hHat_3_3_3 + v*(hHat_3_3_4 + hHat_spin_Symm_3_3_4 + v*(hHat_3_3_5 + hHat_3_3_6*v))))
            Asymm = hHat_spin_Asymm_3_3_3*rhOverM_coeff*v**3
            Modes[11] = Symm + Asymm
            Modes[5] = -conjugate(Symm - Asymm)
            # (ell, m) = (4, +/- 0)
            Symm = hHat_4_0_0*rhOverM_coeff
            Asymm = hHat_spin_Asymm_4_0_4*rhOverM_coeff*v**4
            Modes[16] = Symm + Asymm
            # (ell, m) = (4, +/- 1)
            Symm = rhOverM_coeff*v**3*(hHat_4_1_3 + v*(hHat_spin_Symm_4_1_4 + v*(hHat_4_1_5 + hHat_4_1_6*v)))
            Asymm = 0
            Modes[17] = Symm + Asymm
            Modes[15] = conjugate(Symm - Asymm)
            # (ell, m) = (4, +/- 2)
            Symm = rhOverM_coeff*v**2*(hHat_4_2_2 + v**2*(hHat_4_2_4 + v*(hHat_4_2_5 + hHat_4_2_6*v)))
            Asymm = hHat_spin_Asymm_4_2_4*rhOverM_coeff*v**4
            Modes[18] = Symm + Asymm
            Modes[14] = conjugate(Symm - Asymm)
            # (ell, m) = (4, +/- 3)
            Symm = rhOverM_coeff*v**3*(hHat_4_3_3 + v*(hHat_spin_Symm_4_3_4 + v*(hHat_4_3_5 + hHat_4_3_6*v)))
            Asymm = 0
            Modes[19] = Symm + Asymm
            Modes[13] = conjugate(Symm - Asymm)
            # (ell, m) = (4, +/- 4)
            Symm = rhOverM_coeff*v**2*(hHat_4_4_2 + v**2*(hHat_4_4_4 + v*(hHat_4_4_5 + hHat_4_4_6*v)))
            Asymm = hHat_spin_Asymm_4_4_4*rhOverM_coeff*v**4
            Modes[20] = Symm + Asymm
            Modes[12] = conjugate(Symm - Asymm)
            # (ell, m) = (5, +/- 0)
            Symm = 0
            Asymm = 0
            Modes[26] = Symm + Asymm
            # (ell, m) = (5, +/- 1)
            Symm = rhOverM_coeff*v**3*(hHat_5_1_3 + v**2*(hHat_5_1_5 + hHat_5_1_6*v))
            Asymm = 0
            Modes[27] = Symm + Asymm
            Modes[25] = -conjugate(Symm - Asymm)
            # (ell, m) = (5, +/- 2)
            Symm = rhOverM_coeff*v**4*(hHat_5_2_4 + hHat_5_2_6*v**2)
            Asymm = 0
            Modes[28] = Symm + Asymm
            Modes[24] = -conjugate(Symm - Asymm)
            # (ell, m) = (5, +/- 3)
            Symm = rhOverM_coeff*v**3*(hHat_5_3_3 + v**2*(hHat_5_3_5 + hHat_5_3_6*v))
            Asymm = 0
            Modes[29] = Symm + Asymm
            Modes[23] = -conjugate(Symm - Asymm)
            # (ell, m) = (5, +/- 4)
            Symm = rhOverM_coeff*v**4*(hHat_5_4_4 + hHat_5_4_6*v**2)
            Asymm = 0
            Modes[30] = Symm + Asymm
            Modes[22] = -conjugate(Symm - Asymm)
            # (ell, m) = (5, +/- 5)
            Symm = rhOverM_coeff*v**3*(hHat_5_5_3 + v**2*(hHat_5_5_5 + hHat_5_5_6*v))
            Asymm = 0
            Modes[31] = Symm + Asymm
            Modes[21] = -conjugate(Symm - Asymm)
            # (ell, m) = (6, +/- 0)
            Symm = 0
            Asymm = 0
            Modes[38] = Symm + Asymm
            # (ell, m) = (6, +/- 1)
            Symm = hHat_6_1_5*rhOverM_coeff*v**5
            Asymm = 0
            Modes[39] = Symm + Asymm
            Modes[37] = conjugate(Symm - Asymm)
            # (ell, m) = (6, +/- 2)
            Symm = rhOverM_coeff*v**4*(hHat_6_2_4 + hHat_6_2_6*v**2)
            Asymm = 0
            Modes[40] = Symm + Asymm
            Modes[36] = conjugate(Symm - Asymm)
            # (ell, m) = (6, +/- 3)
            Symm = hHat_6_3_5*rhOverM_coeff*v**5
            Asymm = 0
            Modes[41] = Symm + Asymm
            Modes[35] = conjugate(Symm - Asymm)
            # (ell, m) = (6, +/- 4)
            Symm = rhOverM_coeff*v**4*(hHat_6_4_4 + hHat_6_4_6*v**2)
            Asymm = 0
            Modes[42] = Symm + Asymm
            Modes[34] = conjugate(Symm - Asymm)
            # (ell, m) = (6, +/- 5)
            Symm = hHat_6_5_5*rhOverM_coeff*v**5
            Asymm = 0
            Modes[43] = Symm + Asymm
            Modes[33] = conjugate(Symm - Asymm)
            # (ell, m) = (6, +/- 6)
            Symm = rhOverM_coeff*v**4*(hHat_6_6_4 + hHat_6_6_6*v**2)
            Asymm = 0
            Modes[44] = Symm + Asymm
            Modes[32] = conjugate(Symm - Asymm)
            # (ell, m) = (7, +/- 0)
            Symm = 0
            Asymm = 0
            Modes[52] = Symm + Asymm
            # (ell, m) = (7, +/- 1)
            Symm = hHat_7_1_5*rhOverM_coeff*v**5
            Asymm = 0
            Modes[53] = Symm + Asymm
            Modes[51] = -conjugate(Symm - Asymm)
            # (ell, m) = (7, +/- 2)
            Symm = hHat_7_2_6*rhOverM_coeff*v**6
            Asymm = 0
            Modes[54] = Symm + Asymm
            Modes[50] = -conjugate(Symm - Asymm)
            # (ell, m) = (7, +/- 3)
            Symm = hHat_7_3_5*rhOverM_coeff*v**5
            Asymm = 0
            Modes[55] = Symm + Asymm
            Modes[49] = -conjugate(Symm - Asymm)
            # (ell, m) = (7, +/- 4)
            Symm = hHat_7_4_6*rhOverM_coeff*v**6
            Asymm = 0
            Modes[56] = Symm + Asymm
            Modes[48] = -conjugate(Symm - Asymm)
            # (ell, m) = (7, +/- 5)
            Symm = hHat_7_5_5*rhOverM_coeff*v**5
            Asymm = 0
            Modes[57] = Symm + Asymm
            Modes[47] = -conjugate(Symm - Asymm)
            # (ell, m) = (7, +/- 6)
            Symm = hHat_7_6_6*rhOverM_coeff*v**6
            Asymm = 0
            Modes[58] = Symm + Asymm
            Modes[46] = -conjugate(Symm - Asymm)
            # (ell, m) = (7, +/- 7)
            Symm = hHat_7_7_5*rhOverM_coeff*v**5
            Asymm = 0
            Modes[59] = Symm + Asymm
            Modes[45] = -conjugate(Symm - Asymm)
            # (ell, m) = (8, +/- 0)
            Symm = 0
            Asymm = 0
            Modes[68] = Symm + Asymm
            # (ell, m) = (8, +/- 1)
            Symm = 0
            Asymm = 0
            Modes[69] = Symm + Asymm
            Modes[67] = conjugate(Symm - Asymm)
            # (ell, m) = (8, +/- 2)
            Symm = hHat_8_2_6*rhOverM_coeff*v**6
            Asymm = 0
            Modes[70] = Symm + Asymm
            Modes[66] = conjugate(Symm - Asymm)
            # (ell, m) = (8, +/- 3)
            Symm = 0
            Asymm = 0
            Modes[71] = Symm + Asymm
            Modes[65] = conjugate(Symm - Asymm)
            # (ell, m) = (8, +/- 4)
            Symm = hHat_8_4_6*rhOverM_coeff*v**6
            Asymm = 0
            Modes[72] = Symm + Asymm
            Modes[64] = conjugate(Symm - Asymm)
            # (ell, m) = (8, +/- 5)
            Symm = 0
            Asymm = 0
            Modes[73] = Symm + Asymm
            Modes[63] = conjugate(Symm - Asymm)
            # (ell, m) = (8, +/- 6)
            Symm = hHat_8_6_6*rhOverM_coeff*v**6
            Asymm = 0
            Modes[74] = Symm + Asymm
            Modes[62] = conjugate(Symm - Asymm)
            # (ell, m) = (8, +/- 7)
            Symm = 0
            Asymm = 0
            Modes[75] = Symm + Asymm
            Modes[61] = conjugate(Symm - Asymm)
            # (ell, m) = (8, +/- 8)
            Symm = hHat_8_8_6*rhOverM_coeff*v**6
            Asymm = 0
            Modes[76] = Symm + Asymm
            Modes[60] = conjugate(Symm - Asymm)

            return Modes
            
        for i in range(len(y[0])):
            v = y[0,i]
            rfrak_chi1_x = y[1,i]
            rfrak_chi1_y = y[2,i]
            rfrak_chi2_x = y[3,i]
            rfrak_chi2_y = y[4,i]
            rfrak_frame_x = y[5,i]
            rfrak_frame_y = y[6,i]
            rfrak_frame_z = y[7,i]
            Phi = y[8,i]
            xHat=xHat_i
            yHat=yHat_i
            zHat=zHat_i
            M1=M1_i
            M2=M2_i 
            S_chi1=S_chi1_i
            S_chi2=S_chi2_i
            M=M1 + M2
            delta=(M1 - M2)/M
            nu=M1*M2/M**2
            R=exp(rfrak_frame_x*xHat + rfrak_frame_y*yHat + rfrak_frame_z*zHat)
            nHat=R*xHat*conjugate(R)
            lambdaHat=R*yHat*conjugate(R)
            ellHat=R*zHat*conjugate(R)
            R_S1=exp(rfrak_chi1_x*xHat + rfrak_chi1_y*yHat)
            R_S2=exp(rfrak_chi2_x*xHat + rfrak_chi2_y*yHat)
            chiVec1=S_chi1*R_S1*zHat*conjugate(R_S1)*conjugate(S_chi1)
            chiVec2=S_chi2*R_S2*zHat*conjugate(R_S2)*conjugate(S_chi2)
            v_i=v
            rfrak_frame_x_i=rfrak_frame_x
            rfrak_frame_y_i=rfrak_frame_y
            rfrak_frame_z_i=rfrak_frame_z
            chiVec1_i=quaternion.as_float_array(chiVec1)
            chiVec2_i=quaternion.as_float_array(chiVec2)
            ModeData[i,:]=WaveformModes(M1_i, M2_i, v_i, chiVec1_i, chiVec2_i)
        return ModeData

class WaveformModes_3p5PN :
    def WaveformModes_3p5PN(xHat_i, yHat_i, zHat_i, M1_i, M2_i, v_i, S_chi1_i, S_chi2_i, rfrak_chi1_x_i, rfrak_chi1_y_i, rfrak_chi2_x_i, rfrak_chi2_y_i, rfrak_frame_x_i, rfrak_frame_y_i, rfrak_frame_z_i,y) :            
        ModeData=np.empty((len(y[0]),77), dtype=complex)
        def WaveformModes(M1_i, M2_i, v_i, chiVec1_i, chiVec2_i):
            I=1j
            global chi1_n
            global chi1_lambda
            global chi1_ell
            global chi2_n
            global chi2_lambda
            global chi2_ell
            global S_ell
            global S_n
            global S_lambda
            global Sigma_ell
            global Sigma_n
            global Sigma_lambda
            global S1_ell
            global S1_n
            global S1_lambda
            global S2_ell
            global S2_n
            global S2_lambda
            global logv
            global rhOverM_coeff
            global hHat_spin_Symm_2_2_3
            global hHat_spin_Symm_2_2_4
            global hHat_spin_Symm_2_1_2
            global hHat_spin_Symm_2_1_4
            global hHat_spin_Symm_2_0_4
            global hHat_spin_Symm_3_3_4
            global hHat_spin_Symm_3_2_3
            global hHat_spin_Symm_3_1_4
            global hHat_spin_Symm_4_3_4
            global hHat_spin_Symm_4_1_4
            global hHat_spin_Asymm_2_2_2
            global hHat_spin_Asymm_2_2_4
            global hHat_spin_Asymm_2_1_3
            global hHat_spin_Asymm_2_1_4
            global hHat_spin_Asymm_2_0_2
            global hHat_spin_Asymm_2_0_4
            global hHat_spin_Asymm_3_3_3
            global hHat_spin_Asymm_3_2_4
            global hHat_spin_Asymm_3_1_3
            global hHat_spin_Asymm_3_0_4
            global hHat_spin_Asymm_4_4_4
            global hHat_spin_Asymm_4_2_4
            global hHat_spin_Asymm_4_0_4
            M1=M1_i
            M2=M2_i
            v=v_i
            M=M1 + M2
            delta=(M1 - M2)/M
            nu=M1*M2/M**2
            chiVec1=chiVec1_i
            chiVec2=chiVec2_i
            chi1_n=chiVec1[1]
            chi1_lambda=chiVec1[2]
            chi1_ell=chiVec1[3]
            chi2_n=chiVec2[1]
            chi2_lambda=chiVec2[2]
            chi2_ell=chiVec2[3]
            S_ell=M1**2*chi1_ell + M2**2*chi2_ell
            S_n=M1**2*chi1_n + M2**2*chi2_n
            S_lambda=M1**2*chi1_lambda + M2**2*chi2_lambda
            Sigma_ell=M*(-M1*chi1_ell + M2*chi2_ell)
            Sigma_n=M*(-M1*chi1_n + M2*chi2_n)
            Sigma_lambda=M*(-M1*chi1_lambda + M2*chi2_lambda)
            S1_ell=M1**2*chi1_ell
            S1_n=M1**2*chi1_n
            S1_lambda=M1**2*chi1_lambda
            S2_ell=M2**2*chi2_ell
            S2_n=M2**2*chi2_n
            S2_lambda=M2**2*chi2_lambda
            logv=log(v)
            rhOverM_coeff=8*sqrt(5)*nu*v**2*sqrt(pi)/5
            hHat_2_0_0=-5*sqrt(6)/84
            hHat_2_1_1=I*delta/3
            hHat_2_1_3=I*delta*(20*nu - 17)/84
            hHat_2_1_4=delta*(1 + log(16) + 2*I*pi)/6
            hHat_2_1_5=I*delta*(nu*(237*nu - 2036) - 172)/1512
            hHat_2_1_6=delta*(2*nu*(log(4096) + 353 + 6*I*pi) - 17*log(16) - 17 - 34*I*pi)/168
            hHat_2_2_0=1
            hHat_2_2_2=55*nu/42 - 107/42
            hHat_2_2_3=2*pi
            hHat_2_2_4=nu*(2047*nu - 7483)/1512 - 2173/1512
            hHat_2_2_5=-24*I*nu + pi*(34*nu - 107)/21
            hHat_2_2_6=nu*(nu*(114635*nu - 729396) - 834555)/99792 + 41*nu*pi**2/96 - 428*log(16)/105 - 856*EulerGamma/105 + 27027409/646800 + 2*pi*(35*pi + 214*I)/105
            hHat_2_2_lnv_6=-856/105
            hHat_2_2_7=-I*nu*(24396*nu - 501655)/5670 + pi*(30*nu*(560*nu - 2459) - 32595)/11340
            hHat_3_0_5=-2*sqrt(42)*I*nu/35
            hHat_3_1_1=sqrt(14)*I*delta/168
            hHat_3_1_3=-sqrt(14)*I*delta*(nu + 4)/252
            hHat_3_1_4=sqrt(14)*delta*(log(1024) + 7 + 5*I*pi)/840
            hHat_3_1_5=-sqrt(14)*I*delta*(nu*(247*nu + 272) - 607)/33264
            hHat_3_1_6=sqrt(14)*delta*(2*nu - 5*I*pi*(7*nu + 16) - 2*(35*nu + 80)*log(2) - 112)/5040
            hHat_3_1_lnv_7=-13*sqrt(14)*I*delta/1764
            hHat_3_1_7=sqrt(14)*I*delta*(-17525*nu**3/15444 + 327059*nu**2/30888 + nu*(-1738843/19305 + 41*pi**2/8)/8 - 2*(log(2) + 212/105)*log(2) - 26*EulerGamma/21 + pi**2/6 + 10753397/1513512 - 2*I*pi*(41/105 + log(2)))/168
            hHat_3_2_2=sqrt(35)*(1 - 3*nu)/21
            hHat_3_2_4=sqrt(35)*(nu*(725 - 365*nu) - 193)/1890
            hHat_3_2_5=sqrt(35)*(-30*nu*pi + 66*I*nu + 10*pi - 15*I)/105
            hHat_3_2_6=sqrt(35)*(nu*(nu*(100026 - 16023*nu) - 17387) - 1451)/83160
            hHat_3_3_1=-3*sqrt(210)*I*delta/56
            hHat_3_3_3=-3*sqrt(210)*I*delta*(nu - 2)/28
            hHat_3_3_4=9*sqrt(210)*delta*(-7 + 10*log(3/2) - 5*I*pi)/280
            hHat_3_3_5=-sqrt(210)*I*delta*(nu*(887*nu - 3676) + 369)/6160
            hHat_3_3_6=sqrt(210)*delta*(-96206*nu - 3645*I*pi*(3*nu - 8) + 2*(10935*nu - 29160)*log(3/2) + 40824)/45360
            hHat_3_3_lnv_7=117*sqrt(210)*I*delta*log(4)/196
            hHat_3_3_7=-3*sqrt(210)*I*delta*(8237*nu**3/2860 - 318841*nu**2/17160 + nu*(-7055/429 + 41*pi**2/8)/8 - 78*EulerGamma/7 - 18*log(3/2)**2 + 492*log(3/2)/35 + 3*pi**2/2 + 19388147/280280 + 6*I*pi*(-41/35 + 3*log(3/2)))/56
            hHat_4_0_0=-sqrt(2)/1008
            hHat_4_1_3=sqrt(10)*I*delta*(1 - 2*nu)/840
            hHat_4_1_5=-sqrt(10)*I*delta*(nu*(332*nu - 1011) + 404)/110880
            hHat_4_1_6=sqrt(10)*delta*(-1661*nu + (60 - 120*nu)*log(2) - 30*I*pi*(2*nu - 1) + 64)/25200
            hHat_4_2_2=sqrt(5)*(1 - 3*nu)/63
            hHat_4_2_4=sqrt(5)*(nu*(4025 - 285*nu) - 1311)/20790
            hHat_4_2_5=sqrt(5)*(nu*(-30*pi + 84*I) + 10*pi - 21*I)/315
            hHat_4_2_6=sqrt(5)*(7*nu*(115*nu*(3363*nu + 34822) - 5460759) + 9342351)/113513400
            hHat_4_3_3=9*sqrt(70)*I*delta*(2*nu - 1)/280
            hHat_4_3_5=3*sqrt(70)*I*delta*(nu*(524*nu - 1267) + 468)/12320
            hHat_4_3_6=sqrt(70)*delta*(16301*nu + (4860 - 9720*nu)*log(3/2) + 2430*I*pi*(2*nu - 1) - 5184)/25200
            hHat_4_4_2=8*sqrt(35)*(3*nu - 1)/63
            hHat_4_4_4=sqrt(35)*(20*nu*(525*nu - 1273) + 7116)/10395
            hHat_4_4_5=sqrt(35)*(pi*(480*nu - 160) + I*(-1193*nu + (960*nu - 320)*log(2) + 336))/315
            hHat_4_4_6=sqrt(35)*(7*nu*(5*nu*(678291*nu - 3231338) + 9793071) - 9618039)/14189175
            hHat_5_1_3=sqrt(385)*I*delta*(1 - 2*nu)/110880
            hHat_5_1_5=-sqrt(385)*I*delta*(nu*(4*nu - 352) + 179)/4324320
            hHat_5_1_6=sqrt(385)*delta*(-28*nu*(log(1024) + 313) - 70*I*pi*(2*nu - 1) + 140*log(2) + 181)/7761600
            hHat_5_2_4=sqrt(55)*(2*nu*(5*nu - 5) + 2)/1485
            hHat_5_2_6=sqrt(55)*(7*nu*(35*nu*(33*nu - 118) + 3079) - 3911)/675675
            hHat_5_3_3=9*sqrt(330)*I*delta*(2*nu - 1)/3520
            hHat_5_3_5=3*sqrt(330)*I*delta*(8*nu*(11*nu - 58) + 207)/45760
            hHat_5_3_6=sqrt(330)*delta*(1171828*nu + 153090*I*pi*(2*nu - 1) - (612360*nu - 306180)*log(3/2) - 395847)/19958400
            hHat_5_4_4=sqrt(165)*(-32*nu*(5*nu - 5) - 32)/1485
            hHat_5_4_6=sqrt(165)*(-112*nu*(5*nu*(339*nu - 1042) + 3619) + 71216)/675675
            hHat_5_5_3=-625*sqrt(66)*I*delta*(2*nu - 1)/6336
            hHat_5_5_5=-625*sqrt(66)*I*delta*(16*nu*(16*nu - 43) + 263)/247104
            hHat_5_5_6=sqrt(66)*delta*(-1481676*nu - 218750*I*pi*(2*nu - 1) + (875000*nu - 437500)*log(5/2) + 565625)/443520
            hHat_6_1_5=sqrt(26)*I*delta*(nu - 1)*(3*nu - 1)/216216
            hHat_6_2_4=sqrt(65)*(2*nu*(5*nu - 5) + 2)/19305
            hHat_6_2_6=sqrt(65)*(7*nu*(nu*(7*nu - 64) + 59) - 81)/135135
            hHat_6_3_5=-81*sqrt(65)*I*delta*(nu - 1)*(3*nu - 1)/40040
            hHat_6_4_4=-128*sqrt(78)*(nu*(5*nu - 5) + 1)/19305
            hHat_6_4_6=-64*sqrt(78)*(7*nu*(nu*(19*nu - 88) + 71) - 93)/135135
            hHat_6_5_5=3125*sqrt(429)*I*delta*(nu - 1)*(3*nu - 1)/216216
            hHat_6_6_4=sqrt(143)*(54*nu*(5*nu - 5) + 54)/715
            hHat_6_6_6=sqrt(143)*(189*nu*(nu*(39*nu - 128) + 91) - 3051)/5005
            hHat_7_1_5=sqrt(2)*I*delta*(nu - 1)*(3*nu - 1)/1729728
            hHat_7_2_6=sqrt(3)*(-7*nu*(nu - 1)**2 + 1)/9009
            hHat_7_3_5=-243*sqrt(6)*I*delta*(nu - 1)*(3*nu - 1)/320320
            hHat_7_4_6=128*sqrt(66)*(7*nu*(nu - 1)**2 - 1)/45045
            hHat_7_5_5=15625*sqrt(66)*I*delta*(nu - 1)*(3*nu - 1)/1729728
            hHat_7_6_6=-81*sqrt(429)*(7*nu*(nu - 1)**2 - 1)/5005
            hHat_7_7_5=-16807*sqrt(6006)*I*delta*(nu - 1)*(3*nu - 1)/1235520
            hHat_8_2_6=sqrt(85)*(-7*nu*(nu - 1)**2 + 1)/765765
            hHat_8_4_6=128*sqrt(374)*(7*nu*(nu - 1)**2 - 1)/765765
            hHat_8_6_6=-243*sqrt(51051)*(7*nu*(nu - 1)**2 - 1)/595595
            hHat_8_8_6=16384*sqrt(170170)*(7*nu*(nu - 1)**2 - 1)/5360355
            hHat_spin_Symm_2_2_3=(-6*S_ell - 2*Sigma_ell*delta)/(3*M**2)
            hHat_spin_Symm_2_2_4=(12*S1_ell*S2_ell + 10*S1_lambda*S2_lambda - 15*I*S1_lambda*S2_n - 15*I*S1_n*S2_lambda - 22*S1_n*S2_n)/(6*M**4*nu)
            hHat_spin_Symm_2_1_2=I*Sigma_ell/(2*M**2)
            hHat_spin_Symm_2_1_4=I*(-86*S_ell*delta + Sigma_ell*(139*nu - 79))/(42*M**2)
            hHat_spin_Symm_2_0_4=sqrt(6)*(-S1_lambda*S2_lambda + S1_n*S2_n)/(3*M**4*nu)
            hHat_spin_Symm_3_3_4=3*sqrt(210)*I*(7*S_ell*delta - 3*Sigma_ell*(3*nu - 1))/(112*M**2)
            hHat_spin_Symm_3_2_3=2*sqrt(35)*(S_ell + Sigma_ell*delta)/(21*M**2)
            hHat_spin_Symm_3_1_4=sqrt(14)*I*(S_ell*delta - 5*Sigma_ell*(3*nu - 1))/(336*M**2)
            hHat_spin_Symm_4_3_4=9*sqrt(70)*I*(-S_ell*delta + 3*Sigma_ell*nu - Sigma_ell)/(112*M**2)
            hHat_spin_Symm_4_1_4=sqrt(10)*I*(S_ell*delta - 3*Sigma_ell*nu + Sigma_ell)/(336*M**2)
            hHat_spin_Asymm_2_2_2=(-Sigma_lambda - I*Sigma_n)/(2*M**2)
            hHat_spin_Asymm_2_2_4=(19*S_lambda*delta + 182*I*S_n*delta - 43*Sigma_lambda*nu + 5*Sigma_lambda - 280*I*Sigma_n*nu + 98*I*Sigma_n)/(84*M**2)
            hHat_spin_Asymm_2_1_3=(4*I*S_lambda + 25*S_n + 4*I*Sigma_lambda*delta + 13*Sigma_n*delta)/(6*M**2)
            hHat_spin_Asymm_2_1_4=(-3*S1_ell*S2_n - 3*S1_n*S2_ell)/(2*M**4*nu)
            hHat_spin_Asymm_2_0_2=sqrt(6)*I*Sigma_n/(6*M**2)
            hHat_spin_Asymm_2_0_4=sqrt(6)*I*(255*S_n*delta - Sigma_n*(506*nu - 45))/(126*M**2)
            hHat_spin_Asymm_3_3_3=sqrt(210)*I*(S_lambda + I*S_n + delta*(Sigma_lambda + I*Sigma_n))/(21*M**2)
            hHat_spin_Asymm_3_2_4=sqrt(35)*(-Sigma_lambda*(83*nu - 17) + 4*I*Sigma_n*(55*nu - 13) + 25*delta*(S_lambda - 4*I*S_n))/(168*M**2)
            hHat_spin_Asymm_3_1_3=sqrt(14)*(I*S_lambda + S_n + delta*(I*Sigma_lambda + Sigma_n))/(21*M**2)
            hHat_spin_Asymm_3_0_4=sqrt(42)*(-17*S_lambda*delta + Sigma_lambda*(35*nu - 9))/(168*M**2)
            hHat_spin_Asymm_4_4_4=9*sqrt(35)*(-3*Sigma_lambda*nu + Sigma_lambda - I*Sigma_n*(3*nu - 1) + delta*(S_lambda + I*S_n))/(56*M**2)
            hHat_spin_Asymm_4_2_4=sqrt(5)*(-13*Sigma_lambda*(3*nu - 1) + 14*I*Sigma_n*(3*nu - 1) + delta*(13*S_lambda - 14*I*S_n))/(168*M**2)
            hHat_spin_Asymm_4_0_4=sqrt(2)*I*(S_n*delta - 3*Sigma_n*nu + Sigma_n)/(168*M**2)
            
            chiVec1 = [0., chi1_n, chi1_lambda, chi1_ell]
            chiVec2 = [0., chi2_n, chi2_lambda, chi2_ell]


            chi1_n = chiVec1[1]
            chi1_lambda = chiVec1[2]
            chi1_ell = chiVec1[3]
            chi2_n = chiVec2[1]
            chi2_lambda = chiVec2[2]
            chi2_ell = chiVec2[3]
            S_ell = M1**2*chi1_ell + M2**2*chi2_ell
            S_n = M1**2*chi1_n + M2**2*chi2_n
            S_lambda = M1**2*chi1_lambda + M2**2*chi2_lambda
            Sigma_ell = M*(-M1*chi1_ell + M2*chi2_ell)
            Sigma_n = M*(-M1*chi1_n + M2*chi2_n)
            Sigma_lambda = M*(-M1*chi1_lambda + M2*chi2_lambda)
            S1_ell = M1**2*chi1_ell
            S1_n = M1**2*chi1_n
            S1_lambda = M1**2*chi1_lambda
            S2_ell = M2**2*chi2_ell
            S2_n = M2**2*chi2_n
            S2_lambda = M2**2*chi2_lambda
            logv = log(v)
            rhOverM_coeff = 8*sqrt(5)*nu*v**2*sqrt(pi)/5
            hHat_spin_Symm_2_2_3 = (-6*S_ell - 2*Sigma_ell*delta)/(3*M**2)
            hHat_spin_Symm_2_2_4 = (12*S1_ell*S2_ell + 10*S1_lambda*S2_lambda - 15*I*S1_lambda*S2_n - 15*I*S1_n*S2_lambda - 22*S1_n*S2_n)/(6*M**4*nu)
            hHat_spin_Symm_2_1_2 = I*Sigma_ell/(2*M**2)
            hHat_spin_Symm_2_1_4 = I*(-86*S_ell*delta + Sigma_ell*(139*nu - 79))/(42*M**2)
            hHat_spin_Symm_2_0_4 = sqrt(6)*(-S1_lambda*S2_lambda + S1_n*S2_n)/(3*M**4*nu)
            hHat_spin_Symm_3_3_4 = 3*sqrt(210)*I*(7*S_ell*delta - 3*Sigma_ell*(3*nu - 1))/(112*M**2)
            hHat_spin_Symm_3_2_3 = 2*sqrt(35)*(S_ell + Sigma_ell*delta)/(21*M**2)
            hHat_spin_Symm_3_1_4 = sqrt(14)*I*(S_ell*delta - 5*Sigma_ell*(3*nu - 1))/(336*M**2)
            hHat_spin_Symm_4_3_4 = 9*sqrt(70)*I*(-S_ell*delta + 3*Sigma_ell*nu - Sigma_ell)/(112*M**2)
            hHat_spin_Symm_4_1_4 = sqrt(10)*I*(S_ell*delta - 3*Sigma_ell*nu + Sigma_ell)/(336*M**2)
            hHat_spin_Asymm_2_2_2 = (-Sigma_lambda - I*Sigma_n)/(2*M**2)
            hHat_spin_Asymm_2_2_4 = (19*S_lambda*delta + 182*I*S_n*delta - 43*Sigma_lambda*nu + 5*Sigma_lambda - 280*I*Sigma_n*nu + 98*I*Sigma_n)/(84*M**2)
            hHat_spin_Asymm_2_1_3 = (4*I*S_lambda + 25*S_n + 4*I*Sigma_lambda*delta + 13*Sigma_n*delta)/(6*M**2)
            hHat_spin_Asymm_2_1_4 = (-3*S1_ell*S2_n - 3*S1_n*S2_ell)/(2*M**4*nu)
            hHat_spin_Asymm_2_0_2 = sqrt(6)*I*Sigma_n/(6*M**2)
            hHat_spin_Asymm_2_0_4 = sqrt(6)*I*(255*S_n*delta - Sigma_n*(506*nu - 45))/(126*M**2)
            hHat_spin_Asymm_3_3_3 = sqrt(210)*I*(S_lambda + I*S_n + delta*(Sigma_lambda + I*Sigma_n))/(21*M**2)
            hHat_spin_Asymm_3_2_4 = sqrt(35)*(-Sigma_lambda*(83*nu - 17) + 4*I*Sigma_n*(55*nu - 13) + 25*delta*(S_lambda - 4*I*S_n))/(168*M**2)
            hHat_spin_Asymm_3_1_3 = sqrt(14)*(I*S_lambda + S_n + delta*(I*Sigma_lambda + Sigma_n))/(21*M**2)
            hHat_spin_Asymm_3_0_4 = sqrt(42)*(-17*S_lambda*delta + Sigma_lambda*(35*nu - 9))/(168*M**2)
            hHat_spin_Asymm_4_4_4 = 9*sqrt(35)*(-3*Sigma_lambda*nu + Sigma_lambda - I*Sigma_n*(3*nu - 1) + delta*(S_lambda + I*S_n))/(56*M**2)
            hHat_spin_Asymm_4_2_4 = sqrt(5)*(-13*Sigma_lambda*(3*nu - 1) + 14*I*Sigma_n*(3*nu - 1) + delta*(13*S_lambda - 14*I*S_n))/(168*M**2)
            hHat_spin_Asymm_4_0_4 = sqrt(2)*I*(S_n*delta - 3*Sigma_n*nu + Sigma_n)/(168*M**2)

            Modes=np.empty(77, dtype=complex)
            # (ell, m) = (2, +/- 0)
            Symm = rhOverM_coeff*(hHat_2_0_0 + hHat_spin_Symm_2_0_4*v**4)
            Asymm = rhOverM_coeff*v**2*(hHat_spin_Asymm_2_0_2 + hHat_spin_Asymm_2_0_4*v**2)
            Modes[2] = Symm + Asymm
            # (ell, m) = (2, +/- 1)
            Symm = rhOverM_coeff*v*(hHat_2_1_1 + v*(hHat_spin_Symm_2_1_2 + v*(hHat_2_1_3 + v*(hHat_2_1_4 + hHat_spin_Symm_2_1_4 + v*(hHat_2_1_5 + hHat_2_1_6*v)))))
            Asymm = rhOverM_coeff*v**3*(hHat_spin_Asymm_2_1_3 + hHat_spin_Asymm_2_1_4*v)
            Modes[3] = Symm + Asymm
            Modes[1] = conjugate(Symm - Asymm)
            # (ell, m) = (2, +/- 2)
            Symm = rhOverM_coeff*(hHat_2_2_0 + v**2*(hHat_2_2_2 + v*(hHat_2_2_3 + hHat_spin_Symm_2_2_3 + v*(hHat_2_2_4 + hHat_spin_Symm_2_2_4 + v*(hHat_2_2_5 + v*(hHat_2_2_6 + hHat_2_2_7*v + hHat_2_2_lnv_6*logv))))))
            Asymm = rhOverM_coeff*v**2*(hHat_spin_Asymm_2_2_2 + hHat_spin_Asymm_2_2_4*v**2)
            Modes[4] = Symm + Asymm
            Modes[0] = conjugate(Symm - Asymm)
            # (ell, m) = (3, +/- 0)
            Symm = hHat_3_0_5*rhOverM_coeff*v**5
            Asymm = hHat_spin_Asymm_3_0_4*rhOverM_coeff*v**4
            Modes[8] = Symm + Asymm
            # (ell, m) = (3, +/- 1)
            Symm = rhOverM_coeff*v*(hHat_3_1_1 + v**2*(hHat_3_1_3 + v*(hHat_3_1_4 + hHat_spin_Symm_3_1_4 + v*(hHat_3_1_5 + v*(hHat_3_1_6 + v*(hHat_3_1_7 + hHat_3_1_lnv_7*logv))))))
            Asymm = hHat_spin_Asymm_3_1_3*rhOverM_coeff*v**3
            Modes[9] = Symm + Asymm
            Modes[7] = -conjugate(Symm - Asymm)
            # (ell, m) = (3, +/- 2)
            Symm = rhOverM_coeff*v**2*(hHat_3_2_2 + v*(hHat_spin_Symm_3_2_3 + v*(hHat_3_2_4 + v*(hHat_3_2_5 + hHat_3_2_6*v))))
            Asymm = hHat_spin_Asymm_3_2_4*rhOverM_coeff*v**4
            Modes[10] = Symm + Asymm
            Modes[6] = -conjugate(Symm - Asymm)
            # (ell, m) = (3, +/- 3)
            Symm = rhOverM_coeff*v*(hHat_3_3_1 + v**2*(hHat_3_3_3 + v*(hHat_3_3_4 + hHat_spin_Symm_3_3_4 + v*(hHat_3_3_5 + v*(hHat_3_3_6 + v*(hHat_3_3_7 + hHat_3_3_lnv_7*logv))))))
            Asymm = hHat_spin_Asymm_3_3_3*rhOverM_coeff*v**3
            Modes[11] = Symm + Asymm
            Modes[5] = -conjugate(Symm - Asymm)
            # (ell, m) = (4, +/- 0)
            Symm = hHat_4_0_0*rhOverM_coeff
            Asymm = hHat_spin_Asymm_4_0_4*rhOverM_coeff*v**4
            Modes[16] = Symm + Asymm
            # (ell, m) = (4, +/- 1)
            Symm = rhOverM_coeff*v**3*(hHat_4_1_3 + v*(hHat_spin_Symm_4_1_4 + v*(hHat_4_1_5 + hHat_4_1_6*v)))
            Asymm = 0
            Modes[17] = Symm + Asymm
            Modes[15] = conjugate(Symm - Asymm)
            # (ell, m) = (4, +/- 2)
            Symm = rhOverM_coeff*v**2*(hHat_4_2_2 + v**2*(hHat_4_2_4 + v*(hHat_4_2_5 + hHat_4_2_6*v)))
            Asymm = hHat_spin_Asymm_4_2_4*rhOverM_coeff*v**4
            Modes[18] = Symm + Asymm
            Modes[14] = conjugate(Symm - Asymm)
            # (ell, m) = (4, +/- 3)
            Symm = rhOverM_coeff*v**3*(hHat_4_3_3 + v*(hHat_spin_Symm_4_3_4 + v*(hHat_4_3_5 + hHat_4_3_6*v)))
            Asymm = 0
            Modes[19] = Symm + Asymm
            Modes[13] = conjugate(Symm - Asymm)
            # (ell, m) = (4, +/- 4)
            Symm = rhOverM_coeff*v**2*(hHat_4_4_2 + v**2*(hHat_4_4_4 + v*(hHat_4_4_5 + hHat_4_4_6*v)))
            Asymm = hHat_spin_Asymm_4_4_4*rhOverM_coeff*v**4
            Modes[20] = Symm + Asymm
            Modes[12] = conjugate(Symm - Asymm)
            # (ell, m) = (5, +/- 0)
            Symm = 0
            Asymm = 0
            Modes[26] = Symm + Asymm
            # (ell, m) = (5, +/- 1)
            Symm = rhOverM_coeff*v**3*(hHat_5_1_3 + v**2*(hHat_5_1_5 + hHat_5_1_6*v))
            Asymm = 0
            Modes[27] = Symm + Asymm
            Modes[25] = -conjugate(Symm - Asymm)
            # (ell, m) = (5, +/- 2)
            Symm = rhOverM_coeff*v**4*(hHat_5_2_4 + hHat_5_2_6*v**2)
            Asymm = 0
            Modes[28] = Symm + Asymm
            Modes[24] = -conjugate(Symm - Asymm)
            # (ell, m) = (5, +/- 3)
            Symm = rhOverM_coeff*v**3*(hHat_5_3_3 + v**2*(hHat_5_3_5 + hHat_5_3_6*v))
            Asymm = 0
            Modes[29] = Symm + Asymm
            Modes[23] = -conjugate(Symm - Asymm)
            # (ell, m) = (5, +/- 4)
            Symm = rhOverM_coeff*v**4*(hHat_5_4_4 + hHat_5_4_6*v**2)
            Asymm = 0
            Modes[30] = Symm + Asymm
            Modes[22] = -conjugate(Symm - Asymm)
            # (ell, m) = (5, +/- 5)
            Symm = rhOverM_coeff*v**3*(hHat_5_5_3 + v**2*(hHat_5_5_5 + hHat_5_5_6*v))
            Asymm = 0
            Modes[31] = Symm + Asymm
            Modes[21] = -conjugate(Symm - Asymm)
            # (ell, m) = (6, +/- 0)
            Symm = 0
            Asymm = 0
            Modes[38] = Symm + Asymm
            # (ell, m) = (6, +/- 1)
            Symm = hHat_6_1_5*rhOverM_coeff*v**5
            Asymm = 0
            Modes[39] = Symm + Asymm
            Modes[37] = conjugate(Symm - Asymm)
            # (ell, m) = (6, +/- 2)
            Symm = rhOverM_coeff*v**4*(hHat_6_2_4 + hHat_6_2_6*v**2)
            Asymm = 0
            Modes[40] = Symm + Asymm
            Modes[36] = conjugate(Symm - Asymm)
            # (ell, m) = (6, +/- 3)
            Symm = hHat_6_3_5*rhOverM_coeff*v**5
            Asymm = 0
            Modes[41] = Symm + Asymm
            Modes[35] = conjugate(Symm - Asymm)
            # (ell, m) = (6, +/- 4)
            Symm = rhOverM_coeff*v**4*(hHat_6_4_4 + hHat_6_4_6*v**2)
            Asymm = 0
            Modes[42] = Symm + Asymm
            Modes[34] = conjugate(Symm - Asymm)
            # (ell, m) = (6, +/- 5)
            Symm = hHat_6_5_5*rhOverM_coeff*v**5
            Asymm = 0
            Modes[43] = Symm + Asymm
            Modes[33] = conjugate(Symm - Asymm)
            # (ell, m) = (6, +/- 6)
            Symm = rhOverM_coeff*v**4*(hHat_6_6_4 + hHat_6_6_6*v**2)
            Asymm = 0
            Modes[44] = Symm + Asymm
            Modes[32] = conjugate(Symm - Asymm)
            # (ell, m) = (7, +/- 0)
            Symm = 0
            Asymm = 0
            Modes[52] = Symm + Asymm
            # (ell, m) = (7, +/- 1)
            Symm = hHat_7_1_5*rhOverM_coeff*v**5
            Asymm = 0
            Modes[53] = Symm + Asymm
            Modes[51] = -conjugate(Symm - Asymm)
            # (ell, m) = (7, +/- 2)
            Symm = hHat_7_2_6*rhOverM_coeff*v**6
            Asymm = 0
            Modes[54] = Symm + Asymm
            Modes[50] = -conjugate(Symm - Asymm)
            # (ell, m) = (7, +/- 3)
            Symm = hHat_7_3_5*rhOverM_coeff*v**5
            Asymm = 0
            Modes[55] = Symm + Asymm
            Modes[49] = -conjugate(Symm - Asymm)
            # (ell, m) = (7, +/- 4)
            Symm = hHat_7_4_6*rhOverM_coeff*v**6
            Asymm = 0
            Modes[56] = Symm + Asymm
            Modes[48] = -conjugate(Symm - Asymm)
            # (ell, m) = (7, +/- 5)
            Symm = hHat_7_5_5*rhOverM_coeff*v**5
            Asymm = 0
            Modes[57] = Symm + Asymm
            Modes[47] = -conjugate(Symm - Asymm)
            # (ell, m) = (7, +/- 6)
            Symm = hHat_7_6_6*rhOverM_coeff*v**6
            Asymm = 0
            Modes[58] = Symm + Asymm
            Modes[46] = -conjugate(Symm - Asymm)
            # (ell, m) = (7, +/- 7)
            Symm = hHat_7_7_5*rhOverM_coeff*v**5
            Asymm = 0
            Modes[59] = Symm + Asymm
            Modes[45] = -conjugate(Symm - Asymm)
            # (ell, m) = (8, +/- 0)
            Symm = 0
            Asymm = 0
            Modes[68] = Symm + Asymm
            # (ell, m) = (8, +/- 1)
            Symm = 0
            Asymm = 0
            Modes[69] = Symm + Asymm
            Modes[67] = conjugate(Symm - Asymm)
            # (ell, m) = (8, +/- 2)
            Symm = hHat_8_2_6*rhOverM_coeff*v**6
            Asymm = 0
            Modes[70] = Symm + Asymm
            Modes[66] = conjugate(Symm - Asymm)
            # (ell, m) = (8, +/- 3)
            Symm = 0
            Asymm = 0
            Modes[71] = Symm + Asymm
            Modes[65] = conjugate(Symm - Asymm)
            # (ell, m) = (8, +/- 4)
            Symm = hHat_8_4_6*rhOverM_coeff*v**6
            Asymm = 0
            Modes[72] = Symm + Asymm
            Modes[64] = conjugate(Symm - Asymm)
            # (ell, m) = (8, +/- 5)
            Symm = 0
            Asymm = 0
            Modes[73] = Symm + Asymm
            Modes[63] = conjugate(Symm - Asymm)
            # (ell, m) = (8, +/- 6)
            Symm = hHat_8_6_6*rhOverM_coeff*v**6
            Asymm = 0
            Modes[74] = Symm + Asymm
            Modes[62] = conjugate(Symm - Asymm)
            # (ell, m) = (8, +/- 7)
            Symm = 0
            Asymm = 0
            Modes[75] = Symm + Asymm
            Modes[61] = conjugate(Symm - Asymm)
            # (ell, m) = (8, +/- 8)
            Symm = hHat_8_8_6*rhOverM_coeff*v**6
            Asymm = 0
            Modes[76] = Symm + Asymm
            Modes[60] = conjugate(Symm - Asymm)

            return Modes
            
        for i in range(len(y[0])):
            v = y[0,i]
            rfrak_chi1_x = y[1,i]
            rfrak_chi1_y = y[2,i]
            rfrak_chi2_x = y[3,i]
            rfrak_chi2_y = y[4,i]
            rfrak_frame_x = y[5,i]
            rfrak_frame_y = y[6,i]
            rfrak_frame_z = y[7,i]
            Phi = y[8,i]
            xHat=xHat_i
            yHat=yHat_i
            zHat=zHat_i
            M1=M1_i
            M2=M2_i 
            S_chi1=S_chi1_i
            S_chi2=S_chi2_i
            M=M1 + M2
            delta=(M1 - M2)/M
            nu=M1*M2/M**2
            R=exp(rfrak_frame_x*xHat + rfrak_frame_y*yHat + rfrak_frame_z*zHat)
            nHat=R*xHat*conjugate(R)
            lambdaHat=R*yHat*conjugate(R)
            ellHat=R*zHat*conjugate(R)
            R_S1=exp(rfrak_chi1_x*xHat + rfrak_chi1_y*yHat)
            R_S2=exp(rfrak_chi2_x*xHat + rfrak_chi2_y*yHat)
            chiVec1=S_chi1*R_S1*zHat*conjugate(R_S1)*conjugate(S_chi1)
            chiVec2=S_chi2*R_S2*zHat*conjugate(R_S2)*conjugate(S_chi2)
            v_i=v
            rfrak_frame_x_i=rfrak_frame_x
            rfrak_frame_y_i=rfrak_frame_y
            rfrak_frame_z_i=rfrak_frame_z
            chiVec1_i=quaternion.as_float_array(chiVec1)
            chiVec2_i=quaternion.as_float_array(chiVec2)
            ModeData[i,:]=WaveformModes(M1_i, M2_i, v_i, chiVec1_i, chiVec2_i)
        return ModeData
