# File produced automatically by OrbitalEvolutionCodeGen_Q.ipynb
from scipy.integrate import solve_ivp
import numpy as np
from numpy import conjugate, dot, exp, log, sqrt, pi
from numpy import euler_gamma as EulerGamma
import Quaternions

class TaylorTn_0PN_Q : 
    def TaylorTn_0PN_Q(xHat_i, yHat_i, zHat_i, M1_i, M2_i, v_i, S_chi1_i, S_chi2_i, rfrak_chi1_x_i, rfrak_chi1_y_i, rfrak_chi2_x_i,
           rfrak_chi2_y_i, rfrak_frame_x_i, rfrak_frame_y_i, rfrak_frame_z_i) :
        xHat=xHat_i
        yHat=yHat_i
        zHat=zHat_i
        M1=M1_i
        M2=M2_i
        v=v_i
        S_chi1=S_chi1_i
        S_chi2=S_chi2_i
        rfrak_chi1_x=rfrak_chi1_x_i
        rfrak_chi1_y=rfrak_chi1_y_i
        rfrak_chi2_x=rfrak_chi2_x_i
        rfrak_chi2_y=rfrak_chi2_y_i
        rfrak_frame_x=rfrak_frame_x_i
        rfrak_frame_y=rfrak_frame_y_i
        rfrak_frame_z=rfrak_frame_z_i
        M=M1 + M2
        delta=(M1 - M2)/M
        nu=M1*M2/M**2
        R=exp(rfrak_frame_x*xHat + rfrak_frame_y*yHat + rfrak_frame_z*zHat)
        nHat=R*xHat*conjugate(R)
        ellHat=R*zHat*conjugate(R)
        R_S1=exp(rfrak_chi1_x*xHat + rfrak_chi1_y*yHat)
        R_S2=exp(rfrak_chi2_x*xHat + rfrak_chi2_y*yHat)
        chiVec1=S_chi1*R_S1*zHat*conjugate(R_S1)*conjugate(S_chi1)
        chiVec2=S_chi2*R_S2*zHat*conjugate(R_S2)*conjugate(S_chi2)
        chi1chi1=dot(chiVec1,chiVec1)
        chi1chi2=dot(chiVec1,chiVec2)
        chi2chi2=dot(chiVec2,chiVec2)
        chi1_n=dot(chiVec1,nHat)
        chi1_ell=dot(chiVec1,ellHat)
        chi2_n=dot(chiVec2,nHat)
        chi2_ell=dot(chiVec2,ellHat)
        S_ell=M1**2*chi1_ell + M2**2*chi2_ell
        S_n=M1**2*chi1_n + M2**2*chi2_n
        Sigma_ell=M*(-M1*chi1_ell + M2*chi2_ell)
        Sigma_n=M*(-M1*chi1_n + M2*chi2_n)
        chi_s_ell=chi1_ell/2 + chi2_ell/2
        chi_a_ell=chi1_ell/2 - chi2_ell/2
        Fcal_coeff=32*nu**2*v**10/5
        Fcal_0=1
        E_0=1
        Phi=0.0
        EvolveSpin1=dot(S_chi1,S_chi1)>1e-12
        EvolveSpin2=dot(S_chi2,S_chi2)>1e-12

        def Recalculate(t, y):
            v = y[0];
            rfrak_chi1_x = y[1];
            rfrak_chi1_y = y[2];
            rfrak_chi2_x = y[3];
            rfrak_chi2_y = y[4];
            rfrak_frame_x = y[5];
            rfrak_frame_y = y[6];
            rfrak_frame_z = y[7];
            Phi = y[8];
            R = exp(rfrak_frame_x*xHat + rfrak_frame_y*yHat + rfrak_frame_z*zHat)
            nHat = R*xHat*conjugate(R)
            ellHat = R*zHat*conjugate(R)
            R_S1 = exp(rfrak_chi1_x*xHat + rfrak_chi1_y*yHat)
            R_S2 = exp(rfrak_chi2_x*xHat + rfrak_chi2_y*yHat)
            chiVec1 = S_chi1*R_S1*zHat*conjugate(R_S1)*conjugate(S_chi1)
            chiVec2 = S_chi2*R_S2*zHat*conjugate(R_S2)*conjugate(S_chi2)
            chi1_n = dot(chiVec1,nHat)
            chi1_ell = dot(chiVec1,ellHat)
            chi2_n = dot(chiVec2,nHat)
            chi2_ell = dot(chiVec2,ellHat)
            S_ell = M1**2*chi1_ell + M2**2*chi2_ell
            S_n = M1**2*chi1_n + M2**2*chi2_n
            Sigma_ell = M*(-M1*chi1_ell + M2*chi2_ell)
            Sigma_n = M*(-M1*chi1_n + M2*chi2_n)
            chi_s_ell = chi1_ell/2 + chi2_ell/2
            chi_a_ell = chi1_ell/2 - chi2_ell/2
            Fcal_coeff = 32*nu**2*v**10/5

        def OmegaVec_chiVec_1():
            Omega1_coeff = v**5/M
            return Omega1_coeff*ellHat*(-0.75*delta + 0.5*nu + 0.75)

        def OmegaVec_chiVec_2():
            Omega2_coeff = v**5/M
            return Omega2_coeff*ellHat*(0.75*delta + 0.5*nu + 0.75)

        def OmegaVec():
            a_ell_0 = 7.0*S_n + 3.0*Sigma_n*delta
            gamma_PN_0 = 1.00000000000000
            return ellHat*v**3/M + a_ell_0*gamma_PN_0*nHat*v**6/M**3

        def OrbitalAngularMomentum():
            L_0 = ellHat
            L_coeff = M**2*nu/v
            return L_0*L_coeff


        #Evolve PN
        def TaylorT1(t, y):
            Recalculate(t, y)
            if v>=1.0: 
                return GSL_EDOM # Beyond domain of PN validity
            Flux = Fcal_0*Fcal_coeff
            dEdv = -E_0*M*nu*v
            Absorption = 0
            dvdt_T1 = (-Absorption - Flux)/dEdv
            if dvdt_T1<1.0e-12:
                return GSL_EDIVERGE # v is decreasing
            return CommonRHS(dvdt_T1, y)

        def CommonRHS(dvdt, y):
            dydt=np.zeros(9)
            rfrak_frame=np.zeros(3)
            rfrak_frame[0] = y[5]
            rfrak_frame[1] = y[6]
            rfrak_frame[2] = y[7]
            rfrakdot_frame = Quaternions.FrameFromAngularVelocity_Integrand(rfrak_frame, OmegaVec().vec())
            dydt[0] = dvdt
            if(EvolveSpin1):
                Quaternions.FrameFromAngularVelocity_2D_Integrand(y[1], y[2],\
                    (S_chi1.inverse()*OmegaVec_chiVec_1()*S_chi1).vec(),dydt[1], dydt[2])
            else:
                dydt[1] = 0.0
                dydt[2] = 0.0
            if(EvolveSpin2):
                Quaternions.FrameFromAngularVelocity_2D_Integrand(y[3], y[4],\
                    (S_chi2.inverse()*OmegaVec_chiVec_2()*S_chi2).vec(),dydt[3], dydt[4])
            else:
                dydt[3] = 0.0
                dydt[4] = 0.0
            dydt[5] = rfrakdot_frame[0]
            dydt[6] = rfrakdot_frame[1]
            dydt[7] = rfrakdot_frame[2]
            dydt[8] = v*v*v/M

            return dydt
        
        
        y=solve_ivp(TaylorT1, [0,20000], [0.01**(1/3),rfrak_chi1_x,\
            rfrak_chi1_y,rfrak_chi2_x,rfrak_chi2_y,rfrak_frame_x,\
            rfrak_frame_y,rfrak_frame_z,Phi], method='DOP853',\
            t_eval=np.arange(0,20000,100), dense_output=True)
            

#        def TaylorT1_phi(t,phi):
#            dphidt_T1=v_of_t.y[0][v_of_t.t>=t][0]**3.0/M
#            return dphidt_T1
#        phi_of_t=solve_ivp(TaylorT1_phi, [0, v_of_t.t[-1]], [0.0], method='RK45', t_eval=v_of_t.t)
        
        #Calculate waveform modes

        def WaveformModes():
            I=1j
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

            Modes=np.zeros(77,len(y.t))
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

class TaylorTn_0p50PN_Q : 
    def TaylorTn_0p50PN_Q(xHat_i, yHat_i, zHat_i, M1_i, M2_i, v_i, S_chi1_i, S_chi2_i, rfrak_chi1_x_i, rfrak_chi1_y_i, rfrak_chi2_x_i,
           rfrak_chi2_y_i, rfrak_frame_x_i, rfrak_frame_y_i, rfrak_frame_z_i) :
        xHat=xHat_i
        yHat=yHat_i
        zHat=zHat_i
        M1=M1_i
        M2=M2_i
        v=v_i
        S_chi1=S_chi1_i
        S_chi2=S_chi2_i
        rfrak_chi1_x=rfrak_chi1_x_i
        rfrak_chi1_y=rfrak_chi1_y_i
        rfrak_chi2_x=rfrak_chi2_x_i
        rfrak_chi2_y=rfrak_chi2_y_i
        rfrak_frame_x=rfrak_frame_x_i
        rfrak_frame_y=rfrak_frame_y_i
        rfrak_frame_z=rfrak_frame_z_i
        M=M1 + M2
        delta=(M1 - M2)/M
        nu=M1*M2/M**2
        R=exp(rfrak_frame_x*xHat + rfrak_frame_y*yHat + rfrak_frame_z*zHat)
        nHat=R*xHat*conjugate(R)
        ellHat=R*zHat*conjugate(R)
        R_S1=exp(rfrak_chi1_x*xHat + rfrak_chi1_y*yHat)
        R_S2=exp(rfrak_chi2_x*xHat + rfrak_chi2_y*yHat)
        chiVec1=S_chi1*R_S1*zHat*conjugate(R_S1)*conjugate(S_chi1)
        chiVec2=S_chi2*R_S2*zHat*conjugate(R_S2)*conjugate(S_chi2)
        chi1chi1=dot(chiVec1,chiVec1)
        chi1chi2=dot(chiVec1,chiVec2)
        chi2chi2=dot(chiVec2,chiVec2)
        chi1_n=dot(chiVec1,nHat)
        chi1_ell=dot(chiVec1,ellHat)
        chi2_n=dot(chiVec2,nHat)
        chi2_ell=dot(chiVec2,ellHat)
        S_ell=M1**2*chi1_ell + M2**2*chi2_ell
        S_n=M1**2*chi1_n + M2**2*chi2_n
        Sigma_ell=M*(-M1*chi1_ell + M2*chi2_ell)
        Sigma_n=M*(-M1*chi1_n + M2*chi2_n)
        chi_s_ell=chi1_ell/2 + chi2_ell/2
        chi_a_ell=chi1_ell/2 - chi2_ell/2
        Fcal_coeff=32*nu**2*v**10/5
        Fcal_0=1
        E_0=1
        Phi=0.0
        EvolveSpin1=dot(S_chi1,S_chi1)>1e-12
        EvolveSpin2=dot(S_chi2,S_chi2)>1e-12

        def Recalculate(t, y):
            v = y[0];
            rfrak_chi1_x = y[1];
            rfrak_chi1_y = y[2];
            rfrak_chi2_x = y[3];
            rfrak_chi2_y = y[4];
            rfrak_frame_x = y[5];
            rfrak_frame_y = y[6];
            rfrak_frame_z = y[7];
            Phi = y[8];
            R = exp(rfrak_frame_x*xHat + rfrak_frame_y*yHat + rfrak_frame_z*zHat)
            nHat = R*xHat*conjugate(R)
            ellHat = R*zHat*conjugate(R)
            R_S1 = exp(rfrak_chi1_x*xHat + rfrak_chi1_y*yHat)
            R_S2 = exp(rfrak_chi2_x*xHat + rfrak_chi2_y*yHat)
            chiVec1 = S_chi1*R_S1*zHat*conjugate(R_S1)*conjugate(S_chi1)
            chiVec2 = S_chi2*R_S2*zHat*conjugate(R_S2)*conjugate(S_chi2)
            chi1_n = dot(chiVec1,nHat)
            chi1_ell = dot(chiVec1,ellHat)
            chi2_n = dot(chiVec2,nHat)
            chi2_ell = dot(chiVec2,ellHat)
            S_ell = M1**2*chi1_ell + M2**2*chi2_ell
            S_n = M1**2*chi1_n + M2**2*chi2_n
            Sigma_ell = M*(-M1*chi1_ell + M2*chi2_ell)
            Sigma_n = M*(-M1*chi1_n + M2*chi2_n)
            chi_s_ell = chi1_ell/2 + chi2_ell/2
            chi_a_ell = chi1_ell/2 - chi2_ell/2
            Fcal_coeff = 32*nu**2*v**10/5

        def OmegaVec_chiVec_1():
            Omega1_coeff = v**5/M
            return Omega1_coeff*(ellHat*(-0.75*delta + 0.5*nu + 0.75) + nHat*v*(3.0*chi1_n*nu + 3.0*M2**2*chi2_n/M**2) - M2**2*chiVec2*v/M**2)

        def OmegaVec_chiVec_2():
            Omega2_coeff = v**5/M
            return Omega2_coeff*(ellHat*(0.75*delta + 0.5*nu + 0.75) + nHat*v*(3.0*chi2_n*nu + 3.0*M1**2*chi1_n/M**2) - M1**2*chiVec1*v/M**2)

        def OmegaVec():
            a_ell_0 = 7.0*S_n + 3.0*Sigma_n*delta
            gamma_PN_0 = 1.00000000000000
            return ellHat*v**3/M + a_ell_0*gamma_PN_0*nHat*v**6/M**3

        def OrbitalAngularMomentum():
            L_0 = ellHat
            L_coeff = M**2*nu/v
            return L_0*L_coeff


        #Evolve PN
        def TaylorT1(t, y):
            Recalculate(t, y)
            if v>=1.0: 
                return GSL_EDOM # Beyond domain of PN validity
            Flux = Fcal_0*Fcal_coeff
            dEdv = -E_0*M*nu*v
            Absorption = 0
            dvdt_T1 = (-Absorption - Flux)/dEdv
            if dvdt_T1<1.0e-12:
                return GSL_EDIVERGE # v is decreasing
            return CommonRHS(dvdt_T1, y)

        def CommonRHS(dvdt, y):
            dydt=np.zeros(9)
            rfrak_frame=np.zeros(3)
            rfrak_frame[0] = y[5]
            rfrak_frame[1] = y[6]
            rfrak_frame[2] = y[7]
            rfrakdot_frame = Quaternions.FrameFromAngularVelocity_Integrand(rfrak_frame, OmegaVec().vec())
            dydt[0] = dvdt
            if(EvolveSpin1):
                Quaternions.FrameFromAngularVelocity_2D_Integrand(y[1], y[2],\
                    (S_chi1.inverse()*OmegaVec_chiVec_1()*S_chi1).vec(),dydt[1], dydt[2])
            else:
                dydt[1] = 0.0
                dydt[2] = 0.0
            if(EvolveSpin2):
                Quaternions.FrameFromAngularVelocity_2D_Integrand(y[3], y[4],\
                    (S_chi2.inverse()*OmegaVec_chiVec_2()*S_chi2).vec(),dydt[3], dydt[4])
            else:
                dydt[3] = 0.0
                dydt[4] = 0.0
            dydt[5] = rfrakdot_frame[0]
            dydt[6] = rfrakdot_frame[1]
            dydt[7] = rfrakdot_frame[2]
            dydt[8] = v*v*v/M

            return dydt
        
        
        y=solve_ivp(TaylorT1, [0,20000], [0.01**(1/3),rfrak_chi1_x,\
            rfrak_chi1_y,rfrak_chi2_x,rfrak_chi2_y,rfrak_frame_x,\
            rfrak_frame_y,rfrak_frame_z,Phi], method='DOP853',\
            t_eval=np.arange(0,20000,100), dense_output=True)
            

#        def TaylorT1_phi(t,phi):
#            dphidt_T1=v_of_t.y[0][v_of_t.t>=t][0]**3.0/M
#            return dphidt_T1
#        phi_of_t=solve_ivp(TaylorT1_phi, [0, v_of_t.t[-1]], [0.0], method='RK45', t_eval=v_of_t.t)
        
        #Calculate waveform modes

        def WaveformModes():
            I=1j
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

            Modes=np.zeros(77,len(y.t))
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

class TaylorTn_1p0PN_Q : 
    def TaylorTn_1p0PN_Q(xHat_i, yHat_i, zHat_i, M1_i, M2_i, v_i, S_chi1_i, S_chi2_i, rfrak_chi1_x_i, rfrak_chi1_y_i, rfrak_chi2_x_i,
           rfrak_chi2_y_i, rfrak_frame_x_i, rfrak_frame_y_i, rfrak_frame_z_i) :
        xHat=xHat_i
        yHat=yHat_i
        zHat=zHat_i
        M1=M1_i
        M2=M2_i
        v=v_i
        S_chi1=S_chi1_i
        S_chi2=S_chi2_i
        rfrak_chi1_x=rfrak_chi1_x_i
        rfrak_chi1_y=rfrak_chi1_y_i
        rfrak_chi2_x=rfrak_chi2_x_i
        rfrak_chi2_y=rfrak_chi2_y_i
        rfrak_frame_x=rfrak_frame_x_i
        rfrak_frame_y=rfrak_frame_y_i
        rfrak_frame_z=rfrak_frame_z_i
        M=M1 + M2
        delta=(M1 - M2)/M
        nu=M1*M2/M**2
        R=exp(rfrak_frame_x*xHat + rfrak_frame_y*yHat + rfrak_frame_z*zHat)
        nHat=R*xHat*conjugate(R)
        ellHat=R*zHat*conjugate(R)
        R_S1=exp(rfrak_chi1_x*xHat + rfrak_chi1_y*yHat)
        R_S2=exp(rfrak_chi2_x*xHat + rfrak_chi2_y*yHat)
        chiVec1=S_chi1*R_S1*zHat*conjugate(R_S1)*conjugate(S_chi1)
        chiVec2=S_chi2*R_S2*zHat*conjugate(R_S2)*conjugate(S_chi2)
        chi1chi1=dot(chiVec1,chiVec1)
        chi1chi2=dot(chiVec1,chiVec2)
        chi2chi2=dot(chiVec2,chiVec2)
        chi1_n=dot(chiVec1,nHat)
        chi1_ell=dot(chiVec1,ellHat)
        chi2_n=dot(chiVec2,nHat)
        chi2_ell=dot(chiVec2,ellHat)
        S_ell=M1**2*chi1_ell + M2**2*chi2_ell
        S_n=M1**2*chi1_n + M2**2*chi2_n
        Sigma_ell=M*(-M1*chi1_ell + M2*chi2_ell)
        Sigma_n=M*(-M1*chi1_n + M2*chi2_n)
        chi_s_ell=chi1_ell/2 + chi2_ell/2
        chi_a_ell=chi1_ell/2 - chi2_ell/2
        Fcal_coeff=32*nu**2*v**10/5
        Fcal_0=1
        Fcal_2=-35*nu/12 - 1247/336
        E_0=1
        E_2=-nu/12 - 3/4
        Phi=0.0
        EvolveSpin1=dot(S_chi1,S_chi1)>1e-12
        EvolveSpin2=dot(S_chi2,S_chi2)>1e-12

        def Recalculate(t, y):
            v = y[0];
            rfrak_chi1_x = y[1];
            rfrak_chi1_y = y[2];
            rfrak_chi2_x = y[3];
            rfrak_chi2_y = y[4];
            rfrak_frame_x = y[5];
            rfrak_frame_y = y[6];
            rfrak_frame_z = y[7];
            Phi = y[8];
            R = exp(rfrak_frame_x*xHat + rfrak_frame_y*yHat + rfrak_frame_z*zHat)
            nHat = R*xHat*conjugate(R)
            ellHat = R*zHat*conjugate(R)
            R_S1 = exp(rfrak_chi1_x*xHat + rfrak_chi1_y*yHat)
            R_S2 = exp(rfrak_chi2_x*xHat + rfrak_chi2_y*yHat)
            chiVec1 = S_chi1*R_S1*zHat*conjugate(R_S1)*conjugate(S_chi1)
            chiVec2 = S_chi2*R_S2*zHat*conjugate(R_S2)*conjugate(S_chi2)
            chi1_n = dot(chiVec1,nHat)
            chi1_ell = dot(chiVec1,ellHat)
            chi2_n = dot(chiVec2,nHat)
            chi2_ell = dot(chiVec2,ellHat)
            S_ell = M1**2*chi1_ell + M2**2*chi2_ell
            S_n = M1**2*chi1_n + M2**2*chi2_n
            Sigma_ell = M*(-M1*chi1_ell + M2*chi2_ell)
            Sigma_n = M*(-M1*chi1_n + M2*chi2_n)
            chi_s_ell = chi1_ell/2 + chi2_ell/2
            chi_a_ell = chi1_ell/2 - chi2_ell/2
            Fcal_coeff = 32*nu**2*v**10/5

        def OmegaVec_chiVec_1():
            Omega1_coeff = v**5/M
            return Omega1_coeff*(ellHat*(-0.75*delta + 0.5*nu + v**2*(delta*(0.625*nu - 0.5625) + nu*(1.25 - 0.0416666666666667*nu) + 0.5625) + 0.75) + nHat*v*(3.0*chi1_n*nu + 3.0*M2**2*chi2_n/M**2) - M2**2*chiVec2*v/M**2)

        def OmegaVec_chiVec_2():
            Omega2_coeff = v**5/M
            return Omega2_coeff*(ellHat*(0.75*delta + 0.5*nu + v**2*(delta*(0.5625 - 0.625*nu) + nu*(1.25 - 0.0416666666666667*nu) + 0.5625) + 0.75) + nHat*v*(3.0*chi2_n*nu + 3.0*M1**2*chi1_n/M**2) - M1**2*chiVec1*v/M**2)

        def OmegaVec():
            gamma_PN_2 = 1.0 - 0.333333333333333*nu
            a_ell_2 = S_n*(-9.66666666666667*nu - 10.0) + Sigma_n*delta*(-4.5*nu - 6.0)
            a_ell_0 = 7.0*S_n + 3.0*Sigma_n*delta
            gamma_PN_0 = 1.00000000000000
            return ellHat*v**3/M + nHat*v**6*(a_ell_0 + a_ell_2*v**2)*(gamma_PN_0 + gamma_PN_2*v**2)/M**3

        def OrbitalAngularMomentum():
            L_2 = ellHat*(0.166666666666667*nu + 1.5)
            L_0 = ellHat
            L_coeff = M**2*nu/v
            return L_coeff*(L_0 + L_2*v**2)


        #Evolve PN
        def TaylorT1(t, y):
            Recalculate(t, y)
            if v>=1.0: 
                return GSL_EDOM # Beyond domain of PN validity
            Flux = Fcal_coeff*(Fcal_0 + Fcal_2*v**2)
            dEdv = -M*nu*v*(E_0 + 2.0*E_2*v**2)
            Absorption = 0
            dvdt_T1 = (-Absorption - Flux)/dEdv
            if dvdt_T1<1.0e-12:
                return GSL_EDIVERGE # v is decreasing
            return CommonRHS(dvdt_T1, y)

        def CommonRHS(dvdt, y):
            dydt=np.zeros(9)
            rfrak_frame=np.zeros(3)
            rfrak_frame[0] = y[5]
            rfrak_frame[1] = y[6]
            rfrak_frame[2] = y[7]
            rfrakdot_frame = Quaternions.FrameFromAngularVelocity_Integrand(rfrak_frame, OmegaVec().vec())
            dydt[0] = dvdt
            if(EvolveSpin1):
                Quaternions.FrameFromAngularVelocity_2D_Integrand(y[1], y[2],\
                    (S_chi1.inverse()*OmegaVec_chiVec_1()*S_chi1).vec(),dydt[1], dydt[2])
            else:
                dydt[1] = 0.0
                dydt[2] = 0.0
            if(EvolveSpin2):
                Quaternions.FrameFromAngularVelocity_2D_Integrand(y[3], y[4],\
                    (S_chi2.inverse()*OmegaVec_chiVec_2()*S_chi2).vec(),dydt[3], dydt[4])
            else:
                dydt[3] = 0.0
                dydt[4] = 0.0
            dydt[5] = rfrakdot_frame[0]
            dydt[6] = rfrakdot_frame[1]
            dydt[7] = rfrakdot_frame[2]
            dydt[8] = v*v*v/M

            return dydt
        
        
        y=solve_ivp(TaylorT1, [0,20000], [0.01**(1/3),rfrak_chi1_x,\
            rfrak_chi1_y,rfrak_chi2_x,rfrak_chi2_y,rfrak_frame_x,\
            rfrak_frame_y,rfrak_frame_z,Phi], method='DOP853',\
            t_eval=np.arange(0,20000,100), dense_output=True)
            

#        def TaylorT1_phi(t,phi):
#            dphidt_T1=v_of_t.y[0][v_of_t.t>=t][0]**3.0/M
#            return dphidt_T1
#        phi_of_t=solve_ivp(TaylorT1_phi, [0, v_of_t.t[-1]], [0.0], method='RK45', t_eval=v_of_t.t)
        
        #Calculate waveform modes

        def WaveformModes():
            I=1j
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
            chiVec1 = Quaternions.Quaternion(0., chi1_n_k, chi1_lambda_k, chi1_ell_k);
            chiVec2 = Quaternions.Quaternion(0., chi2_n_k, chi2_lambda_k, chi2_ell_k);

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

            Modes=np.zeros(77,len(y.t))
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

class TaylorTn_1p5PN_Q : 
    def TaylorTn_1p5PN_Q(xHat_i, yHat_i, zHat_i, M1_i, M2_i, v_i, S_chi1_i, S_chi2_i, rfrak_chi1_x_i, rfrak_chi1_y_i, rfrak_chi2_x_i,
           rfrak_chi2_y_i, rfrak_frame_x_i, rfrak_frame_y_i, rfrak_frame_z_i) :
        xHat=xHat_i
        yHat=yHat_i
        zHat=zHat_i
        M1=M1_i
        M2=M2_i
        v=v_i
        S_chi1=S_chi1_i
        S_chi2=S_chi2_i
        rfrak_chi1_x=rfrak_chi1_x_i
        rfrak_chi1_y=rfrak_chi1_y_i
        rfrak_chi2_x=rfrak_chi2_x_i
        rfrak_chi2_y=rfrak_chi2_y_i
        rfrak_frame_x=rfrak_frame_x_i
        rfrak_frame_y=rfrak_frame_y_i
        rfrak_frame_z=rfrak_frame_z_i
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
        chi1chi1=dot(chiVec1,chiVec1)
        chi1chi2=dot(chiVec1,chiVec2)
        chi2chi2=dot(chiVec2,chiVec2)
        chi1_n=dot(chiVec1,nHat)
        chi1_lambda=dot(chiVec1,lambdaHat)
        chi1_ell=dot(chiVec1,ellHat)
        chi2_n=dot(chiVec2,nHat)
        chi2_lambda=dot(chiVec2,lambdaHat)
        chi2_ell=dot(chiVec2,ellHat)
        S_ell=M1**2*chi1_ell + M2**2*chi2_ell
        S_n=M1**2*chi1_n + M2**2*chi2_n
        S_lambda=M1**2*chi1_lambda + M2**2*chi2_lambda
        Sigma_ell=M*(-M1*chi1_ell + M2*chi2_ell)
        Sigma_n=M*(-M1*chi1_n + M2*chi2_n)
        Sigma_lambda=M*(-M1*chi1_lambda + M2*chi2_lambda)
        chi_s_ell=chi1_ell/2 + chi2_ell/2
        chi_a_ell=chi1_ell/2 - chi2_ell/2
        Fcal_coeff=32*nu**2*v**10/5
        Fcal_0=1
        Fcal_2=-35*nu/12 - 1247/336
        Fcal_3=4*pi
        Fcal_SO_3=(-4*S_ell - 5*Sigma_ell*delta/4)/M**2
        E_0=1
        E_2=-nu/12 - 3/4
        E_SO_3=(14*S_ell/3 + 2*Sigma_ell*delta)/M**2
        Phi=0.0
        EvolveSpin1=dot(S_chi1,S_chi1)>1e-12
        EvolveSpin2=dot(S_chi2,S_chi2)>1e-12

        def Recalculate(t, y):
            v = y[0];
            rfrak_chi1_x = y[1];
            rfrak_chi1_y = y[2];
            rfrak_chi2_x = y[3];
            rfrak_chi2_y = y[4];
            rfrak_frame_x = y[5];
            rfrak_frame_y = y[6];
            rfrak_frame_z = y[7];
            Phi = y[8];
            R = exp(rfrak_frame_x*xHat + rfrak_frame_y*yHat + rfrak_frame_z*zHat)
            nHat = R*xHat*conjugate(R)
            lambdaHat = R*yHat*conjugate(R)
            ellHat = R*zHat*conjugate(R)
            R_S1 = exp(rfrak_chi1_x*xHat + rfrak_chi1_y*yHat)
            R_S2 = exp(rfrak_chi2_x*xHat + rfrak_chi2_y*yHat)
            chiVec1 = S_chi1*R_S1*zHat*conjugate(R_S1)*conjugate(S_chi1)
            chiVec2 = S_chi2*R_S2*zHat*conjugate(R_S2)*conjugate(S_chi2)
            chi1_n = dot(chiVec1,nHat)
            chi1_lambda = dot(chiVec1,lambdaHat)
            chi1_ell = dot(chiVec1,ellHat)
            chi2_n = dot(chiVec2,nHat)
            chi2_lambda = dot(chiVec2,lambdaHat)
            chi2_ell = dot(chiVec2,ellHat)
            S_ell = M1**2*chi1_ell + M2**2*chi2_ell
            S_n = M1**2*chi1_n + M2**2*chi2_n
            S_lambda = M1**2*chi1_lambda + M2**2*chi2_lambda
            Sigma_ell = M*(-M1*chi1_ell + M2*chi2_ell)
            Sigma_n = M*(-M1*chi1_n + M2*chi2_n)
            Sigma_lambda = M*(-M1*chi1_lambda + M2*chi2_lambda)
            chi_s_ell = chi1_ell/2 + chi2_ell/2
            chi_a_ell = chi1_ell/2 - chi2_ell/2
            Fcal_coeff = 32*nu**2*v**10/5
            Fcal_SO_3 = (-4*S_ell - 5*Sigma_ell*delta/4)/M**2
            E_SO_3 = (14*S_ell/3 + 2*Sigma_ell*delta)/M**2

        def OmegaVec_chiVec_1():
            Omega1_coeff = v**5/M
            return Omega1_coeff*(ellHat*(-0.75*delta + 0.5*nu + v**2*(delta*(0.625*nu - 0.5625) + nu*(1.25 - 0.0416666666666667*nu) + 0.5625) + 0.75) + nHat*v*(3.0*chi1_n*nu + 3.0*M2**2*chi2_n/M**2) - M2**2*chiVec2*v/M**2)

        def OmegaVec_chiVec_2():
            Omega2_coeff = v**5/M
            return Omega2_coeff*(ellHat*(0.75*delta + 0.5*nu + v**2*(delta*(0.5625 - 0.625*nu) + nu*(1.25 - 0.0416666666666667*nu) + 0.5625) + 0.75) + nHat*v*(3.0*chi2_n*nu + 3.0*M1**2*chi1_n/M**2) - M1**2*chiVec1*v/M**2)

        def OmegaVec():
            gamma_PN_2 = 1.0 - 0.333333333333333*nu
            a_ell_2 = S_n*(-9.66666666666667*nu - 10.0) + Sigma_n*delta*(-4.5*nu - 6.0)
            a_ell_0 = 7.0*S_n + 3.0*Sigma_n*delta
            gamma_PN_3 = (1.66666666666667*S_ell + Sigma_ell*delta)/M**2
            gamma_PN_0 = 1.00000000000000
            return ellHat*v**3/M + nHat*v**6*(a_ell_0 + a_ell_2*v**2)*(gamma_PN_0 + v**2*(gamma_PN_2 + gamma_PN_3*v))/M**3

        def OrbitalAngularMomentum():
            L_2 = ellHat*(0.166666666666667*nu + 1.5)
            L_0 = ellHat
            L_SO_3 = 0.166666666666667*ellHat*(-35.0*S_ell - 15.0*Sigma_ell*delta)/M**2 + lambdaHat*(-3.0*S_lambda - Sigma_lambda*delta)/M**2 + 0.5*nHat*(S_n + Sigma_n*delta)/M**2
            L_coeff = M**2*nu/v
            return L_coeff*(L_0 + v**2*(L_2 + L_SO_3*v))


        #Evolve PN
        def TaylorT1(t, y):
            Recalculate(t, y)
            if v>=1.0: 
                return GSL_EDOM # Beyond domain of PN validity
            Flux = Fcal_coeff*(Fcal_0 + v**2*(Fcal_2 + v*(Fcal_3 + Fcal_SO_3)))
            dEdv = -0.5*M*nu*v*(2.0*E_0 + v**2*(4.0*E_2 + 5.0*E_SO_3*v))
            Absorption = 0
            dvdt_T1 = (-Absorption - Flux)/dEdv
            if dvdt_T1<1.0e-12:
                return GSL_EDIVERGE # v is decreasing
            return CommonRHS(dvdt_T1, y)

        def CommonRHS(dvdt, y):
            dydt=np.zeros(9)
            rfrak_frame=np.zeros(3)
            rfrak_frame[0] = y[5]
            rfrak_frame[1] = y[6]
            rfrak_frame[2] = y[7]
            rfrakdot_frame = Quaternions.FrameFromAngularVelocity_Integrand(rfrak_frame, OmegaVec().vec())
            dydt[0] = dvdt
            if(EvolveSpin1):
                Quaternions.FrameFromAngularVelocity_2D_Integrand(y[1], y[2],\
                    (S_chi1.inverse()*OmegaVec_chiVec_1()*S_chi1).vec(),dydt[1], dydt[2])
            else:
                dydt[1] = 0.0
                dydt[2] = 0.0
            if(EvolveSpin2):
                Quaternions.FrameFromAngularVelocity_2D_Integrand(y[3], y[4],\
                    (S_chi2.inverse()*OmegaVec_chiVec_2()*S_chi2).vec(),dydt[3], dydt[4])
            else:
                dydt[3] = 0.0
                dydt[4] = 0.0
            dydt[5] = rfrakdot_frame[0]
            dydt[6] = rfrakdot_frame[1]
            dydt[7] = rfrakdot_frame[2]
            dydt[8] = v*v*v/M

            return dydt
        
        
        y=solve_ivp(TaylorT1, [0,20000], [0.01**(1/3),rfrak_chi1_x,\
            rfrak_chi1_y,rfrak_chi2_x,rfrak_chi2_y,rfrak_frame_x,\
            rfrak_frame_y,rfrak_frame_z,Phi], method='DOP853',\
            t_eval=np.arange(0,20000,100), dense_output=True)
            

#        def TaylorT1_phi(t,phi):
#            dphidt_T1=v_of_t.y[0][v_of_t.t>=t][0]**3.0/M
#            return dphidt_T1
#        phi_of_t=solve_ivp(TaylorT1_phi, [0, v_of_t.t[-1]], [0.0], method='RK45', t_eval=v_of_t.t)
        
        #Calculate waveform modes

        def WaveformModes():
            I=1j
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
            chiVec1 = Quaternions.Quaternion(0., chi1_n_k, chi1_lambda_k, chi1_ell_k);
            chiVec2 = Quaternions.Quaternion(0., chi2_n_k, chi2_lambda_k, chi2_ell_k);

            R = exp(rfrak_frame_x*xHat + rfrak_frame_y*yHat + rfrak_frame_z*zHat)
            R_S1 = exp(rfrak_chi1_x*xHat + rfrak_chi1_y*yHat)
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

            Modes=np.zeros(77,len(y.t))
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

class TaylorTn_2p0PN_Q : 
    def TaylorTn_2p0PN_Q(xHat_i, yHat_i, zHat_i, M1_i, M2_i, v_i, S_chi1_i, S_chi2_i, rfrak_chi1_x_i, rfrak_chi1_y_i, rfrak_chi2_x_i,
           rfrak_chi2_y_i, rfrak_frame_x_i, rfrak_frame_y_i, rfrak_frame_z_i) :
        xHat=xHat_i
        yHat=yHat_i
        zHat=zHat_i
        M1=M1_i
        M2=M2_i
        v=v_i
        S_chi1=S_chi1_i
        S_chi2=S_chi2_i
        rfrak_chi1_x=rfrak_chi1_x_i
        rfrak_chi1_y=rfrak_chi1_y_i
        rfrak_chi2_x=rfrak_chi2_x_i
        rfrak_chi2_y=rfrak_chi2_y_i
        rfrak_frame_x=rfrak_frame_x_i
        rfrak_frame_y=rfrak_frame_y_i
        rfrak_frame_z=rfrak_frame_z_i
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
        chi1chi1=dot(chiVec1,chiVec1)
        chi1chi2=dot(chiVec1,chiVec2)
        chi2chi2=dot(chiVec2,chiVec2)
        chi1_n=dot(chiVec1,nHat)
        chi1_lambda=dot(chiVec1,lambdaHat)
        chi1_ell=dot(chiVec1,ellHat)
        chi2_n=dot(chiVec2,nHat)
        chi2_lambda=dot(chiVec2,lambdaHat)
        chi2_ell=dot(chiVec2,ellHat)
        S_ell=M1**2*chi1_ell + M2**2*chi2_ell
        S_n=M1**2*chi1_n + M2**2*chi2_n
        S_lambda=M1**2*chi1_lambda + M2**2*chi2_lambda
        Sigma_ell=M*(-M1*chi1_ell + M2*chi2_ell)
        Sigma_n=M*(-M1*chi1_n + M2*chi2_n)
        Sigma_lambda=M*(-M1*chi1_lambda + M2*chi2_lambda)
        chi_s_ell=chi1_ell/2 + chi2_ell/2
        chi_a_ell=chi1_ell/2 - chi2_ell/2
        Fcal_coeff=32*nu**2*v**10/5
        Fcal_0=1
        Fcal_2=-35*nu/12 - 1247/336
        Fcal_3=4*pi
        Fcal_4=65*nu**2/18 + 9271*nu/504 - 44711/9072
        Fcal_SQ_4=chi1chi1*(-89*delta/192 + 89*nu/96 - 89/192) - 103*chi1chi2*nu/48 + chi2chi2*(89*delta/192 + 89*nu/96 - 89/192) + chi_a_ell*(chi_a_ell*(287/96 - 12*nu) + 287*chi_s_ell*delta/48) + chi_s_ell**2*(nu/24 + 287/96)
        Fcal_SO_3=(-4*S_ell - 5*Sigma_ell*delta/4)/M**2
        E_0=1
        E_2=-nu/12 - 3/4
        E_4=-nu**2/24 + 19*nu/8 - 27/8
        E_SQ_4=-3*chi_a_ell**2/2 - 3*chi_s_ell**2/2 - delta*(chi2chi2/2 + 3*chi_a_ell*chi_s_ell) + nu*(chi1chi2 + 6*chi_a_ell**2) + (chi1chi1 + chi2chi2)*(delta - 2*nu + 1)/4
        E_SO_3=(14*S_ell/3 + 2*Sigma_ell*delta)/M**2
        Phi=0.0
        EvolveSpin1=dot(S_chi1,S_chi1)>1e-12
        EvolveSpin2=dot(S_chi2,S_chi2)>1e-12

        def Recalculate(t, y):
            v = y[0];
            rfrak_chi1_x = y[1];
            rfrak_chi1_y = y[2];
            rfrak_chi2_x = y[3];
            rfrak_chi2_y = y[4];
            rfrak_frame_x = y[5];
            rfrak_frame_y = y[6];
            rfrak_frame_z = y[7];
            Phi = y[8];
            R = exp(rfrak_frame_x*xHat + rfrak_frame_y*yHat + rfrak_frame_z*zHat)
            nHat = R*xHat*conjugate(R)
            lambdaHat = R*yHat*conjugate(R)
            ellHat = R*zHat*conjugate(R)
            R_S1 = exp(rfrak_chi1_x*xHat + rfrak_chi1_y*yHat)
            R_S2 = exp(rfrak_chi2_x*xHat + rfrak_chi2_y*yHat)
            chiVec1 = S_chi1*R_S1*zHat*conjugate(R_S1)*conjugate(S_chi1)
            chiVec2 = S_chi2*R_S2*zHat*conjugate(R_S2)*conjugate(S_chi2)
            chi1_n = dot(chiVec1,nHat)
            chi1_lambda = dot(chiVec1,lambdaHat)
            chi1_ell = dot(chiVec1,ellHat)
            chi2_n = dot(chiVec2,nHat)
            chi2_lambda = dot(chiVec2,lambdaHat)
            chi2_ell = dot(chiVec2,ellHat)
            S_ell = M1**2*chi1_ell + M2**2*chi2_ell
            S_n = M1**2*chi1_n + M2**2*chi2_n
            S_lambda = M1**2*chi1_lambda + M2**2*chi2_lambda
            Sigma_ell = M*(-M1*chi1_ell + M2*chi2_ell)
            Sigma_n = M*(-M1*chi1_n + M2*chi2_n)
            Sigma_lambda = M*(-M1*chi1_lambda + M2*chi2_lambda)
            chi_s_ell = chi1_ell/2 + chi2_ell/2
            chi_a_ell = chi1_ell/2 - chi2_ell/2
            Fcal_coeff = 32*nu**2*v**10/5
            Fcal_SQ_4 = chi1chi1*(-89*delta/192 + 89*nu/96 - 89/192) - 103*chi1chi2*nu/48 + chi2chi2*(89*delta/192 + 89*nu/96 - 89/192) + chi_a_ell*(chi_a_ell*(287/96 - 12*nu) + 287*chi_s_ell*delta/48) + chi_s_ell**2*(nu/24 + 287/96)
            Fcal_SO_3 = (-4*S_ell - 5*Sigma_ell*delta/4)/M**2
            E_SQ_4 = -3*chi_a_ell**2/2 - 3*chi_s_ell**2/2 - delta*(chi2chi2/2 + 3*chi_a_ell*chi_s_ell) + nu*(chi1chi2 + 6*chi_a_ell**2) + (chi1chi1 + chi2chi2)*(delta - 2*nu + 1)/4
            E_SO_3 = (14*S_ell/3 + 2*Sigma_ell*delta)/M**2

        def OmegaVec_chiVec_1():
            Omega1_coeff = v**5/M
            return Omega1_coeff*(ellHat*(-0.75*delta + 0.5*nu + v**2*(delta*(0.625*nu - 0.5625) + nu*(1.25 - 0.0416666666666667*nu) + v**2*(delta*(nu*(4.875 - 0.15625*nu) - 0.84375) + nu*(nu*(-0.0208333333333333*nu - 3.28125) + 0.1875) + 0.84375) + 0.5625) + 0.75) + nHat*v*(3.0*chi1_n*nu + 3.0*M2**2*chi2_n/M**2) - M2**2*chiVec2*v/M**2)

        def OmegaVec_chiVec_2():
            Omega2_coeff = v**5/M
            return Omega2_coeff*(ellHat*(0.75*delta + 0.5*nu + v**2*(delta*(0.5625 - 0.625*nu) + nu*(1.25 - 0.0416666666666667*nu) + v**2*(delta*(nu*(0.15625*nu - 4.875) + 0.84375) + nu*(nu*(-0.0208333333333333*nu - 3.28125) + 0.1875) + 0.84375) + 0.5625) + 0.75) + nHat*v*(3.0*chi2_n*nu + 3.0*M1**2*chi1_n/M**2) - M1**2*chiVec1*v/M**2)

        def OmegaVec():
            gamma_PN_2 = 1.0 - 0.333333333333333*nu
            a_ell_2 = S_n*(-9.66666666666667*nu - 10.0) + Sigma_n*delta*(-4.5*nu - 6.0)
            a_ell_0 = 7.0*S_n + 3.0*Sigma_n*delta
            a_ell_4 = S_n*(5.77777777777778*nu**2 + 14.75*nu + 1.5) + Sigma_n*delta*(2.83333333333333*nu**2 + 9.125*nu + 1.5)
            gamma_PN_4 = 1.0 - 5.41666666666667*nu
            gamma_PN_3 = (1.66666666666667*S_ell + Sigma_ell*delta)/M**2
            gamma_PN_0 = 1.00000000000000
            return ellHat*v**3/M + nHat*v**6*(a_ell_0 + v**2*(a_ell_2 + a_ell_4*v**2))*(gamma_PN_0 + v**2*(gamma_PN_2 + v*(gamma_PN_3 + gamma_PN_4*v)))/M**3

        def OrbitalAngularMomentum():
            L_4 = ellHat*(0.0416666666666667*nu**2 - 2.375*nu + 3.375)
            L_2 = ellHat*(0.166666666666667*nu + 1.5)
            L_0 = ellHat
            L_SO_3 = 0.166666666666667*ellHat*(-35.0*S_ell - 15.0*Sigma_ell*delta)/M**2 + lambdaHat*(-3.0*S_lambda - Sigma_lambda*delta)/M**2 + 0.5*nHat*(S_n + Sigma_n*delta)/M**2
            L_coeff = M**2*nu/v
            return L_coeff*(L_0 + v**2*(L_2 + v*(L_4*v + L_SO_3)))


        #Evolve PN
        def TaylorT1(t, y):
            Recalculate(t, y)
            if v>=1.0: 
                return GSL_EDOM # Beyond domain of PN validity
            Flux = Fcal_coeff*(Fcal_0 + v**2*(Fcal_2 + v*(Fcal_3 + Fcal_SO_3 + v*(Fcal_4 + Fcal_SQ_4))))
            dEdv = -0.5*M*nu*v*(2.0*E_0 + v**2*(4.0*E_2 + v*(5.0*E_SO_3 + 6.0*v*(E_4 + E_SQ_4))))
            Absorption = 0
            dvdt_T1 = (-Absorption - Flux)/dEdv
            if dvdt_T1<1.0e-12:
                return GSL_EDIVERGE # v is decreasing
            return CommonRHS(dvdt_T1, y)

        def CommonRHS(dvdt, y):
            dydt=np.zeros(9)
            rfrak_frame=np.zeros(3)
            rfrak_frame[0] = y[5]
            rfrak_frame[1] = y[6]
            rfrak_frame[2] = y[7]
            rfrakdot_frame = Quaternions.FrameFromAngularVelocity_Integrand(rfrak_frame, OmegaVec().vec())
            dydt[0] = dvdt
            if(EvolveSpin1):
                Quaternions.FrameFromAngularVelocity_2D_Integrand(y[1], y[2],\
                    (S_chi1.inverse()*OmegaVec_chiVec_1()*S_chi1).vec(),dydt[1], dydt[2])
            else:
                dydt[1] = 0.0
                dydt[2] = 0.0
            if(EvolveSpin2):
                Quaternions.FrameFromAngularVelocity_2D_Integrand(y[3], y[4],\
                    (S_chi2.inverse()*OmegaVec_chiVec_2()*S_chi2).vec(),dydt[3], dydt[4])
            else:
                dydt[3] = 0.0
                dydt[4] = 0.0
            dydt[5] = rfrakdot_frame[0]
            dydt[6] = rfrakdot_frame[1]
            dydt[7] = rfrakdot_frame[2]
            dydt[8] = v*v*v/M

            return dydt
        
        
        y=solve_ivp(TaylorT1, [0,20000], [0.01**(1/3),rfrak_chi1_x,\
            rfrak_chi1_y,rfrak_chi2_x,rfrak_chi2_y,rfrak_frame_x,\
            rfrak_frame_y,rfrak_frame_z,Phi], method='DOP853',\
            t_eval=np.arange(0,20000,100), dense_output=True)
            

#        def TaylorT1_phi(t,phi):
#            dphidt_T1=v_of_t.y[0][v_of_t.t>=t][0]**3.0/M
#            return dphidt_T1
#        phi_of_t=solve_ivp(TaylorT1_phi, [0, v_of_t.t[-1]], [0.0], method='RK45', t_eval=v_of_t.t)
        
        #Calculate waveform modes

        def WaveformModes():
            I=1j
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
            chiVec1 = Quaternions.Quaternion(0., chi1_n_k, chi1_lambda_k, chi1_ell_k);
            chiVec2 = Quaternions.Quaternion(0., chi2_n_k, chi2_lambda_k, chi2_ell_k);

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

            Modes=np.zeros(77,len(y.t))
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

class TaylorTn_2p5PN_Q : 
    def TaylorTn_2p5PN_Q(xHat_i, yHat_i, zHat_i, M1_i, M2_i, v_i, S_chi1_i, S_chi2_i, rfrak_chi1_x_i, rfrak_chi1_y_i, rfrak_chi2_x_i,
           rfrak_chi2_y_i, rfrak_frame_x_i, rfrak_frame_y_i, rfrak_frame_z_i) :
        xHat=xHat_i
        yHat=yHat_i
        zHat=zHat_i
        M1=M1_i
        M2=M2_i
        v=v_i
        S_chi1=S_chi1_i
        S_chi2=S_chi2_i
        rfrak_chi1_x=rfrak_chi1_x_i
        rfrak_chi1_y=rfrak_chi1_y_i
        rfrak_chi2_x=rfrak_chi2_x_i
        rfrak_chi2_y=rfrak_chi2_y_i
        rfrak_frame_x=rfrak_frame_x_i
        rfrak_frame_y=rfrak_frame_y_i
        rfrak_frame_z=rfrak_frame_z_i
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
        chi1chi1=dot(chiVec1,chiVec1)
        chi1chi2=dot(chiVec1,chiVec2)
        chi2chi2=dot(chiVec2,chiVec2)
        chi1_n=dot(chiVec1,nHat)
        chi1_lambda=dot(chiVec1,lambdaHat)
        chi1_ell=dot(chiVec1,ellHat)
        chi2_n=dot(chiVec2,nHat)
        chi2_lambda=dot(chiVec2,lambdaHat)
        chi2_ell=dot(chiVec2,ellHat)
        S_ell=M1**2*chi1_ell + M2**2*chi2_ell
        S_n=M1**2*chi1_n + M2**2*chi2_n
        S_lambda=M1**2*chi1_lambda + M2**2*chi2_lambda
        Sigma_ell=M*(-M1*chi1_ell + M2*chi2_ell)
        Sigma_n=M*(-M1*chi1_n + M2*chi2_n)
        Sigma_lambda=M*(-M1*chi1_lambda + M2*chi2_lambda)
        chi_s_ell=chi1_ell/2 + chi2_ell/2
        chi_a_ell=chi1_ell/2 - chi2_ell/2
        Fcal_coeff=32*nu**2*v**10/5
        Fcal_0=1
        Fcal_2=-35*nu/12 - 1247/336
        Fcal_3=4*pi
        Fcal_4=65*nu**2/18 + 9271*nu/504 - 44711/9072
        Fcal_5=pi*(-583*nu/24 - 8191/672)
        Fcal_SQ_4=chi1chi1*(-89*delta/192 + 89*nu/96 - 89/192) - 103*chi1chi2*nu/48 + chi2chi2*(89*delta/192 + 89*nu/96 - 89/192) + chi_a_ell*(chi_a_ell*(287/96 - 12*nu) + 287*chi_s_ell*delta/48) + chi_s_ell**2*(nu/24 + 287/96)
        Fcal_SO_3=(-4*S_ell - 5*Sigma_ell*delta/4)/M**2
        Fcal_SO_5=(S_ell*(272*nu/9 - 9/2) + Sigma_ell*delta*(43*nu/4 - 13/16))/M**2
        E_0=1
        E_2=-nu/12 - 3/4
        E_4=-nu**2/24 + 19*nu/8 - 27/8
        E_SQ_4=-3*chi_a_ell**2/2 - 3*chi_s_ell**2/2 - delta*(chi2chi2/2 + 3*chi_a_ell*chi_s_ell) + nu*(chi1chi2 + 6*chi_a_ell**2) + (chi1chi1 + chi2chi2)*(delta - 2*nu + 1)/4
        E_SO_3=(14*S_ell/3 + 2*Sigma_ell*delta)/M**2
        E_SO_5=(S_ell*(11 - 61*nu/9) + Sigma_ell*delta*(3 - 10*nu/3))/M**2
        Phi=0.0
        EvolveSpin1=dot(S_chi1,S_chi1)>1e-12
        EvolveSpin2=dot(S_chi2,S_chi2)>1e-12

        def Recalculate(t, y):
            v = y[0];
            rfrak_chi1_x = y[1];
            rfrak_chi1_y = y[2];
            rfrak_chi2_x = y[3];
            rfrak_chi2_y = y[4];
            rfrak_frame_x = y[5];
            rfrak_frame_y = y[6];
            rfrak_frame_z = y[7];
            Phi = y[8];
            R = exp(rfrak_frame_x*xHat + rfrak_frame_y*yHat + rfrak_frame_z*zHat)
            nHat = R*xHat*conjugate(R)
            lambdaHat = R*yHat*conjugate(R)
            ellHat = R*zHat*conjugate(R)
            R_S1 = exp(rfrak_chi1_x*xHat + rfrak_chi1_y*yHat)
            R_S2 = exp(rfrak_chi2_x*xHat + rfrak_chi2_y*yHat)
            chiVec1 = S_chi1*R_S1*zHat*conjugate(R_S1)*conjugate(S_chi1)
            chiVec2 = S_chi2*R_S2*zHat*conjugate(R_S2)*conjugate(S_chi2)
            chi1_n = dot(chiVec1,nHat)
            chi1_lambda = dot(chiVec1,lambdaHat)
            chi1_ell = dot(chiVec1,ellHat)
            chi2_n = dot(chiVec2,nHat)
            chi2_lambda = dot(chiVec2,lambdaHat)
            chi2_ell = dot(chiVec2,ellHat)
            S_ell = M1**2*chi1_ell + M2**2*chi2_ell
            S_n = M1**2*chi1_n + M2**2*chi2_n
            S_lambda = M1**2*chi1_lambda + M2**2*chi2_lambda
            Sigma_ell = M*(-M1*chi1_ell + M2*chi2_ell)
            Sigma_n = M*(-M1*chi1_n + M2*chi2_n)
            Sigma_lambda = M*(-M1*chi1_lambda + M2*chi2_lambda)
            chi_s_ell = chi1_ell/2 + chi2_ell/2
            chi_a_ell = chi1_ell/2 - chi2_ell/2
            Fcal_coeff = 32*nu**2*v**10/5
            Fcal_SQ_4 = chi1chi1*(-89*delta/192 + 89*nu/96 - 89/192) - 103*chi1chi2*nu/48 + chi2chi2*(89*delta/192 + 89*nu/96 - 89/192) + chi_a_ell*(chi_a_ell*(287/96 - 12*nu) + 287*chi_s_ell*delta/48) + chi_s_ell**2*(nu/24 + 287/96)
            Fcal_SO_3 = (-4*S_ell - 5*Sigma_ell*delta/4)/M**2
            Fcal_SO_5 = (S_ell*(272*nu/9 - 9/2) + Sigma_ell*delta*(43*nu/4 - 13/16))/M**2
            E_SQ_4 = -3*chi_a_ell**2/2 - 3*chi_s_ell**2/2 - delta*(chi2chi2/2 + 3*chi_a_ell*chi_s_ell) + nu*(chi1chi2 + 6*chi_a_ell**2) + (chi1chi1 + chi2chi2)*(delta - 2*nu + 1)/4
            E_SO_3 = (14*S_ell/3 + 2*Sigma_ell*delta)/M**2
            E_SO_5 = (S_ell*(11 - 61*nu/9) + Sigma_ell*delta*(3 - 10*nu/3))/M**2

        def OmegaVec_chiVec_1():
            Omega1_coeff = v**5/M
            return Omega1_coeff*(ellHat*(-0.75*delta + 0.5*nu + v**2*(delta*(0.625*nu - 0.5625) + nu*(1.25 - 0.0416666666666667*nu) + v**2*(delta*(nu*(4.875 - 0.15625*nu) - 0.84375) + nu*(nu*(-0.0208333333333333*nu - 3.28125) + 0.1875) + 0.84375) + 0.5625) + 0.75) + nHat*v*(3.0*chi1_n*nu + 3.0*M2**2*chi2_n/M**2) - M2**2*chiVec2*v/M**2)

        def OmegaVec_chiVec_2():
            Omega2_coeff = v**5/M
            return Omega2_coeff*(ellHat*(0.75*delta + 0.5*nu + v**2*(delta*(0.5625 - 0.625*nu) + nu*(1.25 - 0.0416666666666667*nu) + v**2*(delta*(nu*(0.15625*nu - 4.875) + 0.84375) + nu*(nu*(-0.0208333333333333*nu - 3.28125) + 0.1875) + 0.84375) + 0.5625) + 0.75) + nHat*v*(3.0*chi2_n*nu + 3.0*M1**2*chi1_n/M**2) - M1**2*chiVec1*v/M**2)

        def OmegaVec():
            gamma_PN_2 = 1.0 - 0.333333333333333*nu
            a_ell_2 = S_n*(-9.66666666666667*nu - 10.0) + Sigma_n*delta*(-4.5*nu - 6.0)
            a_ell_0 = 7.0*S_n + 3.0*Sigma_n*delta
            a_ell_4 = S_n*(5.77777777777778*nu**2 + 14.75*nu + 1.5) + Sigma_n*delta*(2.83333333333333*nu**2 + 9.125*nu + 1.5)
            gamma_PN_4 = 1.0 - 5.41666666666667*nu
            gamma_PN_3 = (1.66666666666667*S_ell + Sigma_ell*delta)/M**2
            gamma_PN_5 = (S_ell*(0.888888888888889*nu + 3.33333333333333) + 2.0*Sigma_ell*delta)/M**2
            gamma_PN_0 = 1.00000000000000
            return ellHat*v**3/M + nHat*v**6*(a_ell_0 + v**2*(a_ell_2 + a_ell_4*v**2))*(gamma_PN_0 + v**2*(gamma_PN_2 + v*(gamma_PN_3 + v*(gamma_PN_4 + gamma_PN_5*v))))/M**3

        def OrbitalAngularMomentum():
            L_4 = ellHat*(0.0416666666666667*nu**2 - 2.375*nu + 3.375)
            L_SO_5 = 0.0138888888888889*ellHat*(427.0*S_ell*nu - 693.0*S_ell + 210.0*Sigma_ell*delta*nu - 189.0*Sigma_ell*delta)/M**2 + 0.166666666666667*lambdaHat*(18.0*S_lambda*nu - 21.0*S_lambda + 8.0*Sigma_lambda*delta*nu - 3.0*Sigma_lambda*delta)/M**2 + 0.0416666666666667*nHat*(-19.0*S_n*nu + 33.0*S_n - 10.0*Sigma_n*delta*nu + 33.0*Sigma_n*delta)/M**2
            L_2 = ellHat*(0.166666666666667*nu + 1.5)
            L_0 = ellHat
            L_SO_3 = 0.166666666666667*ellHat*(-35.0*S_ell - 15.0*Sigma_ell*delta)/M**2 + lambdaHat*(-3.0*S_lambda - Sigma_lambda*delta)/M**2 + 0.5*nHat*(S_n + Sigma_n*delta)/M**2
            L_coeff = M**2*nu/v
            return L_coeff*(L_0 + v**2*(L_2 + v*(L_SO_3 + v*(L_4 + L_SO_5*v))))


        #Evolve PN
        def TaylorT1(t, y):
            Recalculate(t, y)
            if v>=1.0: 
                return GSL_EDOM # Beyond domain of PN validity
            Flux = Fcal_coeff*(Fcal_0 + v**2*(Fcal_2 + v*(Fcal_3 + Fcal_SO_3 + v*(Fcal_4 + Fcal_SQ_4 + v*(Fcal_5 + Fcal_SO_5)))))
            dEdv = -0.5*M*nu*v*(2.0*E_0 + v**2*(4.0*E_2 + v*(5.0*E_SO_3 + v*(6.0*E_4 + 7.0*E_SO_5*v + 6.0*E_SQ_4))))
            Absorption = 0
            dvdt_T1 = (-Absorption - Flux)/dEdv
            if dvdt_T1<1.0e-12:
                return GSL_EDIVERGE # v is decreasing
            return CommonRHS(dvdt_T1, y)

        def CommonRHS(dvdt, y):
            dydt=np.zeros(9)
            rfrak_frame=np.zeros(3)
            rfrak_frame[0] = y[5]
            rfrak_frame[1] = y[6]
            rfrak_frame[2] = y[7]
            rfrakdot_frame = Quaternions.FrameFromAngularVelocity_Integrand(rfrak_frame, OmegaVec().vec())
            dydt[0] = dvdt
            if(EvolveSpin1):
                Quaternions.FrameFromAngularVelocity_2D_Integrand(y[1], y[2],\
                    (S_chi1.inverse()*OmegaVec_chiVec_1()*S_chi1).vec(),dydt[1], dydt[2])
            else:
                dydt[1] = 0.0
                dydt[2] = 0.0
            if(EvolveSpin2):
                Quaternions.FrameFromAngularVelocity_2D_Integrand(y[3], y[4],\
                    (S_chi2.inverse()*OmegaVec_chiVec_2()*S_chi2).vec(),dydt[3], dydt[4])
            else:
                dydt[3] = 0.0
                dydt[4] = 0.0
            dydt[5] = rfrakdot_frame[0]
            dydt[6] = rfrakdot_frame[1]
            dydt[7] = rfrakdot_frame[2]
            dydt[8] = v*v*v/M

            return dydt
        
        
        y=solve_ivp(TaylorT1, [0,20000], [0.01**(1/3),rfrak_chi1_x,\
            rfrak_chi1_y,rfrak_chi2_x,rfrak_chi2_y,rfrak_frame_x,\
            rfrak_frame_y,rfrak_frame_z,Phi], method='DOP853',\
            t_eval=np.arange(0,20000,100), dense_output=True)
            

#        def TaylorT1_phi(t,phi):
#            dphidt_T1=v_of_t.y[0][v_of_t.t>=t][0]**3.0/M
#            return dphidt_T1
#        phi_of_t=solve_ivp(TaylorT1_phi, [0, v_of_t.t[-1]], [0.0], method='RK45', t_eval=v_of_t.t)
        
        #Calculate waveform modes

        def WaveformModes():
            I=1j
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
            chiVec1 = Quaternions.Quaternion(0., chi1_n_k, chi1_lambda_k, chi1_ell_k);
            chiVec2 = Quaternions.Quaternion(0., chi2_n_k, chi2_lambda_k, chi2_ell_k);

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

            Modes=np.zeros(77,len(y.t))
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

class TaylorTn_3p0PN_Q : 
    def TaylorTn_3p0PN_Q(xHat_i, yHat_i, zHat_i, M1_i, M2_i, v_i, S_chi1_i, S_chi2_i, rfrak_chi1_x_i, rfrak_chi1_y_i, rfrak_chi2_x_i,
           rfrak_chi2_y_i, rfrak_frame_x_i, rfrak_frame_y_i, rfrak_frame_z_i) :
        xHat=xHat_i
        yHat=yHat_i
        zHat=zHat_i
        M1=M1_i
        M2=M2_i
        v=v_i
        S_chi1=S_chi1_i
        S_chi2=S_chi2_i
        rfrak_chi1_x=rfrak_chi1_x_i
        rfrak_chi1_y=rfrak_chi1_y_i
        rfrak_chi2_x=rfrak_chi2_x_i
        rfrak_chi2_y=rfrak_chi2_y_i
        rfrak_frame_x=rfrak_frame_x_i
        rfrak_frame_y=rfrak_frame_y_i
        rfrak_frame_z=rfrak_frame_z_i
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
        chi1chi1=dot(chiVec1,chiVec1)
        chi1chi2=dot(chiVec1,chiVec2)
        chi2chi2=dot(chiVec2,chiVec2)
        chi1_n=dot(chiVec1,nHat)
        chi1_lambda=dot(chiVec1,lambdaHat)
        chi1_ell=dot(chiVec1,ellHat)
        chi2_n=dot(chiVec2,nHat)
        chi2_lambda=dot(chiVec2,lambdaHat)
        chi2_ell=dot(chiVec2,ellHat)
        S_ell=M1**2*chi1_ell + M2**2*chi2_ell
        S_n=M1**2*chi1_n + M2**2*chi2_n
        S_lambda=M1**2*chi1_lambda + M2**2*chi2_lambda
        Sigma_ell=M*(-M1*chi1_ell + M2*chi2_ell)
        Sigma_n=M*(-M1*chi1_n + M2*chi2_n)
        Sigma_lambda=M*(-M1*chi1_lambda + M2*chi2_lambda)
        chi_s_ell=chi1_ell/2 + chi2_ell/2
        chi_a_ell=chi1_ell/2 - chi2_ell/2
        logv=log(v)
        Fcal_coeff=32*nu**2*v**10/5
        Fcal_0=1
        Fcal_2=-35*nu/12 - 1247/336
        Fcal_3=4*pi
        Fcal_4=65*nu**2/18 + 9271*nu/504 - 44711/9072
        Fcal_5=pi*(-583*nu/24 - 8191/672)
        Fcal_6=-775*nu**3/324 - 94403*nu**2/3024 + nu*(-134543/7776 + 41*pi**2/48) - 1712*log(4)/105 - 1712*EulerGamma/105 + 16*pi**2/3 + 6643739519/69854400
        Fcal_lnv_6=-1712/105
        Fcal_SQ_4=chi1chi1*(-89*delta/192 + 89*nu/96 - 89/192) - 103*chi1chi2*nu/48 + chi2chi2*(89*delta/192 + 89*nu/96 - 89/192) + chi_a_ell*(chi_a_ell*(287/96 - 12*nu) + 287*chi_s_ell*delta/48) + chi_s_ell**2*(nu/24 + 287/96)
        Fcal_SO_3=(-4*S_ell - 5*Sigma_ell*delta/4)/M**2
        Fcal_SO_5=(S_ell*(272*nu/9 - 9/2) + Sigma_ell*delta*(43*nu/4 - 13/16))/M**2
        Fcal_SO_6=(-16*S_ell*pi - 31*Sigma_ell*delta*pi/6)/M**2
        E_0=1
        E_2=-nu/12 - 3/4
        E_4=-nu**2/24 + 19*nu/8 - 27/8
        E_6=-35*nu**3/5184 - 155*nu**2/96 + nu*(34445/576 - 205*pi**2/96) - 675/64
        E_SQ_4=-3*chi_a_ell**2/2 - 3*chi_s_ell**2/2 - delta*(chi2chi2/2 + 3*chi_a_ell*chi_s_ell) + nu*(chi1chi2 + 6*chi_a_ell**2) + (chi1chi1 + chi2chi2)*(delta - 2*nu + 1)/4
        E_SO_3=(14*S_ell/3 + 2*Sigma_ell*delta)/M**2
        E_SO_5=(S_ell*(11 - 61*nu/9) + Sigma_ell*delta*(3 - 10*nu/3))/M**2
        Phi=0.0
        EvolveSpin1=dot(S_chi1,S_chi1)>1e-12
        EvolveSpin2=dot(S_chi2,S_chi2)>1e-12

        def Recalculate(t, y):
            v = y[0];
            rfrak_chi1_x = y[1];
            rfrak_chi1_y = y[2];
            rfrak_chi2_x = y[3];
            rfrak_chi2_y = y[4];
            rfrak_frame_x = y[5];
            rfrak_frame_y = y[6];
            rfrak_frame_z = y[7];
            Phi = y[8];
            R = exp(rfrak_frame_x*xHat + rfrak_frame_y*yHat + rfrak_frame_z*zHat)
            nHat = R*xHat*conjugate(R)
            lambdaHat = R*yHat*conjugate(R)
            ellHat = R*zHat*conjugate(R)
            R_S1 = exp(rfrak_chi1_x*xHat + rfrak_chi1_y*yHat)
            R_S2 = exp(rfrak_chi2_x*xHat + rfrak_chi2_y*yHat)
            chiVec1 = S_chi1*R_S1*zHat*conjugate(R_S1)*conjugate(S_chi1)
            chiVec2 = S_chi2*R_S2*zHat*conjugate(R_S2)*conjugate(S_chi2)
            chi1_n = dot(chiVec1,nHat)
            chi1_lambda = dot(chiVec1,lambdaHat)
            chi1_ell = dot(chiVec1,ellHat)
            chi2_n = dot(chiVec2,nHat)
            chi2_lambda = dot(chiVec2,lambdaHat)
            chi2_ell = dot(chiVec2,ellHat)
            S_ell = M1**2*chi1_ell + M2**2*chi2_ell
            S_n = M1**2*chi1_n + M2**2*chi2_n
            S_lambda = M1**2*chi1_lambda + M2**2*chi2_lambda
            Sigma_ell = M*(-M1*chi1_ell + M2*chi2_ell)
            Sigma_n = M*(-M1*chi1_n + M2*chi2_n)
            Sigma_lambda = M*(-M1*chi1_lambda + M2*chi2_lambda)
            chi_s_ell = chi1_ell/2 + chi2_ell/2
            chi_a_ell = chi1_ell/2 - chi2_ell/2
            logv = log(v)
            Fcal_coeff = 32*nu**2*v**10/5
            Fcal_SQ_4 = chi1chi1*(-89*delta/192 + 89*nu/96 - 89/192) - 103*chi1chi2*nu/48 + chi2chi2*(89*delta/192 + 89*nu/96 - 89/192) + chi_a_ell*(chi_a_ell*(287/96 - 12*nu) + 287*chi_s_ell*delta/48) + chi_s_ell**2*(nu/24 + 287/96)
            Fcal_SO_3 = (-4*S_ell - 5*Sigma_ell*delta/4)/M**2
            Fcal_SO_5 = (S_ell*(272*nu/9 - 9/2) + Sigma_ell*delta*(43*nu/4 - 13/16))/M**2
            Fcal_SO_6 = (-16*S_ell*pi - 31*Sigma_ell*delta*pi/6)/M**2
            E_SQ_4 = -3*chi_a_ell**2/2 - 3*chi_s_ell**2/2 - delta*(chi2chi2/2 + 3*chi_a_ell*chi_s_ell) + nu*(chi1chi2 + 6*chi_a_ell**2) + (chi1chi1 + chi2chi2)*(delta - 2*nu + 1)/4
            E_SO_3 = (14*S_ell/3 + 2*Sigma_ell*delta)/M**2
            E_SO_5 = (S_ell*(11 - 61*nu/9) + Sigma_ell*delta*(3 - 10*nu/3))/M**2

        def OmegaVec_chiVec_1():
            Omega1_coeff = v**5/M
            return Omega1_coeff*(ellHat*(-0.75*delta + 0.5*nu + v**2*(delta*(0.625*nu - 0.5625) + nu*(1.25 - 0.0416666666666667*nu) + v**2*(delta*(nu*(4.875 - 0.15625*nu) - 0.84375) + nu*(nu*(-0.0208333333333333*nu - 3.28125) + 0.1875) + 0.84375) + 0.5625) + 0.75) + nHat*v*(3.0*chi1_n*nu + 3.0*M2**2*chi2_n/M**2) - M2**2*chiVec2*v/M**2)

        def OmegaVec_chiVec_2():
            Omega2_coeff = v**5/M
            return Omega2_coeff*(ellHat*(0.75*delta + 0.5*nu + v**2*(delta*(0.5625 - 0.625*nu) + nu*(1.25 - 0.0416666666666667*nu) + v**2*(delta*(nu*(0.15625*nu - 4.875) + 0.84375) + nu*(nu*(-0.0208333333333333*nu - 3.28125) + 0.1875) + 0.84375) + 0.5625) + 0.75) + nHat*v*(3.0*chi2_n*nu + 3.0*M1**2*chi1_n/M**2) - M1**2*chiVec1*v/M**2)

        def OmegaVec():
            gamma_PN_6 = 0.0123456790123457*nu**3 + 6.36111111111111*nu**2 - 2.98177812235564*nu + 1.0
            gamma_PN_2 = 1.0 - 0.333333333333333*nu
            a_ell_2 = S_n*(-9.66666666666667*nu - 10.0) + Sigma_n*delta*(-4.5*nu - 6.0)
            a_ell_0 = 7.0*S_n + 3.0*Sigma_n*delta
            a_ell_4 = S_n*(5.77777777777778*nu**2 + 14.75*nu + 1.5) + Sigma_n*delta*(2.83333333333333*nu**2 + 9.125*nu + 1.5)
            gamma_PN_4 = 1.0 - 5.41666666666667*nu
            gamma_PN_3 = (1.66666666666667*S_ell + Sigma_ell*delta)/M**2
            gamma_PN_5 = (S_ell*(0.888888888888889*nu + 3.33333333333333) + 2.0*Sigma_ell*delta)/M**2
            gamma_PN_0 = 1.00000000000000
            return ellHat*v**3/M + nHat*v**6*(a_ell_0 + v**2*(a_ell_2 + a_ell_4*v**2))*(gamma_PN_0 + v**2*(gamma_PN_2 + v*(gamma_PN_3 + v*(gamma_PN_4 + v*(gamma_PN_5 + gamma_PN_6*v)))))/M**3

        def OrbitalAngularMomentum():
            L_6 = ellHat*(0.00540123456790123*nu**3 + 1.29166666666667*nu**2 - 30.9797035925835*nu + 8.4375)
            L_4 = ellHat*(0.0416666666666667*nu**2 - 2.375*nu + 3.375)
            L_SO_5 = 0.0138888888888889*ellHat*(427.0*S_ell*nu - 693.0*S_ell + 210.0*Sigma_ell*delta*nu - 189.0*Sigma_ell*delta)/M**2 + 0.166666666666667*lambdaHat*(18.0*S_lambda*nu - 21.0*S_lambda + 8.0*Sigma_lambda*delta*nu - 3.0*Sigma_lambda*delta)/M**2 + 0.0416666666666667*nHat*(-19.0*S_n*nu + 33.0*S_n - 10.0*Sigma_n*delta*nu + 33.0*Sigma_n*delta)/M**2
            L_2 = ellHat*(0.166666666666667*nu + 1.5)
            L_0 = ellHat
            L_SO_3 = 0.166666666666667*ellHat*(-35.0*S_ell - 15.0*Sigma_ell*delta)/M**2 + lambdaHat*(-3.0*S_lambda - Sigma_lambda*delta)/M**2 + 0.5*nHat*(S_n + Sigma_n*delta)/M**2
            L_coeff = M**2*nu/v
            return L_coeff*(L_0 + v**2*(L_2 + v*(L_SO_3 + v*(L_4 + v*(L_6*v + L_SO_5)))))


        #Evolve PN
        def TaylorT1(t, y):
            Recalculate(t, y)
            if v>=1.0: 
                return GSL_EDOM # Beyond domain of PN validity
            Flux = Fcal_coeff*(Fcal_0 + v**2*(Fcal_2 + v*(Fcal_3 + Fcal_SO_3 + v*(Fcal_4 + Fcal_SQ_4 + v*(Fcal_5 + Fcal_SO_5 + v*(Fcal_6 + Fcal_SO_6 + Fcal_lnv_6*logv))))))
            dEdv = -0.5*M*nu*v*(2.0*E_0 + v**2*(4.0*E_2 + v*(5.0*E_SO_3 + v*(6.0*E_4 + 6.0*E_SQ_4 + v*(8.0*E_6*v + 7.0*E_SO_5)))))
            Absorption = 0
            dvdt_T1 = (-Absorption - Flux)/dEdv
            if dvdt_T1<1.0e-12:
                return GSL_EDIVERGE # v is decreasing
            return CommonRHS(dvdt_T1, y)

        def CommonRHS(dvdt, y):
            dydt=np.zeros(9)
            rfrak_frame=np.zeros(3)
            rfrak_frame[0] = y[5]
            rfrak_frame[1] = y[6]
            rfrak_frame[2] = y[7]
            rfrakdot_frame = Quaternions.FrameFromAngularVelocity_Integrand(rfrak_frame, OmegaVec().vec())
            dydt[0] = dvdt
            if(EvolveSpin1):
                Quaternions.FrameFromAngularVelocity_2D_Integrand(y[1], y[2],\
                    (S_chi1.inverse()*OmegaVec_chiVec_1()*S_chi1).vec(),dydt[1], dydt[2])
            else:
                dydt[1] = 0.0
                dydt[2] = 0.0
            if(EvolveSpin2):
                Quaternions.FrameFromAngularVelocity_2D_Integrand(y[3], y[4],\
                    (S_chi2.inverse()*OmegaVec_chiVec_2()*S_chi2).vec(),dydt[3], dydt[4])
            else:
                dydt[3] = 0.0
                dydt[4] = 0.0
            dydt[5] = rfrakdot_frame[0]
            dydt[6] = rfrakdot_frame[1]
            dydt[7] = rfrakdot_frame[2]
            dydt[8] = v*v*v/M

            return dydt
        
        
        y=solve_ivp(TaylorT1, [0,20000], [0.01**(1/3),rfrak_chi1_x,\
            rfrak_chi1_y,rfrak_chi2_x,rfrak_chi2_y,rfrak_frame_x,\
            rfrak_frame_y,rfrak_frame_z,Phi], method='DOP853',\
            t_eval=np.arange(0,20000,100), dense_output=True)
            

#        def TaylorT1_phi(t,phi):
#            dphidt_T1=v_of_t.y[0][v_of_t.t>=t][0]**3.0/M
#            return dphidt_T1
#        phi_of_t=solve_ivp(TaylorT1_phi, [0, v_of_t.t[-1]], [0.0], method='RK45', t_eval=v_of_t.t)
        
        #Calculate waveform modes

        def WaveformModes():
            I=1j
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
            chiVec1 = Quaternions.Quaternion(0., chi1_n_k, chi1_lambda_k, chi1_ell_k);
            chiVec2 = Quaternions.Quaternion(0., chi2_n_k, chi2_lambda_k, chi2_ell_k);

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

            Modes=np.zeros(77,len(y.t))
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

class TaylorTn_3p5PN_Q : 
    def TaylorTn_3p5PN_Q(xHat_i, yHat_i, zHat_i, M1_i, M2_i, v_i, S_chi1_i, S_chi2_i, rfrak_chi1_x_i, rfrak_chi1_y_i, rfrak_chi2_x_i,
           rfrak_chi2_y_i, rfrak_frame_x_i, rfrak_frame_y_i, rfrak_frame_z_i) :
        xHat=xHat_i
        yHat=yHat_i
        zHat=zHat_i
        M1=M1_i
        M2=M2_i
        v=v_i
        S_chi1=S_chi1_i
        S_chi2=S_chi2_i
        rfrak_chi1_x=rfrak_chi1_x_i
        rfrak_chi1_y=rfrak_chi1_y_i
        rfrak_chi2_x=rfrak_chi2_x_i
        rfrak_chi2_y=rfrak_chi2_y_i
        rfrak_frame_x=rfrak_frame_x_i
        rfrak_frame_y=rfrak_frame_y_i
        rfrak_frame_z=rfrak_frame_z_i
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
        chi1chi1=dot(chiVec1,chiVec1)
        chi1chi2=dot(chiVec1,chiVec2)
        chi2chi2=dot(chiVec2,chiVec2)
        chi1_n=dot(chiVec1,nHat)
        chi1_lambda=dot(chiVec1,lambdaHat)
        chi1_ell=dot(chiVec1,ellHat)
        chi2_n=dot(chiVec2,nHat)
        chi2_lambda=dot(chiVec2,lambdaHat)
        chi2_ell=dot(chiVec2,ellHat)
        S_ell=M1**2*chi1_ell + M2**2*chi2_ell
        S_n=M1**2*chi1_n + M2**2*chi2_n
        S_lambda=M1**2*chi1_lambda + M2**2*chi2_lambda
        Sigma_ell=M*(-M1*chi1_ell + M2*chi2_ell)
        Sigma_n=M*(-M1*chi1_n + M2*chi2_n)
        Sigma_lambda=M*(-M1*chi1_lambda + M2*chi2_lambda)
        chi_s_ell=chi1_ell/2 + chi2_ell/2
        chi_a_ell=chi1_ell/2 - chi2_ell/2
        logv=log(v)
        Fcal_coeff=32*nu**2*v**10/5
        Fcal_0=1
        Fcal_2=-35*nu/12 - 1247/336
        Fcal_3=4*pi
        Fcal_4=65*nu**2/18 + 9271*nu/504 - 44711/9072
        Fcal_5=pi*(-583*nu/24 - 8191/672)
        Fcal_6=-775*nu**3/324 - 94403*nu**2/3024 + nu*(-134543/7776 + 41*pi**2/48) - 1712*log(4)/105 - 1712*EulerGamma/105 + 16*pi**2/3 + 6643739519/69854400
        Fcal_lnv_6=-1712/105
        Fcal_7=pi*(193385*nu**2/3024 + 214745*nu/1728 - 16285/504)
        Fcal_SQ_4=chi1chi1*(-89*delta/192 + 89*nu/96 - 89/192) - 103*chi1chi2*nu/48 + chi2chi2*(89*delta/192 + 89*nu/96 - 89/192) + chi_a_ell*(chi_a_ell*(287/96 - 12*nu) + 287*chi_s_ell*delta/48) + chi_s_ell**2*(nu/24 + 287/96)
        Fcal_SO_3=(-4*S_ell - 5*Sigma_ell*delta/4)/M**2
        Fcal_SO_5=(S_ell*(272*nu/9 - 9/2) + Sigma_ell*delta*(43*nu/4 - 13/16))/M**2
        Fcal_SO_6=(-16*S_ell*pi - 31*Sigma_ell*delta*pi/6)/M**2
        Fcal_SO_7=(S_ell*(-2810*nu**2/27 + 6172*nu/189 + 476645/6804) + Sigma_ell*delta*(-1501*nu**2/36 + 1849*nu/126 + 9535/336))/M**2
        E_0=1
        E_2=-nu/12 - 3/4
        E_4=-nu**2/24 + 19*nu/8 - 27/8
        E_6=-35*nu**3/5184 - 155*nu**2/96 + nu*(34445/576 - 205*pi**2/96) - 675/64
        E_SQ_4=-3*chi_a_ell**2/2 - 3*chi_s_ell**2/2 - delta*(chi2chi2/2 + 3*chi_a_ell*chi_s_ell) + nu*(chi1chi2 + 6*chi_a_ell**2) + (chi1chi1 + chi2chi2)*(delta - 2*nu + 1)/4
        E_SO_3=(14*S_ell/3 + 2*Sigma_ell*delta)/M**2
        E_SO_5=(S_ell*(11 - 61*nu/9) + Sigma_ell*delta*(3 - 10*nu/3))/M**2
        E_SO_7=(S_ell*(29*nu**2/12 - 367*nu/4 + 135/4) + Sigma_ell*delta*(5*nu**2/4 - 39*nu + 27/4))/M**2
        Phi=0.0
        EvolveSpin1=dot(S_chi1,S_chi1)>1e-12
        EvolveSpin2=dot(S_chi2,S_chi2)>1e-12

        def Recalculate(t, y):
            v = y[0];
            rfrak_chi1_x = y[1];
            rfrak_chi1_y = y[2];
            rfrak_chi2_x = y[3];
            rfrak_chi2_y = y[4];
            rfrak_frame_x = y[5];
            rfrak_frame_y = y[6];
            rfrak_frame_z = y[7];
            Phi = y[8];
            R = exp(rfrak_frame_x*xHat + rfrak_frame_y*yHat + rfrak_frame_z*zHat)
            nHat = R*xHat*conjugate(R)
            lambdaHat = R*yHat*conjugate(R)
            ellHat = R*zHat*conjugate(R)
            R_S1 = exp(rfrak_chi1_x*xHat + rfrak_chi1_y*yHat)
            R_S2 = exp(rfrak_chi2_x*xHat + rfrak_chi2_y*yHat)
            chiVec1 = S_chi1*R_S1*zHat*conjugate(R_S1)*conjugate(S_chi1)
            chiVec2 = S_chi2*R_S2*zHat*conjugate(R_S2)*conjugate(S_chi2)
            chi1_n = dot(chiVec1,nHat)
            chi1_lambda = dot(chiVec1,lambdaHat)
            chi1_ell = dot(chiVec1,ellHat)
            chi2_n = dot(chiVec2,nHat)
            chi2_lambda = dot(chiVec2,lambdaHat)
            chi2_ell = dot(chiVec2,ellHat)
            S_ell = M1**2*chi1_ell + M2**2*chi2_ell
            S_n = M1**2*chi1_n + M2**2*chi2_n
            S_lambda = M1**2*chi1_lambda + M2**2*chi2_lambda
            Sigma_ell = M*(-M1*chi1_ell + M2*chi2_ell)
            Sigma_n = M*(-M1*chi1_n + M2*chi2_n)
            Sigma_lambda = M*(-M1*chi1_lambda + M2*chi2_lambda)
            chi_s_ell = chi1_ell/2 + chi2_ell/2
            chi_a_ell = chi1_ell/2 - chi2_ell/2
            logv = log(v)
            Fcal_coeff = 32*nu**2*v**10/5
            Fcal_SQ_4 = chi1chi1*(-89*delta/192 + 89*nu/96 - 89/192) - 103*chi1chi2*nu/48 + chi2chi2*(89*delta/192 + 89*nu/96 - 89/192) + chi_a_ell*(chi_a_ell*(287/96 - 12*nu) + 287*chi_s_ell*delta/48) + chi_s_ell**2*(nu/24 + 287/96)
            Fcal_SO_3 = (-4*S_ell - 5*Sigma_ell*delta/4)/M**2
            Fcal_SO_5 = (S_ell*(272*nu/9 - 9/2) + Sigma_ell*delta*(43*nu/4 - 13/16))/M**2
            Fcal_SO_6 = (-16*S_ell*pi - 31*Sigma_ell*delta*pi/6)/M**2
            Fcal_SO_7 = (S_ell*(-2810*nu**2/27 + 6172*nu/189 + 476645/6804) + Sigma_ell*delta*(-1501*nu**2/36 + 1849*nu/126 + 9535/336))/M**2
            E_SQ_4 = -3*chi_a_ell**2/2 - 3*chi_s_ell**2/2 - delta*(chi2chi2/2 + 3*chi_a_ell*chi_s_ell) + nu*(chi1chi2 + 6*chi_a_ell**2) + (chi1chi1 + chi2chi2)*(delta - 2*nu + 1)/4
            E_SO_3 = (14*S_ell/3 + 2*Sigma_ell*delta)/M**2
            E_SO_5 = (S_ell*(11 - 61*nu/9) + Sigma_ell*delta*(3 - 10*nu/3))/M**2
            E_SO_7 = (S_ell*(29*nu**2/12 - 367*nu/4 + 135/4) + Sigma_ell*delta*(5*nu**2/4 - 39*nu + 27/4))/M**2

        def OmegaVec_chiVec_1():
            Omega1_coeff = v**5/M
            return Omega1_coeff*(ellHat*(-0.75*delta + 0.5*nu + v**2*(delta*(0.625*nu - 0.5625) + nu*(1.25 - 0.0416666666666667*nu) + v**2*(delta*(nu*(4.875 - 0.15625*nu) - 0.84375) + nu*(nu*(-0.0208333333333333*nu - 3.28125) + 0.1875) + 0.84375) + 0.5625) + 0.75) + nHat*v*(3.0*chi1_n*nu + 3.0*M2**2*chi2_n/M**2) - M2**2*chiVec2*v/M**2)

        def OmegaVec_chiVec_2():
            Omega2_coeff = v**5/M
            return Omega2_coeff*(ellHat*(0.75*delta + 0.5*nu + v**2*(delta*(0.5625 - 0.625*nu) + nu*(1.25 - 0.0416666666666667*nu) + v**2*(delta*(nu*(0.15625*nu - 4.875) + 0.84375) + nu*(nu*(-0.0208333333333333*nu - 3.28125) + 0.1875) + 0.84375) + 0.5625) + 0.75) + nHat*v*(3.0*chi2_n*nu + 3.0*M1**2*chi1_n/M**2) - M1**2*chiVec1*v/M**2)

        def OmegaVec():
            gamma_PN_7 = (S_ell*(-6.0*nu**2 - 10.5833333333333*nu + 5.0) - 2.66666666666667*Sigma_ell*delta*nu**2 + Sigma_ell*delta*(3.0 - 10.1666666666667*nu))/M**2
            gamma_PN_6 = 0.0123456790123457*nu**3 + 6.36111111111111*nu**2 - 2.98177812235564*nu + 1.0
            gamma_PN_2 = 1.0 - 0.333333333333333*nu
            a_ell_2 = S_n*(-9.66666666666667*nu - 10.0) + Sigma_n*delta*(-4.5*nu - 6.0)
            a_ell_0 = 7.0*S_n + 3.0*Sigma_n*delta
            a_ell_4 = S_n*(5.77777777777778*nu**2 + 14.75*nu + 1.5) + Sigma_n*delta*(2.83333333333333*nu**2 + 9.125*nu + 1.5)
            gamma_PN_4 = 1.0 - 5.41666666666667*nu
            gamma_PN_3 = (1.66666666666667*S_ell + Sigma_ell*delta)/M**2
            gamma_PN_5 = (S_ell*(0.888888888888889*nu + 3.33333333333333) + 2.0*Sigma_ell*delta)/M**2
            gamma_PN_0 = 1.00000000000000
            return ellHat*v**3/M + nHat*v**6*(a_ell_0 + v**2*(a_ell_2 + a_ell_4*v**2))*(gamma_PN_0 + v**2*(gamma_PN_2 + v*(gamma_PN_3 + v*(gamma_PN_4 + v*(gamma_PN_5 + v*(gamma_PN_6 + gamma_PN_7*v))))))/M**3

        def OrbitalAngularMomentum():
            L_6 = ellHat*(0.00540123456790123*nu**3 + 1.29166666666667*nu**2 - 30.9797035925835*nu + 8.4375)
            L_4 = ellHat*(0.0416666666666667*nu**2 - 2.375*nu + 3.375)
            L_SO_5 = 0.0138888888888889*ellHat*(427.0*S_ell*nu - 693.0*S_ell + 210.0*Sigma_ell*delta*nu - 189.0*Sigma_ell*delta)/M**2 + 0.166666666666667*lambdaHat*(18.0*S_lambda*nu - 21.0*S_lambda + 8.0*Sigma_lambda*delta*nu - 3.0*Sigma_lambda*delta)/M**2 + 0.0416666666666667*nHat*(-19.0*S_n*nu + 33.0*S_n - 10.0*Sigma_n*delta*nu + 33.0*Sigma_n*delta)/M**2
            L_2 = ellHat*(0.166666666666667*nu + 1.5)
            L_0 = ellHat
            L_SO_3 = 0.166666666666667*ellHat*(-35.0*S_ell - 15.0*Sigma_ell*delta)/M**2 + lambdaHat*(-3.0*S_lambda - Sigma_lambda*delta)/M**2 + 0.5*nHat*(S_n + Sigma_n*delta)/M**2
            L_coeff = M**2*nu/v
            L_SO_7 = 0.0625*ellHat*(-29.0*S_ell*nu**2 + 1101.0*S_ell*nu - 405.0*S_ell - 15.0*Sigma_ell*delta*nu**2 + 468.0*Sigma_ell*delta*nu - 81.0*Sigma_ell*delta)/M**2 + 0.0416666666666667*lambdaHat*(-32.0*S_lambda*nu**2 + 2.0*S_lambda*nu - 174.0*S_lambda - 16.0*Sigma_lambda*delta*nu**2 - 79.0*Sigma_lambda*delta*nu - 12.0*Sigma_lambda*delta)/M**2 + 0.0208333333333333*nHat*(11.0*S_n*nu**2 - 1331.0*S_n*nu + 183.0*S_n + 5.0*Sigma_n*delta*nu**2 - 734.0*Sigma_n*delta*nu + 183.0*Sigma_n*delta)/M**2
            return L_coeff*(L_0 + v**2*(L_2 + v*(L_SO_3 + v*(L_4 + v*(L_SO_5 + v*(L_6 + L_SO_7*v))))))


        #Evolve PN
        def TaylorT1(t, y):
            Recalculate(t, y)
            if v>=1.0: 
                return GSL_EDOM # Beyond domain of PN validity
            Flux = Fcal_coeff*(Fcal_0 + v**2*(Fcal_2 + v*(Fcal_3 + Fcal_SO_3 + v*(Fcal_4 + Fcal_SQ_4 + v*(Fcal_5 + Fcal_SO_5 + v*(Fcal_6 + Fcal_SO_6 + Fcal_lnv_6*logv + v*(Fcal_7 + Fcal_SO_7)))))))
            dEdv = -0.5*M*nu*v*(2.0*E_0 + v**2*(4.0*E_2 + v*(5.0*E_SO_3 + v*(6.0*E_4 + 6.0*E_SQ_4 + v*(7.0*E_SO_5 + v*(8.0*E_6 + 9.0*E_SO_7*v))))))
            Absorption = 0
            dvdt_T1 = (-Absorption - Flux)/dEdv
            if dvdt_T1<1.0e-12:
                return GSL_EDIVERGE # v is decreasing
            return CommonRHS(dvdt_T1, y)

        def CommonRHS(dvdt, y):
            dydt=np.zeros(9)
            rfrak_frame=np.zeros(3)
            rfrak_frame[0] = y[5]
            rfrak_frame[1] = y[6]
            rfrak_frame[2] = y[7]
            rfrakdot_frame = Quaternions.FrameFromAngularVelocity_Integrand(rfrak_frame, OmegaVec().vec())
            dydt[0] = dvdt
            if(EvolveSpin1):
                Quaternions.FrameFromAngularVelocity_2D_Integrand(y[1], y[2],\
                    (S_chi1.inverse()*OmegaVec_chiVec_1()*S_chi1).vec(),dydt[1], dydt[2])
            else:
                dydt[1] = 0.0
                dydt[2] = 0.0
            if(EvolveSpin2):
                Quaternions.FrameFromAngularVelocity_2D_Integrand(y[3], y[4],\
                    (S_chi2.inverse()*OmegaVec_chiVec_2()*S_chi2).vec(),dydt[3], dydt[4])
            else:
                dydt[3] = 0.0
                dydt[4] = 0.0
            dydt[5] = rfrakdot_frame[0]
            dydt[6] = rfrakdot_frame[1]
            dydt[7] = rfrakdot_frame[2]
            dydt[8] = v*v*v/M

            return dydt
        
        
        y=solve_ivp(TaylorT1, [0,20000], [0.01**(1/3),rfrak_chi1_x,\
            rfrak_chi1_y,rfrak_chi2_x,rfrak_chi2_y,rfrak_frame_x,\
            rfrak_frame_y,rfrak_frame_z,Phi], method='DOP853',\
            t_eval=np.arange(0,20000,100), dense_output=True)
            

#        def TaylorT1_phi(t,phi):
#            dphidt_T1=v_of_t.y[0][v_of_t.t>=t][0]**3.0/M
#            return dphidt_T1
#        phi_of_t=solve_ivp(TaylorT1_phi, [0, v_of_t.t[-1]], [0.0], method='RK45', t_eval=v_of_t.t)
        
        #Calculate waveform modes

        def WaveformModes():
            I=1j
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
            chiVec1 = Quaternions.Quaternion(0., chi1_n_k, chi1_lambda_k, chi1_ell_k);
            chiVec2 = Quaternions.Quaternion(0., chi2_n_k, chi2_lambda_k, chi2_ell_k);

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

            Modes=np.zeros(77,len(y.t))
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

class TaylorTn_4p0PN_Q : 
    def TaylorTn_4p0PN_Q(xHat_i, yHat_i, zHat_i, M1_i, M2_i, v_i, S_chi1_i, S_chi2_i, rfrak_chi1_x_i, rfrak_chi1_y_i, rfrak_chi2_x_i,
           rfrak_chi2_y_i, rfrak_frame_x_i, rfrak_frame_y_i, rfrak_frame_z_i) :
        xHat=xHat_i
        yHat=yHat_i
        zHat=zHat_i
        M1=M1_i
        M2=M2_i
        v=v_i
        S_chi1=S_chi1_i
        S_chi2=S_chi2_i
        rfrak_chi1_x=rfrak_chi1_x_i
        rfrak_chi1_y=rfrak_chi1_y_i
        rfrak_chi2_x=rfrak_chi2_x_i
        rfrak_chi2_y=rfrak_chi2_y_i
        rfrak_frame_x=rfrak_frame_x_i
        rfrak_frame_y=rfrak_frame_y_i
        rfrak_frame_z=rfrak_frame_z_i
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
        chi1chi1=dot(chiVec1,chiVec1)
        chi1chi2=dot(chiVec1,chiVec2)
        chi2chi2=dot(chiVec2,chiVec2)
        chi1_n=dot(chiVec1,nHat)
        chi1_lambda=dot(chiVec1,lambdaHat)
        chi1_ell=dot(chiVec1,ellHat)
        chi2_n=dot(chiVec2,nHat)
        chi2_lambda=dot(chiVec2,lambdaHat)
        chi2_ell=dot(chiVec2,ellHat)
        S_ell=M1**2*chi1_ell + M2**2*chi2_ell
        S_n=M1**2*chi1_n + M2**2*chi2_n
        S_lambda=M1**2*chi1_lambda + M2**2*chi2_lambda
        Sigma_ell=M*(-M1*chi1_ell + M2*chi2_ell)
        Sigma_n=M*(-M1*chi1_n + M2*chi2_n)
        Sigma_lambda=M*(-M1*chi1_lambda + M2*chi2_lambda)
        chi_s_ell=chi1_ell/2 + chi2_ell/2
        chi_a_ell=chi1_ell/2 - chi2_ell/2
        logv=log(v)
        Fcal_coeff=32*nu**2*v**10/5
        Fcal_0=1
        Fcal_2=-35*nu/12 - 1247/336
        Fcal_3=4*pi
        Fcal_4=65*nu**2/18 + 9271*nu/504 - 44711/9072
        Fcal_5=pi*(-583*nu/24 - 8191/672)
        Fcal_6=-775*nu**3/324 - 94403*nu**2/3024 + nu*(-134543/7776 + 41*pi**2/48) - 1712*log(4)/105 - 1712*EulerGamma/105 + 16*pi**2/3 + 6643739519/69854400
        Fcal_lnv_6=-1712/105
        Fcal_7=pi*(193385*nu**2/3024 + 214745*nu/1728 - 16285/504)
        Fcal_8=-1369*pi**2/126 - 323105549467/3178375200 - 47385*log(3)/1568 + 232597*EulerGamma/4410 + 39931*log(2)/294
        Fcal_lnv_8=232597/4410
        Fcal_SQ_4=chi1chi1*(-89*delta/192 + 89*nu/96 - 89/192) - 103*chi1chi2*nu/48 + chi2chi2*(89*delta/192 + 89*nu/96 - 89/192) + chi_a_ell*(chi_a_ell*(287/96 - 12*nu) + 287*chi_s_ell*delta/48) + chi_s_ell**2*(nu/24 + 287/96)
        Fcal_SO_3=(-4*S_ell - 5*Sigma_ell*delta/4)/M**2
        Fcal_SO_5=(S_ell*(272*nu/9 - 9/2) + Sigma_ell*delta*(43*nu/4 - 13/16))/M**2
        Fcal_SO_6=(-16*S_ell*pi - 31*Sigma_ell*delta*pi/6)/M**2
        Fcal_SO_7=(S_ell*(-2810*nu**2/27 + 6172*nu/189 + 476645/6804) + Sigma_ell*delta*(-1501*nu**2/36 + 1849*nu/126 + 9535/336))/M**2
        Fcal_SO_8=(S_ell*pi*(13879*nu/72 - 3485/96) + Sigma_ell*delta*pi*(130583*nu/2016 - 7163/672))/M**2
        E_0=1
        E_2=-nu/12 - 3/4
        E_4=-nu**2/24 + 19*nu/8 - 27/8
        E_6=-35*nu**3/5184 - 155*nu**2/96 + nu*(34445/576 - 205*pi**2/96) - 675/64
        E_8=77*nu**4/31104 + 301*nu**3/1728 + nu**2*(-498449/3456 + 3157*pi**2/576) + nu*(-123671/5760 + 896*EulerGamma/15 + 9037*pi**2/1536 + 1792*log(2)/15) - 3969/128
        E_lnv_8=896*nu/15
        E_SQ_4=-3*chi_a_ell**2/2 - 3*chi_s_ell**2/2 - delta*(chi2chi2/2 + 3*chi_a_ell*chi_s_ell) + nu*(chi1chi2 + 6*chi_a_ell**2) + (chi1chi1 + chi2chi2)*(delta - 2*nu + 1)/4
        E_SO_3=(14*S_ell/3 + 2*Sigma_ell*delta)/M**2
        E_SO_5=(S_ell*(11 - 61*nu/9) + Sigma_ell*delta*(3 - 10*nu/3))/M**2
        E_SO_7=(S_ell*(29*nu**2/12 - 367*nu/4 + 135/4) + Sigma_ell*delta*(5*nu**2/4 - 39*nu + 27/4))/M**2
        Phi=0.0
        EvolveSpin1=dot(S_chi1,S_chi1)>1e-12
        EvolveSpin2=dot(S_chi2,S_chi2)>1e-12

        def Recalculate(t, y):
            v = y[0];
            rfrak_chi1_x = y[1];
            rfrak_chi1_y = y[2];
            rfrak_chi2_x = y[3];
            rfrak_chi2_y = y[4];
            rfrak_frame_x = y[5];
            rfrak_frame_y = y[6];
            rfrak_frame_z = y[7];
            Phi = y[8];
            R = exp(rfrak_frame_x*xHat + rfrak_frame_y*yHat + rfrak_frame_z*zHat)
            nHat = R*xHat*conjugate(R)
            lambdaHat = R*yHat*conjugate(R)
            ellHat = R*zHat*conjugate(R)
            R_S1 = exp(rfrak_chi1_x*xHat + rfrak_chi1_y*yHat)
            R_S2 = exp(rfrak_chi2_x*xHat + rfrak_chi2_y*yHat)
            chiVec1 = S_chi1*R_S1*zHat*conjugate(R_S1)*conjugate(S_chi1)
            chiVec2 = S_chi2*R_S2*zHat*conjugate(R_S2)*conjugate(S_chi2)
            chi1_n = dot(chiVec1,nHat)
            chi1_lambda = dot(chiVec1,lambdaHat)
            chi1_ell = dot(chiVec1,ellHat)
            chi2_n = dot(chiVec2,nHat)
            chi2_lambda = dot(chiVec2,lambdaHat)
            chi2_ell = dot(chiVec2,ellHat)
            S_ell = M1**2*chi1_ell + M2**2*chi2_ell
            S_n = M1**2*chi1_n + M2**2*chi2_n
            S_lambda = M1**2*chi1_lambda + M2**2*chi2_lambda
            Sigma_ell = M*(-M1*chi1_ell + M2*chi2_ell)
            Sigma_n = M*(-M1*chi1_n + M2*chi2_n)
            Sigma_lambda = M*(-M1*chi1_lambda + M2*chi2_lambda)
            chi_s_ell = chi1_ell/2 + chi2_ell/2
            chi_a_ell = chi1_ell/2 - chi2_ell/2
            logv = log(v)
            Fcal_coeff = 32*nu**2*v**10/5
            Fcal_SQ_4 = chi1chi1*(-89*delta/192 + 89*nu/96 - 89/192) - 103*chi1chi2*nu/48 + chi2chi2*(89*delta/192 + 89*nu/96 - 89/192) + chi_a_ell*(chi_a_ell*(287/96 - 12*nu) + 287*chi_s_ell*delta/48) + chi_s_ell**2*(nu/24 + 287/96)
            Fcal_SO_3 = (-4*S_ell - 5*Sigma_ell*delta/4)/M**2
            Fcal_SO_5 = (S_ell*(272*nu/9 - 9/2) + Sigma_ell*delta*(43*nu/4 - 13/16))/M**2
            Fcal_SO_6 = (-16*S_ell*pi - 31*Sigma_ell*delta*pi/6)/M**2
            Fcal_SO_7 = (S_ell*(-2810*nu**2/27 + 6172*nu/189 + 476645/6804) + Sigma_ell*delta*(-1501*nu**2/36 + 1849*nu/126 + 9535/336))/M**2
            Fcal_SO_8 = (S_ell*pi*(13879*nu/72 - 3485/96) + Sigma_ell*delta*pi*(130583*nu/2016 - 7163/672))/M**2
            E_SQ_4 = -3*chi_a_ell**2/2 - 3*chi_s_ell**2/2 - delta*(chi2chi2/2 + 3*chi_a_ell*chi_s_ell) + nu*(chi1chi2 + 6*chi_a_ell**2) + (chi1chi1 + chi2chi2)*(delta - 2*nu + 1)/4
            E_SO_3 = (14*S_ell/3 + 2*Sigma_ell*delta)/M**2
            E_SO_5 = (S_ell*(11 - 61*nu/9) + Sigma_ell*delta*(3 - 10*nu/3))/M**2
            E_SO_7 = (S_ell*(29*nu**2/12 - 367*nu/4 + 135/4) + Sigma_ell*delta*(5*nu**2/4 - 39*nu + 27/4))/M**2

        def OmegaVec_chiVec_1():
            Omega1_coeff = v**5/M
            return Omega1_coeff*(ellHat*(-0.75*delta + 0.5*nu + v**2*(delta*(0.625*nu - 0.5625) + nu*(1.25 - 0.0416666666666667*nu) + v**2*(delta*(nu*(4.875 - 0.15625*nu) - 0.84375) + nu*(nu*(-0.0208333333333333*nu - 3.28125) + 0.1875) + 0.84375) + 0.5625) + 0.75) + nHat*v*(3.0*chi1_n*nu + 3.0*M2**2*chi2_n/M**2) - M2**2*chiVec2*v/M**2)

        def OmegaVec_chiVec_2():
            Omega2_coeff = v**5/M
            return Omega2_coeff*(ellHat*(0.75*delta + 0.5*nu + v**2*(delta*(0.5625 - 0.625*nu) + nu*(1.25 - 0.0416666666666667*nu) + v**2*(delta*(nu*(0.15625*nu - 4.875) + 0.84375) + nu*(nu*(-0.0208333333333333*nu - 3.28125) + 0.1875) + 0.84375) + 0.5625) + 0.75) + nHat*v*(3.0*chi2_n*nu + 3.0*M1**2*chi1_n/M**2) - M1**2*chiVec1*v/M**2)

        def OmegaVec():
            gamma_PN_7 = (S_ell*(-6.0*nu**2 - 10.5833333333333*nu + 5.0) - 2.66666666666667*Sigma_ell*delta*nu**2 + Sigma_ell*delta*(3.0 - 10.1666666666667*nu))/M**2
            gamma_PN_6 = 0.0123456790123457*nu**3 + 6.36111111111111*nu**2 - 2.98177812235564*nu + 1.0
            gamma_PN_2 = 1.0 - 0.333333333333333*nu
            a_ell_2 = S_n*(-9.66666666666667*nu - 10.0) + Sigma_n*delta*(-4.5*nu - 6.0)
            a_ell_0 = 7.0*S_n + 3.0*Sigma_n*delta
            a_ell_4 = S_n*(5.77777777777778*nu**2 + 14.75*nu + 1.5) + Sigma_n*delta*(2.83333333333333*nu**2 + 9.125*nu + 1.5)
            gamma_PN_4 = 1.0 - 5.41666666666667*nu
            gamma_PN_3 = (1.66666666666667*S_ell + Sigma_ell*delta)/M**2
            gamma_PN_5 = (S_ell*(0.888888888888889*nu + 3.33333333333333) + 2.0*Sigma_ell*delta)/M**2
            gamma_PN_0 = 1.00000000000000
            return ellHat*v**3/M + nHat*v**6*(a_ell_0 + v**2*(a_ell_2 + a_ell_4*v**2))*(gamma_PN_0 + v**2*(gamma_PN_2 + v*(gamma_PN_3 + v*(gamma_PN_4 + v*(gamma_PN_5 + v*(gamma_PN_6 + gamma_PN_7*v))))))/M**3

        def OrbitalAngularMomentum():
            L_6 = ellHat*(0.00540123456790123*nu**3 + 1.29166666666667*nu**2 - 30.9797035925835*nu + 8.4375)
            L_4 = ellHat*(0.0416666666666667*nu**2 - 2.375*nu + 3.375)
            L_SO_5 = 0.0138888888888889*ellHat*(427.0*S_ell*nu - 693.0*S_ell + 210.0*Sigma_ell*delta*nu - 189.0*Sigma_ell*delta)/M**2 + 0.166666666666667*lambdaHat*(18.0*S_lambda*nu - 21.0*S_lambda + 8.0*Sigma_lambda*delta*nu - 3.0*Sigma_lambda*delta)/M**2 + 0.0416666666666667*nHat*(-19.0*S_n*nu + 33.0*S_n - 10.0*Sigma_n*delta*nu + 33.0*Sigma_n*delta)/M**2
            L_2 = ellHat*(0.166666666666667*nu + 1.5)
            L_0 = ellHat
            L_SO_3 = 0.166666666666667*ellHat*(-35.0*S_ell - 15.0*Sigma_ell*delta)/M**2 + lambdaHat*(-3.0*S_lambda - Sigma_lambda*delta)/M**2 + 0.5*nHat*(S_n + Sigma_n*delta)/M**2
            L_lnv_8 = -42.6666666666667*ellHat*nu
            L_SO_7 = 0.0625*ellHat*(-29.0*S_ell*nu**2 + 1101.0*S_ell*nu - 405.0*S_ell - 15.0*Sigma_ell*delta*nu**2 + 468.0*Sigma_ell*delta*nu - 81.0*Sigma_ell*delta)/M**2 + 0.0416666666666667*lambdaHat*(-32.0*S_lambda*nu**2 + 2.0*S_lambda*nu - 174.0*S_lambda - 16.0*Sigma_lambda*delta*nu**2 - 79.0*Sigma_lambda*delta*nu - 12.0*Sigma_lambda*delta)/M**2 + 0.0208333333333333*nHat*(11.0*S_n*nu**2 - 1331.0*S_n*nu + 183.0*S_n + 5.0*Sigma_n*delta*nu**2 - 734.0*Sigma_n*delta*nu + 183.0*Sigma_n*delta)/M**2
            L_coeff = M**2*nu/v
            L_8 = ellHat*(-0.714285714285714*nu*(0.0024755658436214*nu**3 + 0.174189814814815*nu**2 - 90.1327990262052*nu + 153.88379682994) + 1.82857142857143*nu + 22.1484375)
            return L_coeff*(L_0 + v**2*(L_2 + v*(L_SO_3 + v*(L_4 + v*(L_SO_5 + v*(L_6 + v*(L_SO_7 + v*(L_8 + L_lnv_8*math.log(v)))))))))


        #Evolve PN
        def TaylorT1(t, y):
            Recalculate(t, y)
            if v>=1.0: 
                return GSL_EDOM # Beyond domain of PN validity
            Flux = Fcal_coeff*(Fcal_0 + v**2*(Fcal_2 + v*(Fcal_3 + Fcal_SO_3 + v*(Fcal_4 + Fcal_SQ_4 + v*(Fcal_5 + Fcal_SO_5 + v*(Fcal_6 + Fcal_SO_6 + Fcal_lnv_6*logv + v*(Fcal_7 + Fcal_SO_7 + v*(Fcal_8 + Fcal_SO_8 + Fcal_lnv_8*logv))))))))
            dEdv = -0.5*M*nu*v*(2.0*E_0 + v**2*(4.0*E_2 + v*(5.0*E_SO_3 + v*(6.0*E_4 + 6.0*E_SQ_4 + v*(7.0*E_SO_5 + v*(8.0*E_6 + v*(9.0*E_SO_7 + v*(10.0*E_8 + E_lnv_8*(10.0*logv + 1.0)))))))))
            Absorption = 0
            dvdt_T1 = (-Absorption - Flux)/dEdv
            if dvdt_T1<1.0e-12:
                return GSL_EDIVERGE # v is decreasing
            return CommonRHS(dvdt_T1, y)

        def CommonRHS(dvdt, y):
            dydt=np.zeros(9)
            rfrak_frame=np.zeros(3)
            rfrak_frame[0] = y[5]
            rfrak_frame[1] = y[6]
            rfrak_frame[2] = y[7]
            rfrakdot_frame = Quaternions.FrameFromAngularVelocity_Integrand(rfrak_frame, OmegaVec().vec())
            dydt[0] = dvdt
            if(EvolveSpin1):
                Quaternions.FrameFromAngularVelocity_2D_Integrand(y[1], y[2],\
                    (S_chi1.inverse()*OmegaVec_chiVec_1()*S_chi1).vec(),dydt[1], dydt[2])
            else:
                dydt[1] = 0.0
                dydt[2] = 0.0
            if(EvolveSpin2):
                Quaternions.FrameFromAngularVelocity_2D_Integrand(y[3], y[4],\
                    (S_chi2.inverse()*OmegaVec_chiVec_2()*S_chi2).vec(),dydt[3], dydt[4])
            else:
                dydt[3] = 0.0
                dydt[4] = 0.0
            dydt[5] = rfrakdot_frame[0]
            dydt[6] = rfrakdot_frame[1]
            dydt[7] = rfrakdot_frame[2]
            dydt[8] = v*v*v/M

            return dydt
        
        
        y=solve_ivp(TaylorT1, [0,20000], [0.01**(1/3),rfrak_chi1_x,\
            rfrak_chi1_y,rfrak_chi2_x,rfrak_chi2_y,rfrak_frame_x,\
            rfrak_frame_y,rfrak_frame_z,Phi], method='DOP853',\
            t_eval=np.arange(0,20000,100), dense_output=True)
            

#        def TaylorT1_phi(t,phi):
#            dphidt_T1=v_of_t.y[0][v_of_t.t>=t][0]**3.0/M
#            return dphidt_T1
#        phi_of_t=solve_ivp(TaylorT1_phi, [0, v_of_t.t[-1]], [0.0], method='RK45', t_eval=v_of_t.t)
        
        #Calculate waveform modes

        def WaveformModes():
            I=1j
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
            chiVec1 = Quaternions.Quaternion(0., chi1_n_k, chi1_lambda_k, chi1_ell_k);
            chiVec2 = Quaternions.Quaternion(0., chi2_n_k, chi2_lambda_k, chi2_ell_k);

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

            Modes=np.zeros(77,len(y.t))
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

class TaylorTn_4p5PN_Q : 
    def TaylorTn_4p5PN_Q(xHat_i, yHat_i, zHat_i, M1_i, M2_i, v_i, S_chi1_i, S_chi2_i, rfrak_chi1_x_i, rfrak_chi1_y_i, rfrak_chi2_x_i,
           rfrak_chi2_y_i, rfrak_frame_x_i, rfrak_frame_y_i, rfrak_frame_z_i) :
        xHat=xHat_i
        yHat=yHat_i
        zHat=zHat_i
        M1=M1_i
        M2=M2_i
        v=v_i
        S_chi1=S_chi1_i
        S_chi2=S_chi2_i
        rfrak_chi1_x=rfrak_chi1_x_i
        rfrak_chi1_y=rfrak_chi1_y_i
        rfrak_chi2_x=rfrak_chi2_x_i
        rfrak_chi2_y=rfrak_chi2_y_i
        rfrak_frame_x=rfrak_frame_x_i
        rfrak_frame_y=rfrak_frame_y_i
        rfrak_frame_z=rfrak_frame_z_i
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
        chi1chi1=dot(chiVec1,chiVec1)
        chi1chi2=dot(chiVec1,chiVec2)
        chi2chi2=dot(chiVec2,chiVec2)
        chi1_n=dot(chiVec1,nHat)
        chi1_lambda=dot(chiVec1,lambdaHat)
        chi1_ell=dot(chiVec1,ellHat)
        chi2_n=dot(chiVec2,nHat)
        chi2_lambda=dot(chiVec2,lambdaHat)
        chi2_ell=dot(chiVec2,ellHat)
        S_ell=M1**2*chi1_ell + M2**2*chi2_ell
        S_n=M1**2*chi1_n + M2**2*chi2_n
        S_lambda=M1**2*chi1_lambda + M2**2*chi2_lambda
        Sigma_ell=M*(-M1*chi1_ell + M2*chi2_ell)
        Sigma_n=M*(-M1*chi1_n + M2*chi2_n)
        Sigma_lambda=M*(-M1*chi1_lambda + M2*chi2_lambda)
        chi_s_ell=chi1_ell/2 + chi2_ell/2
        chi_a_ell=chi1_ell/2 - chi2_ell/2
        logv=log(v)
        Fcal_coeff=32*nu**2*v**10/5
        Fcal_0=1
        Fcal_2=-35*nu/12 - 1247/336
        Fcal_3=4*pi
        Fcal_4=65*nu**2/18 + 9271*nu/504 - 44711/9072
        Fcal_5=pi*(-583*nu/24 - 8191/672)
        Fcal_6=-775*nu**3/324 - 94403*nu**2/3024 + nu*(-134543/7776 + 41*pi**2/48) - 1712*log(4)/105 - 1712*EulerGamma/105 + 16*pi**2/3 + 6643739519/69854400
        Fcal_lnv_6=-1712/105
        Fcal_7=pi*(193385*nu**2/3024 + 214745*nu/1728 - 16285/504)
        Fcal_8=-1369*pi**2/126 - 323105549467/3178375200 - 47385*log(3)/1568 + 232597*EulerGamma/4410 + 39931*log(2)/294
        Fcal_lnv_8=232597/4410
        Fcal_9=-13696*pi*log(2)/105 - 6848*EulerGamma*pi/105 + 265978667519*pi/745113600
        Fcal_lnv_9=-6848*pi/105
        Fcal_SQ_4=chi1chi1*(-89*delta/192 + 89*nu/96 - 89/192) - 103*chi1chi2*nu/48 + chi2chi2*(89*delta/192 + 89*nu/96 - 89/192) + chi_a_ell*(chi_a_ell*(287/96 - 12*nu) + 287*chi_s_ell*delta/48) + chi_s_ell**2*(nu/24 + 287/96)
        Fcal_SO_3=(-4*S_ell - 5*Sigma_ell*delta/4)/M**2
        Fcal_SO_5=(S_ell*(272*nu/9 - 9/2) + Sigma_ell*delta*(43*nu/4 - 13/16))/M**2
        Fcal_SO_6=(-16*S_ell*pi - 31*Sigma_ell*delta*pi/6)/M**2
        Fcal_SO_7=(S_ell*(-2810*nu**2/27 + 6172*nu/189 + 476645/6804) + Sigma_ell*delta*(-1501*nu**2/36 + 1849*nu/126 + 9535/336))/M**2
        Fcal_SO_8=(S_ell*pi*(13879*nu/72 - 3485/96) + Sigma_ell*delta*pi*(130583*nu/2016 - 7163/672))/M**2
        E_0=1
        E_2=-nu/12 - 3/4
        E_4=-nu**2/24 + 19*nu/8 - 27/8
        E_6=-35*nu**3/5184 - 155*nu**2/96 + nu*(34445/576 - 205*pi**2/96) - 675/64
        E_8=77*nu**4/31104 + 301*nu**3/1728 + nu**2*(-498449/3456 + 3157*pi**2/576) + nu*(-123671/5760 + 896*EulerGamma/15 + 9037*pi**2/1536 + 1792*log(2)/15) - 3969/128
        E_lnv_8=896*nu/15
        E_SQ_4=-3*chi_a_ell**2/2 - 3*chi_s_ell**2/2 - delta*(chi2chi2/2 + 3*chi_a_ell*chi_s_ell) + nu*(chi1chi2 + 6*chi_a_ell**2) + (chi1chi1 + chi2chi2)*(delta - 2*nu + 1)/4
        E_SO_3=(14*S_ell/3 + 2*Sigma_ell*delta)/M**2
        E_SO_5=(S_ell*(11 - 61*nu/9) + Sigma_ell*delta*(3 - 10*nu/3))/M**2
        E_SO_7=(S_ell*(29*nu**2/12 - 367*nu/4 + 135/4) + Sigma_ell*delta*(5*nu**2/4 - 39*nu + 27/4))/M**2
        Phi=0.0
        EvolveSpin1=dot(S_chi1,S_chi1)>1e-12
        EvolveSpin2=dot(S_chi2,S_chi2)>1e-12

        def Recalculate(t, y):
            v = y[0];
            rfrak_chi1_x = y[1];
            rfrak_chi1_y = y[2];
            rfrak_chi2_x = y[3];
            rfrak_chi2_y = y[4];
            rfrak_frame_x = y[5];
            rfrak_frame_y = y[6];
            rfrak_frame_z = y[7];
            Phi = y[8];
            R = exp(rfrak_frame_x*xHat + rfrak_frame_y*yHat + rfrak_frame_z*zHat)
            nHat = R*xHat*conjugate(R)
            lambdaHat = R*yHat*conjugate(R)
            ellHat = R*zHat*conjugate(R)
            R_S1 = exp(rfrak_chi1_x*xHat + rfrak_chi1_y*yHat)
            R_S2 = exp(rfrak_chi2_x*xHat + rfrak_chi2_y*yHat)
            chiVec1 = S_chi1*R_S1*zHat*conjugate(R_S1)*conjugate(S_chi1)
            chiVec2 = S_chi2*R_S2*zHat*conjugate(R_S2)*conjugate(S_chi2)
            chi1_n = dot(chiVec1,nHat)
            chi1_lambda = dot(chiVec1,lambdaHat)
            chi1_ell = dot(chiVec1,ellHat)
            chi2_n = dot(chiVec2,nHat)
            chi2_lambda = dot(chiVec2,lambdaHat)
            chi2_ell = dot(chiVec2,ellHat)
            S_ell = M1**2*chi1_ell + M2**2*chi2_ell
            S_n = M1**2*chi1_n + M2**2*chi2_n
            S_lambda = M1**2*chi1_lambda + M2**2*chi2_lambda
            Sigma_ell = M*(-M1*chi1_ell + M2*chi2_ell)
            Sigma_n = M*(-M1*chi1_n + M2*chi2_n)
            Sigma_lambda = M*(-M1*chi1_lambda + M2*chi2_lambda)
            chi_s_ell = chi1_ell/2 + chi2_ell/2
            chi_a_ell = chi1_ell/2 - chi2_ell/2
            logv = log(v)
            Fcal_coeff = 32*nu**2*v**10/5
            Fcal_SQ_4 = chi1chi1*(-89*delta/192 + 89*nu/96 - 89/192) - 103*chi1chi2*nu/48 + chi2chi2*(89*delta/192 + 89*nu/96 - 89/192) + chi_a_ell*(chi_a_ell*(287/96 - 12*nu) + 287*chi_s_ell*delta/48) + chi_s_ell**2*(nu/24 + 287/96)
            Fcal_SO_3 = (-4*S_ell - 5*Sigma_ell*delta/4)/M**2
            Fcal_SO_5 = (S_ell*(272*nu/9 - 9/2) + Sigma_ell*delta*(43*nu/4 - 13/16))/M**2
            Fcal_SO_6 = (-16*S_ell*pi - 31*Sigma_ell*delta*pi/6)/M**2
            Fcal_SO_7 = (S_ell*(-2810*nu**2/27 + 6172*nu/189 + 476645/6804) + Sigma_ell*delta*(-1501*nu**2/36 + 1849*nu/126 + 9535/336))/M**2
            Fcal_SO_8 = (S_ell*pi*(13879*nu/72 - 3485/96) + Sigma_ell*delta*pi*(130583*nu/2016 - 7163/672))/M**2
            E_SQ_4 = -3*chi_a_ell**2/2 - 3*chi_s_ell**2/2 - delta*(chi2chi2/2 + 3*chi_a_ell*chi_s_ell) + nu*(chi1chi2 + 6*chi_a_ell**2) + (chi1chi1 + chi2chi2)*(delta - 2*nu + 1)/4
            E_SO_3 = (14*S_ell/3 + 2*Sigma_ell*delta)/M**2
            E_SO_5 = (S_ell*(11 - 61*nu/9) + Sigma_ell*delta*(3 - 10*nu/3))/M**2
            E_SO_7 = (S_ell*(29*nu**2/12 - 367*nu/4 + 135/4) + Sigma_ell*delta*(5*nu**2/4 - 39*nu + 27/4))/M**2

        def OmegaVec_chiVec_1():
            Omega1_coeff = v**5/M
            return Omega1_coeff*(ellHat*(-0.75*delta + 0.5*nu + v**2*(delta*(0.625*nu - 0.5625) + nu*(1.25 - 0.0416666666666667*nu) + v**2*(delta*(nu*(4.875 - 0.15625*nu) - 0.84375) + nu*(nu*(-0.0208333333333333*nu - 3.28125) + 0.1875) + 0.84375) + 0.5625) + 0.75) + nHat*v*(3.0*chi1_n*nu + 3.0*M2**2*chi2_n/M**2) - M2**2*chiVec2*v/M**2)

        def OmegaVec_chiVec_2():
            Omega2_coeff = v**5/M
            return Omega2_coeff*(ellHat*(0.75*delta + 0.5*nu + v**2*(delta*(0.5625 - 0.625*nu) + nu*(1.25 - 0.0416666666666667*nu) + v**2*(delta*(nu*(0.15625*nu - 4.875) + 0.84375) + nu*(nu*(-0.0208333333333333*nu - 3.28125) + 0.1875) + 0.84375) + 0.5625) + 0.75) + nHat*v*(3.0*chi2_n*nu + 3.0*M1**2*chi1_n/M**2) - M1**2*chiVec1*v/M**2)

        def OmegaVec():
            gamma_PN_7 = (S_ell*(-6.0*nu**2 - 10.5833333333333*nu + 5.0) - 2.66666666666667*Sigma_ell*delta*nu**2 + Sigma_ell*delta*(3.0 - 10.1666666666667*nu))/M**2
            gamma_PN_6 = 0.0123456790123457*nu**3 + 6.36111111111111*nu**2 - 2.98177812235564*nu + 1.0
            gamma_PN_2 = 1.0 - 0.333333333333333*nu
            a_ell_2 = S_n*(-9.66666666666667*nu - 10.0) + Sigma_n*delta*(-4.5*nu - 6.0)
            a_ell_0 = 7.0*S_n + 3.0*Sigma_n*delta
            a_ell_4 = S_n*(5.77777777777778*nu**2 + 14.75*nu + 1.5) + Sigma_n*delta*(2.83333333333333*nu**2 + 9.125*nu + 1.5)
            gamma_PN_4 = 1.0 - 5.41666666666667*nu
            gamma_PN_3 = (1.66666666666667*S_ell + Sigma_ell*delta)/M**2
            gamma_PN_5 = (S_ell*(0.888888888888889*nu + 3.33333333333333) + 2.0*Sigma_ell*delta)/M**2
            gamma_PN_0 = 1.00000000000000
            return ellHat*v**3/M + nHat*v**6*(a_ell_0 + v**2*(a_ell_2 + a_ell_4*v**2))*(gamma_PN_0 + v**2*(gamma_PN_2 + v*(gamma_PN_3 + v*(gamma_PN_4 + v*(gamma_PN_5 + v*(gamma_PN_6 + gamma_PN_7*v))))))/M**3

        def OrbitalAngularMomentum():
            L_6 = ellHat*(0.00540123456790123*nu**3 + 1.29166666666667*nu**2 - 30.9797035925835*nu + 8.4375)
            L_4 = ellHat*(0.0416666666666667*nu**2 - 2.375*nu + 3.375)
            L_SO_5 = 0.0138888888888889*ellHat*(427.0*S_ell*nu - 693.0*S_ell + 210.0*Sigma_ell*delta*nu - 189.0*Sigma_ell*delta)/M**2 + 0.166666666666667*lambdaHat*(18.0*S_lambda*nu - 21.0*S_lambda + 8.0*Sigma_lambda*delta*nu - 3.0*Sigma_lambda*delta)/M**2 + 0.0416666666666667*nHat*(-19.0*S_n*nu + 33.0*S_n - 10.0*Sigma_n*delta*nu + 33.0*Sigma_n*delta)/M**2
            L_2 = ellHat*(0.166666666666667*nu + 1.5)
            L_0 = ellHat
            L_SO_3 = 0.166666666666667*ellHat*(-35.0*S_ell - 15.0*Sigma_ell*delta)/M**2 + lambdaHat*(-3.0*S_lambda - Sigma_lambda*delta)/M**2 + 0.5*nHat*(S_n + Sigma_n*delta)/M**2
            L_lnv_8 = -42.6666666666667*ellHat*nu
            L_SO_7 = 0.0625*ellHat*(-29.0*S_ell*nu**2 + 1101.0*S_ell*nu - 405.0*S_ell - 15.0*Sigma_ell*delta*nu**2 + 468.0*Sigma_ell*delta*nu - 81.0*Sigma_ell*delta)/M**2 + 0.0416666666666667*lambdaHat*(-32.0*S_lambda*nu**2 + 2.0*S_lambda*nu - 174.0*S_lambda - 16.0*Sigma_lambda*delta*nu**2 - 79.0*Sigma_lambda*delta*nu - 12.0*Sigma_lambda*delta)/M**2 + 0.0208333333333333*nHat*(11.0*S_n*nu**2 - 1331.0*S_n*nu + 183.0*S_n + 5.0*Sigma_n*delta*nu**2 - 734.0*Sigma_n*delta*nu + 183.0*Sigma_n*delta)/M**2
            L_coeff = M**2*nu/v
            L_8 = ellHat*(-0.714285714285714*nu*(0.0024755658436214*nu**3 + 0.174189814814815*nu**2 - 90.1327990262052*nu + 153.88379682994) + 1.82857142857143*nu + 22.1484375)
            return L_coeff*(L_0 + v**2*(L_2 + v*(L_SO_3 + v*(L_4 + v*(L_SO_5 + v*(L_6 + v*(L_SO_7 + v*(L_8 + L_lnv_8*math.log(v)))))))))


        #Evolve PN
        def TaylorT1(t, y):
            Recalculate(t, y)
            if v>=1.0: 
                return GSL_EDOM # Beyond domain of PN validity
            Flux = Fcal_coeff*(Fcal_0 + v**2*(Fcal_2 + v*(Fcal_3 + Fcal_SO_3 + v*(Fcal_4 + Fcal_SQ_4 + v*(Fcal_5 + Fcal_SO_5 + v*(Fcal_6 + Fcal_SO_6 + Fcal_lnv_6*logv + v*(Fcal_7 + Fcal_SO_7 + v*(Fcal_8 + Fcal_SO_8 + Fcal_lnv_8*logv + v*(Fcal_9 + Fcal_lnv_9*logv)))))))))
            dEdv = -0.5*M*nu*v*(2.0*E_0 + v**2*(4.0*E_2 + v*(5.0*E_SO_3 + v*(6.0*E_4 + 6.0*E_SQ_4 + v*(7.0*E_SO_5 + v*(8.0*E_6 + v*(9.0*E_SO_7 + v*(10.0*E_8 + E_lnv_8*(10.0*logv + 1.0)))))))))
            Absorption = 0
            dvdt_T1 = (-Absorption - Flux)/dEdv
            if dvdt_T1<1.0e-12:
                return GSL_EDIVERGE # v is decreasing
            return CommonRHS(dvdt_T1, y)

        def CommonRHS(dvdt, y):
            dydt=np.zeros(9)
            rfrak_frame=np.zeros(3)
            rfrak_frame[0] = y[5]
            rfrak_frame[1] = y[6]
            rfrak_frame[2] = y[7]
            rfrakdot_frame = Quaternions.FrameFromAngularVelocity_Integrand(rfrak_frame, OmegaVec().vec())
            dydt[0] = dvdt
            if(EvolveSpin1):
                Quaternions.FrameFromAngularVelocity_2D_Integrand(y[1], y[2],\
                    (S_chi1.inverse()*OmegaVec_chiVec_1()*S_chi1).vec(),dydt[1], dydt[2])
            else:
                dydt[1] = 0.0
                dydt[2] = 0.0
            if(EvolveSpin2):
                Quaternions.FrameFromAngularVelocity_2D_Integrand(y[3], y[4],\
                    (S_chi2.inverse()*OmegaVec_chiVec_2()*S_chi2).vec(),dydt[3], dydt[4])
            else:
                dydt[3] = 0.0
                dydt[4] = 0.0
            dydt[5] = rfrakdot_frame[0]
            dydt[6] = rfrakdot_frame[1]
            dydt[7] = rfrakdot_frame[2]
            dydt[8] = v*v*v/M

            return dydt
        
        
        y=solve_ivp(TaylorT1, [0,20000], [0.01**(1/3),rfrak_chi1_x,\
            rfrak_chi1_y,rfrak_chi2_x,rfrak_chi2_y,rfrak_frame_x,\
            rfrak_frame_y,rfrak_frame_z,Phi], method='DOP853',\
            t_eval=np.arange(0,20000,100), dense_output=True)
            

#        def TaylorT1_phi(t,phi):
#            dphidt_T1=v_of_t.y[0][v_of_t.t>=t][0]**3.0/M
#            return dphidt_T1
#        phi_of_t=solve_ivp(TaylorT1_phi, [0, v_of_t.t[-1]], [0.0], method='RK45', t_eval=v_of_t.t)
        
        #Calculate waveform modes

        def WaveformModes():
            I=1j
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
            chiVec1 = Quaternions.Quaternion(0., chi1_n_k, chi1_lambda_k, chi1_ell_k);
            chiVec2 = Quaternions.Quaternion(0., chi2_n_k, chi2_lambda_k, chi2_ell_k);

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

            Modes=np.zeros(77,len(y.t))
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

class TaylorTn_5p0PN_Q : 
    def TaylorTn_5p0PN_Q(xHat_i, yHat_i, zHat_i, M1_i, M2_i, v_i, S_chi1_i, S_chi2_i, rfrak_chi1_x_i, rfrak_chi1_y_i, rfrak_chi2_x_i,
           rfrak_chi2_y_i, rfrak_frame_x_i, rfrak_frame_y_i, rfrak_frame_z_i) :
        xHat=xHat_i
        yHat=yHat_i
        zHat=zHat_i
        M1=M1_i
        M2=M2_i
        v=v_i
        S_chi1=S_chi1_i
        S_chi2=S_chi2_i
        rfrak_chi1_x=rfrak_chi1_x_i
        rfrak_chi1_y=rfrak_chi1_y_i
        rfrak_chi2_x=rfrak_chi2_x_i
        rfrak_chi2_y=rfrak_chi2_y_i
        rfrak_frame_x=rfrak_frame_x_i
        rfrak_frame_y=rfrak_frame_y_i
        rfrak_frame_z=rfrak_frame_z_i
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
        chi1chi1=dot(chiVec1,chiVec1)
        chi1chi2=dot(chiVec1,chiVec2)
        chi2chi2=dot(chiVec2,chiVec2)
        chi1_n=dot(chiVec1,nHat)
        chi1_lambda=dot(chiVec1,lambdaHat)
        chi1_ell=dot(chiVec1,ellHat)
        chi2_n=dot(chiVec2,nHat)
        chi2_lambda=dot(chiVec2,lambdaHat)
        chi2_ell=dot(chiVec2,ellHat)
        S_ell=M1**2*chi1_ell + M2**2*chi2_ell
        S_n=M1**2*chi1_n + M2**2*chi2_n
        S_lambda=M1**2*chi1_lambda + M2**2*chi2_lambda
        Sigma_ell=M*(-M1*chi1_ell + M2*chi2_ell)
        Sigma_n=M*(-M1*chi1_n + M2*chi2_n)
        Sigma_lambda=M*(-M1*chi1_lambda + M2*chi2_lambda)
        chi_s_ell=chi1_ell/2 + chi2_ell/2
        chi_a_ell=chi1_ell/2 - chi2_ell/2
        logv=log(v)
        Fcal_coeff=32*nu**2*v**10/5
        Fcal_0=1
        Fcal_2=-35*nu/12 - 1247/336
        Fcal_3=4*pi
        Fcal_4=65*nu**2/18 + 9271*nu/504 - 44711/9072
        Fcal_5=pi*(-583*nu/24 - 8191/672)
        Fcal_6=-775*nu**3/324 - 94403*nu**2/3024 + nu*(-134543/7776 + 41*pi**2/48) - 1712*log(4)/105 - 1712*EulerGamma/105 + 16*pi**2/3 + 6643739519/69854400
        Fcal_lnv_6=-1712/105
        Fcal_7=pi*(193385*nu**2/3024 + 214745*nu/1728 - 16285/504)
        Fcal_8=-1369*pi**2/126 - 323105549467/3178375200 - 47385*log(3)/1568 + 232597*EulerGamma/4410 + 39931*log(2)/294
        Fcal_lnv_8=232597/4410
        Fcal_9=-13696*pi*log(2)/105 - 6848*EulerGamma*pi/105 + 265978667519*pi/745113600
        Fcal_lnv_9=-6848*pi/105
        Fcal_10=-2500861660823683/2831932303200 - 424223*pi**2/6804 - 83217611*log(2)/1122660 + 916628467*EulerGamma/7858620 + 47385*log(3)/196
        Fcal_lnv_10=916628467/7858620
        Fcal_SQ_4=chi1chi1*(-89*delta/192 + 89*nu/96 - 89/192) - 103*chi1chi2*nu/48 + chi2chi2*(89*delta/192 + 89*nu/96 - 89/192) + chi_a_ell*(chi_a_ell*(287/96 - 12*nu) + 287*chi_s_ell*delta/48) + chi_s_ell**2*(nu/24 + 287/96)
        Fcal_SO_3=(-4*S_ell - 5*Sigma_ell*delta/4)/M**2
        Fcal_SO_5=(S_ell*(272*nu/9 - 9/2) + Sigma_ell*delta*(43*nu/4 - 13/16))/M**2
        Fcal_SO_6=(-16*S_ell*pi - 31*Sigma_ell*delta*pi/6)/M**2
        Fcal_SO_7=(S_ell*(-2810*nu**2/27 + 6172*nu/189 + 476645/6804) + Sigma_ell*delta*(-1501*nu**2/36 + 1849*nu/126 + 9535/336))/M**2
        Fcal_SO_8=(S_ell*pi*(13879*nu/72 - 3485/96) + Sigma_ell*delta*pi*(130583*nu/2016 - 7163/672))/M**2
        E_0=1
        E_2=-nu/12 - 3/4
        E_4=-nu**2/24 + 19*nu/8 - 27/8
        E_6=-35*nu**3/5184 - 155*nu**2/96 + nu*(34445/576 - 205*pi**2/96) - 675/64
        E_8=77*nu**4/31104 + 301*nu**3/1728 + nu**2*(-498449/3456 + 3157*pi**2/576) + nu*(-123671/5760 + 896*EulerGamma/15 + 9037*pi**2/1536 + 1792*log(2)/15) - 3969/128
        E_lnv_8=896*nu/15
        E_10=nu**5/512 + 55*nu**4/512 + nu**3*(69423/512 - 1353*pi**2/256) + nu**2*(-21337*pi**2/1024 - 896*log(2)/5 - 448*EulerGamma/5 + 893429/2880) + nu*(-228916843/115200 - 23672*log(2)/35 - 9976*EulerGamma/35 + 729*log(3)/7 + 126779*pi**2/512) - 45927/512
        E_lnv_10=-1312*nu**2/5 - 9976*nu/35
        E_SQ_4=-3*chi_a_ell**2/2 - 3*chi_s_ell**2/2 - delta*(chi2chi2/2 + 3*chi_a_ell*chi_s_ell) + nu*(chi1chi2 + 6*chi_a_ell**2) + (chi1chi1 + chi2chi2)*(delta - 2*nu + 1)/4
        E_SO_3=(14*S_ell/3 + 2*Sigma_ell*delta)/M**2
        E_SO_5=(S_ell*(11 - 61*nu/9) + Sigma_ell*delta*(3 - 10*nu/3))/M**2
        E_SO_7=(S_ell*(29*nu**2/12 - 367*nu/4 + 135/4) + Sigma_ell*delta*(5*nu**2/4 - 39*nu + 27/4))/M**2
        Phi=0.0
        EvolveSpin1=dot(S_chi1,S_chi1)>1e-12
        EvolveSpin2=dot(S_chi2,S_chi2)>1e-12

        def Recalculate(t, y):
            v = y[0];
            rfrak_chi1_x = y[1];
            rfrak_chi1_y = y[2];
            rfrak_chi2_x = y[3];
            rfrak_chi2_y = y[4];
            rfrak_frame_x = y[5];
            rfrak_frame_y = y[6];
            rfrak_frame_z = y[7];
            Phi = y[8];
            R = exp(rfrak_frame_x*xHat + rfrak_frame_y*yHat + rfrak_frame_z*zHat)
            nHat = R*xHat*conjugate(R)
            lambdaHat = R*yHat*conjugate(R)
            ellHat = R*zHat*conjugate(R)
            R_S1 = exp(rfrak_chi1_x*xHat + rfrak_chi1_y*yHat)
            R_S2 = exp(rfrak_chi2_x*xHat + rfrak_chi2_y*yHat)
            chiVec1 = S_chi1*R_S1*zHat*conjugate(R_S1)*conjugate(S_chi1)
            chiVec2 = S_chi2*R_S2*zHat*conjugate(R_S2)*conjugate(S_chi2)
            chi1_n = dot(chiVec1,nHat)
            chi1_lambda = dot(chiVec1,lambdaHat)
            chi1_ell = dot(chiVec1,ellHat)
            chi2_n = dot(chiVec2,nHat)
            chi2_lambda = dot(chiVec2,lambdaHat)
            chi2_ell = dot(chiVec2,ellHat)
            S_ell = M1**2*chi1_ell + M2**2*chi2_ell
            S_n = M1**2*chi1_n + M2**2*chi2_n
            S_lambda = M1**2*chi1_lambda + M2**2*chi2_lambda
            Sigma_ell = M*(-M1*chi1_ell + M2*chi2_ell)
            Sigma_n = M*(-M1*chi1_n + M2*chi2_n)
            Sigma_lambda = M*(-M1*chi1_lambda + M2*chi2_lambda)
            chi_s_ell = chi1_ell/2 + chi2_ell/2
            chi_a_ell = chi1_ell/2 - chi2_ell/2
            logv = log(v)
            Fcal_coeff = 32*nu**2*v**10/5
            Fcal_SQ_4 = chi1chi1*(-89*delta/192 + 89*nu/96 - 89/192) - 103*chi1chi2*nu/48 + chi2chi2*(89*delta/192 + 89*nu/96 - 89/192) + chi_a_ell*(chi_a_ell*(287/96 - 12*nu) + 287*chi_s_ell*delta/48) + chi_s_ell**2*(nu/24 + 287/96)
            Fcal_SO_3 = (-4*S_ell - 5*Sigma_ell*delta/4)/M**2
            Fcal_SO_5 = (S_ell*(272*nu/9 - 9/2) + Sigma_ell*delta*(43*nu/4 - 13/16))/M**2
            Fcal_SO_6 = (-16*S_ell*pi - 31*Sigma_ell*delta*pi/6)/M**2
            Fcal_SO_7 = (S_ell*(-2810*nu**2/27 + 6172*nu/189 + 476645/6804) + Sigma_ell*delta*(-1501*nu**2/36 + 1849*nu/126 + 9535/336))/M**2
            Fcal_SO_8 = (S_ell*pi*(13879*nu/72 - 3485/96) + Sigma_ell*delta*pi*(130583*nu/2016 - 7163/672))/M**2
            E_SQ_4 = -3*chi_a_ell**2/2 - 3*chi_s_ell**2/2 - delta*(chi2chi2/2 + 3*chi_a_ell*chi_s_ell) + nu*(chi1chi2 + 6*chi_a_ell**2) + (chi1chi1 + chi2chi2)*(delta - 2*nu + 1)/4
            E_SO_3 = (14*S_ell/3 + 2*Sigma_ell*delta)/M**2
            E_SO_5 = (S_ell*(11 - 61*nu/9) + Sigma_ell*delta*(3 - 10*nu/3))/M**2
            E_SO_7 = (S_ell*(29*nu**2/12 - 367*nu/4 + 135/4) + Sigma_ell*delta*(5*nu**2/4 - 39*nu + 27/4))/M**2

        def OmegaVec_chiVec_1():
            Omega1_coeff = v**5/M
            return Omega1_coeff*(ellHat*(-0.75*delta + 0.5*nu + v**2*(delta*(0.625*nu - 0.5625) + nu*(1.25 - 0.0416666666666667*nu) + v**2*(delta*(nu*(4.875 - 0.15625*nu) - 0.84375) + nu*(nu*(-0.0208333333333333*nu - 3.28125) + 0.1875) + 0.84375) + 0.5625) + 0.75) + nHat*v*(3.0*chi1_n*nu + 3.0*M2**2*chi2_n/M**2) - M2**2*chiVec2*v/M**2)

        def OmegaVec_chiVec_2():
            Omega2_coeff = v**5/M
            return Omega2_coeff*(ellHat*(0.75*delta + 0.5*nu + v**2*(delta*(0.5625 - 0.625*nu) + nu*(1.25 - 0.0416666666666667*nu) + v**2*(delta*(nu*(0.15625*nu - 4.875) + 0.84375) + nu*(nu*(-0.0208333333333333*nu - 3.28125) + 0.1875) + 0.84375) + 0.5625) + 0.75) + nHat*v*(3.0*chi2_n*nu + 3.0*M1**2*chi1_n/M**2) - M1**2*chiVec1*v/M**2)

        def OmegaVec():
            gamma_PN_7 = (S_ell*(-6.0*nu**2 - 10.5833333333333*nu + 5.0) - 2.66666666666667*Sigma_ell*delta*nu**2 + Sigma_ell*delta*(3.0 - 10.1666666666667*nu))/M**2
            gamma_PN_6 = 0.0123456790123457*nu**3 + 6.36111111111111*nu**2 - 2.98177812235564*nu + 1.0
            gamma_PN_2 = 1.0 - 0.333333333333333*nu
            a_ell_2 = S_n*(-9.66666666666667*nu - 10.0) + Sigma_n*delta*(-4.5*nu - 6.0)
            a_ell_0 = 7.0*S_n + 3.0*Sigma_n*delta
            a_ell_4 = S_n*(5.77777777777778*nu**2 + 14.75*nu + 1.5) + Sigma_n*delta*(2.83333333333333*nu**2 + 9.125*nu + 1.5)
            gamma_PN_4 = 1.0 - 5.41666666666667*nu
            gamma_PN_3 = (1.66666666666667*S_ell + Sigma_ell*delta)/M**2
            gamma_PN_5 = (S_ell*(0.888888888888889*nu + 3.33333333333333) + 2.0*Sigma_ell*delta)/M**2
            gamma_PN_0 = 1.00000000000000
            return ellHat*v**3/M + nHat*v**6*(a_ell_0 + v**2*(a_ell_2 + a_ell_4*v**2))*(gamma_PN_0 + v**2*(gamma_PN_2 + v*(gamma_PN_3 + v*(gamma_PN_4 + v*(gamma_PN_5 + v*(gamma_PN_6 + gamma_PN_7*v))))))/M**3

        def OrbitalAngularMomentum():
            L_6 = ellHat*(0.00540123456790123*nu**3 + 1.29166666666667*nu**2 - 30.9797035925835*nu + 8.4375)
            L_4 = ellHat*(0.0416666666666667*nu**2 - 2.375*nu + 3.375)
            L_SO_5 = 0.0138888888888889*ellHat*(427.0*S_ell*nu - 693.0*S_ell + 210.0*Sigma_ell*delta*nu - 189.0*Sigma_ell*delta)/M**2 + 0.166666666666667*lambdaHat*(18.0*S_lambda*nu - 21.0*S_lambda + 8.0*Sigma_lambda*delta*nu - 3.0*Sigma_lambda*delta)/M**2 + 0.0416666666666667*nHat*(-19.0*S_n*nu + 33.0*S_n - 10.0*Sigma_n*delta*nu + 33.0*Sigma_n*delta)/M**2
            L_2 = ellHat*(0.166666666666667*nu + 1.5)
            L_10 = ellHat*(nu*(-4.85925925925926*nu - 5.27830687830688) + 59.80078125)
            L_0 = ellHat
            L_SO_3 = 0.166666666666667*ellHat*(-35.0*S_ell - 15.0*Sigma_ell*delta)/M**2 + lambdaHat*(-3.0*S_lambda - Sigma_lambda*delta)/M**2 + 0.5*nHat*(S_n + Sigma_n*delta)/M**2
            L_lnv_8 = -42.6666666666667*ellHat*nu
            L_lnv_10 = 2.0*ellHat*nu*(87.4666666666667*nu + 95.0095238095238)
            L_coeff = M**2*nu/v
            L_SO_7 = 0.0625*ellHat*(-29.0*S_ell*nu**2 + 1101.0*S_ell*nu - 405.0*S_ell - 15.0*Sigma_ell*delta*nu**2 + 468.0*Sigma_ell*delta*nu - 81.0*Sigma_ell*delta)/M**2 + 0.0416666666666667*lambdaHat*(-32.0*S_lambda*nu**2 + 2.0*S_lambda*nu - 174.0*S_lambda - 16.0*Sigma_lambda*delta*nu**2 - 79.0*Sigma_lambda*delta*nu - 12.0*Sigma_lambda*delta)/M**2 + 0.0208333333333333*nHat*(11.0*S_n*nu**2 - 1331.0*S_n*nu + 183.0*S_n + 5.0*Sigma_n*delta*nu**2 - 734.0*Sigma_n*delta*nu + 183.0*Sigma_n*delta)/M**2
            L_8 = ellHat*(-0.714285714285714*nu*(0.0024755658436214*nu**3 + 0.174189814814815*nu**2 - 90.1327990262052*nu + 153.88379682994) + 1.82857142857143*nu + 22.1484375)
            return L_coeff*(L_0 + v**2*(L_2 + v*(L_SO_3 + v*(L_4 + v*(L_SO_5 + v*(L_6 + v*(L_SO_7 + v*(L_8 + L_lnv_8*math.log(v) + v**2*(L_10 + L_lnv_10*math.log(v))))))))))


        #Evolve PN
        def TaylorT1(t, y):
            Recalculate(t, y)
            if v>=1.0: 
                return GSL_EDOM # Beyond domain of PN validity
            Flux = Fcal_coeff*(Fcal_0 + v**2*(Fcal_2 + v*(Fcal_3 + Fcal_SO_3 + v*(Fcal_4 + Fcal_SQ_4 + v*(Fcal_5 + Fcal_SO_5 + v*(Fcal_6 + Fcal_SO_6 + Fcal_lnv_6*logv + v*(Fcal_7 + Fcal_SO_7 + v*(Fcal_8 + Fcal_SO_8 + Fcal_lnv_8*logv + v*(Fcal_9 + Fcal_lnv_9*logv + v*(Fcal_10 + Fcal_lnv_10*logv))))))))))
            dEdv = -0.5*M*nu*v*(2.0*E_0 + v**2*(4.0*E_2 + v*(5.0*E_SO_3 + v*(6.0*E_4 + 6.0*E_SQ_4 + v*(7.0*E_SO_5 + v*(8.0*E_6 + v*(9.0*E_SO_7 + v*(10.0*E_8 + E_lnv_8*(10.0*logv + 1.0) + v**2*(12.0*E_10 + E_lnv_10*(12.0*logv + 1.0))))))))))
            Absorption = 0
            dvdt_T1 = (-Absorption - Flux)/dEdv
            if dvdt_T1<1.0e-12:
                return GSL_EDIVERGE # v is decreasing
            return CommonRHS(dvdt_T1, y)

        def CommonRHS(dvdt, y):
            dydt=np.zeros(9)
            rfrak_frame=np.zeros(3)
            rfrak_frame[0] = y[5]
            rfrak_frame[1] = y[6]
            rfrak_frame[2] = y[7]
            rfrakdot_frame = Quaternions.FrameFromAngularVelocity_Integrand(rfrak_frame, OmegaVec().vec())
            dydt[0] = dvdt
            if(EvolveSpin1):
                Quaternions.FrameFromAngularVelocity_2D_Integrand(y[1], y[2],\
                    (S_chi1.inverse()*OmegaVec_chiVec_1()*S_chi1).vec(),dydt[1], dydt[2])
            else:
                dydt[1] = 0.0
                dydt[2] = 0.0
            if(EvolveSpin2):
                Quaternions.FrameFromAngularVelocity_2D_Integrand(y[3], y[4],\
                    (S_chi2.inverse()*OmegaVec_chiVec_2()*S_chi2).vec(),dydt[3], dydt[4])
            else:
                dydt[3] = 0.0
                dydt[4] = 0.0
            dydt[5] = rfrakdot_frame[0]
            dydt[6] = rfrakdot_frame[1]
            dydt[7] = rfrakdot_frame[2]
            dydt[8] = v*v*v/M

            return dydt
        
        
        y=solve_ivp(TaylorT1, [0,20000], [0.01**(1/3),rfrak_chi1_x,\
            rfrak_chi1_y,rfrak_chi2_x,rfrak_chi2_y,rfrak_frame_x,\
            rfrak_frame_y,rfrak_frame_z,Phi], method='DOP853',\
            t_eval=np.arange(0,20000,100), dense_output=True)
            

#        def TaylorT1_phi(t,phi):
#            dphidt_T1=v_of_t.y[0][v_of_t.t>=t][0]**3.0/M
#            return dphidt_T1
#        phi_of_t=solve_ivp(TaylorT1_phi, [0, v_of_t.t[-1]], [0.0], method='RK45', t_eval=v_of_t.t)
        
        #Calculate waveform modes

        def WaveformModes():
            I=1j
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
            chiVec1 = Quaternions.Quaternion(0., chi1_n_k, chi1_lambda_k, chi1_ell_k);
            chiVec2 = Quaternions.Quaternion(0., chi2_n_k, chi2_lambda_k, chi2_ell_k);

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

            Modes=np.zeros(77,len(y.t))
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

class TaylorTn_5p5PN_Q : 
    def TaylorTn_5p5PN_Q(xHat_i, yHat_i, zHat_i, M1_i, M2_i, v_i, S_chi1_i, S_chi2_i, rfrak_chi1_x_i, rfrak_chi1_y_i, rfrak_chi2_x_i,
           rfrak_chi2_y_i, rfrak_frame_x_i, rfrak_frame_y_i, rfrak_frame_z_i) :
        xHat=xHat_i
        yHat=yHat_i
        zHat=zHat_i
        M1=M1_i
        M2=M2_i
        v=v_i
        S_chi1=S_chi1_i
        S_chi2=S_chi2_i
        rfrak_chi1_x=rfrak_chi1_x_i
        rfrak_chi1_y=rfrak_chi1_y_i
        rfrak_chi2_x=rfrak_chi2_x_i
        rfrak_chi2_y=rfrak_chi2_y_i
        rfrak_frame_x=rfrak_frame_x_i
        rfrak_frame_y=rfrak_frame_y_i
        rfrak_frame_z=rfrak_frame_z_i
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
        chi1chi1=dot(chiVec1,chiVec1)
        chi1chi2=dot(chiVec1,chiVec2)
        chi2chi2=dot(chiVec2,chiVec2)
        chi1_n=dot(chiVec1,nHat)
        chi1_lambda=dot(chiVec1,lambdaHat)
        chi1_ell=dot(chiVec1,ellHat)
        chi2_n=dot(chiVec2,nHat)
        chi2_lambda=dot(chiVec2,lambdaHat)
        chi2_ell=dot(chiVec2,ellHat)
        S_ell=M1**2*chi1_ell + M2**2*chi2_ell
        S_n=M1**2*chi1_n + M2**2*chi2_n
        S_lambda=M1**2*chi1_lambda + M2**2*chi2_lambda
        Sigma_ell=M*(-M1*chi1_ell + M2*chi2_ell)
        Sigma_n=M*(-M1*chi1_n + M2*chi2_n)
        Sigma_lambda=M*(-M1*chi1_lambda + M2*chi2_lambda)
        chi_s_ell=chi1_ell/2 + chi2_ell/2
        chi_a_ell=chi1_ell/2 - chi2_ell/2
        logv=log(v)
        Fcal_coeff=32*nu**2*v**10/5
        Fcal_0=1
        Fcal_2=-35*nu/12 - 1247/336
        Fcal_3=4*pi
        Fcal_4=65*nu**2/18 + 9271*nu/504 - 44711/9072
        Fcal_5=pi*(-583*nu/24 - 8191/672)
        Fcal_6=-775*nu**3/324 - 94403*nu**2/3024 + nu*(-134543/7776 + 41*pi**2/48) - 1712*log(4)/105 - 1712*EulerGamma/105 + 16*pi**2/3 + 6643739519/69854400
        Fcal_lnv_6=-1712/105
        Fcal_7=pi*(193385*nu**2/3024 + 214745*nu/1728 - 16285/504)
        Fcal_8=-1369*pi**2/126 - 323105549467/3178375200 - 47385*log(3)/1568 + 232597*EulerGamma/4410 + 39931*log(2)/294
        Fcal_lnv_8=232597/4410
        Fcal_9=-13696*pi*log(2)/105 - 6848*EulerGamma*pi/105 + 265978667519*pi/745113600
        Fcal_lnv_9=-6848*pi/105
        Fcal_10=-2500861660823683/2831932303200 - 424223*pi**2/6804 - 83217611*log(2)/1122660 + 916628467*EulerGamma/7858620 + 47385*log(3)/196
        Fcal_lnv_10=916628467/7858620
        Fcal_11=-142155*pi*log(3)/784 + 8399309750401*pi/101708006400 + 177293*EulerGamma*pi/1176 + 8521283*pi*log(2)/17640
        Fcal_lnv_11=177293*pi/1176
        Fcal_SQ_4=chi1chi1*(-89*delta/192 + 89*nu/96 - 89/192) - 103*chi1chi2*nu/48 + chi2chi2*(89*delta/192 + 89*nu/96 - 89/192) + chi_a_ell*(chi_a_ell*(287/96 - 12*nu) + 287*chi_s_ell*delta/48) + chi_s_ell**2*(nu/24 + 287/96)
        Fcal_SO_3=(-4*S_ell - 5*Sigma_ell*delta/4)/M**2
        Fcal_SO_5=(S_ell*(272*nu/9 - 9/2) + Sigma_ell*delta*(43*nu/4 - 13/16))/M**2
        Fcal_SO_6=(-16*S_ell*pi - 31*Sigma_ell*delta*pi/6)/M**2
        Fcal_SO_7=(S_ell*(-2810*nu**2/27 + 6172*nu/189 + 476645/6804) + Sigma_ell*delta*(-1501*nu**2/36 + 1849*nu/126 + 9535/336))/M**2
        Fcal_SO_8=(S_ell*pi*(13879*nu/72 - 3485/96) + Sigma_ell*delta*pi*(130583*nu/2016 - 7163/672))/M**2
        E_0=1
        E_2=-nu/12 - 3/4
        E_4=-nu**2/24 + 19*nu/8 - 27/8
        E_6=-35*nu**3/5184 - 155*nu**2/96 + nu*(34445/576 - 205*pi**2/96) - 675/64
        E_8=77*nu**4/31104 + 301*nu**3/1728 + nu**2*(-498449/3456 + 3157*pi**2/576) + nu*(-123671/5760 + 896*EulerGamma/15 + 9037*pi**2/1536 + 1792*log(2)/15) - 3969/128
        E_lnv_8=896*nu/15
        E_10=nu**5/512 + 55*nu**4/512 + nu**3*(69423/512 - 1353*pi**2/256) + nu**2*(-21337*pi**2/1024 - 896*log(2)/5 - 448*EulerGamma/5 + 893429/2880) + nu*(-228916843/115200 - 23672*log(2)/35 - 9976*EulerGamma/35 + 729*log(3)/7 + 126779*pi**2/512) - 45927/512
        E_lnv_10=-1312*nu**2/5 - 9976*nu/35
        E_11=27392*nu*pi/315
        E_SQ_4=-3*chi_a_ell**2/2 - 3*chi_s_ell**2/2 - delta*(chi2chi2/2 + 3*chi_a_ell*chi_s_ell) + nu*(chi1chi2 + 6*chi_a_ell**2) + (chi1chi1 + chi2chi2)*(delta - 2*nu + 1)/4
        E_SO_3=(14*S_ell/3 + 2*Sigma_ell*delta)/M**2
        E_SO_5=(S_ell*(11 - 61*nu/9) + Sigma_ell*delta*(3 - 10*nu/3))/M**2
        E_SO_7=(S_ell*(29*nu**2/12 - 367*nu/4 + 135/4) + Sigma_ell*delta*(5*nu**2/4 - 39*nu + 27/4))/M**2
        Phi=0.0
        EvolveSpin1=dot(S_chi1,S_chi1)>1e-12
        EvolveSpin2=dot(S_chi2,S_chi2)>1e-12

        def Recalculate(t, y):
            v = y[0];
            rfrak_chi1_x = y[1];
            rfrak_chi1_y = y[2];
            rfrak_chi2_x = y[3];
            rfrak_chi2_y = y[4];
            rfrak_frame_x = y[5];
            rfrak_frame_y = y[6];
            rfrak_frame_z = y[7];
            Phi = y[8];
            R = exp(rfrak_frame_x*xHat + rfrak_frame_y*yHat + rfrak_frame_z*zHat)
            nHat = R*xHat*conjugate(R)
            lambdaHat = R*yHat*conjugate(R)
            ellHat = R*zHat*conjugate(R)
            R_S1 = exp(rfrak_chi1_x*xHat + rfrak_chi1_y*yHat)
            R_S2 = exp(rfrak_chi2_x*xHat + rfrak_chi2_y*yHat)
            chiVec1 = S_chi1*R_S1*zHat*conjugate(R_S1)*conjugate(S_chi1)
            chiVec2 = S_chi2*R_S2*zHat*conjugate(R_S2)*conjugate(S_chi2)
            chi1_n = dot(chiVec1,nHat)
            chi1_lambda = dot(chiVec1,lambdaHat)
            chi1_ell = dot(chiVec1,ellHat)
            chi2_n = dot(chiVec2,nHat)
            chi2_lambda = dot(chiVec2,lambdaHat)
            chi2_ell = dot(chiVec2,ellHat)
            S_ell = M1**2*chi1_ell + M2**2*chi2_ell
            S_n = M1**2*chi1_n + M2**2*chi2_n
            S_lambda = M1**2*chi1_lambda + M2**2*chi2_lambda
            Sigma_ell = M*(-M1*chi1_ell + M2*chi2_ell)
            Sigma_n = M*(-M1*chi1_n + M2*chi2_n)
            Sigma_lambda = M*(-M1*chi1_lambda + M2*chi2_lambda)
            chi_s_ell = chi1_ell/2 + chi2_ell/2
            chi_a_ell = chi1_ell/2 - chi2_ell/2
            logv = log(v)
            Fcal_coeff = 32*nu**2*v**10/5
            Fcal_SQ_4 = chi1chi1*(-89*delta/192 + 89*nu/96 - 89/192) - 103*chi1chi2*nu/48 + chi2chi2*(89*delta/192 + 89*nu/96 - 89/192) + chi_a_ell*(chi_a_ell*(287/96 - 12*nu) + 287*chi_s_ell*delta/48) + chi_s_ell**2*(nu/24 + 287/96)
            Fcal_SO_3 = (-4*S_ell - 5*Sigma_ell*delta/4)/M**2
            Fcal_SO_5 = (S_ell*(272*nu/9 - 9/2) + Sigma_ell*delta*(43*nu/4 - 13/16))/M**2
            Fcal_SO_6 = (-16*S_ell*pi - 31*Sigma_ell*delta*pi/6)/M**2
            Fcal_SO_7 = (S_ell*(-2810*nu**2/27 + 6172*nu/189 + 476645/6804) + Sigma_ell*delta*(-1501*nu**2/36 + 1849*nu/126 + 9535/336))/M**2
            Fcal_SO_8 = (S_ell*pi*(13879*nu/72 - 3485/96) + Sigma_ell*delta*pi*(130583*nu/2016 - 7163/672))/M**2
            E_SQ_4 = -3*chi_a_ell**2/2 - 3*chi_s_ell**2/2 - delta*(chi2chi2/2 + 3*chi_a_ell*chi_s_ell) + nu*(chi1chi2 + 6*chi_a_ell**2) + (chi1chi1 + chi2chi2)*(delta - 2*nu + 1)/4
            E_SO_3 = (14*S_ell/3 + 2*Sigma_ell*delta)/M**2
            E_SO_5 = (S_ell*(11 - 61*nu/9) + Sigma_ell*delta*(3 - 10*nu/3))/M**2
            E_SO_7 = (S_ell*(29*nu**2/12 - 367*nu/4 + 135/4) + Sigma_ell*delta*(5*nu**2/4 - 39*nu + 27/4))/M**2

        def OmegaVec_chiVec_1():
            Omega1_coeff = v**5/M
            return Omega1_coeff*(ellHat*(-0.75*delta + 0.5*nu + v**2*(delta*(0.625*nu - 0.5625) + nu*(1.25 - 0.0416666666666667*nu) + v**2*(delta*(nu*(4.875 - 0.15625*nu) - 0.84375) + nu*(nu*(-0.0208333333333333*nu - 3.28125) + 0.1875) + 0.84375) + 0.5625) + 0.75) + nHat*v*(3.0*chi1_n*nu + 3.0*M2**2*chi2_n/M**2) - M2**2*chiVec2*v/M**2)

        def OmegaVec_chiVec_2():
            Omega2_coeff = v**5/M
            return Omega2_coeff*(ellHat*(0.75*delta + 0.5*nu + v**2*(delta*(0.5625 - 0.625*nu) + nu*(1.25 - 0.0416666666666667*nu) + v**2*(delta*(nu*(0.15625*nu - 4.875) + 0.84375) + nu*(nu*(-0.0208333333333333*nu - 3.28125) + 0.1875) + 0.84375) + 0.5625) + 0.75) + nHat*v*(3.0*chi2_n*nu + 3.0*M1**2*chi1_n/M**2) - M1**2*chiVec1*v/M**2)

        def OmegaVec():
            gamma_PN_7 = (S_ell*(-6.0*nu**2 - 10.5833333333333*nu + 5.0) - 2.66666666666667*Sigma_ell*delta*nu**2 + Sigma_ell*delta*(3.0 - 10.1666666666667*nu))/M**2
            gamma_PN_6 = 0.0123456790123457*nu**3 + 6.36111111111111*nu**2 - 2.98177812235564*nu + 1.0
            gamma_PN_2 = 1.0 - 0.333333333333333*nu
            a_ell_2 = S_n*(-9.66666666666667*nu - 10.0) + Sigma_n*delta*(-4.5*nu - 6.0)
            a_ell_0 = 7.0*S_n + 3.0*Sigma_n*delta
            a_ell_4 = S_n*(5.77777777777778*nu**2 + 14.75*nu + 1.5) + Sigma_n*delta*(2.83333333333333*nu**2 + 9.125*nu + 1.5)
            gamma_PN_4 = 1.0 - 5.41666666666667*nu
            gamma_PN_3 = (1.66666666666667*S_ell + Sigma_ell*delta)/M**2
            gamma_PN_5 = (S_ell*(0.888888888888889*nu + 3.33333333333333) + 2.0*Sigma_ell*delta)/M**2
            gamma_PN_0 = 1.00000000000000
            return ellHat*v**3/M + nHat*v**6*(a_ell_0 + v**2*(a_ell_2 + a_ell_4*v**2))*(gamma_PN_0 + v**2*(gamma_PN_2 + v*(gamma_PN_3 + v*(gamma_PN_4 + v*(gamma_PN_5 + v*(gamma_PN_6 + gamma_PN_7*v))))))/M**3

        def OrbitalAngularMomentum():
            L_6 = ellHat*(0.00540123456790123*nu**3 + 1.29166666666667*nu**2 - 30.9797035925835*nu + 8.4375)
            L_4 = ellHat*(0.0416666666666667*nu**2 - 2.375*nu + 3.375)
            L_SO_5 = 0.0138888888888889*ellHat*(427.0*S_ell*nu - 693.0*S_ell + 210.0*Sigma_ell*delta*nu - 189.0*Sigma_ell*delta)/M**2 + 0.166666666666667*lambdaHat*(18.0*S_lambda*nu - 21.0*S_lambda + 8.0*Sigma_lambda*delta*nu - 3.0*Sigma_lambda*delta)/M**2 + 0.0416666666666667*nHat*(-19.0*S_n*nu + 33.0*S_n - 10.0*Sigma_n*delta*nu + 33.0*Sigma_n*delta)/M**2
            L_2 = ellHat*(0.166666666666667*nu + 1.5)
            L_10 = ellHat*(nu*(-4.85925925925926*nu - 5.27830687830688) + 59.80078125)
            L_0 = ellHat
            L_SO_3 = 0.166666666666667*ellHat*(-35.0*S_ell - 15.0*Sigma_ell*delta)/M**2 + lambdaHat*(-3.0*S_lambda - Sigma_lambda*delta)/M**2 + 0.5*nHat*(S_n + Sigma_n*delta)/M**2
            L_lnv_8 = -42.6666666666667*ellHat*nu
            L_lnv_10 = 2.0*ellHat*nu*(87.4666666666667*nu + 95.0095238095238)
            L_coeff = M**2*nu/v
            L_SO_7 = 0.0625*ellHat*(-29.0*S_ell*nu**2 + 1101.0*S_ell*nu - 405.0*S_ell - 15.0*Sigma_ell*delta*nu**2 + 468.0*Sigma_ell*delta*nu - 81.0*Sigma_ell*delta)/M**2 + 0.0416666666666667*lambdaHat*(-32.0*S_lambda*nu**2 + 2.0*S_lambda*nu - 174.0*S_lambda - 16.0*Sigma_lambda*delta*nu**2 - 79.0*Sigma_lambda*delta*nu - 12.0*Sigma_lambda*delta)/M**2 + 0.0208333333333333*nHat*(11.0*S_n*nu**2 - 1331.0*S_n*nu + 183.0*S_n + 5.0*Sigma_n*delta*nu**2 - 734.0*Sigma_n*delta*nu + 183.0*Sigma_n*delta)/M**2
            L_8 = ellHat*(-0.714285714285714*nu*(0.0024755658436214*nu**3 + 0.174189814814815*nu**2 - 90.1327990262052*nu + 153.88379682994) + 1.82857142857143*nu + 22.1484375)
            return L_coeff*(L_0 + v**2*(L_2 + v*(L_SO_3 + v*(L_4 + v*(L_SO_5 + v*(L_6 + v*(L_SO_7 + v*(L_8 + L_lnv_8*math.log(v) + v**2*(L_10 + L_lnv_10*math.log(v))))))))))


        #Evolve PN
        def TaylorT1(t, y):
            Recalculate(t, y)
            if v>=1.0: 
                return GSL_EDOM # Beyond domain of PN validity
            Flux = Fcal_coeff*(Fcal_0 + v**2*(Fcal_2 + v*(Fcal_3 + Fcal_SO_3 + v*(Fcal_4 + Fcal_SQ_4 + v*(Fcal_5 + Fcal_SO_5 + v*(Fcal_6 + Fcal_SO_6 + Fcal_lnv_6*logv + v*(Fcal_7 + Fcal_SO_7 + v*(Fcal_8 + Fcal_SO_8 + Fcal_lnv_8*logv + v*(Fcal_9 + Fcal_lnv_9*logv + v*(Fcal_10 + Fcal_lnv_10*logv + v*(Fcal_11 + Fcal_lnv_11*logv)))))))))))
            dEdv = -0.5*M*nu*v*(2.0*E_0 + v**2*(4.0*E_2 + v*(5.0*E_SO_3 + v*(6.0*E_4 + 6.0*E_SQ_4 + v*(7.0*E_SO_5 + v*(8.0*E_6 + v*(9.0*E_SO_7 + v*(10.0*E_8 + E_lnv_8*(10.0*logv + 1.0) + v**2*(12.0*E_10 + 13.0*E_11*v + E_lnv_10*(12.0*logv + 1.0))))))))))
            Absorption = 0
            dvdt_T1 = (-Absorption - Flux)/dEdv
            if dvdt_T1<1.0e-12:
                return GSL_EDIVERGE # v is decreasing
            return CommonRHS(dvdt_T1, y)

        def CommonRHS(dvdt, y):
            dydt=np.zeros(9)
            rfrak_frame=np.zeros(3)
            rfrak_frame[0] = y[5]
            rfrak_frame[1] = y[6]
            rfrak_frame[2] = y[7]
            rfrakdot_frame = Quaternions.FrameFromAngularVelocity_Integrand(rfrak_frame, OmegaVec().vec())
            dydt[0] = dvdt
            if(EvolveSpin1):
                Quaternions.FrameFromAngularVelocity_2D_Integrand(y[1], y[2],\
                    (S_chi1.inverse()*OmegaVec_chiVec_1()*S_chi1).vec(),dydt[1], dydt[2])
            else:
                dydt[1] = 0.0
                dydt[2] = 0.0
            if(EvolveSpin2):
                Quaternions.FrameFromAngularVelocity_2D_Integrand(y[3], y[4],\
                    (S_chi2.inverse()*OmegaVec_chiVec_2()*S_chi2).vec(),dydt[3], dydt[4])
            else:
                dydt[3] = 0.0
                dydt[4] = 0.0
            dydt[5] = rfrakdot_frame[0]
            dydt[6] = rfrakdot_frame[1]
            dydt[7] = rfrakdot_frame[2]
            dydt[8] = v*v*v/M

            return dydt
        
        
        y=solve_ivp(TaylorT1, [0,20000], [0.01**(1/3),rfrak_chi1_x,\
            rfrak_chi1_y,rfrak_chi2_x,rfrak_chi2_y,rfrak_frame_x,\
            rfrak_frame_y,rfrak_frame_z,Phi], method='DOP853',\
            t_eval=np.arange(0,20000,100), dense_output=True)
            

#        def TaylorT1_phi(t,phi):
#            dphidt_T1=v_of_t.y[0][v_of_t.t>=t][0]**3.0/M
#            return dphidt_T1
#        phi_of_t=solve_ivp(TaylorT1_phi, [0, v_of_t.t[-1]], [0.0], method='RK45', t_eval=v_of_t.t)
        
        #Calculate waveform modes

        def WaveformModes():
            I=1j
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
            chiVec1 = Quaternions.Quaternion(0., chi1_n_k, chi1_lambda_k, chi1_ell_k);
            chiVec2 = Quaternions.Quaternion(0., chi2_n_k, chi2_lambda_k, chi2_ell_k);

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

            Modes=np.zeros(77,len(y.t))
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

class TaylorTn_6p0PN_Q : 
    def TaylorTn_6p0PN_Q(xHat_i, yHat_i, zHat_i, M1_i, M2_i, v_i, S_chi1_i, S_chi2_i, rfrak_chi1_x_i, rfrak_chi1_y_i, rfrak_chi2_x_i,
           rfrak_chi2_y_i, rfrak_frame_x_i, rfrak_frame_y_i, rfrak_frame_z_i) :
        xHat=xHat_i
        yHat=yHat_i
        zHat=zHat_i
        M1=M1_i
        M2=M2_i
        v=v_i
        S_chi1=S_chi1_i
        S_chi2=S_chi2_i
        rfrak_chi1_x=rfrak_chi1_x_i
        rfrak_chi1_y=rfrak_chi1_y_i
        rfrak_chi2_x=rfrak_chi2_x_i
        rfrak_chi2_y=rfrak_chi2_y_i
        rfrak_frame_x=rfrak_frame_x_i
        rfrak_frame_y=rfrak_frame_y_i
        rfrak_frame_z=rfrak_frame_z_i
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
        chi1chi1=dot(chiVec1,chiVec1)
        chi1chi2=dot(chiVec1,chiVec2)
        chi2chi2=dot(chiVec2,chiVec2)
        chi1_n=dot(chiVec1,nHat)
        chi1_lambda=dot(chiVec1,lambdaHat)
        chi1_ell=dot(chiVec1,ellHat)
        chi2_n=dot(chiVec2,nHat)
        chi2_lambda=dot(chiVec2,lambdaHat)
        chi2_ell=dot(chiVec2,ellHat)
        S_ell=M1**2*chi1_ell + M2**2*chi2_ell
        S_n=M1**2*chi1_n + M2**2*chi2_n
        S_lambda=M1**2*chi1_lambda + M2**2*chi2_lambda
        Sigma_ell=M*(-M1*chi1_ell + M2*chi2_ell)
        Sigma_n=M*(-M1*chi1_n + M2*chi2_n)
        Sigma_lambda=M*(-M1*chi1_lambda + M2*chi2_lambda)
        chi_s_ell=chi1_ell/2 + chi2_ell/2
        chi_a_ell=chi1_ell/2 - chi2_ell/2
        logv=log(v)
        Fcal_coeff=32*nu**2*v**10/5
        Fcal_0=1
        Fcal_2=-35*nu/12 - 1247/336
        Fcal_3=4*pi
        Fcal_4=65*nu**2/18 + 9271*nu/504 - 44711/9072
        Fcal_5=pi*(-583*nu/24 - 8191/672)
        Fcal_6=-775*nu**3/324 - 94403*nu**2/3024 + nu*(-134543/7776 + 41*pi**2/48) - 1712*log(4)/105 - 1712*EulerGamma/105 + 16*pi**2/3 + 6643739519/69854400
        Fcal_lnv_6=-1712/105
        Fcal_7=pi*(193385*nu**2/3024 + 214745*nu/1728 - 16285/504)
        Fcal_8=-1369*pi**2/126 - 323105549467/3178375200 - 47385*log(3)/1568 + 232597*EulerGamma/4410 + 39931*log(2)/294
        Fcal_lnv_8=232597/4410
        Fcal_9=-13696*pi*log(2)/105 - 6848*EulerGamma*pi/105 + 265978667519*pi/745113600
        Fcal_lnv_9=-6848*pi/105
        Fcal_10=-2500861660823683/2831932303200 - 424223*pi**2/6804 - 83217611*log(2)/1122660 + 916628467*EulerGamma/7858620 + 47385*log(3)/196
        Fcal_lnv_10=916628467/7858620
        Fcal_11=-142155*pi*log(3)/784 + 8399309750401*pi/101708006400 + 177293*EulerGamma*pi/1176 + 8521283*pi*log(2)/17640
        Fcal_lnv_11=177293*pi/1176
        Fcal_12=-271272899815409*log(2)/157329572400 - 54784*pi**2*log(2)/315 - 246137536815857*EulerGamma/157329572400 - 437114506833*log(3)/789268480 - 256*pi**4/45 - 27392*EulerGamma*pi**2/315 - 27392*zeta(3)/105 - 37744140625*log(5)/260941824 + 1465472*EulerGamma**2/11025 + 5861888*EulerGamma*log(2)/11025 + 5861888*log(2)**2/11025 + 2067586193789233570693/602387400044430000 + 3803225263*pi**2/10478160
        Fcal_lnv_12=-246137536815857/157329572400 - 27392*pi**2/315 + 2930944*EulerGamma/11025 + 5861888*log(2)/11025
        Fcal_lnv2_12=1465472/11025
        Fcal_SQ_4=chi1chi1*(-89*delta/192 + 89*nu/96 - 89/192) - 103*chi1chi2*nu/48 + chi2chi2*(89*delta/192 + 89*nu/96 - 89/192) + chi_a_ell*(chi_a_ell*(287/96 - 12*nu) + 287*chi_s_ell*delta/48) + chi_s_ell**2*(nu/24 + 287/96)
        Fcal_SO_3=(-4*S_ell - 5*Sigma_ell*delta/4)/M**2
        Fcal_SO_5=(S_ell*(272*nu/9 - 9/2) + Sigma_ell*delta*(43*nu/4 - 13/16))/M**2
        Fcal_SO_6=(-16*S_ell*pi - 31*Sigma_ell*delta*pi/6)/M**2
        Fcal_SO_7=(S_ell*(-2810*nu**2/27 + 6172*nu/189 + 476645/6804) + Sigma_ell*delta*(-1501*nu**2/36 + 1849*nu/126 + 9535/336))/M**2
        Fcal_SO_8=(S_ell*pi*(13879*nu/72 - 3485/96) + Sigma_ell*delta*pi*(130583*nu/2016 - 7163/672))/M**2
        E_0=1
        E_2=-nu/12 - 3/4
        E_4=-nu**2/24 + 19*nu/8 - 27/8
        E_6=-35*nu**3/5184 - 155*nu**2/96 + nu*(34445/576 - 205*pi**2/96) - 675/64
        E_8=77*nu**4/31104 + 301*nu**3/1728 + nu**2*(-498449/3456 + 3157*pi**2/576) + nu*(-123671/5760 + 896*EulerGamma/15 + 9037*pi**2/1536 + 1792*log(2)/15) - 3969/128
        E_lnv_8=896*nu/15
        E_10=nu**5/512 + 55*nu**4/512 + nu**3*(69423/512 - 1353*pi**2/256) + nu**2*(-21337*pi**2/1024 - 896*log(2)/5 - 448*EulerGamma/5 + 893429/2880) + nu*(-228916843/115200 - 23672*log(2)/35 - 9976*EulerGamma/35 + 729*log(3)/7 + 126779*pi**2/512) - 45927/512
        E_lnv_10=-1312*nu**2/5 - 9976*nu/35
        E_11=27392*nu*pi/315
        E_12=2717*nu**6/6718464 + 5159*nu**5/248832 + nu**4*(-20543435/373248 + 272855*pi**2/124416) + nu**3*(-71700787/51840 + 1232*EulerGamma/27 + 2464*log(2)/27 + 6634243*pi**2/110592) + nu**2*(-86017789*pi**2/110592 - 2673*log(3)/14 + 112772*EulerGamma/105 + 18491*pi**4/2304 + 246004*log(2)/105 + 113176680983/14515200) + nu*(-389727504721/43545600 - 30809603*pi**4/786432 - 7128*log(3)/7 - 3934568*EulerGamma/8505 + 74888*log(2)/243 + 9118627045*pi**2/5308416) - 264627/1024
        E_lnv_12=48928*nu**3/135 + 79508*nu**2/105 - 3934568*nu/8505
        E_SQ_4=-3*chi_a_ell**2/2 - 3*chi_s_ell**2/2 - delta*(chi2chi2/2 + 3*chi_a_ell*chi_s_ell) + nu*(chi1chi2 + 6*chi_a_ell**2) + (chi1chi1 + chi2chi2)*(delta - 2*nu + 1)/4
        E_SO_3=(14*S_ell/3 + 2*Sigma_ell*delta)/M**2
        E_SO_5=(S_ell*(11 - 61*nu/9) + Sigma_ell*delta*(3 - 10*nu/3))/M**2
        E_SO_7=(S_ell*(29*nu**2/12 - 367*nu/4 + 135/4) + Sigma_ell*delta*(5*nu**2/4 - 39*nu + 27/4))/M**2
        Phi=0.0
        EvolveSpin1=dot(S_chi1,S_chi1)>1e-12
        EvolveSpin2=dot(S_chi2,S_chi2)>1e-12

        def Recalculate(t, y):
            v = y[0];
            rfrak_chi1_x = y[1];
            rfrak_chi1_y = y[2];
            rfrak_chi2_x = y[3];
            rfrak_chi2_y = y[4];
            rfrak_frame_x = y[5];
            rfrak_frame_y = y[6];
            rfrak_frame_z = y[7];
            Phi = y[8];
            R = exp(rfrak_frame_x*xHat + rfrak_frame_y*yHat + rfrak_frame_z*zHat)
            nHat = R*xHat*conjugate(R)
            lambdaHat = R*yHat*conjugate(R)
            ellHat = R*zHat*conjugate(R)
            R_S1 = exp(rfrak_chi1_x*xHat + rfrak_chi1_y*yHat)
            R_S2 = exp(rfrak_chi2_x*xHat + rfrak_chi2_y*yHat)
            chiVec1 = S_chi1*R_S1*zHat*conjugate(R_S1)*conjugate(S_chi1)
            chiVec2 = S_chi2*R_S2*zHat*conjugate(R_S2)*conjugate(S_chi2)
            chi1_n = dot(chiVec1,nHat)
            chi1_lambda = dot(chiVec1,lambdaHat)
            chi1_ell = dot(chiVec1,ellHat)
            chi2_n = dot(chiVec2,nHat)
            chi2_lambda = dot(chiVec2,lambdaHat)
            chi2_ell = dot(chiVec2,ellHat)
            S_ell = M1**2*chi1_ell + M2**2*chi2_ell
            S_n = M1**2*chi1_n + M2**2*chi2_n
            S_lambda = M1**2*chi1_lambda + M2**2*chi2_lambda
            Sigma_ell = M*(-M1*chi1_ell + M2*chi2_ell)
            Sigma_n = M*(-M1*chi1_n + M2*chi2_n)
            Sigma_lambda = M*(-M1*chi1_lambda + M2*chi2_lambda)
            chi_s_ell = chi1_ell/2 + chi2_ell/2
            chi_a_ell = chi1_ell/2 - chi2_ell/2
            logv = log(v)
            Fcal_coeff = 32*nu**2*v**10/5
            Fcal_SQ_4 = chi1chi1*(-89*delta/192 + 89*nu/96 - 89/192) - 103*chi1chi2*nu/48 + chi2chi2*(89*delta/192 + 89*nu/96 - 89/192) + chi_a_ell*(chi_a_ell*(287/96 - 12*nu) + 287*chi_s_ell*delta/48) + chi_s_ell**2*(nu/24 + 287/96)
            Fcal_SO_3 = (-4*S_ell - 5*Sigma_ell*delta/4)/M**2
            Fcal_SO_5 = (S_ell*(272*nu/9 - 9/2) + Sigma_ell*delta*(43*nu/4 - 13/16))/M**2
            Fcal_SO_6 = (-16*S_ell*pi - 31*Sigma_ell*delta*pi/6)/M**2
            Fcal_SO_7 = (S_ell*(-2810*nu**2/27 + 6172*nu/189 + 476645/6804) + Sigma_ell*delta*(-1501*nu**2/36 + 1849*nu/126 + 9535/336))/M**2
            Fcal_SO_8 = (S_ell*pi*(13879*nu/72 - 3485/96) + Sigma_ell*delta*pi*(130583*nu/2016 - 7163/672))/M**2
            E_SQ_4 = -3*chi_a_ell**2/2 - 3*chi_s_ell**2/2 - delta*(chi2chi2/2 + 3*chi_a_ell*chi_s_ell) + nu*(chi1chi2 + 6*chi_a_ell**2) + (chi1chi1 + chi2chi2)*(delta - 2*nu + 1)/4
            E_SO_3 = (14*S_ell/3 + 2*Sigma_ell*delta)/M**2
            E_SO_5 = (S_ell*(11 - 61*nu/9) + Sigma_ell*delta*(3 - 10*nu/3))/M**2
            E_SO_7 = (S_ell*(29*nu**2/12 - 367*nu/4 + 135/4) + Sigma_ell*delta*(5*nu**2/4 - 39*nu + 27/4))/M**2

        def OmegaVec_chiVec_1():
            Omega1_coeff = v**5/M
            return Omega1_coeff*(ellHat*(-0.75*delta + 0.5*nu + v**2*(delta*(0.625*nu - 0.5625) + nu*(1.25 - 0.0416666666666667*nu) + v**2*(delta*(nu*(4.875 - 0.15625*nu) - 0.84375) + nu*(nu*(-0.0208333333333333*nu - 3.28125) + 0.1875) + 0.84375) + 0.5625) + 0.75) + nHat*v*(3.0*chi1_n*nu + 3.0*M2**2*chi2_n/M**2) - M2**2*chiVec2*v/M**2)

        def OmegaVec_chiVec_2():
            Omega2_coeff = v**5/M
            return Omega2_coeff*(ellHat*(0.75*delta + 0.5*nu + v**2*(delta*(0.5625 - 0.625*nu) + nu*(1.25 - 0.0416666666666667*nu) + v**2*(delta*(nu*(0.15625*nu - 4.875) + 0.84375) + nu*(nu*(-0.0208333333333333*nu - 3.28125) + 0.1875) + 0.84375) + 0.5625) + 0.75) + nHat*v*(3.0*chi2_n*nu + 3.0*M1**2*chi1_n/M**2) - M1**2*chiVec1*v/M**2)

        def OmegaVec():
            gamma_PN_7 = (S_ell*(-6.0*nu**2 - 10.5833333333333*nu + 5.0) - 2.66666666666667*Sigma_ell*delta*nu**2 + Sigma_ell*delta*(3.0 - 10.1666666666667*nu))/M**2
            gamma_PN_6 = 0.0123456790123457*nu**3 + 6.36111111111111*nu**2 - 2.98177812235564*nu + 1.0
            gamma_PN_2 = 1.0 - 0.333333333333333*nu
            a_ell_2 = S_n*(-9.66666666666667*nu - 10.0) + Sigma_n*delta*(-4.5*nu - 6.0)
            a_ell_0 = 7.0*S_n + 3.0*Sigma_n*delta
            a_ell_4 = S_n*(5.77777777777778*nu**2 + 14.75*nu + 1.5) + Sigma_n*delta*(2.83333333333333*nu**2 + 9.125*nu + 1.5)
            gamma_PN_4 = 1.0 - 5.41666666666667*nu
            gamma_PN_3 = (1.66666666666667*S_ell + Sigma_ell*delta)/M**2
            gamma_PN_5 = (S_ell*(0.888888888888889*nu + 3.33333333333333) + 2.0*Sigma_ell*delta)/M**2
            gamma_PN_0 = 1.00000000000000
            return ellHat*v**3/M + nHat*v**6*(a_ell_0 + v**2*(a_ell_2 + a_ell_4*v**2))*(gamma_PN_0 + v**2*(gamma_PN_2 + v*(gamma_PN_3 + v*(gamma_PN_4 + v*(gamma_PN_5 + v*(gamma_PN_6 + gamma_PN_7*v))))))/M**3

        def OrbitalAngularMomentum():
            L_6 = ellHat*(0.00540123456790123*nu**3 + 1.29166666666667*nu**2 - 30.9797035925835*nu + 8.4375)
            L_4 = ellHat*(0.0416666666666667*nu**2 - 2.375*nu + 3.375)
            L_SO_5 = 0.0138888888888889*ellHat*(427.0*S_ell*nu - 693.0*S_ell + 210.0*Sigma_ell*delta*nu - 189.0*Sigma_ell*delta)/M**2 + 0.166666666666667*lambdaHat*(18.0*S_lambda*nu - 21.0*S_lambda + 8.0*Sigma_lambda*delta*nu - 3.0*Sigma_lambda*delta)/M**2 + 0.0416666666666667*nHat*(-19.0*S_n*nu + 33.0*S_n - 10.0*Sigma_n*delta*nu + 33.0*Sigma_n*delta)/M**2
            L_2 = ellHat*(0.166666666666667*nu + 1.5)
            L_10 = ellHat*(nu*(-4.85925925925926*nu - 5.27830687830688) + 59.80078125)
            L_0 = ellHat
            L_SO_3 = 0.166666666666667*ellHat*(-35.0*S_ell - 15.0*Sigma_ell*delta)/M**2 + lambdaHat*(-3.0*S_lambda - Sigma_lambda*delta)/M**2 + 0.5*nHat*(S_n + Sigma_n*delta)/M**2
            L_lnv_8 = -42.6666666666667*ellHat*nu
            L_lnv_10 = 2.0*ellHat*nu*(87.4666666666667*nu + 95.0095238095238)
            L_coeff = M**2*nu/v
            L_SO_7 = 0.0625*ellHat*(-29.0*S_ell*nu**2 + 1101.0*S_ell*nu - 405.0*S_ell - 15.0*Sigma_ell*delta*nu**2 + 468.0*Sigma_ell*delta*nu - 81.0*Sigma_ell*delta)/M**2 + 0.0416666666666667*lambdaHat*(-32.0*S_lambda*nu**2 + 2.0*S_lambda*nu - 174.0*S_lambda - 16.0*Sigma_lambda*delta*nu**2 - 79.0*Sigma_lambda*delta*nu - 12.0*Sigma_lambda*delta)/M**2 + 0.0208333333333333*nHat*(11.0*S_n*nu**2 - 1331.0*S_n*nu + 183.0*S_n + 5.0*Sigma_n*delta*nu**2 - 734.0*Sigma_n*delta*nu + 183.0*Sigma_n*delta)/M**2
            L_8 = ellHat*(-0.714285714285714*nu*(0.0024755658436214*nu**3 + 0.174189814814815*nu**2 - 90.1327990262052*nu + 153.88379682994) + 1.82857142857143*nu + 22.1484375)
            return L_coeff*(L_0 + v**2*(L_2 + v*(L_SO_3 + v*(L_4 + v*(L_SO_5 + v*(L_6 + v*(L_SO_7 + v*(L_8 + L_lnv_8*math.log(v) + v**2*(L_10 + L_lnv_10*math.log(v))))))))))


        #Evolve PN
        def TaylorT1(t, y):
            Recalculate(t, y)
            if v>=1.0: 
                return GSL_EDOM # Beyond domain of PN validity
            Flux = Fcal_coeff*(Fcal_0 + v**2*(Fcal_2 + v*(Fcal_3 + Fcal_SO_3 + v*(Fcal_4 + Fcal_SQ_4 + v*(Fcal_5 + Fcal_SO_5 + v*(Fcal_6 + Fcal_SO_6 + Fcal_lnv_6*logv + v*(Fcal_7 + Fcal_SO_7 + v*(Fcal_8 + Fcal_SO_8 + Fcal_lnv_8*logv + v*(Fcal_9 + Fcal_lnv_9*logv + v*(Fcal_10 + Fcal_lnv_10*logv + v*(Fcal_11 + Fcal_lnv_11*logv + v*(Fcal_12 + Fcal_lnv2_12*logv**2 + Fcal_lnv_12*logv))))))))))))
            dEdv = -0.5*M*nu*v*(2.0*E_0 + v**2*(4.0*E_2 + v*(5.0*E_SO_3 + v*(6.0*E_4 + 6.0*E_SQ_4 + v*(7.0*E_SO_5 + v*(8.0*E_6 + v*(9.0*E_SO_7 + v*(10.0*E_8 + E_lnv_8*(10.0*logv + 1.0) + v**2*(12.0*E_10 + E_lnv_10*(12.0*logv + 1.0) + v*(13.0*E_11 + v*(14.0*E_12 + E_lnv_12*(14.0*logv + 1.0))))))))))))
            Absorption = 0
            dvdt_T1 = (-Absorption - Flux)/dEdv
            if dvdt_T1<1.0e-12:
                return GSL_EDIVERGE # v is decreasing
            return CommonRHS(dvdt_T1, y)

        def CommonRHS(dvdt, y):
            dydt=np.zeros(9)
            rfrak_frame=np.zeros(3)
            rfrak_frame[0] = y[5]
            rfrak_frame[1] = y[6]
            rfrak_frame[2] = y[7]
            rfrakdot_frame = Quaternions.FrameFromAngularVelocity_Integrand(rfrak_frame, OmegaVec().vec())
            dydt[0] = dvdt
            if(EvolveSpin1):
                Quaternions.FrameFromAngularVelocity_2D_Integrand(y[1], y[2],\
                    (S_chi1.inverse()*OmegaVec_chiVec_1()*S_chi1).vec(),dydt[1], dydt[2])
            else:
                dydt[1] = 0.0
                dydt[2] = 0.0
            if(EvolveSpin2):
                Quaternions.FrameFromAngularVelocity_2D_Integrand(y[3], y[4],\
                    (S_chi2.inverse()*OmegaVec_chiVec_2()*S_chi2).vec(),dydt[3], dydt[4])
            else:
                dydt[3] = 0.0
                dydt[4] = 0.0
            dydt[5] = rfrakdot_frame[0]
            dydt[6] = rfrakdot_frame[1]
            dydt[7] = rfrakdot_frame[2]
            dydt[8] = v*v*v/M

            return dydt
        
        
        y=solve_ivp(TaylorT1, [0,20000], [0.01**(1/3),rfrak_chi1_x,\
            rfrak_chi1_y,rfrak_chi2_x,rfrak_chi2_y,rfrak_frame_x,\
            rfrak_frame_y,rfrak_frame_z,Phi], method='DOP853',\
            t_eval=np.arange(0,20000,100), dense_output=True)
            

#        def TaylorT1_phi(t,phi):
#            dphidt_T1=v_of_t.y[0][v_of_t.t>=t][0]**3.0/M
#            return dphidt_T1
#        phi_of_t=solve_ivp(TaylorT1_phi, [0, v_of_t.t[-1]], [0.0], method='RK45', t_eval=v_of_t.t)
        
        #Calculate waveform modes

        def WaveformModes():
            I=1j
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
            chiVec1 = Quaternions.Quaternion(0., chi1_n_k, chi1_lambda_k, chi1_ell_k);
            chiVec2 = Quaternions.Quaternion(0., chi2_n_k, chi2_lambda_k, chi2_ell_k);

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

            Modes=np.zeros(77,len(y.t))
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
