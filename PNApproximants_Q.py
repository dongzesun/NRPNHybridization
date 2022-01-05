# File produced automatically by OrbitalEvolutionCodeGen_Q.ipynb
from scipy.integrate import solve_ivp
import numpy as np
from numpy import conjugate, dot, exp, log, sqrt, pi
from numpy import euler_gamma as EulerGamma
import quaternion

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
        EvolveSpin1=dot(S_chi1,S_chi1).norm()>1e-12
        EvolveSpin2=dot(S_chi2,S_chi2).norm()>1e-12

        def Recalculate(t, y):
            v = y[0]
            rfrak_chi1_x = y[1]
            rfrak_chi1_y = y[2]
            rfrak_chi2_x = y[3]
            rfrak_chi2_y = y[4]
            rfrak_frame_x = y[5]
            rfrak_frame_y = y[6]
            rfrak_frame_z = y[7]
            Phi = y[8]
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

        def FrameFromAngularVelocity_2D_Integrand(rfrak_x, rfrak_y, Omega):
            rfrakMag = np.sqrt(rfrak_x*rfrak_x+rfrak_y*rfrak_y)
            rfrakDot_x = Omega[0]/2.0
            rfrakDot_y = Omega[1]/2.0
            if abs(np.sin(rfrakMag)) > 1e-12 and abs(np.cos(rfrakMag)) > 1e-12:
                omega_v = (Omega[0]*(-rfrak_y/rfrakMag)+Omega[1]*(rfrak_x/rfrakMag))*np.tan(rfrakMag)-Omega[2]
                Omega[0] += -omega_v*np.sin(2*rfrakMag)*(-rfrak_y/rfrakMag)
                Omega[1] += -omega_v*np.sin(2*rfrakMag)*(rfrak_x/rfrakMag)
                Omega[2] +=  omega_v*np.cos(2*rfrakMag)
                dotTerm = (rfrak_x*Omega[0]+rfrak_y*Omega[1])/(rfrakMag*rfrakMag)
                cotTerm = rfrakMag/(2*np.tan(rfrakMag))
                rfrakDot_x = (Omega[0] - rfrak_x*dotTerm)*cotTerm + rfrak_x*dotTerm/2. - 0.5*Omega[2]*rfrak_y
                rfrakDot_y = (Omega[1] - rfrak_y*dotTerm)*cotTerm + rfrak_y*dotTerm/2. + 0.5*Omega[2]*rfrak_x
            return rfrakDot_x, rfrakDot_y
            

        #Evolve PN
        def TaylorT1(t, y):
            v = y[0]
            rfrak_chi1_x = y[1]
            rfrak_chi1_y = y[2]
            rfrak_chi2_x = y[3]
            rfrak_chi2_y = y[4]
            rfrak_frame_x = y[5]
            rfrak_frame_y = y[6]
            rfrak_frame_z = y[7]
            Phi = y[8]
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

            if v>=1.0: 
                raise ValueError("Beyond domain of PN validity")
            Flux = Fcal_0*Fcal_coeff
            dEdv = -E_0*M*nu*v
            Absorption = 0
            dvdt_T1 = (-Absorption - Flux)/dEdv
            try:
                dvdt_T1=quaternion.as_float_array(dvdt_T1)[0]
            except:
                pass
            if dvdt_T1<1.0e-12:
                raise ValueError("v is decreasing")
            return CommonRHS(dvdt_T1, y)

        def CommonRHS(dvdt, y):
            v = y[0]
            rfrak_chi1_x = y[1]
            rfrak_chi1_y = y[2]
            rfrak_chi2_x = y[3]
            rfrak_chi2_y = y[4]
            rfrak_frame_x = y[5]
            rfrak_frame_y = y[6]
            rfrak_frame_z = y[7]
            Phi = y[8]
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

            dydt=np.zeros(9)
            rfrak_frame=np.zeros(3)
            rfrak_frame[0] = y[5]
            rfrak_frame[1] = y[6]
            rfrak_frame[2] = y[7]
            rfrakdot_frame = quaternion.quaternion_time_series.frame_from_angular_velocity_integrand(rfrak_frame,\
                quaternion.as_float_array(OmegaVec())[1:])
            dydt[0] = dvdt
            if(EvolveSpin1):
                dydt[1], dydt[2]= FrameFromAngularVelocity_2D_Integrand(y[1], y[2],\
                    quaternion.as_float_array(S_chi1.inverse()*OmegaVec_chiVec_1()*S_chi1))
            else:
                dydt[1] = 0.0
                dydt[2] = 0.0
            if(EvolveSpin2):
                dydt[3], dydt[4] = FrameFromAngularVelocity_2D_Integrand(y[3], y[4],\
                    quaternion.as_float_array(S_chi2.inverse()*OmegaVec_chiVec_2()*S_chi2))
            else:
                dydt[3] = 0.0
                dydt[4] = 0.0
            dydt[5] = rfrakdot_frame[0]
            dydt[6] = rfrakdot_frame[1]
            dydt[7] = rfrakdot_frame[2]
            dydt[8] = v*v*v/M

            return dydt
        
        
        yy=solve_ivp(TaylorT1, [0,114000], [v,rfrak_chi1_x,\
            rfrak_chi1_y,rfrak_chi2_x,rfrak_chi2_y,rfrak_frame_x,\
            rfrak_frame_y,rfrak_frame_z,Phi], method='DOP853',\
            t_eval=np.arange(0,114000,100), dense_output=True)
        return yy    

#        def TaylorT1_phi(t,phi):
#            dphidt_T1=v_of_t.y[0][v_of_t.t>=t][0]**3.0/M
#            return dphidt_T1
#        phi_of_t=solve_ivp(TaylorT1_phi, [0, v_of_t.t[-1]], [0.0], method='RK45', t_eval=v_of_t.t)
        
        #Calculate waveform modes

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
        EvolveSpin1=dot(S_chi1,S_chi1).norm()>1e-12
        EvolveSpin2=dot(S_chi2,S_chi2).norm()>1e-12

        def Recalculate(t, y):
            v = y[0]
            rfrak_chi1_x = y[1]
            rfrak_chi1_y = y[2]
            rfrak_chi2_x = y[3]
            rfrak_chi2_y = y[4]
            rfrak_frame_x = y[5]
            rfrak_frame_y = y[6]
            rfrak_frame_z = y[7]
            Phi = y[8]
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

        def FrameFromAngularVelocity_2D_Integrand(rfrak_x, rfrak_y, Omega):
            rfrakMag = np.sqrt(rfrak_x*rfrak_x+rfrak_y*rfrak_y)
            rfrakDot_x = Omega[0]/2.0
            rfrakDot_y = Omega[1]/2.0
            if abs(np.sin(rfrakMag)) > 1e-12 and abs(np.cos(rfrakMag)) > 1e-12:
                omega_v = (Omega[0]*(-rfrak_y/rfrakMag)+Omega[1]*(rfrak_x/rfrakMag))*np.tan(rfrakMag)-Omega[2]
                Omega[0] += -omega_v*np.sin(2*rfrakMag)*(-rfrak_y/rfrakMag)
                Omega[1] += -omega_v*np.sin(2*rfrakMag)*(rfrak_x/rfrakMag)
                Omega[2] +=  omega_v*np.cos(2*rfrakMag)
                dotTerm = (rfrak_x*Omega[0]+rfrak_y*Omega[1])/(rfrakMag*rfrakMag)
                cotTerm = rfrakMag/(2*np.tan(rfrakMag))
                rfrakDot_x = (Omega[0] - rfrak_x*dotTerm)*cotTerm + rfrak_x*dotTerm/2. - 0.5*Omega[2]*rfrak_y
                rfrakDot_y = (Omega[1] - rfrak_y*dotTerm)*cotTerm + rfrak_y*dotTerm/2. + 0.5*Omega[2]*rfrak_x
            return rfrakDot_x, rfrakDot_y
            

        #Evolve PN
        def TaylorT1(t, y):
            v = y[0]
            rfrak_chi1_x = y[1]
            rfrak_chi1_y = y[2]
            rfrak_chi2_x = y[3]
            rfrak_chi2_y = y[4]
            rfrak_frame_x = y[5]
            rfrak_frame_y = y[6]
            rfrak_frame_z = y[7]
            Phi = y[8]
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

            if v>=1.0: 
                raise ValueError("Beyond domain of PN validity")
            Flux = Fcal_0*Fcal_coeff
            dEdv = -E_0*M*nu*v
            Absorption = 0
            dvdt_T1 = (-Absorption - Flux)/dEdv
            try:
                dvdt_T1=quaternion.as_float_array(dvdt_T1)[0]
            except:
                pass
            if dvdt_T1<1.0e-12:
                raise ValueError("v is decreasing")
            return CommonRHS(dvdt_T1, y)

        def CommonRHS(dvdt, y):
            v = y[0]
            rfrak_chi1_x = y[1]
            rfrak_chi1_y = y[2]
            rfrak_chi2_x = y[3]
            rfrak_chi2_y = y[4]
            rfrak_frame_x = y[5]
            rfrak_frame_y = y[6]
            rfrak_frame_z = y[7]
            Phi = y[8]
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

            dydt=np.zeros(9)
            rfrak_frame=np.zeros(3)
            rfrak_frame[0] = y[5]
            rfrak_frame[1] = y[6]
            rfrak_frame[2] = y[7]
            rfrakdot_frame = quaternion.quaternion_time_series.frame_from_angular_velocity_integrand(rfrak_frame,\
                quaternion.as_float_array(OmegaVec())[1:])
            dydt[0] = dvdt
            if(EvolveSpin1):
                dydt[1], dydt[2]= FrameFromAngularVelocity_2D_Integrand(y[1], y[2],\
                    quaternion.as_float_array(S_chi1.inverse()*OmegaVec_chiVec_1()*S_chi1))
            else:
                dydt[1] = 0.0
                dydt[2] = 0.0
            if(EvolveSpin2):
                dydt[3], dydt[4] = FrameFromAngularVelocity_2D_Integrand(y[3], y[4],\
                    quaternion.as_float_array(S_chi2.inverse()*OmegaVec_chiVec_2()*S_chi2))
            else:
                dydt[3] = 0.0
                dydt[4] = 0.0
            dydt[5] = rfrakdot_frame[0]
            dydt[6] = rfrakdot_frame[1]
            dydt[7] = rfrakdot_frame[2]
            dydt[8] = v*v*v/M

            return dydt
        
        
        yy=solve_ivp(TaylorT1, [0,114000], [v,rfrak_chi1_x,\
            rfrak_chi1_y,rfrak_chi2_x,rfrak_chi2_y,rfrak_frame_x,\
            rfrak_frame_y,rfrak_frame_z,Phi], method='DOP853',\
            t_eval=np.arange(0,114000,100), dense_output=True)
        return yy    

#        def TaylorT1_phi(t,phi):
#            dphidt_T1=v_of_t.y[0][v_of_t.t>=t][0]**3.0/M
#            return dphidt_T1
#        phi_of_t=solve_ivp(TaylorT1_phi, [0, v_of_t.t[-1]], [0.0], method='RK45', t_eval=v_of_t.t)
        
        #Calculate waveform modes

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
        EvolveSpin1=dot(S_chi1,S_chi1).norm()>1e-12
        EvolveSpin2=dot(S_chi2,S_chi2).norm()>1e-12

        def Recalculate(t, y):
            v = y[0]
            rfrak_chi1_x = y[1]
            rfrak_chi1_y = y[2]
            rfrak_chi2_x = y[3]
            rfrak_chi2_y = y[4]
            rfrak_frame_x = y[5]
            rfrak_frame_y = y[6]
            rfrak_frame_z = y[7]
            Phi = y[8]
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

        def FrameFromAngularVelocity_2D_Integrand(rfrak_x, rfrak_y, Omega):
            rfrakMag = np.sqrt(rfrak_x*rfrak_x+rfrak_y*rfrak_y)
            rfrakDot_x = Omega[0]/2.0
            rfrakDot_y = Omega[1]/2.0
            if abs(np.sin(rfrakMag)) > 1e-12 and abs(np.cos(rfrakMag)) > 1e-12:
                omega_v = (Omega[0]*(-rfrak_y/rfrakMag)+Omega[1]*(rfrak_x/rfrakMag))*np.tan(rfrakMag)-Omega[2]
                Omega[0] += -omega_v*np.sin(2*rfrakMag)*(-rfrak_y/rfrakMag)
                Omega[1] += -omega_v*np.sin(2*rfrakMag)*(rfrak_x/rfrakMag)
                Omega[2] +=  omega_v*np.cos(2*rfrakMag)
                dotTerm = (rfrak_x*Omega[0]+rfrak_y*Omega[1])/(rfrakMag*rfrakMag)
                cotTerm = rfrakMag/(2*np.tan(rfrakMag))
                rfrakDot_x = (Omega[0] - rfrak_x*dotTerm)*cotTerm + rfrak_x*dotTerm/2. - 0.5*Omega[2]*rfrak_y
                rfrakDot_y = (Omega[1] - rfrak_y*dotTerm)*cotTerm + rfrak_y*dotTerm/2. + 0.5*Omega[2]*rfrak_x
            return rfrakDot_x, rfrakDot_y
            

        #Evolve PN
        def TaylorT1(t, y):
            v = y[0]
            rfrak_chi1_x = y[1]
            rfrak_chi1_y = y[2]
            rfrak_chi2_x = y[3]
            rfrak_chi2_y = y[4]
            rfrak_frame_x = y[5]
            rfrak_frame_y = y[6]
            rfrak_frame_z = y[7]
            Phi = y[8]
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

            if v>=1.0: 
                raise ValueError("Beyond domain of PN validity")
            Flux = Fcal_coeff*(Fcal_0 + Fcal_2*v**2)
            dEdv = -M*nu*v*(E_0 + 2.0*E_2*v**2)
            Absorption = 0
            dvdt_T1 = (-Absorption - Flux)/dEdv
            try:
                dvdt_T1=quaternion.as_float_array(dvdt_T1)[0]
            except:
                pass
            if dvdt_T1<1.0e-12:
                raise ValueError("v is decreasing")
            return CommonRHS(dvdt_T1, y)

        def CommonRHS(dvdt, y):
            v = y[0]
            rfrak_chi1_x = y[1]
            rfrak_chi1_y = y[2]
            rfrak_chi2_x = y[3]
            rfrak_chi2_y = y[4]
            rfrak_frame_x = y[5]
            rfrak_frame_y = y[6]
            rfrak_frame_z = y[7]
            Phi = y[8]
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
                a_ell_2 = S_n*(-9.66666666666667*nu - 10.0) + Sigma_n*delta*(-4.5*nu - 6.0)
                a_ell_0 = 7.0*S_n + 3.0*Sigma_n*delta
                gamma_PN_0 = 1.00000000000000
                gamma_PN_2 = 1.0 - 0.333333333333333*nu
                return ellHat*v**3/M + nHat*v**6*(a_ell_0 + a_ell_2*v**2)*(gamma_PN_0 + gamma_PN_2*v**2)/M**3

            def OrbitalAngularMomentum():
                L_0 = ellHat
                L_coeff = M**2*nu/v
                L_2 = ellHat*(0.166666666666667*nu + 1.5)
                return L_coeff*(L_0 + L_2*v**2)

            dydt=np.zeros(9)
            rfrak_frame=np.zeros(3)
            rfrak_frame[0] = y[5]
            rfrak_frame[1] = y[6]
            rfrak_frame[2] = y[7]
            rfrakdot_frame = quaternion.quaternion_time_series.frame_from_angular_velocity_integrand(rfrak_frame,\
                quaternion.as_float_array(OmegaVec())[1:])
            dydt[0] = dvdt
            if(EvolveSpin1):
                dydt[1], dydt[2]= FrameFromAngularVelocity_2D_Integrand(y[1], y[2],\
                    quaternion.as_float_array(S_chi1.inverse()*OmegaVec_chiVec_1()*S_chi1))
            else:
                dydt[1] = 0.0
                dydt[2] = 0.0
            if(EvolveSpin2):
                dydt[3], dydt[4] = FrameFromAngularVelocity_2D_Integrand(y[3], y[4],\
                    quaternion.as_float_array(S_chi2.inverse()*OmegaVec_chiVec_2()*S_chi2))
            else:
                dydt[3] = 0.0
                dydt[4] = 0.0
            dydt[5] = rfrakdot_frame[0]
            dydt[6] = rfrakdot_frame[1]
            dydt[7] = rfrakdot_frame[2]
            dydt[8] = v*v*v/M

            return dydt
        
        
        yy=solve_ivp(TaylorT1, [0,114000], [v,rfrak_chi1_x,\
            rfrak_chi1_y,rfrak_chi2_x,rfrak_chi2_y,rfrak_frame_x,\
            rfrak_frame_y,rfrak_frame_z,Phi], method='DOP853',\
            t_eval=np.arange(0,114000,100), dense_output=True)
        return yy    

#        def TaylorT1_phi(t,phi):
#            dphidt_T1=v_of_t.y[0][v_of_t.t>=t][0]**3.0/M
#            return dphidt_T1
#        phi_of_t=solve_ivp(TaylorT1_phi, [0, v_of_t.t[-1]], [0.0], method='RK45', t_eval=v_of_t.t)
        
        #Calculate waveform modes

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
        EvolveSpin1=dot(S_chi1,S_chi1).norm()>1e-12
        EvolveSpin2=dot(S_chi2,S_chi2).norm()>1e-12

        def Recalculate(t, y):
            v = y[0]
            rfrak_chi1_x = y[1]
            rfrak_chi1_y = y[2]
            rfrak_chi2_x = y[3]
            rfrak_chi2_y = y[4]
            rfrak_frame_x = y[5]
            rfrak_frame_y = y[6]
            rfrak_frame_z = y[7]
            Phi = y[8]
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

        def FrameFromAngularVelocity_2D_Integrand(rfrak_x, rfrak_y, Omega):
            rfrakMag = np.sqrt(rfrak_x*rfrak_x+rfrak_y*rfrak_y)
            rfrakDot_x = Omega[0]/2.0
            rfrakDot_y = Omega[1]/2.0
            if abs(np.sin(rfrakMag)) > 1e-12 and abs(np.cos(rfrakMag)) > 1e-12:
                omega_v = (Omega[0]*(-rfrak_y/rfrakMag)+Omega[1]*(rfrak_x/rfrakMag))*np.tan(rfrakMag)-Omega[2]
                Omega[0] += -omega_v*np.sin(2*rfrakMag)*(-rfrak_y/rfrakMag)
                Omega[1] += -omega_v*np.sin(2*rfrakMag)*(rfrak_x/rfrakMag)
                Omega[2] +=  omega_v*np.cos(2*rfrakMag)
                dotTerm = (rfrak_x*Omega[0]+rfrak_y*Omega[1])/(rfrakMag*rfrakMag)
                cotTerm = rfrakMag/(2*np.tan(rfrakMag))
                rfrakDot_x = (Omega[0] - rfrak_x*dotTerm)*cotTerm + rfrak_x*dotTerm/2. - 0.5*Omega[2]*rfrak_y
                rfrakDot_y = (Omega[1] - rfrak_y*dotTerm)*cotTerm + rfrak_y*dotTerm/2. + 0.5*Omega[2]*rfrak_x
            return rfrakDot_x, rfrakDot_y
            

        #Evolve PN
        def TaylorT1(t, y):
            v = y[0]
            rfrak_chi1_x = y[1]
            rfrak_chi1_y = y[2]
            rfrak_chi2_x = y[3]
            rfrak_chi2_y = y[4]
            rfrak_frame_x = y[5]
            rfrak_frame_y = y[6]
            rfrak_frame_z = y[7]
            Phi = y[8]
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

            if v>=1.0: 
                raise ValueError("Beyond domain of PN validity")
            Flux = Fcal_coeff*(Fcal_0 + v**2*(Fcal_2 + v*(Fcal_3 + Fcal_SO_3)))
            dEdv = -0.5*M*nu*v*(2.0*E_0 + v**2*(4.0*E_2 + 5.0*E_SO_3*v))
            Absorption = 0
            dvdt_T1 = (-Absorption - Flux)/dEdv
            try:
                dvdt_T1=quaternion.as_float_array(dvdt_T1)[0]
            except:
                pass
            if dvdt_T1<1.0e-12:
                raise ValueError("v is decreasing")
            return CommonRHS(dvdt_T1, y)

        def CommonRHS(dvdt, y):
            v = y[0]
            rfrak_chi1_x = y[1]
            rfrak_chi1_y = y[2]
            rfrak_chi2_x = y[3]
            rfrak_chi2_y = y[4]
            rfrak_frame_x = y[5]
            rfrak_frame_y = y[6]
            rfrak_frame_z = y[7]
            Phi = y[8]
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
                a_ell_2 = S_n*(-9.66666666666667*nu - 10.0) + Sigma_n*delta*(-4.5*nu - 6.0)
                a_ell_0 = 7.0*S_n + 3.0*Sigma_n*delta
                gamma_PN_3 = (1.66666666666667*S_ell + Sigma_ell*delta)/M**2
                gamma_PN_0 = 1.00000000000000
                gamma_PN_2 = 1.0 - 0.333333333333333*nu
                return ellHat*v**3/M + nHat*v**6*(a_ell_0 + a_ell_2*v**2)*(gamma_PN_0 + v**2*(gamma_PN_2 + gamma_PN_3*v))/M**3

            def OrbitalAngularMomentum():
                L_SO_3 = 0.166666666666667*ellHat*(-35.0*S_ell - 15.0*Sigma_ell*delta)/M**2 + lambdaHat*(-3.0*S_lambda - Sigma_lambda*delta)/M**2 + 0.5*nHat*(S_n + Sigma_n*delta)/M**2
                L_0 = ellHat
                L_coeff = M**2*nu/v
                L_2 = ellHat*(0.166666666666667*nu + 1.5)
                return L_coeff*(L_0 + v**2*(L_2 + L_SO_3*v))

            dydt=np.zeros(9)
            rfrak_frame=np.zeros(3)
            rfrak_frame[0] = y[5]
            rfrak_frame[1] = y[6]
            rfrak_frame[2] = y[7]
            rfrakdot_frame = quaternion.quaternion_time_series.frame_from_angular_velocity_integrand(rfrak_frame,\
                quaternion.as_float_array(OmegaVec())[1:])
            dydt[0] = dvdt
            if(EvolveSpin1):
                dydt[1], dydt[2]= FrameFromAngularVelocity_2D_Integrand(y[1], y[2],\
                    quaternion.as_float_array(S_chi1.inverse()*OmegaVec_chiVec_1()*S_chi1))
            else:
                dydt[1] = 0.0
                dydt[2] = 0.0
            if(EvolveSpin2):
                dydt[3], dydt[4] = FrameFromAngularVelocity_2D_Integrand(y[3], y[4],\
                    quaternion.as_float_array(S_chi2.inverse()*OmegaVec_chiVec_2()*S_chi2))
            else:
                dydt[3] = 0.0
                dydt[4] = 0.0
            dydt[5] = rfrakdot_frame[0]
            dydt[6] = rfrakdot_frame[1]
            dydt[7] = rfrakdot_frame[2]
            dydt[8] = v*v*v/M

            return dydt
        
        
        yy=solve_ivp(TaylorT1, [0,114000], [v,rfrak_chi1_x,\
            rfrak_chi1_y,rfrak_chi2_x,rfrak_chi2_y,rfrak_frame_x,\
            rfrak_frame_y,rfrak_frame_z,Phi], method='DOP853',\
            t_eval=np.arange(0,114000,100), dense_output=True)
        return yy    

#        def TaylorT1_phi(t,phi):
#            dphidt_T1=v_of_t.y[0][v_of_t.t>=t][0]**3.0/M
#            return dphidt_T1
#        phi_of_t=solve_ivp(TaylorT1_phi, [0, v_of_t.t[-1]], [0.0], method='RK45', t_eval=v_of_t.t)
        
        #Calculate waveform modes

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
        EvolveSpin1=dot(S_chi1,S_chi1).norm()>1e-12
        EvolveSpin2=dot(S_chi2,S_chi2).norm()>1e-12

        def Recalculate(t, y):
            v = y[0]
            rfrak_chi1_x = y[1]
            rfrak_chi1_y = y[2]
            rfrak_chi2_x = y[3]
            rfrak_chi2_y = y[4]
            rfrak_frame_x = y[5]
            rfrak_frame_y = y[6]
            rfrak_frame_z = y[7]
            Phi = y[8]
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

        def FrameFromAngularVelocity_2D_Integrand(rfrak_x, rfrak_y, Omega):
            rfrakMag = np.sqrt(rfrak_x*rfrak_x+rfrak_y*rfrak_y)
            rfrakDot_x = Omega[0]/2.0
            rfrakDot_y = Omega[1]/2.0
            if abs(np.sin(rfrakMag)) > 1e-12 and abs(np.cos(rfrakMag)) > 1e-12:
                omega_v = (Omega[0]*(-rfrak_y/rfrakMag)+Omega[1]*(rfrak_x/rfrakMag))*np.tan(rfrakMag)-Omega[2]
                Omega[0] += -omega_v*np.sin(2*rfrakMag)*(-rfrak_y/rfrakMag)
                Omega[1] += -omega_v*np.sin(2*rfrakMag)*(rfrak_x/rfrakMag)
                Omega[2] +=  omega_v*np.cos(2*rfrakMag)
                dotTerm = (rfrak_x*Omega[0]+rfrak_y*Omega[1])/(rfrakMag*rfrakMag)
                cotTerm = rfrakMag/(2*np.tan(rfrakMag))
                rfrakDot_x = (Omega[0] - rfrak_x*dotTerm)*cotTerm + rfrak_x*dotTerm/2. - 0.5*Omega[2]*rfrak_y
                rfrakDot_y = (Omega[1] - rfrak_y*dotTerm)*cotTerm + rfrak_y*dotTerm/2. + 0.5*Omega[2]*rfrak_x
            return rfrakDot_x, rfrakDot_y
            

        #Evolve PN
        def TaylorT1(t, y):
            v = y[0]
            rfrak_chi1_x = y[1]
            rfrak_chi1_y = y[2]
            rfrak_chi2_x = y[3]
            rfrak_chi2_y = y[4]
            rfrak_frame_x = y[5]
            rfrak_frame_y = y[6]
            rfrak_frame_z = y[7]
            Phi = y[8]
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

            if v>=1.0: 
                raise ValueError("Beyond domain of PN validity")
            Flux = Fcal_coeff*(Fcal_0 + v**2*(Fcal_2 + v*(Fcal_3 + Fcal_SO_3 + v*(Fcal_4 + Fcal_SQ_4))))
            dEdv = -0.5*M*nu*v*(2.0*E_0 + v**2*(4.0*E_2 + v*(5.0*E_SO_3 + 6.0*v*(E_4 + E_SQ_4))))
            Absorption = 0
            dvdt_T1 = (-Absorption - Flux)/dEdv
            try:
                dvdt_T1=quaternion.as_float_array(dvdt_T1)[0]
            except:
                pass
            if dvdt_T1<1.0e-12:
                raise ValueError("v is decreasing")
            return CommonRHS(dvdt_T1, y)

        def CommonRHS(dvdt, y):
            v = y[0]
            rfrak_chi1_x = y[1]
            rfrak_chi1_y = y[2]
            rfrak_chi2_x = y[3]
            rfrak_chi2_y = y[4]
            rfrak_frame_x = y[5]
            rfrak_frame_y = y[6]
            rfrak_frame_z = y[7]
            Phi = y[8]
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
                a_ell_4 = S_n*(5.77777777777778*nu**2 + 14.75*nu + 1.5) + Sigma_n*delta*(2.83333333333333*nu**2 + 9.125*nu + 1.5)
                a_ell_2 = S_n*(-9.66666666666667*nu - 10.0) + Sigma_n*delta*(-4.5*nu - 6.0)
                a_ell_0 = 7.0*S_n + 3.0*Sigma_n*delta
                gamma_PN_4 = 1.0 - 5.41666666666667*nu
                gamma_PN_3 = (1.66666666666667*S_ell + Sigma_ell*delta)/M**2
                gamma_PN_0 = 1.00000000000000
                gamma_PN_2 = 1.0 - 0.333333333333333*nu
                return ellHat*v**3/M + nHat*v**6*(a_ell_0 + v**2*(a_ell_2 + a_ell_4*v**2))*(gamma_PN_0 + v**2*(gamma_PN_2 + v*(gamma_PN_3 + gamma_PN_4*v)))/M**3

            def OrbitalAngularMomentum():
                L_SO_3 = 0.166666666666667*ellHat*(-35.0*S_ell - 15.0*Sigma_ell*delta)/M**2 + lambdaHat*(-3.0*S_lambda - Sigma_lambda*delta)/M**2 + 0.5*nHat*(S_n + Sigma_n*delta)/M**2
                L_0 = ellHat
                L_coeff = M**2*nu/v
                L_4 = ellHat*(0.0416666666666667*nu**2 - 2.375*nu + 3.375)
                L_2 = ellHat*(0.166666666666667*nu + 1.5)
                return L_coeff*(L_0 + v**2*(L_2 + v*(L_4*v + L_SO_3)))

            dydt=np.zeros(9)
            rfrak_frame=np.zeros(3)
            rfrak_frame[0] = y[5]
            rfrak_frame[1] = y[6]
            rfrak_frame[2] = y[7]
            rfrakdot_frame = quaternion.quaternion_time_series.frame_from_angular_velocity_integrand(rfrak_frame,\
                quaternion.as_float_array(OmegaVec())[1:])
            dydt[0] = dvdt
            if(EvolveSpin1):
                dydt[1], dydt[2]= FrameFromAngularVelocity_2D_Integrand(y[1], y[2],\
                    quaternion.as_float_array(S_chi1.inverse()*OmegaVec_chiVec_1()*S_chi1))
            else:
                dydt[1] = 0.0
                dydt[2] = 0.0
            if(EvolveSpin2):
                dydt[3], dydt[4] = FrameFromAngularVelocity_2D_Integrand(y[3], y[4],\
                    quaternion.as_float_array(S_chi2.inverse()*OmegaVec_chiVec_2()*S_chi2))
            else:
                dydt[3] = 0.0
                dydt[4] = 0.0
            dydt[5] = rfrakdot_frame[0]
            dydt[6] = rfrakdot_frame[1]
            dydt[7] = rfrakdot_frame[2]
            dydt[8] = v*v*v/M

            return dydt
        
        
        yy=solve_ivp(TaylorT1, [0,114000], [v,rfrak_chi1_x,\
            rfrak_chi1_y,rfrak_chi2_x,rfrak_chi2_y,rfrak_frame_x,\
            rfrak_frame_y,rfrak_frame_z,Phi], method='DOP853',\
            t_eval=np.arange(0,114000,100), dense_output=True)
        return yy    

#        def TaylorT1_phi(t,phi):
#            dphidt_T1=v_of_t.y[0][v_of_t.t>=t][0]**3.0/M
#            return dphidt_T1
#        phi_of_t=solve_ivp(TaylorT1_phi, [0, v_of_t.t[-1]], [0.0], method='RK45', t_eval=v_of_t.t)
        
        #Calculate waveform modes

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
        EvolveSpin1=dot(S_chi1,S_chi1).norm()>1e-12
        EvolveSpin2=dot(S_chi2,S_chi2).norm()>1e-12

        def Recalculate(t, y):
            v = y[0]
            rfrak_chi1_x = y[1]
            rfrak_chi1_y = y[2]
            rfrak_chi2_x = y[3]
            rfrak_chi2_y = y[4]
            rfrak_frame_x = y[5]
            rfrak_frame_y = y[6]
            rfrak_frame_z = y[7]
            Phi = y[8]
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

        def FrameFromAngularVelocity_2D_Integrand(rfrak_x, rfrak_y, Omega):
            rfrakMag = np.sqrt(rfrak_x*rfrak_x+rfrak_y*rfrak_y)
            rfrakDot_x = Omega[0]/2.0
            rfrakDot_y = Omega[1]/2.0
            if abs(np.sin(rfrakMag)) > 1e-12 and abs(np.cos(rfrakMag)) > 1e-12:
                omega_v = (Omega[0]*(-rfrak_y/rfrakMag)+Omega[1]*(rfrak_x/rfrakMag))*np.tan(rfrakMag)-Omega[2]
                Omega[0] += -omega_v*np.sin(2*rfrakMag)*(-rfrak_y/rfrakMag)
                Omega[1] += -omega_v*np.sin(2*rfrakMag)*(rfrak_x/rfrakMag)
                Omega[2] +=  omega_v*np.cos(2*rfrakMag)
                dotTerm = (rfrak_x*Omega[0]+rfrak_y*Omega[1])/(rfrakMag*rfrakMag)
                cotTerm = rfrakMag/(2*np.tan(rfrakMag))
                rfrakDot_x = (Omega[0] - rfrak_x*dotTerm)*cotTerm + rfrak_x*dotTerm/2. - 0.5*Omega[2]*rfrak_y
                rfrakDot_y = (Omega[1] - rfrak_y*dotTerm)*cotTerm + rfrak_y*dotTerm/2. + 0.5*Omega[2]*rfrak_x
            return rfrakDot_x, rfrakDot_y
            

        #Evolve PN
        def TaylorT1(t, y):
            v = y[0]
            rfrak_chi1_x = y[1]
            rfrak_chi1_y = y[2]
            rfrak_chi2_x = y[3]
            rfrak_chi2_y = y[4]
            rfrak_frame_x = y[5]
            rfrak_frame_y = y[6]
            rfrak_frame_z = y[7]
            Phi = y[8]
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

            if v>=1.0: 
                raise ValueError("Beyond domain of PN validity")
            Flux = Fcal_coeff*(Fcal_0 + v**2*(Fcal_2 + v*(Fcal_3 + Fcal_SO_3 + v*(Fcal_4 + Fcal_SQ_4 + v*(Fcal_5 + Fcal_SO_5)))))
            dEdv = -0.5*M*nu*v*(2.0*E_0 + v**2*(4.0*E_2 + v*(5.0*E_SO_3 + v*(6.0*E_4 + 7.0*E_SO_5*v + 6.0*E_SQ_4))))
            Absorption = 0
            dvdt_T1 = (-Absorption - Flux)/dEdv
            try:
                dvdt_T1=quaternion.as_float_array(dvdt_T1)[0]
            except:
                pass
            if dvdt_T1<1.0e-12:
                raise ValueError("v is decreasing")
            return CommonRHS(dvdt_T1, y)

        def CommonRHS(dvdt, y):
            v = y[0]
            rfrak_chi1_x = y[1]
            rfrak_chi1_y = y[2]
            rfrak_chi2_x = y[3]
            rfrak_chi2_y = y[4]
            rfrak_frame_x = y[5]
            rfrak_frame_y = y[6]
            rfrak_frame_z = y[7]
            Phi = y[8]
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
                a_ell_4 = S_n*(5.77777777777778*nu**2 + 14.75*nu + 1.5) + Sigma_n*delta*(2.83333333333333*nu**2 + 9.125*nu + 1.5)
                gamma_PN_5 = (S_ell*(0.888888888888889*nu + 3.33333333333333) + 2.0*Sigma_ell*delta)/M**2
                a_ell_2 = S_n*(-9.66666666666667*nu - 10.0) + Sigma_n*delta*(-4.5*nu - 6.0)
                a_ell_0 = 7.0*S_n + 3.0*Sigma_n*delta
                gamma_PN_4 = 1.0 - 5.41666666666667*nu
                gamma_PN_3 = (1.66666666666667*S_ell + Sigma_ell*delta)/M**2
                gamma_PN_0 = 1.00000000000000
                gamma_PN_2 = 1.0 - 0.333333333333333*nu
                return ellHat*v**3/M + nHat*v**6*(a_ell_0 + v**2*(a_ell_2 + a_ell_4*v**2))*(gamma_PN_0 + v**2*(gamma_PN_2 + v*(gamma_PN_3 + v*(gamma_PN_4 + gamma_PN_5*v))))/M**3

            def OrbitalAngularMomentum():
                L_SO_3 = 0.166666666666667*ellHat*(-35.0*S_ell - 15.0*Sigma_ell*delta)/M**2 + lambdaHat*(-3.0*S_lambda - Sigma_lambda*delta)/M**2 + 0.5*nHat*(S_n + Sigma_n*delta)/M**2
                L_0 = ellHat
                L_coeff = M**2*nu/v
                L_4 = ellHat*(0.0416666666666667*nu**2 - 2.375*nu + 3.375)
                L_2 = ellHat*(0.166666666666667*nu + 1.5)
                L_SO_5 = 0.0138888888888889*ellHat*(427.0*S_ell*nu - 693.0*S_ell + 210.0*Sigma_ell*delta*nu - 189.0*Sigma_ell*delta)/M**2 + 0.166666666666667*lambdaHat*(18.0*S_lambda*nu - 21.0*S_lambda + 8.0*Sigma_lambda*delta*nu - 3.0*Sigma_lambda*delta)/M**2 + 0.0416666666666667*nHat*(-19.0*S_n*nu + 33.0*S_n - 10.0*Sigma_n*delta*nu + 33.0*Sigma_n*delta)/M**2
                return L_coeff*(L_0 + v**2*(L_2 + v*(L_SO_3 + v*(L_4 + L_SO_5*v))))

            dydt=np.zeros(9)
            rfrak_frame=np.zeros(3)
            rfrak_frame[0] = y[5]
            rfrak_frame[1] = y[6]
            rfrak_frame[2] = y[7]
            rfrakdot_frame = quaternion.quaternion_time_series.frame_from_angular_velocity_integrand(rfrak_frame,\
                quaternion.as_float_array(OmegaVec())[1:])
            dydt[0] = dvdt
            if(EvolveSpin1):
                dydt[1], dydt[2]= FrameFromAngularVelocity_2D_Integrand(y[1], y[2],\
                    quaternion.as_float_array(S_chi1.inverse()*OmegaVec_chiVec_1()*S_chi1))
            else:
                dydt[1] = 0.0
                dydt[2] = 0.0
            if(EvolveSpin2):
                dydt[3], dydt[4] = FrameFromAngularVelocity_2D_Integrand(y[3], y[4],\
                    quaternion.as_float_array(S_chi2.inverse()*OmegaVec_chiVec_2()*S_chi2))
            else:
                dydt[3] = 0.0
                dydt[4] = 0.0
            dydt[5] = rfrakdot_frame[0]
            dydt[6] = rfrakdot_frame[1]
            dydt[7] = rfrakdot_frame[2]
            dydt[8] = v*v*v/M

            return dydt
        
        
        yy=solve_ivp(TaylorT1, [0,114000], [v,rfrak_chi1_x,\
            rfrak_chi1_y,rfrak_chi2_x,rfrak_chi2_y,rfrak_frame_x,\
            rfrak_frame_y,rfrak_frame_z,Phi], method='DOP853',\
            t_eval=np.arange(0,114000,100), dense_output=True)
        return yy    

#        def TaylorT1_phi(t,phi):
#            dphidt_T1=v_of_t.y[0][v_of_t.t>=t][0]**3.0/M
#            return dphidt_T1
#        phi_of_t=solve_ivp(TaylorT1_phi, [0, v_of_t.t[-1]], [0.0], method='RK45', t_eval=v_of_t.t)
        
        #Calculate waveform modes

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
        EvolveSpin1=dot(S_chi1,S_chi1).norm()>1e-12
        EvolveSpin2=dot(S_chi2,S_chi2).norm()>1e-12

        def Recalculate(t, y):
            v = y[0]
            rfrak_chi1_x = y[1]
            rfrak_chi1_y = y[2]
            rfrak_chi2_x = y[3]
            rfrak_chi2_y = y[4]
            rfrak_frame_x = y[5]
            rfrak_frame_y = y[6]
            rfrak_frame_z = y[7]
            Phi = y[8]
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

        def FrameFromAngularVelocity_2D_Integrand(rfrak_x, rfrak_y, Omega):
            rfrakMag = np.sqrt(rfrak_x*rfrak_x+rfrak_y*rfrak_y)
            rfrakDot_x = Omega[0]/2.0
            rfrakDot_y = Omega[1]/2.0
            if abs(np.sin(rfrakMag)) > 1e-12 and abs(np.cos(rfrakMag)) > 1e-12:
                omega_v = (Omega[0]*(-rfrak_y/rfrakMag)+Omega[1]*(rfrak_x/rfrakMag))*np.tan(rfrakMag)-Omega[2]
                Omega[0] += -omega_v*np.sin(2*rfrakMag)*(-rfrak_y/rfrakMag)
                Omega[1] += -omega_v*np.sin(2*rfrakMag)*(rfrak_x/rfrakMag)
                Omega[2] +=  omega_v*np.cos(2*rfrakMag)
                dotTerm = (rfrak_x*Omega[0]+rfrak_y*Omega[1])/(rfrakMag*rfrakMag)
                cotTerm = rfrakMag/(2*np.tan(rfrakMag))
                rfrakDot_x = (Omega[0] - rfrak_x*dotTerm)*cotTerm + rfrak_x*dotTerm/2. - 0.5*Omega[2]*rfrak_y
                rfrakDot_y = (Omega[1] - rfrak_y*dotTerm)*cotTerm + rfrak_y*dotTerm/2. + 0.5*Omega[2]*rfrak_x
            return rfrakDot_x, rfrakDot_y
            

        #Evolve PN
        def TaylorT1(t, y):
            v = y[0]
            rfrak_chi1_x = y[1]
            rfrak_chi1_y = y[2]
            rfrak_chi2_x = y[3]
            rfrak_chi2_y = y[4]
            rfrak_frame_x = y[5]
            rfrak_frame_y = y[6]
            rfrak_frame_z = y[7]
            Phi = y[8]
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

            if v>=1.0: 
                raise ValueError("Beyond domain of PN validity")
            Flux = Fcal_coeff*(Fcal_0 + v**2*(Fcal_2 + v*(Fcal_3 + Fcal_SO_3 + v*(Fcal_4 + Fcal_SQ_4 + v*(Fcal_5 + Fcal_SO_5 + v*(Fcal_6 + Fcal_SO_6 + Fcal_lnv_6*logv))))))
            dEdv = -0.5*M*nu*v*(2.0*E_0 + v**2*(4.0*E_2 + v*(5.0*E_SO_3 + v*(6.0*E_4 + 6.0*E_SQ_4 + v*(8.0*E_6*v + 7.0*E_SO_5)))))
            Absorption = 0
            dvdt_T1 = (-Absorption - Flux)/dEdv
            try:
                dvdt_T1=quaternion.as_float_array(dvdt_T1)[0]
            except:
                pass
            if dvdt_T1<1.0e-12:
                raise ValueError("v is decreasing")
            return CommonRHS(dvdt_T1, y)

        def CommonRHS(dvdt, y):
            v = y[0]
            rfrak_chi1_x = y[1]
            rfrak_chi1_y = y[2]
            rfrak_chi2_x = y[3]
            rfrak_chi2_y = y[4]
            rfrak_frame_x = y[5]
            rfrak_frame_y = y[6]
            rfrak_frame_z = y[7]
            Phi = y[8]
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
                a_ell_4 = S_n*(5.77777777777778*nu**2 + 14.75*nu + 1.5) + Sigma_n*delta*(2.83333333333333*nu**2 + 9.125*nu + 1.5)
                gamma_PN_5 = (S_ell*(0.888888888888889*nu + 3.33333333333333) + 2.0*Sigma_ell*delta)/M**2
                a_ell_2 = S_n*(-9.66666666666667*nu - 10.0) + Sigma_n*delta*(-4.5*nu - 6.0)
                a_ell_0 = 7.0*S_n + 3.0*Sigma_n*delta
                gamma_PN_6 = 0.0123456790123457*nu**3 + 6.36111111111111*nu**2 - 2.98177812235564*nu + 1.0
                gamma_PN_4 = 1.0 - 5.41666666666667*nu
                gamma_PN_3 = (1.66666666666667*S_ell + Sigma_ell*delta)/M**2
                gamma_PN_0 = 1.00000000000000
                gamma_PN_2 = 1.0 - 0.333333333333333*nu
                return ellHat*v**3/M + nHat*v**6*(a_ell_0 + v**2*(a_ell_2 + a_ell_4*v**2))*(gamma_PN_0 + v**2*(gamma_PN_2 + v*(gamma_PN_3 + v*(gamma_PN_4 + v*(gamma_PN_5 + gamma_PN_6*v)))))/M**3

            def OrbitalAngularMomentum():
                L_SO_3 = 0.166666666666667*ellHat*(-35.0*S_ell - 15.0*Sigma_ell*delta)/M**2 + lambdaHat*(-3.0*S_lambda - Sigma_lambda*delta)/M**2 + 0.5*nHat*(S_n + Sigma_n*delta)/M**2
                L_0 = ellHat
                L_coeff = M**2*nu/v
                L_4 = ellHat*(0.0416666666666667*nu**2 - 2.375*nu + 3.375)
                L_2 = ellHat*(0.166666666666667*nu + 1.5)
                L_6 = ellHat*(0.00540123456790123*nu**3 + 1.29166666666667*nu**2 - 30.9797035925835*nu + 8.4375)
                L_SO_5 = 0.0138888888888889*ellHat*(427.0*S_ell*nu - 693.0*S_ell + 210.0*Sigma_ell*delta*nu - 189.0*Sigma_ell*delta)/M**2 + 0.166666666666667*lambdaHat*(18.0*S_lambda*nu - 21.0*S_lambda + 8.0*Sigma_lambda*delta*nu - 3.0*Sigma_lambda*delta)/M**2 + 0.0416666666666667*nHat*(-19.0*S_n*nu + 33.0*S_n - 10.0*Sigma_n*delta*nu + 33.0*Sigma_n*delta)/M**2
                return L_coeff*(L_0 + v**2*(L_2 + v*(L_SO_3 + v*(L_4 + v*(L_6*v + L_SO_5)))))

            dydt=np.zeros(9)
            rfrak_frame=np.zeros(3)
            rfrak_frame[0] = y[5]
            rfrak_frame[1] = y[6]
            rfrak_frame[2] = y[7]
            rfrakdot_frame = quaternion.quaternion_time_series.frame_from_angular_velocity_integrand(rfrak_frame,\
                quaternion.as_float_array(OmegaVec())[1:])
            dydt[0] = dvdt
            if(EvolveSpin1):
                dydt[1], dydt[2]= FrameFromAngularVelocity_2D_Integrand(y[1], y[2],\
                    quaternion.as_float_array(S_chi1.inverse()*OmegaVec_chiVec_1()*S_chi1))
            else:
                dydt[1] = 0.0
                dydt[2] = 0.0
            if(EvolveSpin2):
                dydt[3], dydt[4] = FrameFromAngularVelocity_2D_Integrand(y[3], y[4],\
                    quaternion.as_float_array(S_chi2.inverse()*OmegaVec_chiVec_2()*S_chi2))
            else:
                dydt[3] = 0.0
                dydt[4] = 0.0
            dydt[5] = rfrakdot_frame[0]
            dydt[6] = rfrakdot_frame[1]
            dydt[7] = rfrakdot_frame[2]
            dydt[8] = v*v*v/M

            return dydt
        
        
        yy=solve_ivp(TaylorT1, [0,114000], [v,rfrak_chi1_x,\
            rfrak_chi1_y,rfrak_chi2_x,rfrak_chi2_y,rfrak_frame_x,\
            rfrak_frame_y,rfrak_frame_z,Phi], method='DOP853',\
            t_eval=np.arange(0,114000,100), dense_output=True)
        return yy    

#        def TaylorT1_phi(t,phi):
#            dphidt_T1=v_of_t.y[0][v_of_t.t>=t][0]**3.0/M
#            return dphidt_T1
#        phi_of_t=solve_ivp(TaylorT1_phi, [0, v_of_t.t[-1]], [0.0], method='RK45', t_eval=v_of_t.t)
        
        #Calculate waveform modes

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
        EvolveSpin1=dot(S_chi1,S_chi1).norm()>1e-12
        EvolveSpin2=dot(S_chi2,S_chi2).norm()>1e-12

        def Recalculate(t, y):
            v = y[0]
            rfrak_chi1_x = y[1]
            rfrak_chi1_y = y[2]
            rfrak_chi2_x = y[3]
            rfrak_chi2_y = y[4]
            rfrak_frame_x = y[5]
            rfrak_frame_y = y[6]
            rfrak_frame_z = y[7]
            Phi = y[8]
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

        def FrameFromAngularVelocity_2D_Integrand(rfrak_x, rfrak_y, Omega):
            rfrakMag = np.sqrt(rfrak_x*rfrak_x+rfrak_y*rfrak_y)
            rfrakDot_x = Omega[0]/2.0
            rfrakDot_y = Omega[1]/2.0
            if abs(np.sin(rfrakMag)) > 1e-12 and abs(np.cos(rfrakMag)) > 1e-12:
                omega_v = (Omega[0]*(-rfrak_y/rfrakMag)+Omega[1]*(rfrak_x/rfrakMag))*np.tan(rfrakMag)-Omega[2]
                Omega[0] += -omega_v*np.sin(2*rfrakMag)*(-rfrak_y/rfrakMag)
                Omega[1] += -omega_v*np.sin(2*rfrakMag)*(rfrak_x/rfrakMag)
                Omega[2] +=  omega_v*np.cos(2*rfrakMag)
                dotTerm = (rfrak_x*Omega[0]+rfrak_y*Omega[1])/(rfrakMag*rfrakMag)
                cotTerm = rfrakMag/(2*np.tan(rfrakMag))
                rfrakDot_x = (Omega[0] - rfrak_x*dotTerm)*cotTerm + rfrak_x*dotTerm/2. - 0.5*Omega[2]*rfrak_y
                rfrakDot_y = (Omega[1] - rfrak_y*dotTerm)*cotTerm + rfrak_y*dotTerm/2. + 0.5*Omega[2]*rfrak_x
            return rfrakDot_x, rfrakDot_y
            

        #Evolve PN
        def TaylorT1(t, y):
            v = y[0]
            rfrak_chi1_x = y[1]
            rfrak_chi1_y = y[2]
            rfrak_chi2_x = y[3]
            rfrak_chi2_y = y[4]
            rfrak_frame_x = y[5]
            rfrak_frame_y = y[6]
            rfrak_frame_z = y[7]
            Phi = y[8]
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

            if v>=1.0: 
                raise ValueError("Beyond domain of PN validity")
            Flux = Fcal_coeff*(Fcal_0 + v**2*(Fcal_2 + v*(Fcal_3 + Fcal_SO_3 + v*(Fcal_4 + Fcal_SQ_4 + v*(Fcal_5 + Fcal_SO_5 + v*(Fcal_6 + Fcal_SO_6 + Fcal_lnv_6*logv + v*(Fcal_7 + Fcal_SO_7)))))))
            dEdv = -0.5*M*nu*v*(2.0*E_0 + v**2*(4.0*E_2 + v*(5.0*E_SO_3 + v*(6.0*E_4 + 6.0*E_SQ_4 + v*(7.0*E_SO_5 + v*(8.0*E_6 + 9.0*E_SO_7*v))))))
            Absorption = 0
            dvdt_T1 = (-Absorption - Flux)/dEdv
            try:
                dvdt_T1=quaternion.as_float_array(dvdt_T1)[0]
            except:
                pass
            if dvdt_T1<1.0e-12:
                raise ValueError("v is decreasing")
            return CommonRHS(dvdt_T1, y)

        def CommonRHS(dvdt, y):
            v = y[0]
            rfrak_chi1_x = y[1]
            rfrak_chi1_y = y[2]
            rfrak_chi2_x = y[3]
            rfrak_chi2_y = y[4]
            rfrak_frame_x = y[5]
            rfrak_frame_y = y[6]
            rfrak_frame_z = y[7]
            Phi = y[8]
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
                a_ell_4 = S_n*(5.77777777777778*nu**2 + 14.75*nu + 1.5) + Sigma_n*delta*(2.83333333333333*nu**2 + 9.125*nu + 1.5)
                gamma_PN_5 = (S_ell*(0.888888888888889*nu + 3.33333333333333) + 2.0*Sigma_ell*delta)/M**2
                a_ell_2 = S_n*(-9.66666666666667*nu - 10.0) + Sigma_n*delta*(-4.5*nu - 6.0)
                a_ell_0 = 7.0*S_n + 3.0*Sigma_n*delta
                gamma_PN_6 = 0.0123456790123457*nu**3 + 6.36111111111111*nu**2 - 2.98177812235564*nu + 1.0
                gamma_PN_4 = 1.0 - 5.41666666666667*nu
                gamma_PN_3 = (1.66666666666667*S_ell + Sigma_ell*delta)/M**2
                gamma_PN_0 = 1.00000000000000
                gamma_PN_7 = (S_ell*(-6.0*nu**2 - 10.5833333333333*nu + 5.0) - 2.66666666666667*Sigma_ell*delta*nu**2 + Sigma_ell*delta*(3.0 - 10.1666666666667*nu))/M**2
                gamma_PN_2 = 1.0 - 0.333333333333333*nu
                return ellHat*v**3/M + nHat*v**6*(a_ell_0 + v**2*(a_ell_2 + a_ell_4*v**2))*(gamma_PN_0 + v**2*(gamma_PN_2 + v*(gamma_PN_3 + v*(gamma_PN_4 + v*(gamma_PN_5 + v*(gamma_PN_6 + gamma_PN_7*v))))))/M**3

            def OrbitalAngularMomentum():
                L_SO_3 = 0.166666666666667*ellHat*(-35.0*S_ell - 15.0*Sigma_ell*delta)/M**2 + lambdaHat*(-3.0*S_lambda - Sigma_lambda*delta)/M**2 + 0.5*nHat*(S_n + Sigma_n*delta)/M**2
                L_0 = ellHat
                L_SO_7 = 0.0625*ellHat*(-29.0*S_ell*nu**2 + 1101.0*S_ell*nu - 405.0*S_ell - 15.0*Sigma_ell*delta*nu**2 + 468.0*Sigma_ell*delta*nu - 81.0*Sigma_ell*delta)/M**2 + 0.0416666666666667*lambdaHat*(-32.0*S_lambda*nu**2 + 2.0*S_lambda*nu - 174.0*S_lambda - 16.0*Sigma_lambda*delta*nu**2 - 79.0*Sigma_lambda*delta*nu - 12.0*Sigma_lambda*delta)/M**2 + 0.0208333333333333*nHat*(11.0*S_n*nu**2 - 1331.0*S_n*nu + 183.0*S_n + 5.0*Sigma_n*delta*nu**2 - 734.0*Sigma_n*delta*nu + 183.0*Sigma_n*delta)/M**2
                L_coeff = M**2*nu/v
                L_4 = ellHat*(0.0416666666666667*nu**2 - 2.375*nu + 3.375)
                L_2 = ellHat*(0.166666666666667*nu + 1.5)
                L_6 = ellHat*(0.00540123456790123*nu**3 + 1.29166666666667*nu**2 - 30.9797035925835*nu + 8.4375)
                L_SO_5 = 0.0138888888888889*ellHat*(427.0*S_ell*nu - 693.0*S_ell + 210.0*Sigma_ell*delta*nu - 189.0*Sigma_ell*delta)/M**2 + 0.166666666666667*lambdaHat*(18.0*S_lambda*nu - 21.0*S_lambda + 8.0*Sigma_lambda*delta*nu - 3.0*Sigma_lambda*delta)/M**2 + 0.0416666666666667*nHat*(-19.0*S_n*nu + 33.0*S_n - 10.0*Sigma_n*delta*nu + 33.0*Sigma_n*delta)/M**2
                return L_coeff*(L_0 + v**2*(L_2 + v*(L_SO_3 + v*(L_4 + v*(L_SO_5 + v*(L_6 + L_SO_7*v))))))

            dydt=np.zeros(9)
            rfrak_frame=np.zeros(3)
            rfrak_frame[0] = y[5]
            rfrak_frame[1] = y[6]
            rfrak_frame[2] = y[7]
            rfrakdot_frame = quaternion.quaternion_time_series.frame_from_angular_velocity_integrand(rfrak_frame,\
                quaternion.as_float_array(OmegaVec())[1:])
            dydt[0] = dvdt
            if(EvolveSpin1):
                dydt[1], dydt[2]= FrameFromAngularVelocity_2D_Integrand(y[1], y[2],\
                    quaternion.as_float_array(S_chi1.inverse()*OmegaVec_chiVec_1()*S_chi1))
            else:
                dydt[1] = 0.0
                dydt[2] = 0.0
            if(EvolveSpin2):
                dydt[3], dydt[4] = FrameFromAngularVelocity_2D_Integrand(y[3], y[4],\
                    quaternion.as_float_array(S_chi2.inverse()*OmegaVec_chiVec_2()*S_chi2))
            else:
                dydt[3] = 0.0
                dydt[4] = 0.0
            dydt[5] = rfrakdot_frame[0]
            dydt[6] = rfrakdot_frame[1]
            dydt[7] = rfrakdot_frame[2]
            dydt[8] = v*v*v/M

            return dydt
        
        
        yy=solve_ivp(TaylorT1, [0,114000], [v,rfrak_chi1_x,\
            rfrak_chi1_y,rfrak_chi2_x,rfrak_chi2_y,rfrak_frame_x,\
            rfrak_frame_y,rfrak_frame_z,Phi], method='DOP853',\
            t_eval=np.arange(0,114000,100), dense_output=True)
        return yy    

#        def TaylorT1_phi(t,phi):
#            dphidt_T1=v_of_t.y[0][v_of_t.t>=t][0]**3.0/M
#            return dphidt_T1
#        phi_of_t=solve_ivp(TaylorT1_phi, [0, v_of_t.t[-1]], [0.0], method='RK45', t_eval=v_of_t.t)
        
        #Calculate waveform modes
