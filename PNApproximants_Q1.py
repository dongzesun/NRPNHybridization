# File produced automatically by OrbitalEvolutionCodeGen_Q.ipynb
from scipy.integrate import solve_ivp
import numpy as np
from numpy import conjugate, dot, exp, log, pi
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
            a_ell_4 = S_n*(5.77777777777778*nu**2 + 14.75*nu + 1.5) + Sigma_n*delta*(2.83333333333333*nu**2 + 9.125*nu + 1.5)
            gamma_PN_0 = 1.00000000000000
            a_ell_2 = S_n*(-9.66666666666667*nu - 10.0) + Sigma_n*delta*(-4.5*nu - 6.0)
            gamma_PN_2 = 1.0 - 0.333333333333333*nu
            gamma_PN_3 = (1.66666666666667*S_ell + Sigma_ell*delta)/M**2
            gamma_PN_4 = 1.0 - 5.41666666666667*nu
            gamma_PN_5 = (S_ell*(0.888888888888889*nu + 3.33333333333333) + 2.0*Sigma_ell*delta)/M**2
            a_ell_0 = 7.0*S_n + 3.0*Sigma_n*delta
            gamma_PN_6 = 0.0123456790123457*nu**3 + 6.36111111111111*nu**2 - 2.98177812235564*nu + 1.0
            return ellHat*v**3/M + nHat*v**6*(a_ell_0 + v**2*(a_ell_2 + a_ell_4*v**2))*(gamma_PN_0 + v**2*(gamma_PN_2 + v*(gamma_PN_3 + v*(gamma_PN_4 + v*(gamma_PN_5 + v*(gamma_PN_6 + gamma_PN_7*v))))))/M**3

        def OrbitalAngularMomentum():
            L_coeff = M**2*nu/v
            L_0 = ellHat
            return L_0*L_coeff


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


class WaveformModes_0PN : public WaveformModes_Base {
private:
  const double M1, M2;
  double v;
  const double M, nu;

public:
  WaveformModes_0PN(M1_i, M2_i, v_i) :
        M1=M1_i
        M2=M2_i
        v=v_i
        M=M1 + M2
        nu=M1*M2/M**2
  { }

  using WaveformModes_Base::operator();

  std::vector<std::complex<double> > operator()(
    const double v_k,
    const double chi1_n_k, const double chi1_lambda_k, const double chi1_ell_k,
    const double chi2_n_k, const double chi2_lambda_k, const double chi2_ell_k)
  {
    v = v_k;



    std::vector<std::complex<double> > Modes(77);
    std::complex<double> Symm, Asymm;
    // (ell, m) = (2, +/- 0)
    Symm = hHat_2_0_0*rhOverM_coeff;
    Asymm = 0;
    Modes[2] = Symm + Asymm;
    // (ell, m) = (2, +/- 1)
    Symm = 0;
    Asymm = 0;
    Modes[3] = Symm + Asymm;
    Modes[1] = std::conj(Symm - Asymm);
    // (ell, m) = (2, +/- 2)
    Symm = hHat_2_2_0*rhOverM_coeff;
    Asymm = 0;
    Modes[4] = Symm + Asymm;
    Modes[0] = std::conj(Symm - Asymm);
    // (ell, m) = (3, +/- 0)
    Symm = 0;
    Asymm = 0;
    Modes[8] = Symm + Asymm;
    // (ell, m) = (3, +/- 1)
    Symm = 0;
    Asymm = 0;
    Modes[9] = Symm + Asymm;
    Modes[7] = -std::conj(Symm - Asymm);
    // (ell, m) = (3, +/- 2)
    Symm = 0;
    Asymm = 0;
    Modes[10] = Symm + Asymm;
    Modes[6] = -std::conj(Symm - Asymm);
    // (ell, m) = (3, +/- 3)
    Symm = 0;
    Asymm = 0;
    Modes[11] = Symm + Asymm;
    Modes[5] = -std::conj(Symm - Asymm);
    // (ell, m) = (4, +/- 0)
    Symm = hHat_4_0_0*rhOverM_coeff;
    Asymm = 0;
    Modes[16] = Symm + Asymm;
    // (ell, m) = (4, +/- 1)
    Symm = 0;
    Asymm = 0;
    Modes[17] = Symm + Asymm;
    Modes[15] = std::conj(Symm - Asymm);
    // (ell, m) = (4, +/- 2)
    Symm = 0;
    Asymm = 0;
    Modes[18] = Symm + Asymm;
    Modes[14] = std::conj(Symm - Asymm);
    // (ell, m) = (4, +/- 3)
    Symm = 0;
    Asymm = 0;
    Modes[19] = Symm + Asymm;
    Modes[13] = std::conj(Symm - Asymm);
    // (ell, m) = (4, +/- 4)
    Symm = 0;
    Asymm = 0;
    Modes[20] = Symm + Asymm;
    Modes[12] = std::conj(Symm - Asymm);
    // (ell, m) = (5, +/- 0)
    Symm = 0;
    Asymm = 0;
    Modes[26] = Symm + Asymm;
    // (ell, m) = (5, +/- 1)
    Symm = 0;
    Asymm = 0;
    Modes[27] = Symm + Asymm;
    Modes[25] = -std::conj(Symm - Asymm);
    // (ell, m) = (5, +/- 2)
    Symm = 0;
    Asymm = 0;
    Modes[28] = Symm + Asymm;
    Modes[24] = -std::conj(Symm - Asymm);
    // (ell, m) = (5, +/- 3)
    Symm = 0;
    Asymm = 0;
    Modes[29] = Symm + Asymm;
    Modes[23] = -std::conj(Symm - Asymm);
    // (ell, m) = (5, +/- 4)
    Symm = 0;
    Asymm = 0;
    Modes[30] = Symm + Asymm;
    Modes[22] = -std::conj(Symm - Asymm);
    // (ell, m) = (5, +/- 5)
    Symm = 0;
    Asymm = 0;
    Modes[31] = Symm + Asymm;
    Modes[21] = -std::conj(Symm - Asymm);
    // (ell, m) = (6, +/- 0)
    Symm = 0;
    Asymm = 0;
    Modes[38] = Symm + Asymm;
    // (ell, m) = (6, +/- 1)
    Symm = 0;
    Asymm = 0;
    Modes[39] = Symm + Asymm;
    Modes[37] = std::conj(Symm - Asymm);
    // (ell, m) = (6, +/- 2)
    Symm = 0;
    Asymm = 0;
    Modes[40] = Symm + Asymm;
    Modes[36] = std::conj(Symm - Asymm);
    // (ell, m) = (6, +/- 3)
    Symm = 0;
    Asymm = 0;
    Modes[41] = Symm + Asymm;
    Modes[35] = std::conj(Symm - Asymm);
    // (ell, m) = (6, +/- 4)
    Symm = 0;
    Asymm = 0;
    Modes[42] = Symm + Asymm;
    Modes[34] = std::conj(Symm - Asymm);
    // (ell, m) = (6, +/- 5)
    Symm = 0;
    Asymm = 0;
    Modes[43] = Symm + Asymm;
    Modes[33] = std::conj(Symm - Asymm);
    // (ell, m) = (6, +/- 6)
    Symm = 0;
    Asymm = 0;
    Modes[44] = Symm + Asymm;
    Modes[32] = std::conj(Symm - Asymm);
    // (ell, m) = (7, +/- 0)
    Symm = 0;
    Asymm = 0;
    Modes[52] = Symm + Asymm;
    // (ell, m) = (7, +/- 1)
    Symm = 0;
    Asymm = 0;
    Modes[53] = Symm + Asymm;
    Modes[51] = -std::conj(Symm - Asymm);
    // (ell, m) = (7, +/- 2)
    Symm = 0;
    Asymm = 0;
    Modes[54] = Symm + Asymm;
    Modes[50] = -std::conj(Symm - Asymm);
    // (ell, m) = (7, +/- 3)
    Symm = 0;
    Asymm = 0;
    Modes[55] = Symm + Asymm;
    Modes[49] = -std::conj(Symm - Asymm);
    // (ell, m) = (7, +/- 4)
    Symm = 0;
    Asymm = 0;
    Modes[56] = Symm + Asymm;
    Modes[48] = -std::conj(Symm - Asymm);
    // (ell, m) = (7, +/- 5)
    Symm = 0;
    Asymm = 0;
    Modes[57] = Symm + Asymm;
    Modes[47] = -std::conj(Symm - Asymm);
    // (ell, m) = (7, +/- 6)
    Symm = 0;
    Asymm = 0;
    Modes[58] = Symm + Asymm;
    Modes[46] = -std::conj(Symm - Asymm);
    // (ell, m) = (7, +/- 7)
    Symm = 0;
    Asymm = 0;
    Modes[59] = Symm + Asymm;
    Modes[45] = -std::conj(Symm - Asymm);
    // (ell, m) = (8, +/- 0)
    Symm = 0;
    Asymm = 0;
    Modes[68] = Symm + Asymm;
    // (ell, m) = (8, +/- 1)
    Symm = 0;
    Asymm = 0;
    Modes[69] = Symm + Asymm;
    Modes[67] = std::conj(Symm - Asymm);
    // (ell, m) = (8, +/- 2)
    Symm = 0;
    Asymm = 0;
    Modes[70] = Symm + Asymm;
    Modes[66] = std::conj(Symm - Asymm);
    // (ell, m) = (8, +/- 3)
    Symm = 0;
    Asymm = 0;
    Modes[71] = Symm + Asymm;
    Modes[65] = std::conj(Symm - Asymm);
    // (ell, m) = (8, +/- 4)
    Symm = 0;
    Asymm = 0;
    Modes[72] = Symm + Asymm;
    Modes[64] = std::conj(Symm - Asymm);
    // (ell, m) = (8, +/- 5)
    Symm = 0;
    Asymm = 0;
    Modes[73] = Symm + Asymm;
    Modes[63] = std::conj(Symm - Asymm);
    // (ell, m) = (8, +/- 6)
    Symm = 0;
    Asymm = 0;
    Modes[74] = Symm + Asymm;
    Modes[62] = std::conj(Symm - Asymm);
    // (ell, m) = (8, +/- 7)
    Symm = 0;
    Asymm = 0;
    Modes[75] = Symm + Asymm;
    Modes[61] = std::conj(Symm - Asymm);
    // (ell, m) = (8, +/- 8)
    Symm = 0;
    Asymm = 0;
    Modes[76] = Symm + Asymm;
    Modes[60] = std::conj(Symm - Asymm);

    return Modes;
  }

}; // class WaveformModes_0PN : public WaveformModes_Base

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
            a_ell_4 = S_n*(5.77777777777778*nu**2 + 14.75*nu + 1.5) + Sigma_n*delta*(2.83333333333333*nu**2 + 9.125*nu + 1.5)
            gamma_PN_0 = 1.00000000000000
            a_ell_2 = S_n*(-9.66666666666667*nu - 10.0) + Sigma_n*delta*(-4.5*nu - 6.0)
            gamma_PN_2 = 1.0 - 0.333333333333333*nu
            gamma_PN_3 = (1.66666666666667*S_ell + Sigma_ell*delta)/M**2
            gamma_PN_4 = 1.0 - 5.41666666666667*nu
            gamma_PN_5 = (S_ell*(0.888888888888889*nu + 3.33333333333333) + 2.0*Sigma_ell*delta)/M**2
            a_ell_0 = 7.0*S_n + 3.0*Sigma_n*delta
            gamma_PN_6 = 0.0123456790123457*nu**3 + 6.36111111111111*nu**2 - 2.98177812235564*nu + 1.0
            return ellHat*v**3/M + nHat*v**6*(a_ell_0 + v**2*(a_ell_2 + a_ell_4*v**2))*(gamma_PN_0 + v**2*(gamma_PN_2 + v*(gamma_PN_3 + v*(gamma_PN_4 + v*(gamma_PN_5 + v*(gamma_PN_6 + gamma_PN_7*v))))))/M**3

        def OrbitalAngularMomentum():
            L_coeff = M**2*nu/v
            L_0 = ellHat
            return L_0*L_coeff


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


class WaveformModes_0p50PN : public WaveformModes_Base {
private:
  const double M1, M2;
  double v;
  const double M, delta, nu;

public:
  WaveformModes_0p50PN(M1_i, M2_i, v_i) :
        M1=M1_i
        M2=M2_i
        v=v_i
        M=M1 + M2
        delta=(M1 - M2)/M
        nu=M1*M2/M**2
  { }

  using WaveformModes_Base::operator();

  std::vector<std::complex<double> > operator()(
    const double v_k,
    const double chi1_n_k, const double chi1_lambda_k, const double chi1_ell_k,
    const double chi2_n_k, const double chi2_lambda_k, const double chi2_ell_k)
  {
    v = v_k;



    std::vector<std::complex<double> > Modes(77);
    std::complex<double> Symm, Asymm;
    // (ell, m) = (2, +/- 0)
    Symm = hHat_2_0_0*rhOverM_coeff;
    Asymm = 0;
    Modes[2] = Symm + Asymm;
    // (ell, m) = (2, +/- 1)
    Symm = hHat_2_1_1*rhOverM_coeff*v;
    Asymm = 0;
    Modes[3] = Symm + Asymm;
    Modes[1] = std::conj(Symm - Asymm);
    // (ell, m) = (2, +/- 2)
    Symm = hHat_2_2_0*rhOverM_coeff;
    Asymm = 0;
    Modes[4] = Symm + Asymm;
    Modes[0] = std::conj(Symm - Asymm);
    // (ell, m) = (3, +/- 0)
    Symm = 0;
    Asymm = 0;
    Modes[8] = Symm + Asymm;
    // (ell, m) = (3, +/- 1)
    Symm = hHat_3_1_1*rhOverM_coeff*v;
    Asymm = 0;
    Modes[9] = Symm + Asymm;
    Modes[7] = -std::conj(Symm - Asymm);
    // (ell, m) = (3, +/- 2)
    Symm = 0;
    Asymm = 0;
    Modes[10] = Symm + Asymm;
    Modes[6] = -std::conj(Symm - Asymm);
    // (ell, m) = (3, +/- 3)
    Symm = hHat_3_3_1*rhOverM_coeff*v;
    Asymm = 0;
    Modes[11] = Symm + Asymm;
    Modes[5] = -std::conj(Symm - Asymm);
    // (ell, m) = (4, +/- 0)
    Symm = hHat_4_0_0*rhOverM_coeff;
    Asymm = 0;
    Modes[16] = Symm + Asymm;
    // (ell, m) = (4, +/- 1)
    Symm = 0;
    Asymm = 0;
    Modes[17] = Symm + Asymm;
    Modes[15] = std::conj(Symm - Asymm);
    // (ell, m) = (4, +/- 2)
    Symm = 0;
    Asymm = 0;
    Modes[18] = Symm + Asymm;
    Modes[14] = std::conj(Symm - Asymm);
    // (ell, m) = (4, +/- 3)
    Symm = 0;
    Asymm = 0;
    Modes[19] = Symm + Asymm;
    Modes[13] = std::conj(Symm - Asymm);
    // (ell, m) = (4, +/- 4)
    Symm = 0;
    Asymm = 0;
    Modes[20] = Symm + Asymm;
    Modes[12] = std::conj(Symm - Asymm);
    // (ell, m) = (5, +/- 0)
    Symm = 0;
    Asymm = 0;
    Modes[26] = Symm + Asymm;
    // (ell, m) = (5, +/- 1)
    Symm = 0;
    Asymm = 0;
    Modes[27] = Symm + Asymm;
    Modes[25] = -std::conj(Symm - Asymm);
    // (ell, m) = (5, +/- 2)
    Symm = 0;
    Asymm = 0;
    Modes[28] = Symm + Asymm;
    Modes[24] = -std::conj(Symm - Asymm);
    // (ell, m) = (5, +/- 3)
    Symm = 0;
    Asymm = 0;
    Modes[29] = Symm + Asymm;
    Modes[23] = -std::conj(Symm - Asymm);
    // (ell, m) = (5, +/- 4)
    Symm = 0;
    Asymm = 0;
    Modes[30] = Symm + Asymm;
    Modes[22] = -std::conj(Symm - Asymm);
    // (ell, m) = (5, +/- 5)
    Symm = 0;
    Asymm = 0;
    Modes[31] = Symm + Asymm;
    Modes[21] = -std::conj(Symm - Asymm);
    // (ell, m) = (6, +/- 0)
    Symm = 0;
    Asymm = 0;
    Modes[38] = Symm + Asymm;
    // (ell, m) = (6, +/- 1)
    Symm = 0;
    Asymm = 0;
    Modes[39] = Symm + Asymm;
    Modes[37] = std::conj(Symm - Asymm);
    // (ell, m) = (6, +/- 2)
    Symm = 0;
    Asymm = 0;
    Modes[40] = Symm + Asymm;
    Modes[36] = std::conj(Symm - Asymm);
    // (ell, m) = (6, +/- 3)
    Symm = 0;
    Asymm = 0;
    Modes[41] = Symm + Asymm;
    Modes[35] = std::conj(Symm - Asymm);
    // (ell, m) = (6, +/- 4)
    Symm = 0;
    Asymm = 0;
    Modes[42] = Symm + Asymm;
    Modes[34] = std::conj(Symm - Asymm);
    // (ell, m) = (6, +/- 5)
    Symm = 0;
    Asymm = 0;
    Modes[43] = Symm + Asymm;
    Modes[33] = std::conj(Symm - Asymm);
    // (ell, m) = (6, +/- 6)
    Symm = 0;
    Asymm = 0;
    Modes[44] = Symm + Asymm;
    Modes[32] = std::conj(Symm - Asymm);
    // (ell, m) = (7, +/- 0)
    Symm = 0;
    Asymm = 0;
    Modes[52] = Symm + Asymm;
    // (ell, m) = (7, +/- 1)
    Symm = 0;
    Asymm = 0;
    Modes[53] = Symm + Asymm;
    Modes[51] = -std::conj(Symm - Asymm);
    // (ell, m) = (7, +/- 2)
    Symm = 0;
    Asymm = 0;
    Modes[54] = Symm + Asymm;
    Modes[50] = -std::conj(Symm - Asymm);
    // (ell, m) = (7, +/- 3)
    Symm = 0;
    Asymm = 0;
    Modes[55] = Symm + Asymm;
    Modes[49] = -std::conj(Symm - Asymm);
    // (ell, m) = (7, +/- 4)
    Symm = 0;
    Asymm = 0;
    Modes[56] = Symm + Asymm;
    Modes[48] = -std::conj(Symm - Asymm);
    // (ell, m) = (7, +/- 5)
    Symm = 0;
    Asymm = 0;
    Modes[57] = Symm + Asymm;
    Modes[47] = -std::conj(Symm - Asymm);
    // (ell, m) = (7, +/- 6)
    Symm = 0;
    Asymm = 0;
    Modes[58] = Symm + Asymm;
    Modes[46] = -std::conj(Symm - Asymm);
    // (ell, m) = (7, +/- 7)
    Symm = 0;
    Asymm = 0;
    Modes[59] = Symm + Asymm;
    Modes[45] = -std::conj(Symm - Asymm);
    // (ell, m) = (8, +/- 0)
    Symm = 0;
    Asymm = 0;
    Modes[68] = Symm + Asymm;
    // (ell, m) = (8, +/- 1)
    Symm = 0;
    Asymm = 0;
    Modes[69] = Symm + Asymm;
    Modes[67] = std::conj(Symm - Asymm);
    // (ell, m) = (8, +/- 2)
    Symm = 0;
    Asymm = 0;
    Modes[70] = Symm + Asymm;
    Modes[66] = std::conj(Symm - Asymm);
    // (ell, m) = (8, +/- 3)
    Symm = 0;
    Asymm = 0;
    Modes[71] = Symm + Asymm;
    Modes[65] = std::conj(Symm - Asymm);
    // (ell, m) = (8, +/- 4)
    Symm = 0;
    Asymm = 0;
    Modes[72] = Symm + Asymm;
    Modes[64] = std::conj(Symm - Asymm);
    // (ell, m) = (8, +/- 5)
    Symm = 0;
    Asymm = 0;
    Modes[73] = Symm + Asymm;
    Modes[63] = std::conj(Symm - Asymm);
    // (ell, m) = (8, +/- 6)
    Symm = 0;
    Asymm = 0;
    Modes[74] = Symm + Asymm;
    Modes[62] = std::conj(Symm - Asymm);
    // (ell, m) = (8, +/- 7)
    Symm = 0;
    Asymm = 0;
    Modes[75] = Symm + Asymm;
    Modes[61] = std::conj(Symm - Asymm);
    // (ell, m) = (8, +/- 8)
    Symm = 0;
    Asymm = 0;
    Modes[76] = Symm + Asymm;
    Modes[60] = std::conj(Symm - Asymm);

    return Modes;
  }

}; // class WaveformModes_0p50PN : public WaveformModes_Base

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
            a_ell_4 = S_n*(5.77777777777778*nu**2 + 14.75*nu + 1.5) + Sigma_n*delta*(2.83333333333333*nu**2 + 9.125*nu + 1.5)
            gamma_PN_0 = 1.00000000000000
            a_ell_2 = S_n*(-9.66666666666667*nu - 10.0) + Sigma_n*delta*(-4.5*nu - 6.0)
            gamma_PN_2 = 1.0 - 0.333333333333333*nu
            gamma_PN_3 = (1.66666666666667*S_ell + Sigma_ell*delta)/M**2
            gamma_PN_4 = 1.0 - 5.41666666666667*nu
            gamma_PN_5 = (S_ell*(0.888888888888889*nu + 3.33333333333333) + 2.0*Sigma_ell*delta)/M**2
            a_ell_0 = 7.0*S_n + 3.0*Sigma_n*delta
            gamma_PN_6 = 0.0123456790123457*nu**3 + 6.36111111111111*nu**2 - 2.98177812235564*nu + 1.0
            return ellHat*v**3/M + nHat*v**6*(a_ell_0 + v**2*(a_ell_2 + a_ell_4*v**2))*(gamma_PN_0 + v**2*(gamma_PN_2 + v*(gamma_PN_3 + v*(gamma_PN_4 + v*(gamma_PN_5 + v*(gamma_PN_6 + gamma_PN_7*v))))))/M**3

        def OrbitalAngularMomentum():
            L_coeff = M**2*nu/v
            L_2 = ellHat*(0.166666666666667*nu + 1.5)
            L_0 = ellHat
            return L_coeff*(L_0 + L_2*v**2)


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


class WaveformModes_1p0PN : public WaveformModes_Base {
private:
  const Quaternions.Quaternion xHat, yHat, zHat;
  const double M1, M2;
  double v;
  const Quaternions.Quaternion S_chi1, S_chi2;
  double rfrak_chi1_x, rfrak_chi1_y, rfrak_chi2_x, rfrak_chi2_y, rfrak_frame_x, rfrak_frame_y, rfrak_frame_z;
  const double M, delta, nu;
  Quaternions.Quaternion R, nHat, lambdaHat, ellHat, R_S1, R_S2, chiVec1, chiVec2;
  double chi1_n, chi1_lambda, chi1_ell, chi2_n, chi2_lambda, chi2_ell, Sigma_ell, Sigma_n, Sigma_lambda;

public:
  WaveformModes_1p0PN(xHat_i, yHat_i, zHat_i, M1_i, M2_i, v_i, S_chi1_i, S_chi2_i, rfrak_chi1_x_i, rfrak_chi1_y_i,
                      rfrak_chi2_x_i, rfrak_chi2_y_i, rfrak_frame_x_i, rfrak_frame_y_i, rfrak_frame_z_i) :
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
        chi1_n=dot(chiVec1,nHat)
        chi1_lambda=dot(chiVec1,lambdaHat)
        chi1_ell=dot(chiVec1,ellHat)
        chi2_n=dot(chiVec2,nHat)
        chi2_lambda=dot(chiVec2,lambdaHat)
        chi2_ell=dot(chiVec2,ellHat)
        Sigma_ell=M*(-M1*chi1_ell + M2*chi2_ell)
        Sigma_n=M*(-M1*chi1_n + M2*chi2_n)
        Sigma_lambda=M*(-M1*chi1_lambda + M2*chi2_lambda)
  { }

  using WaveformModes_Base::operator();

  std::vector<std::complex<double> > operator()(
    const double v_k,
    const double chi1_n_k, const double chi1_lambda_k, const double chi1_ell_k,
    const double chi2_n_k, const double chi2_lambda_k, const double chi2_ell_k)
  {
    v = v_k;
    chiVec1 = Quaternions::Quaternion(0., chi1_n_k, chi1_lambda_k, chi1_ell_k);
    chiVec2 = Quaternions::Quaternion(0., chi2_n_k, chi2_lambda_k, chi2_ell_k);

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
            Sigma_ell = M*(-M1*chi1_ell + M2*chi2_ell)
            Sigma_n = M*(-M1*chi1_n + M2*chi2_n)
            Sigma_lambda = M*(-M1*chi1_lambda + M2*chi2_lambda)

    std::vector<std::complex<double> > Modes(77);
    std::complex<double> Symm, Asymm;
    // (ell, m) = (2, +/- 0)
    Symm = hHat_2_0_0*rhOverM_coeff;
    Asymm = hHat_spin_Asymm_2_0_2*rhOverM_coeff*pow(v, 2);
    Modes[2] = Symm + Asymm;
    // (ell, m) = (2, +/- 1)
    Symm = rhOverM_coeff*v*(hHat_2_1_1 + hHat_spin_Symm_2_1_2*v);
    Asymm = 0;
    Modes[3] = Symm + Asymm;
    Modes[1] = std::conj(Symm - Asymm);
    // (ell, m) = (2, +/- 2)
    Symm = rhOverM_coeff*(hHat_2_2_0 + hHat_2_2_2*pow(v, 2));
    Asymm = hHat_spin_Asymm_2_2_2*rhOverM_coeff*pow(v, 2);
    Modes[4] = Symm + Asymm;
    Modes[0] = std::conj(Symm - Asymm);
    // (ell, m) = (3, +/- 0)
    Symm = 0;
    Asymm = 0;
    Modes[8] = Symm + Asymm;
    // (ell, m) = (3, +/- 1)
    Symm = hHat_3_1_1*rhOverM_coeff*v;
    Asymm = 0;
    Modes[9] = Symm + Asymm;
    Modes[7] = -std::conj(Symm - Asymm);
    // (ell, m) = (3, +/- 2)
    Symm = hHat_3_2_2*rhOverM_coeff*pow(v, 2);
    Asymm = 0;
    Modes[10] = Symm + Asymm;
    Modes[6] = -std::conj(Symm - Asymm);
    // (ell, m) = (3, +/- 3)
    Symm = hHat_3_3_1*rhOverM_coeff*v;
    Asymm = 0;
    Modes[11] = Symm + Asymm;
    Modes[5] = -std::conj(Symm - Asymm);
    // (ell, m) = (4, +/- 0)
    Symm = hHat_4_0_0*rhOverM_coeff;
    Asymm = 0;
    Modes[16] = Symm + Asymm;
    // (ell, m) = (4, +/- 1)
    Symm = 0;
    Asymm = 0;
    Modes[17] = Symm + Asymm;
    Modes[15] = std::conj(Symm - Asymm);
    // (ell, m) = (4, +/- 2)
    Symm = hHat_4_2_2*rhOverM_coeff*pow(v, 2);
    Asymm = 0;
    Modes[18] = Symm + Asymm;
    Modes[14] = std::conj(Symm - Asymm);
    // (ell, m) = (4, +/- 3)
    Symm = 0;
    Asymm = 0;
    Modes[19] = Symm + Asymm;
    Modes[13] = std::conj(Symm - Asymm);
    // (ell, m) = (4, +/- 4)
    Symm = hHat_4_4_2*rhOverM_coeff*pow(v, 2);
    Asymm = 0;
    Modes[20] = Symm + Asymm;
    Modes[12] = std::conj(Symm - Asymm);
    // (ell, m) = (5, +/- 0)
    Symm = 0;
    Asymm = 0;
    Modes[26] = Symm + Asymm;
    // (ell, m) = (5, +/- 1)
    Symm = 0;
    Asymm = 0;
    Modes[27] = Symm + Asymm;
    Modes[25] = -std::conj(Symm - Asymm);
    // (ell, m) = (5, +/- 2)
    Symm = 0;
    Asymm = 0;
    Modes[28] = Symm + Asymm;
    Modes[24] = -std::conj(Symm - Asymm);
    // (ell, m) = (5, +/- 3)
    Symm = 0;
    Asymm = 0;
    Modes[29] = Symm + Asymm;
    Modes[23] = -std::conj(Symm - Asymm);
    // (ell, m) = (5, +/- 4)
    Symm = 0;
    Asymm = 0;
    Modes[30] = Symm + Asymm;
    Modes[22] = -std::conj(Symm - Asymm);
    // (ell, m) = (5, +/- 5)
    Symm = 0;
    Asymm = 0;
    Modes[31] = Symm + Asymm;
    Modes[21] = -std::conj(Symm - Asymm);
    // (ell, m) = (6, +/- 0)
    Symm = 0;
    Asymm = 0;
    Modes[38] = Symm + Asymm;
    // (ell, m) = (6, +/- 1)
    Symm = 0;
    Asymm = 0;
    Modes[39] = Symm + Asymm;
    Modes[37] = std::conj(Symm - Asymm);
    // (ell, m) = (6, +/- 2)
    Symm = 0;
    Asymm = 0;
    Modes[40] = Symm + Asymm;
    Modes[36] = std::conj(Symm - Asymm);
    // (ell, m) = (6, +/- 3)
    Symm = 0;
    Asymm = 0;
    Modes[41] = Symm + Asymm;
    Modes[35] = std::conj(Symm - Asymm);
    // (ell, m) = (6, +/- 4)
    Symm = 0;
    Asymm = 0;
    Modes[42] = Symm + Asymm;
    Modes[34] = std::conj(Symm - Asymm);
    // (ell, m) = (6, +/- 5)
    Symm = 0;
    Asymm = 0;
    Modes[43] = Symm + Asymm;
    Modes[33] = std::conj(Symm - Asymm);
    // (ell, m) = (6, +/- 6)
    Symm = 0;
    Asymm = 0;
    Modes[44] = Symm + Asymm;
    Modes[32] = std::conj(Symm - Asymm);
    // (ell, m) = (7, +/- 0)
    Symm = 0;
    Asymm = 0;
    Modes[52] = Symm + Asymm;
    // (ell, m) = (7, +/- 1)
    Symm = 0;
    Asymm = 0;
    Modes[53] = Symm + Asymm;
    Modes[51] = -std::conj(Symm - Asymm);
    // (ell, m) = (7, +/- 2)
    Symm = 0;
    Asymm = 0;
    Modes[54] = Symm + Asymm;
    Modes[50] = -std::conj(Symm - Asymm);
    // (ell, m) = (7, +/- 3)
    Symm = 0;
    Asymm = 0;
    Modes[55] = Symm + Asymm;
    Modes[49] = -std::conj(Symm - Asymm);
    // (ell, m) = (7, +/- 4)
    Symm = 0;
    Asymm = 0;
    Modes[56] = Symm + Asymm;
    Modes[48] = -std::conj(Symm - Asymm);
    // (ell, m) = (7, +/- 5)
    Symm = 0;
    Asymm = 0;
    Modes[57] = Symm + Asymm;
    Modes[47] = -std::conj(Symm - Asymm);
    // (ell, m) = (7, +/- 6)
    Symm = 0;
    Asymm = 0;
    Modes[58] = Symm + Asymm;
    Modes[46] = -std::conj(Symm - Asymm);
    // (ell, m) = (7, +/- 7)
    Symm = 0;
    Asymm = 0;
    Modes[59] = Symm + Asymm;
    Modes[45] = -std::conj(Symm - Asymm);
    // (ell, m) = (8, +/- 0)
    Symm = 0;
    Asymm = 0;
    Modes[68] = Symm + Asymm;
    // (ell, m) = (8, +/- 1)
    Symm = 0;
    Asymm = 0;
    Modes[69] = Symm + Asymm;
    Modes[67] = std::conj(Symm - Asymm);
    // (ell, m) = (8, +/- 2)
    Symm = 0;
    Asymm = 0;
    Modes[70] = Symm + Asymm;
    Modes[66] = std::conj(Symm - Asymm);
    // (ell, m) = (8, +/- 3)
    Symm = 0;
    Asymm = 0;
    Modes[71] = Symm + Asymm;
    Modes[65] = std::conj(Symm - Asymm);
    // (ell, m) = (8, +/- 4)
    Symm = 0;
    Asymm = 0;
    Modes[72] = Symm + Asymm;
    Modes[64] = std::conj(Symm - Asymm);
    // (ell, m) = (8, +/- 5)
    Symm = 0;
    Asymm = 0;
    Modes[73] = Symm + Asymm;
    Modes[63] = std::conj(Symm - Asymm);
    // (ell, m) = (8, +/- 6)
    Symm = 0;
    Asymm = 0;
    Modes[74] = Symm + Asymm;
    Modes[62] = std::conj(Symm - Asymm);
    // (ell, m) = (8, +/- 7)
    Symm = 0;
    Asymm = 0;
    Modes[75] = Symm + Asymm;
    Modes[61] = std::conj(Symm - Asymm);
    // (ell, m) = (8, +/- 8)
    Symm = 0;
    Asymm = 0;
    Modes[76] = Symm + Asymm;
    Modes[60] = std::conj(Symm - Asymm);

    return Modes;
  }

}; // class WaveformModes_1p0PN : public WaveformModes_Base

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
            a_ell_4 = S_n*(5.77777777777778*nu**2 + 14.75*nu + 1.5) + Sigma_n*delta*(2.83333333333333*nu**2 + 9.125*nu + 1.5)
            gamma_PN_0 = 1.00000000000000
            a_ell_2 = S_n*(-9.66666666666667*nu - 10.0) + Sigma_n*delta*(-4.5*nu - 6.0)
            gamma_PN_2 = 1.0 - 0.333333333333333*nu
            gamma_PN_3 = (1.66666666666667*S_ell + Sigma_ell*delta)/M**2
            gamma_PN_4 = 1.0 - 5.41666666666667*nu
            gamma_PN_5 = (S_ell*(0.888888888888889*nu + 3.33333333333333) + 2.0*Sigma_ell*delta)/M**2
            a_ell_0 = 7.0*S_n + 3.0*Sigma_n*delta
            gamma_PN_6 = 0.0123456790123457*nu**3 + 6.36111111111111*nu**2 - 2.98177812235564*nu + 1.0
            return ellHat*v**3/M + nHat*v**6*(a_ell_0 + v**2*(a_ell_2 + a_ell_4*v**2))*(gamma_PN_0 + v**2*(gamma_PN_2 + v*(gamma_PN_3 + v*(gamma_PN_4 + v*(gamma_PN_5 + v*(gamma_PN_6 + gamma_PN_7*v))))))/M**3

        def OrbitalAngularMomentum():
            L_0 = ellHat
            L_2 = ellHat*(0.166666666666667*nu + 1.5)
            L_coeff = M**2*nu/v
            L_SO_3 = 0.166666666666667*ellHat*(-35.0*S_ell - 15.0*Sigma_ell*delta)/M**2 + lambdaHat*(-3.0*S_lambda - Sigma_lambda*delta)/M**2 + 0.5*nHat*(S_n + Sigma_n*delta)/M**2
            return L_coeff*(L_0 + v**2*(L_2 + L_SO_3*v))


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


class WaveformModes_1p5PN : public WaveformModes_Base {
private:
  const Quaternions.Quaternion xHat, yHat, zHat;
  const double M1, M2;
  double v;
  const Quaternions.Quaternion S_chi1, S_chi2;
  double rfrak_chi1_x, rfrak_chi1_y, rfrak_chi2_x, rfrak_chi2_y, rfrak_frame_x, rfrak_frame_y, rfrak_frame_z;
  const double M, delta, nu;
  Quaternions.Quaternion R, nHat, lambdaHat, ellHat, R_S1, R_S2, chiVec1, chiVec2;
  double chi1_n, chi1_lambda, chi1_ell, chi2_n, chi2_lambda, chi2_ell, S_ell, S_n, S_lambda, Sigma_ell, Sigma_n,
         Sigma_lambda;

public:
  WaveformModes_1p5PN(xHat_i, yHat_i, zHat_i, M1_i, M2_i, v_i, S_chi1_i, S_chi2_i, rfrak_chi1_x_i, rfrak_chi1_y_i,
                      rfrak_chi2_x_i, rfrak_chi2_y_i, rfrak_frame_x_i, rfrak_frame_y_i, rfrak_frame_z_i) :
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
  { }

  using WaveformModes_Base::operator();

  std::vector<std::complex<double> > operator()(
    const double v_k,
    const double chi1_n_k, const double chi1_lambda_k, const double chi1_ell_k,
    const double chi2_n_k, const double chi2_lambda_k, const double chi2_ell_k)
  {
    v = v_k;
    chiVec1 = Quaternions::Quaternion(0., chi1_n_k, chi1_lambda_k, chi1_ell_k);
    chiVec2 = Quaternions::Quaternion(0., chi2_n_k, chi2_lambda_k, chi2_ell_k);

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

    std::vector<std::complex<double> > Modes(77);
    std::complex<double> Symm, Asymm;
    // (ell, m) = (2, +/- 0)
    Symm = hHat_2_0_0*rhOverM_coeff;
    Asymm = hHat_spin_Asymm_2_0_2*rhOverM_coeff*pow(v, 2);
    Modes[2] = Symm + Asymm;
    // (ell, m) = (2, +/- 1)
    Symm = rhOverM_coeff*v*(hHat_2_1_1 + v*(hHat_2_1_3*v + hHat_spin_Symm_2_1_2));
    Asymm = hHat_spin_Asymm_2_1_3*rhOverM_coeff*pow(v, 3);
    Modes[3] = Symm + Asymm;
    Modes[1] = std::conj(Symm - Asymm);
    // (ell, m) = (2, +/- 2)
    Symm = rhOverM_coeff*(hHat_2_2_0 + pow(v, 2)*(hHat_2_2_2 + v*(hHat_2_2_3 + hHat_spin_Symm_2_2_3)));
    Asymm = hHat_spin_Asymm_2_2_2*rhOverM_coeff*pow(v, 2);
    Modes[4] = Symm + Asymm;
    Modes[0] = std::conj(Symm - Asymm);
    // (ell, m) = (3, +/- 0)
    Symm = 0;
    Asymm = 0;
    Modes[8] = Symm + Asymm;
    // (ell, m) = (3, +/- 1)
    Symm = rhOverM_coeff*v*(hHat_3_1_1 + hHat_3_1_3*pow(v, 2));
    Asymm = hHat_spin_Asymm_3_1_3*rhOverM_coeff*pow(v, 3);
    Modes[9] = Symm + Asymm;
    Modes[7] = -std::conj(Symm - Asymm);
    // (ell, m) = (3, +/- 2)
    Symm = rhOverM_coeff*pow(v, 2)*(hHat_3_2_2 + hHat_spin_Symm_3_2_3*v);
    Asymm = 0;
    Modes[10] = Symm + Asymm;
    Modes[6] = -std::conj(Symm - Asymm);
    // (ell, m) = (3, +/- 3)
    Symm = rhOverM_coeff*v*(hHat_3_3_1 + hHat_3_3_3*pow(v, 2));
    Asymm = hHat_spin_Asymm_3_3_3*rhOverM_coeff*pow(v, 3);
    Modes[11] = Symm + Asymm;
    Modes[5] = -std::conj(Symm - Asymm);
    // (ell, m) = (4, +/- 0)
    Symm = hHat_4_0_0*rhOverM_coeff;
    Asymm = 0;
    Modes[16] = Symm + Asymm;
    // (ell, m) = (4, +/- 1)
    Symm = hHat_4_1_3*rhOverM_coeff*pow(v, 3);
    Asymm = 0;
    Modes[17] = Symm + Asymm;
    Modes[15] = std::conj(Symm - Asymm);
    // (ell, m) = (4, +/- 2)
    Symm = hHat_4_2_2*rhOverM_coeff*pow(v, 2);
    Asymm = 0;
    Modes[18] = Symm + Asymm;
    Modes[14] = std::conj(Symm - Asymm);
    // (ell, m) = (4, +/- 3)
    Symm = hHat_4_3_3*rhOverM_coeff*pow(v, 3);
    Asymm = 0;
    Modes[19] = Symm + Asymm;
    Modes[13] = std::conj(Symm - Asymm);
    // (ell, m) = (4, +/- 4)
    Symm = hHat_4_4_2*rhOverM_coeff*pow(v, 2);
    Asymm = 0;
    Modes[20] = Symm + Asymm;
    Modes[12] = std::conj(Symm - Asymm);
    // (ell, m) = (5, +/- 0)
    Symm = 0;
    Asymm = 0;
    Modes[26] = Symm + Asymm;
    // (ell, m) = (5, +/- 1)
    Symm = hHat_5_1_3*rhOverM_coeff*pow(v, 3);
    Asymm = 0;
    Modes[27] = Symm + Asymm;
    Modes[25] = -std::conj(Symm - Asymm);
    // (ell, m) = (5, +/- 2)
    Symm = 0;
    Asymm = 0;
    Modes[28] = Symm + Asymm;
    Modes[24] = -std::conj(Symm - Asymm);
    // (ell, m) = (5, +/- 3)
    Symm = hHat_5_3_3*rhOverM_coeff*pow(v, 3);
    Asymm = 0;
    Modes[29] = Symm + Asymm;
    Modes[23] = -std::conj(Symm - Asymm);
    // (ell, m) = (5, +/- 4)
    Symm = 0;
    Asymm = 0;
    Modes[30] = Symm + Asymm;
    Modes[22] = -std::conj(Symm - Asymm);
    // (ell, m) = (5, +/- 5)
    Symm = hHat_5_5_3*rhOverM_coeff*pow(v, 3);
    Asymm = 0;
    Modes[31] = Symm + Asymm;
    Modes[21] = -std::conj(Symm - Asymm);
    // (ell, m) = (6, +/- 0)
    Symm = 0;
    Asymm = 0;
    Modes[38] = Symm + Asymm;
    // (ell, m) = (6, +/- 1)
    Symm = 0;
    Asymm = 0;
    Modes[39] = Symm + Asymm;
    Modes[37] = std::conj(Symm - Asymm);
    // (ell, m) = (6, +/- 2)
    Symm = 0;
    Asymm = 0;
    Modes[40] = Symm + Asymm;
    Modes[36] = std::conj(Symm - Asymm);
    // (ell, m) = (6, +/- 3)
    Symm = 0;
    Asymm = 0;
    Modes[41] = Symm + Asymm;
    Modes[35] = std::conj(Symm - Asymm);
    // (ell, m) = (6, +/- 4)
    Symm = 0;
    Asymm = 0;
    Modes[42] = Symm + Asymm;
    Modes[34] = std::conj(Symm - Asymm);
    // (ell, m) = (6, +/- 5)
    Symm = 0;
    Asymm = 0;
    Modes[43] = Symm + Asymm;
    Modes[33] = std::conj(Symm - Asymm);
    // (ell, m) = (6, +/- 6)
    Symm = 0;
    Asymm = 0;
    Modes[44] = Symm + Asymm;
    Modes[32] = std::conj(Symm - Asymm);
    // (ell, m) = (7, +/- 0)
    Symm = 0;
    Asymm = 0;
    Modes[52] = Symm + Asymm;
    // (ell, m) = (7, +/- 1)
    Symm = 0;
    Asymm = 0;
    Modes[53] = Symm + Asymm;
    Modes[51] = -std::conj(Symm - Asymm);
    // (ell, m) = (7, +/- 2)
    Symm = 0;
    Asymm = 0;
    Modes[54] = Symm + Asymm;
    Modes[50] = -std::conj(Symm - Asymm);
    // (ell, m) = (7, +/- 3)
    Symm = 0;
    Asymm = 0;
    Modes[55] = Symm + Asymm;
    Modes[49] = -std::conj(Symm - Asymm);
    // (ell, m) = (7, +/- 4)
    Symm = 0;
    Asymm = 0;
    Modes[56] = Symm + Asymm;
    Modes[48] = -std::conj(Symm - Asymm);
    // (ell, m) = (7, +/- 5)
    Symm = 0;
    Asymm = 0;
    Modes[57] = Symm + Asymm;
    Modes[47] = -std::conj(Symm - Asymm);
    // (ell, m) = (7, +/- 6)
    Symm = 0;
    Asymm = 0;
    Modes[58] = Symm + Asymm;
    Modes[46] = -std::conj(Symm - Asymm);
    // (ell, m) = (7, +/- 7)
    Symm = 0;
    Asymm = 0;
    Modes[59] = Symm + Asymm;
    Modes[45] = -std::conj(Symm - Asymm);
    // (ell, m) = (8, +/- 0)
    Symm = 0;
    Asymm = 0;
    Modes[68] = Symm + Asymm;
    // (ell, m) = (8, +/- 1)
    Symm = 0;
    Asymm = 0;
    Modes[69] = Symm + Asymm;
    Modes[67] = std::conj(Symm - Asymm);
    // (ell, m) = (8, +/- 2)
    Symm = 0;
    Asymm = 0;
    Modes[70] = Symm + Asymm;
    Modes[66] = std::conj(Symm - Asymm);
    // (ell, m) = (8, +/- 3)
    Symm = 0;
    Asymm = 0;
    Modes[71] = Symm + Asymm;
    Modes[65] = std::conj(Symm - Asymm);
    // (ell, m) = (8, +/- 4)
    Symm = 0;
    Asymm = 0;
    Modes[72] = Symm + Asymm;
    Modes[64] = std::conj(Symm - Asymm);
    // (ell, m) = (8, +/- 5)
    Symm = 0;
    Asymm = 0;
    Modes[73] = Symm + Asymm;
    Modes[63] = std::conj(Symm - Asymm);
    // (ell, m) = (8, +/- 6)
    Symm = 0;
    Asymm = 0;
    Modes[74] = Symm + Asymm;
    Modes[62] = std::conj(Symm - Asymm);
    // (ell, m) = (8, +/- 7)
    Symm = 0;
    Asymm = 0;
    Modes[75] = Symm + Asymm;
    Modes[61] = std::conj(Symm - Asymm);
    // (ell, m) = (8, +/- 8)
    Symm = 0;
    Asymm = 0;
    Modes[76] = Symm + Asymm;
    Modes[60] = std::conj(Symm - Asymm);

    return Modes;
  }

}; // class WaveformModes_1p5PN : public WaveformModes_Base

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
            a_ell_4 = S_n*(5.77777777777778*nu**2 + 14.75*nu + 1.5) + Sigma_n*delta*(2.83333333333333*nu**2 + 9.125*nu + 1.5)
            gamma_PN_0 = 1.00000000000000
            a_ell_2 = S_n*(-9.66666666666667*nu - 10.0) + Sigma_n*delta*(-4.5*nu - 6.0)
            gamma_PN_2 = 1.0 - 0.333333333333333*nu
            gamma_PN_3 = (1.66666666666667*S_ell + Sigma_ell*delta)/M**2
            gamma_PN_4 = 1.0 - 5.41666666666667*nu
            gamma_PN_5 = (S_ell*(0.888888888888889*nu + 3.33333333333333) + 2.0*Sigma_ell*delta)/M**2
            a_ell_0 = 7.0*S_n + 3.0*Sigma_n*delta
            gamma_PN_6 = 0.0123456790123457*nu**3 + 6.36111111111111*nu**2 - 2.98177812235564*nu + 1.0
            return ellHat*v**3/M + nHat*v**6*(a_ell_0 + v**2*(a_ell_2 + a_ell_4*v**2))*(gamma_PN_0 + v**2*(gamma_PN_2 + v*(gamma_PN_3 + v*(gamma_PN_4 + v*(gamma_PN_5 + v*(gamma_PN_6 + gamma_PN_7*v))))))/M**3

        def OrbitalAngularMomentum():
            L_4 = ellHat*(0.0416666666666667*nu**2 - 2.375*nu + 3.375)
            L_0 = ellHat
            L_2 = ellHat*(0.166666666666667*nu + 1.5)
            L_coeff = M**2*nu/v
            L_SO_3 = 0.166666666666667*ellHat*(-35.0*S_ell - 15.0*Sigma_ell*delta)/M**2 + lambdaHat*(-3.0*S_lambda - Sigma_lambda*delta)/M**2 + 0.5*nHat*(S_n + Sigma_n*delta)/M**2
            return L_coeff*(L_0 + v**2*(L_2 + v*(L_4*v + L_SO_3)))


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


class WaveformModes_2p0PN : public WaveformModes_Base {
private:
  const Quaternions.Quaternion xHat, yHat, zHat;
  const double M1, M2;
  double v;
  const Quaternions.Quaternion S_chi1, S_chi2;
  double rfrak_chi1_x, rfrak_chi1_y, rfrak_chi2_x, rfrak_chi2_y, rfrak_frame_x, rfrak_frame_y, rfrak_frame_z;
  const double M, delta, nu;
  Quaternions.Quaternion R, nHat, lambdaHat, ellHat, R_S1, R_S2, chiVec1, chiVec2;
  double chi1_n, chi1_lambda, chi1_ell, chi2_n, chi2_lambda, chi2_ell, S_ell, S_n, S_lambda, Sigma_ell, Sigma_n,
         Sigma_lambda, S1_ell, S1_n, S1_lambda, S2_ell, S2_n, S2_lambda;

public:
  WaveformModes_2p0PN(xHat_i, yHat_i, zHat_i, M1_i, M2_i, v_i, S_chi1_i, S_chi2_i, rfrak_chi1_x_i, rfrak_chi1_y_i,
                      rfrak_chi2_x_i, rfrak_chi2_y_i, rfrak_frame_x_i, rfrak_frame_y_i, rfrak_frame_z_i) :
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
        S1_ell=M1**2*chi1_ell
        S1_n=M1**2*chi1_n
        S1_lambda=M1**2*chi1_lambda
        S2_ell=M2**2*chi2_ell
        S2_n=M2**2*chi2_n
        S2_lambda=M2**2*chi2_lambda
  { }

  using WaveformModes_Base::operator();

  std::vector<std::complex<double> > operator()(
    const double v_k,
    const double chi1_n_k, const double chi1_lambda_k, const double chi1_ell_k,
    const double chi2_n_k, const double chi2_lambda_k, const double chi2_ell_k)
  {
    v = v_k;
    chiVec1 = Quaternions::Quaternion(0., chi1_n_k, chi1_lambda_k, chi1_ell_k);
    chiVec2 = Quaternions::Quaternion(0., chi2_n_k, chi2_lambda_k, chi2_ell_k);

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
            S1_ell = M1**2*chi1_ell
            S1_n = M1**2*chi1_n
            S1_lambda = M1**2*chi1_lambda
            S2_ell = M2**2*chi2_ell
            S2_n = M2**2*chi2_n
            S2_lambda = M2**2*chi2_lambda

    std::vector<std::complex<double> > Modes(77);
    std::complex<double> Symm, Asymm;
    // (ell, m) = (2, +/- 0)
    Symm = rhOverM_coeff*(hHat_2_0_0 + hHat_spin_Symm_2_0_4*pow(v, 4));
    Asymm = rhOverM_coeff*pow(v, 2)*(hHat_spin_Asymm_2_0_2 + hHat_spin_Asymm_2_0_4*pow(v, 2));
    Modes[2] = Symm + Asymm;
    // (ell, m) = (2, +/- 1)
    Symm = rhOverM_coeff*v*(hHat_2_1_1 + v*(hHat_spin_Symm_2_1_2 + v*(hHat_2_1_3 + v*(hHat_2_1_4 +
      hHat_spin_Symm_2_1_4))));
    Asymm = rhOverM_coeff*pow(v, 3)*(hHat_spin_Asymm_2_1_3 + hHat_spin_Asymm_2_1_4*v);
    Modes[3] = Symm + Asymm;
    Modes[1] = std::conj(Symm - Asymm);
    // (ell, m) = (2, +/- 2)
    Symm = rhOverM_coeff*(hHat_2_2_0 + pow(v, 2)*(hHat_2_2_2 + v*(hHat_2_2_3 + hHat_spin_Symm_2_2_3 + v*(hHat_2_2_4 +
      hHat_spin_Symm_2_2_4))));
    Asymm = rhOverM_coeff*pow(v, 2)*(hHat_spin_Asymm_2_2_2 + hHat_spin_Asymm_2_2_4*pow(v, 2));
    Modes[4] = Symm + Asymm;
    Modes[0] = std::conj(Symm - Asymm);
    // (ell, m) = (3, +/- 0)
    Symm = 0;
    Asymm = hHat_spin_Asymm_3_0_4*rhOverM_coeff*pow(v, 4);
    Modes[8] = Symm + Asymm;
    // (ell, m) = (3, +/- 1)
    Symm = rhOverM_coeff*v*(hHat_3_1_1 + pow(v, 2)*(hHat_3_1_3 + v*(hHat_3_1_4 + hHat_spin_Symm_3_1_4)));
    Asymm = hHat_spin_Asymm_3_1_3*rhOverM_coeff*pow(v, 3);
    Modes[9] = Symm + Asymm;
    Modes[7] = -std::conj(Symm - Asymm);
    // (ell, m) = (3, +/- 2)
    Symm = rhOverM_coeff*pow(v, 2)*(hHat_3_2_2 + v*(hHat_3_2_4*v + hHat_spin_Symm_3_2_3));
    Asymm = hHat_spin_Asymm_3_2_4*rhOverM_coeff*pow(v, 4);
    Modes[10] = Symm + Asymm;
    Modes[6] = -std::conj(Symm - Asymm);
    // (ell, m) = (3, +/- 3)
    Symm = rhOverM_coeff*v*(hHat_3_3_1 + pow(v, 2)*(hHat_3_3_3 + v*(hHat_3_3_4 + hHat_spin_Symm_3_3_4)));
    Asymm = hHat_spin_Asymm_3_3_3*rhOverM_coeff*pow(v, 3);
    Modes[11] = Symm + Asymm;
    Modes[5] = -std::conj(Symm - Asymm);
    // (ell, m) = (4, +/- 0)
    Symm = hHat_4_0_0*rhOverM_coeff;
    Asymm = hHat_spin_Asymm_4_0_4*rhOverM_coeff*pow(v, 4);
    Modes[16] = Symm + Asymm;
    // (ell, m) = (4, +/- 1)
    Symm = rhOverM_coeff*pow(v, 3)*(hHat_4_1_3 + hHat_spin_Symm_4_1_4*v);
    Asymm = 0;
    Modes[17] = Symm + Asymm;
    Modes[15] = std::conj(Symm - Asymm);
    // (ell, m) = (4, +/- 2)
    Symm = rhOverM_coeff*pow(v, 2)*(hHat_4_2_2 + hHat_4_2_4*pow(v, 2));
    Asymm = hHat_spin_Asymm_4_2_4*rhOverM_coeff*pow(v, 4);
    Modes[18] = Symm + Asymm;
    Modes[14] = std::conj(Symm - Asymm);
    // (ell, m) = (4, +/- 3)
    Symm = rhOverM_coeff*pow(v, 3)*(hHat_4_3_3 + hHat_spin_Symm_4_3_4*v);
    Asymm = 0;
    Modes[19] = Symm + Asymm;
    Modes[13] = std::conj(Symm - Asymm);
    // (ell, m) = (4, +/- 4)
    Symm = rhOverM_coeff*pow(v, 2)*(hHat_4_4_2 + hHat_4_4_4*pow(v, 2));
    Asymm = hHat_spin_Asymm_4_4_4*rhOverM_coeff*pow(v, 4);
    Modes[20] = Symm + Asymm;
    Modes[12] = std::conj(Symm - Asymm);
    // (ell, m) = (5, +/- 0)
    Symm = 0;
    Asymm = 0;
    Modes[26] = Symm + Asymm;
    // (ell, m) = (5, +/- 1)
    Symm = hHat_5_1_3*rhOverM_coeff*pow(v, 3);
    Asymm = 0;
    Modes[27] = Symm + Asymm;
    Modes[25] = -std::conj(Symm - Asymm);
    // (ell, m) = (5, +/- 2)
    Symm = hHat_5_2_4*rhOverM_coeff*pow(v, 4);
    Asymm = 0;
    Modes[28] = Symm + Asymm;
    Modes[24] = -std::conj(Symm - Asymm);
    // (ell, m) = (5, +/- 3)
    Symm = hHat_5_3_3*rhOverM_coeff*pow(v, 3);
    Asymm = 0;
    Modes[29] = Symm + Asymm;
    Modes[23] = -std::conj(Symm - Asymm);
    // (ell, m) = (5, +/- 4)
    Symm = hHat_5_4_4*rhOverM_coeff*pow(v, 4);
    Asymm = 0;
    Modes[30] = Symm + Asymm;
    Modes[22] = -std::conj(Symm - Asymm);
    // (ell, m) = (5, +/- 5)
    Symm = hHat_5_5_3*rhOverM_coeff*pow(v, 3);
    Asymm = 0;
    Modes[31] = Symm + Asymm;
    Modes[21] = -std::conj(Symm - Asymm);
    // (ell, m) = (6, +/- 0)
    Symm = 0;
    Asymm = 0;
    Modes[38] = Symm + Asymm;
    // (ell, m) = (6, +/- 1)
    Symm = 0;
    Asymm = 0;
    Modes[39] = Symm + Asymm;
    Modes[37] = std::conj(Symm - Asymm);
    // (ell, m) = (6, +/- 2)
    Symm = hHat_6_2_4*rhOverM_coeff*pow(v, 4);
    Asymm = 0;
    Modes[40] = Symm + Asymm;
    Modes[36] = std::conj(Symm - Asymm);
    // (ell, m) = (6, +/- 3)
    Symm = 0;
    Asymm = 0;
    Modes[41] = Symm + Asymm;
    Modes[35] = std::conj(Symm - Asymm);
    // (ell, m) = (6, +/- 4)
    Symm = hHat_6_4_4*rhOverM_coeff*pow(v, 4);
    Asymm = 0;
    Modes[42] = Symm + Asymm;
    Modes[34] = std::conj(Symm - Asymm);
    // (ell, m) = (6, +/- 5)
    Symm = 0;
    Asymm = 0;
    Modes[43] = Symm + Asymm;
    Modes[33] = std::conj(Symm - Asymm);
    // (ell, m) = (6, +/- 6)
    Symm = hHat_6_6_4*rhOverM_coeff*pow(v, 4);
    Asymm = 0;
    Modes[44] = Symm + Asymm;
    Modes[32] = std::conj(Symm - Asymm);
    // (ell, m) = (7, +/- 0)
    Symm = 0;
    Asymm = 0;
    Modes[52] = Symm + Asymm;
    // (ell, m) = (7, +/- 1)
    Symm = 0;
    Asymm = 0;
    Modes[53] = Symm + Asymm;
    Modes[51] = -std::conj(Symm - Asymm);
    // (ell, m) = (7, +/- 2)
    Symm = 0;
    Asymm = 0;
    Modes[54] = Symm + Asymm;
    Modes[50] = -std::conj(Symm - Asymm);
    // (ell, m) = (7, +/- 3)
    Symm = 0;
    Asymm = 0;
    Modes[55] = Symm + Asymm;
    Modes[49] = -std::conj(Symm - Asymm);
    // (ell, m) = (7, +/- 4)
    Symm = 0;
    Asymm = 0;
    Modes[56] = Symm + Asymm;
    Modes[48] = -std::conj(Symm - Asymm);
    // (ell, m) = (7, +/- 5)
    Symm = 0;
    Asymm = 0;
    Modes[57] = Symm + Asymm;
    Modes[47] = -std::conj(Symm - Asymm);
    // (ell, m) = (7, +/- 6)
    Symm = 0;
    Asymm = 0;
    Modes[58] = Symm + Asymm;
    Modes[46] = -std::conj(Symm - Asymm);
    // (ell, m) = (7, +/- 7)
    Symm = 0;
    Asymm = 0;
    Modes[59] = Symm + Asymm;
    Modes[45] = -std::conj(Symm - Asymm);
    // (ell, m) = (8, +/- 0)
    Symm = 0;
    Asymm = 0;
    Modes[68] = Symm + Asymm;
    // (ell, m) = (8, +/- 1)
    Symm = 0;
    Asymm = 0;
    Modes[69] = Symm + Asymm;
    Modes[67] = std::conj(Symm - Asymm);
    // (ell, m) = (8, +/- 2)
    Symm = 0;
    Asymm = 0;
    Modes[70] = Symm + Asymm;
    Modes[66] = std::conj(Symm - Asymm);
    // (ell, m) = (8, +/- 3)
    Symm = 0;
    Asymm = 0;
    Modes[71] = Symm + Asymm;
    Modes[65] = std::conj(Symm - Asymm);
    // (ell, m) = (8, +/- 4)
    Symm = 0;
    Asymm = 0;
    Modes[72] = Symm + Asymm;
    Modes[64] = std::conj(Symm - Asymm);
    // (ell, m) = (8, +/- 5)
    Symm = 0;
    Asymm = 0;
    Modes[73] = Symm + Asymm;
    Modes[63] = std::conj(Symm - Asymm);
    // (ell, m) = (8, +/- 6)
    Symm = 0;
    Asymm = 0;
    Modes[74] = Symm + Asymm;
    Modes[62] = std::conj(Symm - Asymm);
    // (ell, m) = (8, +/- 7)
    Symm = 0;
    Asymm = 0;
    Modes[75] = Symm + Asymm;
    Modes[61] = std::conj(Symm - Asymm);
    // (ell, m) = (8, +/- 8)
    Symm = 0;
    Asymm = 0;
    Modes[76] = Symm + Asymm;
    Modes[60] = std::conj(Symm - Asymm);

    return Modes;
  }

}; // class WaveformModes_2p0PN : public WaveformModes_Base

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
            a_ell_4 = S_n*(5.77777777777778*nu**2 + 14.75*nu + 1.5) + Sigma_n*delta*(2.83333333333333*nu**2 + 9.125*nu + 1.5)
            gamma_PN_0 = 1.00000000000000
            a_ell_2 = S_n*(-9.66666666666667*nu - 10.0) + Sigma_n*delta*(-4.5*nu - 6.0)
            gamma_PN_2 = 1.0 - 0.333333333333333*nu
            gamma_PN_3 = (1.66666666666667*S_ell + Sigma_ell*delta)/M**2
            gamma_PN_4 = 1.0 - 5.41666666666667*nu
            gamma_PN_5 = (S_ell*(0.888888888888889*nu + 3.33333333333333) + 2.0*Sigma_ell*delta)/M**2
            a_ell_0 = 7.0*S_n + 3.0*Sigma_n*delta
            gamma_PN_6 = 0.0123456790123457*nu**3 + 6.36111111111111*nu**2 - 2.98177812235564*nu + 1.0
            return ellHat*v**3/M + nHat*v**6*(a_ell_0 + v**2*(a_ell_2 + a_ell_4*v**2))*(gamma_PN_0 + v**2*(gamma_PN_2 + v*(gamma_PN_3 + v*(gamma_PN_4 + v*(gamma_PN_5 + v*(gamma_PN_6 + gamma_PN_7*v))))))/M**3

        def OrbitalAngularMomentum():
            L_4 = ellHat*(0.0416666666666667*nu**2 - 2.375*nu + 3.375)
            L_SO_5 = 0.0138888888888889*ellHat*(427.0*S_ell*nu - 693.0*S_ell + 210.0*Sigma_ell*delta*nu - 189.0*Sigma_ell*delta)/M**2 + 0.166666666666667*lambdaHat*(18.0*S_lambda*nu - 21.0*S_lambda + 8.0*Sigma_lambda*delta*nu - 3.0*Sigma_lambda*delta)/M**2 + 0.0416666666666667*nHat*(-19.0*S_n*nu + 33.0*S_n - 10.0*Sigma_n*delta*nu + 33.0*Sigma_n*delta)/M**2
            L_0 = ellHat
            L_2 = ellHat*(0.166666666666667*nu + 1.5)
            L_coeff = M**2*nu/v
            L_SO_3 = 0.166666666666667*ellHat*(-35.0*S_ell - 15.0*Sigma_ell*delta)/M**2 + lambdaHat*(-3.0*S_lambda - Sigma_lambda*delta)/M**2 + 0.5*nHat*(S_n + Sigma_n*delta)/M**2
            return L_coeff*(L_0 + v**2*(L_2 + v*(L_SO_3 + v*(L_4 + L_SO_5*v))))


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


class WaveformModes_2p5PN : public WaveformModes_Base {
private:
  const Quaternions.Quaternion xHat, yHat, zHat;
  const double M1, M2;
  double v;
  const Quaternions.Quaternion S_chi1, S_chi2;
  double rfrak_chi1_x, rfrak_chi1_y, rfrak_chi2_x, rfrak_chi2_y, rfrak_frame_x, rfrak_frame_y, rfrak_frame_z;
  const double M, delta, nu;
  Quaternions.Quaternion R, nHat, lambdaHat, ellHat, R_S1, R_S2, chiVec1, chiVec2;
  double chi1_n, chi1_lambda, chi1_ell, chi2_n, chi2_lambda, chi2_ell, S_ell, S_n, S_lambda, Sigma_ell, Sigma_n,
         Sigma_lambda, S1_ell, S1_n, S1_lambda, S2_ell, S2_n, S2_lambda;

public:
  WaveformModes_2p5PN(xHat_i, yHat_i, zHat_i, M1_i, M2_i, v_i, S_chi1_i, S_chi2_i, rfrak_chi1_x_i, rfrak_chi1_y_i,
                      rfrak_chi2_x_i, rfrak_chi2_y_i, rfrak_frame_x_i, rfrak_frame_y_i, rfrak_frame_z_i) :
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
        S1_ell=M1**2*chi1_ell
        S1_n=M1**2*chi1_n
        S1_lambda=M1**2*chi1_lambda
        S2_ell=M2**2*chi2_ell
        S2_n=M2**2*chi2_n
        S2_lambda=M2**2*chi2_lambda
  { }

  using WaveformModes_Base::operator();

  std::vector<std::complex<double> > operator()(
    const double v_k,
    const double chi1_n_k, const double chi1_lambda_k, const double chi1_ell_k,
    const double chi2_n_k, const double chi2_lambda_k, const double chi2_ell_k)
  {
    v = v_k;
    chiVec1 = Quaternions::Quaternion(0., chi1_n_k, chi1_lambda_k, chi1_ell_k);
    chiVec2 = Quaternions::Quaternion(0., chi2_n_k, chi2_lambda_k, chi2_ell_k);

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
            S1_ell = M1**2*chi1_ell
            S1_n = M1**2*chi1_n
            S1_lambda = M1**2*chi1_lambda
            S2_ell = M2**2*chi2_ell
            S2_n = M2**2*chi2_n
            S2_lambda = M2**2*chi2_lambda

    std::vector<std::complex<double> > Modes(77);
    std::complex<double> Symm, Asymm;
    // (ell, m) = (2, +/- 0)
    Symm = rhOverM_coeff*(hHat_2_0_0 + hHat_spin_Symm_2_0_4*pow(v, 4));
    Asymm = rhOverM_coeff*pow(v, 2)*(hHat_spin_Asymm_2_0_2 + hHat_spin_Asymm_2_0_4*pow(v, 2));
    Modes[2] = Symm + Asymm;
    // (ell, m) = (2, +/- 1)
    Symm = rhOverM_coeff*v*(hHat_2_1_1 + v*(hHat_spin_Symm_2_1_2 + v*(hHat_2_1_3 + v*(hHat_2_1_4 + hHat_2_1_5*v +
      hHat_spin_Symm_2_1_4))));
    Asymm = rhOverM_coeff*pow(v, 3)*(hHat_spin_Asymm_2_1_3 + hHat_spin_Asymm_2_1_4*v);
    Modes[3] = Symm + Asymm;
    Modes[1] = std::conj(Symm - Asymm);
    // (ell, m) = (2, +/- 2)
    Symm = rhOverM_coeff*(hHat_2_2_0 + pow(v, 2)*(hHat_2_2_2 + v*(hHat_2_2_3 + hHat_spin_Symm_2_2_3 + v*(hHat_2_2_4 +
      hHat_2_2_5*v + hHat_spin_Symm_2_2_4))));
    Asymm = rhOverM_coeff*pow(v, 2)*(hHat_spin_Asymm_2_2_2 + hHat_spin_Asymm_2_2_4*pow(v, 2));
    Modes[4] = Symm + Asymm;
    Modes[0] = std::conj(Symm - Asymm);
    // (ell, m) = (3, +/- 0)
    Symm = hHat_3_0_5*rhOverM_coeff*pow(v, 5);
    Asymm = hHat_spin_Asymm_3_0_4*rhOverM_coeff*pow(v, 4);
    Modes[8] = Symm + Asymm;
    // (ell, m) = (3, +/- 1)
    Symm = rhOverM_coeff*v*(hHat_3_1_1 + pow(v, 2)*(hHat_3_1_3 + v*(hHat_3_1_4 + hHat_3_1_5*v + hHat_spin_Symm_3_1_4)));
    Asymm = hHat_spin_Asymm_3_1_3*rhOverM_coeff*pow(v, 3);
    Modes[9] = Symm + Asymm;
    Modes[7] = -std::conj(Symm - Asymm);
    // (ell, m) = (3, +/- 2)
    Symm = rhOverM_coeff*pow(v, 2)*(hHat_3_2_2 + v*(hHat_spin_Symm_3_2_3 + v*(hHat_3_2_4 + hHat_3_2_5*v)));
    Asymm = hHat_spin_Asymm_3_2_4*rhOverM_coeff*pow(v, 4);
    Modes[10] = Symm + Asymm;
    Modes[6] = -std::conj(Symm - Asymm);
    // (ell, m) = (3, +/- 3)
    Symm = rhOverM_coeff*v*(hHat_3_3_1 + pow(v, 2)*(hHat_3_3_3 + v*(hHat_3_3_4 + hHat_3_3_5*v + hHat_spin_Symm_3_3_4)));
    Asymm = hHat_spin_Asymm_3_3_3*rhOverM_coeff*pow(v, 3);
    Modes[11] = Symm + Asymm;
    Modes[5] = -std::conj(Symm - Asymm);
    // (ell, m) = (4, +/- 0)
    Symm = hHat_4_0_0*rhOverM_coeff;
    Asymm = hHat_spin_Asymm_4_0_4*rhOverM_coeff*pow(v, 4);
    Modes[16] = Symm + Asymm;
    // (ell, m) = (4, +/- 1)
    Symm = rhOverM_coeff*pow(v, 3)*(hHat_4_1_3 + v*(hHat_4_1_5*v + hHat_spin_Symm_4_1_4));
    Asymm = 0;
    Modes[17] = Symm + Asymm;
    Modes[15] = std::conj(Symm - Asymm);
    // (ell, m) = (4, +/- 2)
    Symm = rhOverM_coeff*pow(v, 2)*(hHat_4_2_2 + pow(v, 2)*(hHat_4_2_4 + hHat_4_2_5*v));
    Asymm = hHat_spin_Asymm_4_2_4*rhOverM_coeff*pow(v, 4);
    Modes[18] = Symm + Asymm;
    Modes[14] = std::conj(Symm - Asymm);
    // (ell, m) = (4, +/- 3)
    Symm = rhOverM_coeff*pow(v, 3)*(hHat_4_3_3 + v*(hHat_4_3_5*v + hHat_spin_Symm_4_3_4));
    Asymm = 0;
    Modes[19] = Symm + Asymm;
    Modes[13] = std::conj(Symm - Asymm);
    // (ell, m) = (4, +/- 4)
    Symm = rhOverM_coeff*pow(v, 2)*(hHat_4_4_2 + pow(v, 2)*(hHat_4_4_4 + hHat_4_4_5*v));
    Asymm = hHat_spin_Asymm_4_4_4*rhOverM_coeff*pow(v, 4);
    Modes[20] = Symm + Asymm;
    Modes[12] = std::conj(Symm - Asymm);
    // (ell, m) = (5, +/- 0)
    Symm = 0;
    Asymm = 0;
    Modes[26] = Symm + Asymm;
    // (ell, m) = (5, +/- 1)
    Symm = rhOverM_coeff*pow(v, 3)*(hHat_5_1_3 + hHat_5_1_5*pow(v, 2));
    Asymm = 0;
    Modes[27] = Symm + Asymm;
    Modes[25] = -std::conj(Symm - Asymm);
    // (ell, m) = (5, +/- 2)
    Symm = hHat_5_2_4*rhOverM_coeff*pow(v, 4);
    Asymm = 0;
    Modes[28] = Symm + Asymm;
    Modes[24] = -std::conj(Symm - Asymm);
    // (ell, m) = (5, +/- 3)
    Symm = rhOverM_coeff*pow(v, 3)*(hHat_5_3_3 + hHat_5_3_5*pow(v, 2));
    Asymm = 0;
    Modes[29] = Symm + Asymm;
    Modes[23] = -std::conj(Symm - Asymm);
    // (ell, m) = (5, +/- 4)
    Symm = hHat_5_4_4*rhOverM_coeff*pow(v, 4);
    Asymm = 0;
    Modes[30] = Symm + Asymm;
    Modes[22] = -std::conj(Symm - Asymm);
    // (ell, m) = (5, +/- 5)
    Symm = rhOverM_coeff*pow(v, 3)*(hHat_5_5_3 + hHat_5_5_5*pow(v, 2));
    Asymm = 0;
    Modes[31] = Symm + Asymm;
    Modes[21] = -std::conj(Symm - Asymm);
    // (ell, m) = (6, +/- 0)
    Symm = 0;
    Asymm = 0;
    Modes[38] = Symm + Asymm;
    // (ell, m) = (6, +/- 1)
    Symm = hHat_6_1_5*rhOverM_coeff*pow(v, 5);
    Asymm = 0;
    Modes[39] = Symm + Asymm;
    Modes[37] = std::conj(Symm - Asymm);
    // (ell, m) = (6, +/- 2)
    Symm = hHat_6_2_4*rhOverM_coeff*pow(v, 4);
    Asymm = 0;
    Modes[40] = Symm + Asymm;
    Modes[36] = std::conj(Symm - Asymm);
    // (ell, m) = (6, +/- 3)
    Symm = hHat_6_3_5*rhOverM_coeff*pow(v, 5);
    Asymm = 0;
    Modes[41] = Symm + Asymm;
    Modes[35] = std::conj(Symm - Asymm);
    // (ell, m) = (6, +/- 4)
    Symm = hHat_6_4_4*rhOverM_coeff*pow(v, 4);
    Asymm = 0;
    Modes[42] = Symm + Asymm;
    Modes[34] = std::conj(Symm - Asymm);
    // (ell, m) = (6, +/- 5)
    Symm = hHat_6_5_5*rhOverM_coeff*pow(v, 5);
    Asymm = 0;
    Modes[43] = Symm + Asymm;
    Modes[33] = std::conj(Symm - Asymm);
    // (ell, m) = (6, +/- 6)
    Symm = hHat_6_6_4*rhOverM_coeff*pow(v, 4);
    Asymm = 0;
    Modes[44] = Symm + Asymm;
    Modes[32] = std::conj(Symm - Asymm);
    // (ell, m) = (7, +/- 0)
    Symm = 0;
    Asymm = 0;
    Modes[52] = Symm + Asymm;
    // (ell, m) = (7, +/- 1)
    Symm = hHat_7_1_5*rhOverM_coeff*pow(v, 5);
    Asymm = 0;
    Modes[53] = Symm + Asymm;
    Modes[51] = -std::conj(Symm - Asymm);
    // (ell, m) = (7, +/- 2)
    Symm = 0;
    Asymm = 0;
    Modes[54] = Symm + Asymm;
    Modes[50] = -std::conj(Symm - Asymm);
    // (ell, m) = (7, +/- 3)
    Symm = hHat_7_3_5*rhOverM_coeff*pow(v, 5);
    Asymm = 0;
    Modes[55] = Symm + Asymm;
    Modes[49] = -std::conj(Symm - Asymm);
    // (ell, m) = (7, +/- 4)
    Symm = 0;
    Asymm = 0;
    Modes[56] = Symm + Asymm;
    Modes[48] = -std::conj(Symm - Asymm);
    // (ell, m) = (7, +/- 5)
    Symm = hHat_7_5_5*rhOverM_coeff*pow(v, 5);
    Asymm = 0;
    Modes[57] = Symm + Asymm;
    Modes[47] = -std::conj(Symm - Asymm);
    // (ell, m) = (7, +/- 6)
    Symm = 0;
    Asymm = 0;
    Modes[58] = Symm + Asymm;
    Modes[46] = -std::conj(Symm - Asymm);
    // (ell, m) = (7, +/- 7)
    Symm = hHat_7_7_5*rhOverM_coeff*pow(v, 5);
    Asymm = 0;
    Modes[59] = Symm + Asymm;
    Modes[45] = -std::conj(Symm - Asymm);
    // (ell, m) = (8, +/- 0)
    Symm = 0;
    Asymm = 0;
    Modes[68] = Symm + Asymm;
    // (ell, m) = (8, +/- 1)
    Symm = 0;
    Asymm = 0;
    Modes[69] = Symm + Asymm;
    Modes[67] = std::conj(Symm - Asymm);
    // (ell, m) = (8, +/- 2)
    Symm = 0;
    Asymm = 0;
    Modes[70] = Symm + Asymm;
    Modes[66] = std::conj(Symm - Asymm);
    // (ell, m) = (8, +/- 3)
    Symm = 0;
    Asymm = 0;
    Modes[71] = Symm + Asymm;
    Modes[65] = std::conj(Symm - Asymm);
    // (ell, m) = (8, +/- 4)
    Symm = 0;
    Asymm = 0;
    Modes[72] = Symm + Asymm;
    Modes[64] = std::conj(Symm - Asymm);
    // (ell, m) = (8, +/- 5)
    Symm = 0;
    Asymm = 0;
    Modes[73] = Symm + Asymm;
    Modes[63] = std::conj(Symm - Asymm);
    // (ell, m) = (8, +/- 6)
    Symm = 0;
    Asymm = 0;
    Modes[74] = Symm + Asymm;
    Modes[62] = std::conj(Symm - Asymm);
    // (ell, m) = (8, +/- 7)
    Symm = 0;
    Asymm = 0;
    Modes[75] = Symm + Asymm;
    Modes[61] = std::conj(Symm - Asymm);
    // (ell, m) = (8, +/- 8)
    Symm = 0;
    Asymm = 0;
    Modes[76] = Symm + Asymm;
    Modes[60] = std::conj(Symm - Asymm);

    return Modes;
  }

}; // class WaveformModes_2p5PN : public WaveformModes_Base

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
            a_ell_4 = S_n*(5.77777777777778*nu**2 + 14.75*nu + 1.5) + Sigma_n*delta*(2.83333333333333*nu**2 + 9.125*nu + 1.5)
            gamma_PN_0 = 1.00000000000000
            a_ell_2 = S_n*(-9.66666666666667*nu - 10.0) + Sigma_n*delta*(-4.5*nu - 6.0)
            gamma_PN_2 = 1.0 - 0.333333333333333*nu
            gamma_PN_3 = (1.66666666666667*S_ell + Sigma_ell*delta)/M**2
            gamma_PN_4 = 1.0 - 5.41666666666667*nu
            gamma_PN_5 = (S_ell*(0.888888888888889*nu + 3.33333333333333) + 2.0*Sigma_ell*delta)/M**2
            a_ell_0 = 7.0*S_n + 3.0*Sigma_n*delta
            gamma_PN_6 = 0.0123456790123457*nu**3 + 6.36111111111111*nu**2 - 2.98177812235564*nu + 1.0
            return ellHat*v**3/M + nHat*v**6*(a_ell_0 + v**2*(a_ell_2 + a_ell_4*v**2))*(gamma_PN_0 + v**2*(gamma_PN_2 + v*(gamma_PN_3 + v*(gamma_PN_4 + v*(gamma_PN_5 + v*(gamma_PN_6 + gamma_PN_7*v))))))/M**3

        def OrbitalAngularMomentum():
            L_6 = ellHat*(0.00540123456790123*nu**3 + 1.29166666666667*nu**2 - 30.9797035925835*nu + 8.4375)
            L_4 = ellHat*(0.0416666666666667*nu**2 - 2.375*nu + 3.375)
            L_SO_5 = 0.0138888888888889*ellHat*(427.0*S_ell*nu - 693.0*S_ell + 210.0*Sigma_ell*delta*nu - 189.0*Sigma_ell*delta)/M**2 + 0.166666666666667*lambdaHat*(18.0*S_lambda*nu - 21.0*S_lambda + 8.0*Sigma_lambda*delta*nu - 3.0*Sigma_lambda*delta)/M**2 + 0.0416666666666667*nHat*(-19.0*S_n*nu + 33.0*S_n - 10.0*Sigma_n*delta*nu + 33.0*Sigma_n*delta)/M**2
            L_0 = ellHat
            L_2 = ellHat*(0.166666666666667*nu + 1.5)
            L_coeff = M**2*nu/v
            L_SO_3 = 0.166666666666667*ellHat*(-35.0*S_ell - 15.0*Sigma_ell*delta)/M**2 + lambdaHat*(-3.0*S_lambda - Sigma_lambda*delta)/M**2 + 0.5*nHat*(S_n + Sigma_n*delta)/M**2
            return L_coeff*(L_0 + v**2*(L_2 + v*(L_SO_3 + v*(L_4 + v*(L_6*v + L_SO_5)))))


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


class WaveformModes_3p0PN : public WaveformModes_Base {
private:
  const Quaternions.Quaternion xHat, yHat, zHat;
  const double M1, M2;
  double v;
  const Quaternions.Quaternion S_chi1, S_chi2;
  double rfrak_chi1_x, rfrak_chi1_y, rfrak_chi2_x, rfrak_chi2_y, rfrak_frame_x, rfrak_frame_y, rfrak_frame_z;
  const double M, delta, nu;
  Quaternions.Quaternion R, nHat, lambdaHat, ellHat, R_S1, R_S2, chiVec1, chiVec2;
  double chi1_n, chi1_lambda, chi1_ell, chi2_n, chi2_lambda, chi2_ell, S_ell, S_n, S_lambda, Sigma_ell, Sigma_n,
         Sigma_lambda, S1_ell, S1_n, S1_lambda, S2_ell, S2_n, S2_lambda, logv;

public:
  WaveformModes_3p0PN(xHat_i, yHat_i, zHat_i, M1_i, M2_i, v_i, S_chi1_i, S_chi2_i, rfrak_chi1_x_i, rfrak_chi1_y_i,
                      rfrak_chi2_x_i, rfrak_chi2_y_i, rfrak_frame_x_i, rfrak_frame_y_i, rfrak_frame_z_i) :
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
        S1_ell=M1**2*chi1_ell
        S1_n=M1**2*chi1_n
        S1_lambda=M1**2*chi1_lambda
        S2_ell=M2**2*chi2_ell
        S2_n=M2**2*chi2_n
        S2_lambda=M2**2*chi2_lambda
        logv=log(v)
  { }

  using WaveformModes_Base::operator();

  std::vector<std::complex<double> > operator()(
    const double v_k,
    const double chi1_n_k, const double chi1_lambda_k, const double chi1_ell_k,
    const double chi2_n_k, const double chi2_lambda_k, const double chi2_ell_k)
  {
    v = v_k;
    chiVec1 = Quaternions::Quaternion(0., chi1_n_k, chi1_lambda_k, chi1_ell_k);
    chiVec2 = Quaternions::Quaternion(0., chi2_n_k, chi2_lambda_k, chi2_ell_k);

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
            S1_ell = M1**2*chi1_ell
            S1_n = M1**2*chi1_n
            S1_lambda = M1**2*chi1_lambda
            S2_ell = M2**2*chi2_ell
            S2_n = M2**2*chi2_n
            S2_lambda = M2**2*chi2_lambda
            logv = log(v)

    std::vector<std::complex<double> > Modes(77);
    std::complex<double> Symm, Asymm;
    // (ell, m) = (2, +/- 0)
    Symm = rhOverM_coeff*(hHat_2_0_0 + hHat_spin_Symm_2_0_4*pow(v, 4));
    Asymm = rhOverM_coeff*pow(v, 2)*(hHat_spin_Asymm_2_0_2 + hHat_spin_Asymm_2_0_4*pow(v, 2));
    Modes[2] = Symm + Asymm;
    // (ell, m) = (2, +/- 1)
    Symm = rhOverM_coeff*v*(hHat_2_1_1 + v*(hHat_spin_Symm_2_1_2 + v*(hHat_2_1_3 + v*(hHat_2_1_4 + hHat_spin_Symm_2_1_4
      + v*(hHat_2_1_5 + hHat_2_1_6*v)))));
    Asymm = rhOverM_coeff*pow(v, 3)*(hHat_spin_Asymm_2_1_3 + hHat_spin_Asymm_2_1_4*v);
    Modes[3] = Symm + Asymm;
    Modes[1] = std::conj(Symm - Asymm);
    // (ell, m) = (2, +/- 2)
    Symm = rhOverM_coeff*(hHat_2_2_0 + pow(v, 2)*(hHat_2_2_2 + v*(hHat_2_2_3 + hHat_spin_Symm_2_2_3 + v*(hHat_2_2_4 +
      hHat_spin_Symm_2_2_4 + v*(hHat_2_2_5 + v*(hHat_2_2_6 + hHat_2_2_lnv_6*logv))))));
    Asymm = rhOverM_coeff*pow(v, 2)*(hHat_spin_Asymm_2_2_2 + hHat_spin_Asymm_2_2_4*pow(v, 2));
    Modes[4] = Symm + Asymm;
    Modes[0] = std::conj(Symm - Asymm);
    // (ell, m) = (3, +/- 0)
    Symm = hHat_3_0_5*rhOverM_coeff*pow(v, 5);
    Asymm = hHat_spin_Asymm_3_0_4*rhOverM_coeff*pow(v, 4);
    Modes[8] = Symm + Asymm;
    // (ell, m) = (3, +/- 1)
    Symm = rhOverM_coeff*v*(hHat_3_1_1 + pow(v, 2)*(hHat_3_1_3 + v*(hHat_3_1_4 + hHat_spin_Symm_3_1_4 + v*(hHat_3_1_5 +
      hHat_3_1_6*v))));
    Asymm = hHat_spin_Asymm_3_1_3*rhOverM_coeff*pow(v, 3);
    Modes[9] = Symm + Asymm;
    Modes[7] = -std::conj(Symm - Asymm);
    // (ell, m) = (3, +/- 2)
    Symm = rhOverM_coeff*pow(v, 2)*(hHat_3_2_2 + v*(hHat_spin_Symm_3_2_3 + v*(hHat_3_2_4 + v*(hHat_3_2_5 +
      hHat_3_2_6*v))));
    Asymm = hHat_spin_Asymm_3_2_4*rhOverM_coeff*pow(v, 4);
    Modes[10] = Symm + Asymm;
    Modes[6] = -std::conj(Symm - Asymm);
    // (ell, m) = (3, +/- 3)
    Symm = rhOverM_coeff*v*(hHat_3_3_1 + pow(v, 2)*(hHat_3_3_3 + v*(hHat_3_3_4 + hHat_spin_Symm_3_3_4 + v*(hHat_3_3_5 +
      hHat_3_3_6*v))));
    Asymm = hHat_spin_Asymm_3_3_3*rhOverM_coeff*pow(v, 3);
    Modes[11] = Symm + Asymm;
    Modes[5] = -std::conj(Symm - Asymm);
    // (ell, m) = (4, +/- 0)
    Symm = hHat_4_0_0*rhOverM_coeff;
    Asymm = hHat_spin_Asymm_4_0_4*rhOverM_coeff*pow(v, 4);
    Modes[16] = Symm + Asymm;
    // (ell, m) = (4, +/- 1)
    Symm = rhOverM_coeff*pow(v, 3)*(hHat_4_1_3 + v*(hHat_spin_Symm_4_1_4 + v*(hHat_4_1_5 + hHat_4_1_6*v)));
    Asymm = 0;
    Modes[17] = Symm + Asymm;
    Modes[15] = std::conj(Symm - Asymm);
    // (ell, m) = (4, +/- 2)
    Symm = rhOverM_coeff*pow(v, 2)*(hHat_4_2_2 + pow(v, 2)*(hHat_4_2_4 + v*(hHat_4_2_5 + hHat_4_2_6*v)));
    Asymm = hHat_spin_Asymm_4_2_4*rhOverM_coeff*pow(v, 4);
    Modes[18] = Symm + Asymm;
    Modes[14] = std::conj(Symm - Asymm);
    // (ell, m) = (4, +/- 3)
    Symm = rhOverM_coeff*pow(v, 3)*(hHat_4_3_3 + v*(hHat_spin_Symm_4_3_4 + v*(hHat_4_3_5 + hHat_4_3_6*v)));
    Asymm = 0;
    Modes[19] = Symm + Asymm;
    Modes[13] = std::conj(Symm - Asymm);
    // (ell, m) = (4, +/- 4)
    Symm = rhOverM_coeff*pow(v, 2)*(hHat_4_4_2 + pow(v, 2)*(hHat_4_4_4 + v*(hHat_4_4_5 + hHat_4_4_6*v)));
    Asymm = hHat_spin_Asymm_4_4_4*rhOverM_coeff*pow(v, 4);
    Modes[20] = Symm + Asymm;
    Modes[12] = std::conj(Symm - Asymm);
    // (ell, m) = (5, +/- 0)
    Symm = 0;
    Asymm = 0;
    Modes[26] = Symm + Asymm;
    // (ell, m) = (5, +/- 1)
    Symm = rhOverM_coeff*pow(v, 3)*(hHat_5_1_3 + pow(v, 2)*(hHat_5_1_5 + hHat_5_1_6*v));
    Asymm = 0;
    Modes[27] = Symm + Asymm;
    Modes[25] = -std::conj(Symm - Asymm);
    // (ell, m) = (5, +/- 2)
    Symm = rhOverM_coeff*pow(v, 4)*(hHat_5_2_4 + hHat_5_2_6*pow(v, 2));
    Asymm = 0;
    Modes[28] = Symm + Asymm;
    Modes[24] = -std::conj(Symm - Asymm);
    // (ell, m) = (5, +/- 3)
    Symm = rhOverM_coeff*pow(v, 3)*(hHat_5_3_3 + pow(v, 2)*(hHat_5_3_5 + hHat_5_3_6*v));
    Asymm = 0;
    Modes[29] = Symm + Asymm;
    Modes[23] = -std::conj(Symm - Asymm);
    // (ell, m) = (5, +/- 4)
    Symm = rhOverM_coeff*pow(v, 4)*(hHat_5_4_4 + hHat_5_4_6*pow(v, 2));
    Asymm = 0;
    Modes[30] = Symm + Asymm;
    Modes[22] = -std::conj(Symm - Asymm);
    // (ell, m) = (5, +/- 5)
    Symm = rhOverM_coeff*pow(v, 3)*(hHat_5_5_3 + pow(v, 2)*(hHat_5_5_5 + hHat_5_5_6*v));
    Asymm = 0;
    Modes[31] = Symm + Asymm;
    Modes[21] = -std::conj(Symm - Asymm);
    // (ell, m) = (6, +/- 0)
    Symm = 0;
    Asymm = 0;
    Modes[38] = Symm + Asymm;
    // (ell, m) = (6, +/- 1)
    Symm = hHat_6_1_5*rhOverM_coeff*pow(v, 5);
    Asymm = 0;
    Modes[39] = Symm + Asymm;
    Modes[37] = std::conj(Symm - Asymm);
    // (ell, m) = (6, +/- 2)
    Symm = rhOverM_coeff*pow(v, 4)*(hHat_6_2_4 + hHat_6_2_6*pow(v, 2));
    Asymm = 0;
    Modes[40] = Symm + Asymm;
    Modes[36] = std::conj(Symm - Asymm);
    // (ell, m) = (6, +/- 3)
    Symm = hHat_6_3_5*rhOverM_coeff*pow(v, 5);
    Asymm = 0;
    Modes[41] = Symm + Asymm;
    Modes[35] = std::conj(Symm - Asymm);
    // (ell, m) = (6, +/- 4)
    Symm = rhOverM_coeff*pow(v, 4)*(hHat_6_4_4 + hHat_6_4_6*pow(v, 2));
    Asymm = 0;
    Modes[42] = Symm + Asymm;
    Modes[34] = std::conj(Symm - Asymm);
    // (ell, m) = (6, +/- 5)
    Symm = hHat_6_5_5*rhOverM_coeff*pow(v, 5);
    Asymm = 0;
    Modes[43] = Symm + Asymm;
    Modes[33] = std::conj(Symm - Asymm);
    // (ell, m) = (6, +/- 6)
    Symm = rhOverM_coeff*pow(v, 4)*(hHat_6_6_4 + hHat_6_6_6*pow(v, 2));
    Asymm = 0;
    Modes[44] = Symm + Asymm;
    Modes[32] = std::conj(Symm - Asymm);
    // (ell, m) = (7, +/- 0)
    Symm = 0;
    Asymm = 0;
    Modes[52] = Symm + Asymm;
    // (ell, m) = (7, +/- 1)
    Symm = hHat_7_1_5*rhOverM_coeff*pow(v, 5);
    Asymm = 0;
    Modes[53] = Symm + Asymm;
    Modes[51] = -std::conj(Symm - Asymm);
    // (ell, m) = (7, +/- 2)
    Symm = hHat_7_2_6*rhOverM_coeff*pow(v, 6);
    Asymm = 0;
    Modes[54] = Symm + Asymm;
    Modes[50] = -std::conj(Symm - Asymm);
    // (ell, m) = (7, +/- 3)
    Symm = hHat_7_3_5*rhOverM_coeff*pow(v, 5);
    Asymm = 0;
    Modes[55] = Symm + Asymm;
    Modes[49] = -std::conj(Symm - Asymm);
    // (ell, m) = (7, +/- 4)
    Symm = hHat_7_4_6*rhOverM_coeff*pow(v, 6);
    Asymm = 0;
    Modes[56] = Symm + Asymm;
    Modes[48] = -std::conj(Symm - Asymm);
    // (ell, m) = (7, +/- 5)
    Symm = hHat_7_5_5*rhOverM_coeff*pow(v, 5);
    Asymm = 0;
    Modes[57] = Symm + Asymm;
    Modes[47] = -std::conj(Symm - Asymm);
    // (ell, m) = (7, +/- 6)
    Symm = hHat_7_6_6*rhOverM_coeff*pow(v, 6);
    Asymm = 0;
    Modes[58] = Symm + Asymm;
    Modes[46] = -std::conj(Symm - Asymm);
    // (ell, m) = (7, +/- 7)
    Symm = hHat_7_7_5*rhOverM_coeff*pow(v, 5);
    Asymm = 0;
    Modes[59] = Symm + Asymm;
    Modes[45] = -std::conj(Symm - Asymm);
    // (ell, m) = (8, +/- 0)
    Symm = 0;
    Asymm = 0;
    Modes[68] = Symm + Asymm;
    // (ell, m) = (8, +/- 1)
    Symm = 0;
    Asymm = 0;
    Modes[69] = Symm + Asymm;
    Modes[67] = std::conj(Symm - Asymm);
    // (ell, m) = (8, +/- 2)
    Symm = hHat_8_2_6*rhOverM_coeff*pow(v, 6);
    Asymm = 0;
    Modes[70] = Symm + Asymm;
    Modes[66] = std::conj(Symm - Asymm);
    // (ell, m) = (8, +/- 3)
    Symm = 0;
    Asymm = 0;
    Modes[71] = Symm + Asymm;
    Modes[65] = std::conj(Symm - Asymm);
    // (ell, m) = (8, +/- 4)
    Symm = hHat_8_4_6*rhOverM_coeff*pow(v, 6);
    Asymm = 0;
    Modes[72] = Symm + Asymm;
    Modes[64] = std::conj(Symm - Asymm);
    // (ell, m) = (8, +/- 5)
    Symm = 0;
    Asymm = 0;
    Modes[73] = Symm + Asymm;
    Modes[63] = std::conj(Symm - Asymm);
    // (ell, m) = (8, +/- 6)
    Symm = hHat_8_6_6*rhOverM_coeff*pow(v, 6);
    Asymm = 0;
    Modes[74] = Symm + Asymm;
    Modes[62] = std::conj(Symm - Asymm);
    // (ell, m) = (8, +/- 7)
    Symm = 0;
    Asymm = 0;
    Modes[75] = Symm + Asymm;
    Modes[61] = std::conj(Symm - Asymm);
    // (ell, m) = (8, +/- 8)
    Symm = hHat_8_8_6*rhOverM_coeff*pow(v, 6);
    Asymm = 0;
    Modes[76] = Symm + Asymm;
    Modes[60] = std::conj(Symm - Asymm);

    return Modes;
  }

}; // class WaveformModes_3p0PN : public WaveformModes_Base

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
            a_ell_4 = S_n*(5.77777777777778*nu**2 + 14.75*nu + 1.5) + Sigma_n*delta*(2.83333333333333*nu**2 + 9.125*nu + 1.5)
            gamma_PN_0 = 1.00000000000000
            a_ell_2 = S_n*(-9.66666666666667*nu - 10.0) + Sigma_n*delta*(-4.5*nu - 6.0)
            gamma_PN_2 = 1.0 - 0.333333333333333*nu
            gamma_PN_3 = (1.66666666666667*S_ell + Sigma_ell*delta)/M**2
            gamma_PN_4 = 1.0 - 5.41666666666667*nu
            gamma_PN_5 = (S_ell*(0.888888888888889*nu + 3.33333333333333) + 2.0*Sigma_ell*delta)/M**2
            a_ell_0 = 7.0*S_n + 3.0*Sigma_n*delta
            gamma_PN_6 = 0.0123456790123457*nu**3 + 6.36111111111111*nu**2 - 2.98177812235564*nu + 1.0
            return ellHat*v**3/M + nHat*v**6*(a_ell_0 + v**2*(a_ell_2 + a_ell_4*v**2))*(gamma_PN_0 + v**2*(gamma_PN_2 + v*(gamma_PN_3 + v*(gamma_PN_4 + v*(gamma_PN_5 + v*(gamma_PN_6 + gamma_PN_7*v))))))/M**3

        def OrbitalAngularMomentum():
            L_6 = ellHat*(0.00540123456790123*nu**3 + 1.29166666666667*nu**2 - 30.9797035925835*nu + 8.4375)
            L_4 = ellHat*(0.0416666666666667*nu**2 - 2.375*nu + 3.375)
            L_SO_5 = 0.0138888888888889*ellHat*(427.0*S_ell*nu - 693.0*S_ell + 210.0*Sigma_ell*delta*nu - 189.0*Sigma_ell*delta)/M**2 + 0.166666666666667*lambdaHat*(18.0*S_lambda*nu - 21.0*S_lambda + 8.0*Sigma_lambda*delta*nu - 3.0*Sigma_lambda*delta)/M**2 + 0.0416666666666667*nHat*(-19.0*S_n*nu + 33.0*S_n - 10.0*Sigma_n*delta*nu + 33.0*Sigma_n*delta)/M**2
            L_SO_7 = 0.0625*ellHat*(-29.0*S_ell*nu**2 + 1101.0*S_ell*nu - 405.0*S_ell - 15.0*Sigma_ell*delta*nu**2 + 468.0*Sigma_ell*delta*nu - 81.0*Sigma_ell*delta)/M**2 + 0.0416666666666667*lambdaHat*(-32.0*S_lambda*nu**2 + 2.0*S_lambda*nu - 174.0*S_lambda - 16.0*Sigma_lambda*delta*nu**2 - 79.0*Sigma_lambda*delta*nu - 12.0*Sigma_lambda*delta)/M**2 + 0.0208333333333333*nHat*(11.0*S_n*nu**2 - 1331.0*S_n*nu + 183.0*S_n + 5.0*Sigma_n*delta*nu**2 - 734.0*Sigma_n*delta*nu + 183.0*Sigma_n*delta)/M**2
            L_0 = ellHat
            L_2 = ellHat*(0.166666666666667*nu + 1.5)
            L_coeff = M**2*nu/v
            L_SO_3 = 0.166666666666667*ellHat*(-35.0*S_ell - 15.0*Sigma_ell*delta)/M**2 + lambdaHat*(-3.0*S_lambda - Sigma_lambda*delta)/M**2 + 0.5*nHat*(S_n + Sigma_n*delta)/M**2
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


class WaveformModes_3p5PN : public WaveformModes_Base {
private:
  const Quaternions.Quaternion xHat, yHat, zHat;
  const double M1, M2;
  double v;
  const Quaternions.Quaternion S_chi1, S_chi2;
  double rfrak_chi1_x, rfrak_chi1_y, rfrak_chi2_x, rfrak_chi2_y, rfrak_frame_x, rfrak_frame_y, rfrak_frame_z;
  const double M, delta, nu;
  Quaternions.Quaternion R, nHat, lambdaHat, ellHat, R_S1, R_S2, chiVec1, chiVec2;
  double chi1_n, chi1_lambda, chi1_ell, chi2_n, chi2_lambda, chi2_ell, S_ell, S_n, S_lambda, Sigma_ell, Sigma_n,
         Sigma_lambda, S1_ell, S1_n, S1_lambda, S2_ell, S2_n, S2_lambda, logv;

public:
  WaveformModes_3p5PN(xHat_i, yHat_i, zHat_i, M1_i, M2_i, v_i, S_chi1_i, S_chi2_i, rfrak_chi1_x_i, rfrak_chi1_y_i,
                      rfrak_chi2_x_i, rfrak_chi2_y_i, rfrak_frame_x_i, rfrak_frame_y_i, rfrak_frame_z_i) :
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
        S1_ell=M1**2*chi1_ell
        S1_n=M1**2*chi1_n
        S1_lambda=M1**2*chi1_lambda
        S2_ell=M2**2*chi2_ell
        S2_n=M2**2*chi2_n
        S2_lambda=M2**2*chi2_lambda
        logv=log(v)
  { }

  using WaveformModes_Base::operator();

  std::vector<std::complex<double> > operator()(
    const double v_k,
    const double chi1_n_k, const double chi1_lambda_k, const double chi1_ell_k,
    const double chi2_n_k, const double chi2_lambda_k, const double chi2_ell_k)
  {
    v = v_k;
    chiVec1 = Quaternions::Quaternion(0., chi1_n_k, chi1_lambda_k, chi1_ell_k);
    chiVec2 = Quaternions::Quaternion(0., chi2_n_k, chi2_lambda_k, chi2_ell_k);

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
            S1_ell = M1**2*chi1_ell
            S1_n = M1**2*chi1_n
            S1_lambda = M1**2*chi1_lambda
            S2_ell = M2**2*chi2_ell
            S2_n = M2**2*chi2_n
            S2_lambda = M2**2*chi2_lambda
            logv = log(v)

    std::vector<std::complex<double> > Modes(77);
    std::complex<double> Symm, Asymm;
    // (ell, m) = (2, +/- 0)
    Symm = rhOverM_coeff*(hHat_2_0_0 + hHat_spin_Symm_2_0_4*pow(v, 4));
    Asymm = rhOverM_coeff*pow(v, 2)*(hHat_spin_Asymm_2_0_2 + hHat_spin_Asymm_2_0_4*pow(v, 2));
    Modes[2] = Symm + Asymm;
    // (ell, m) = (2, +/- 1)
    Symm = rhOverM_coeff*v*(hHat_2_1_1 + v*(hHat_spin_Symm_2_1_2 + v*(hHat_2_1_3 + v*(hHat_2_1_4 + hHat_spin_Symm_2_1_4
      + v*(hHat_2_1_5 + hHat_2_1_6*v)))));
    Asymm = rhOverM_coeff*pow(v, 3)*(hHat_spin_Asymm_2_1_3 + hHat_spin_Asymm_2_1_4*v);
    Modes[3] = Symm + Asymm;
    Modes[1] = std::conj(Symm - Asymm);
    // (ell, m) = (2, +/- 2)
    Symm = rhOverM_coeff*(hHat_2_2_0 + pow(v, 2)*(hHat_2_2_2 + v*(hHat_2_2_3 + hHat_spin_Symm_2_2_3 + v*(hHat_2_2_4 +
      hHat_spin_Symm_2_2_4 + v*(hHat_2_2_5 + v*(hHat_2_2_6 + hHat_2_2_7*v + hHat_2_2_lnv_6*logv))))));
    Asymm = rhOverM_coeff*pow(v, 2)*(hHat_spin_Asymm_2_2_2 + hHat_spin_Asymm_2_2_4*pow(v, 2));
    Modes[4] = Symm + Asymm;
    Modes[0] = std::conj(Symm - Asymm);
    // (ell, m) = (3, +/- 0)
    Symm = hHat_3_0_5*rhOverM_coeff*pow(v, 5);
    Asymm = hHat_spin_Asymm_3_0_4*rhOverM_coeff*pow(v, 4);
    Modes[8] = Symm + Asymm;
    // (ell, m) = (3, +/- 1)
    Symm = rhOverM_coeff*v*(hHat_3_1_1 + pow(v, 2)*(hHat_3_1_3 + v*(hHat_3_1_4 + hHat_spin_Symm_3_1_4 + v*(hHat_3_1_5 +
      v*(hHat_3_1_6 + v*(hHat_3_1_7 + hHat_3_1_lnv_7*logv))))));
    Asymm = hHat_spin_Asymm_3_1_3*rhOverM_coeff*pow(v, 3);
    Modes[9] = Symm + Asymm;
    Modes[7] = -std::conj(Symm - Asymm);
    // (ell, m) = (3, +/- 2)
    Symm = rhOverM_coeff*pow(v, 2)*(hHat_3_2_2 + v*(hHat_spin_Symm_3_2_3 + v*(hHat_3_2_4 + v*(hHat_3_2_5 +
      hHat_3_2_6*v))));
    Asymm = hHat_spin_Asymm_3_2_4*rhOverM_coeff*pow(v, 4);
    Modes[10] = Symm + Asymm;
    Modes[6] = -std::conj(Symm - Asymm);
    // (ell, m) = (3, +/- 3)
    Symm = rhOverM_coeff*v*(hHat_3_3_1 + pow(v, 2)*(hHat_3_3_3 + v*(hHat_3_3_4 + hHat_spin_Symm_3_3_4 + v*(hHat_3_3_5 +
      v*(hHat_3_3_6 + v*(hHat_3_3_7 + hHat_3_3_lnv_7*logv))))));
    Asymm = hHat_spin_Asymm_3_3_3*rhOverM_coeff*pow(v, 3);
    Modes[11] = Symm + Asymm;
    Modes[5] = -std::conj(Symm - Asymm);
    // (ell, m) = (4, +/- 0)
    Symm = hHat_4_0_0*rhOverM_coeff;
    Asymm = hHat_spin_Asymm_4_0_4*rhOverM_coeff*pow(v, 4);
    Modes[16] = Symm + Asymm;
    // (ell, m) = (4, +/- 1)
    Symm = rhOverM_coeff*pow(v, 3)*(hHat_4_1_3 + v*(hHat_spin_Symm_4_1_4 + v*(hHat_4_1_5 + hHat_4_1_6*v)));
    Asymm = 0;
    Modes[17] = Symm + Asymm;
    Modes[15] = std::conj(Symm - Asymm);
    // (ell, m) = (4, +/- 2)
    Symm = rhOverM_coeff*pow(v, 2)*(hHat_4_2_2 + pow(v, 2)*(hHat_4_2_4 + v*(hHat_4_2_5 + hHat_4_2_6*v)));
    Asymm = hHat_spin_Asymm_4_2_4*rhOverM_coeff*pow(v, 4);
    Modes[18] = Symm + Asymm;
    Modes[14] = std::conj(Symm - Asymm);
    // (ell, m) = (4, +/- 3)
    Symm = rhOverM_coeff*pow(v, 3)*(hHat_4_3_3 + v*(hHat_spin_Symm_4_3_4 + v*(hHat_4_3_5 + hHat_4_3_6*v)));
    Asymm = 0;
    Modes[19] = Symm + Asymm;
    Modes[13] = std::conj(Symm - Asymm);
    // (ell, m) = (4, +/- 4)
    Symm = rhOverM_coeff*pow(v, 2)*(hHat_4_4_2 + pow(v, 2)*(hHat_4_4_4 + v*(hHat_4_4_5 + hHat_4_4_6*v)));
    Asymm = hHat_spin_Asymm_4_4_4*rhOverM_coeff*pow(v, 4);
    Modes[20] = Symm + Asymm;
    Modes[12] = std::conj(Symm - Asymm);
    // (ell, m) = (5, +/- 0)
    Symm = 0;
    Asymm = 0;
    Modes[26] = Symm + Asymm;
    // (ell, m) = (5, +/- 1)
    Symm = rhOverM_coeff*pow(v, 3)*(hHat_5_1_3 + pow(v, 2)*(hHat_5_1_5 + hHat_5_1_6*v));
    Asymm = 0;
    Modes[27] = Symm + Asymm;
    Modes[25] = -std::conj(Symm - Asymm);
    // (ell, m) = (5, +/- 2)
    Symm = rhOverM_coeff*pow(v, 4)*(hHat_5_2_4 + hHat_5_2_6*pow(v, 2));
    Asymm = 0;
    Modes[28] = Symm + Asymm;
    Modes[24] = -std::conj(Symm - Asymm);
    // (ell, m) = (5, +/- 3)
    Symm = rhOverM_coeff*pow(v, 3)*(hHat_5_3_3 + pow(v, 2)*(hHat_5_3_5 + hHat_5_3_6*v));
    Asymm = 0;
    Modes[29] = Symm + Asymm;
    Modes[23] = -std::conj(Symm - Asymm);
    // (ell, m) = (5, +/- 4)
    Symm = rhOverM_coeff*pow(v, 4)*(hHat_5_4_4 + hHat_5_4_6*pow(v, 2));
    Asymm = 0;
    Modes[30] = Symm + Asymm;
    Modes[22] = -std::conj(Symm - Asymm);
    // (ell, m) = (5, +/- 5)
    Symm = rhOverM_coeff*pow(v, 3)*(hHat_5_5_3 + pow(v, 2)*(hHat_5_5_5 + hHat_5_5_6*v));
    Asymm = 0;
    Modes[31] = Symm + Asymm;
    Modes[21] = -std::conj(Symm - Asymm);
    // (ell, m) = (6, +/- 0)
    Symm = 0;
    Asymm = 0;
    Modes[38] = Symm + Asymm;
    // (ell, m) = (6, +/- 1)
    Symm = hHat_6_1_5*rhOverM_coeff*pow(v, 5);
    Asymm = 0;
    Modes[39] = Symm + Asymm;
    Modes[37] = std::conj(Symm - Asymm);
    // (ell, m) = (6, +/- 2)
    Symm = rhOverM_coeff*pow(v, 4)*(hHat_6_2_4 + hHat_6_2_6*pow(v, 2));
    Asymm = 0;
    Modes[40] = Symm + Asymm;
    Modes[36] = std::conj(Symm - Asymm);
    // (ell, m) = (6, +/- 3)
    Symm = hHat_6_3_5*rhOverM_coeff*pow(v, 5);
    Asymm = 0;
    Modes[41] = Symm + Asymm;
    Modes[35] = std::conj(Symm - Asymm);
    // (ell, m) = (6, +/- 4)
    Symm = rhOverM_coeff*pow(v, 4)*(hHat_6_4_4 + hHat_6_4_6*pow(v, 2));
    Asymm = 0;
    Modes[42] = Symm + Asymm;
    Modes[34] = std::conj(Symm - Asymm);
    // (ell, m) = (6, +/- 5)
    Symm = hHat_6_5_5*rhOverM_coeff*pow(v, 5);
    Asymm = 0;
    Modes[43] = Symm + Asymm;
    Modes[33] = std::conj(Symm - Asymm);
    // (ell, m) = (6, +/- 6)
    Symm = rhOverM_coeff*pow(v, 4)*(hHat_6_6_4 + hHat_6_6_6*pow(v, 2));
    Asymm = 0;
    Modes[44] = Symm + Asymm;
    Modes[32] = std::conj(Symm - Asymm);
    // (ell, m) = (7, +/- 0)
    Symm = 0;
    Asymm = 0;
    Modes[52] = Symm + Asymm;
    // (ell, m) = (7, +/- 1)
    Symm = hHat_7_1_5*rhOverM_coeff*pow(v, 5);
    Asymm = 0;
    Modes[53] = Symm + Asymm;
    Modes[51] = -std::conj(Symm - Asymm);
    // (ell, m) = (7, +/- 2)
    Symm = hHat_7_2_6*rhOverM_coeff*pow(v, 6);
    Asymm = 0;
    Modes[54] = Symm + Asymm;
    Modes[50] = -std::conj(Symm - Asymm);
    // (ell, m) = (7, +/- 3)
    Symm = hHat_7_3_5*rhOverM_coeff*pow(v, 5);
    Asymm = 0;
    Modes[55] = Symm + Asymm;
    Modes[49] = -std::conj(Symm - Asymm);
    // (ell, m) = (7, +/- 4)
    Symm = hHat_7_4_6*rhOverM_coeff*pow(v, 6);
    Asymm = 0;
    Modes[56] = Symm + Asymm;
    Modes[48] = -std::conj(Symm - Asymm);
    // (ell, m) = (7, +/- 5)
    Symm = hHat_7_5_5*rhOverM_coeff*pow(v, 5);
    Asymm = 0;
    Modes[57] = Symm + Asymm;
    Modes[47] = -std::conj(Symm - Asymm);
    // (ell, m) = (7, +/- 6)
    Symm = hHat_7_6_6*rhOverM_coeff*pow(v, 6);
    Asymm = 0;
    Modes[58] = Symm + Asymm;
    Modes[46] = -std::conj(Symm - Asymm);
    // (ell, m) = (7, +/- 7)
    Symm = hHat_7_7_5*rhOverM_coeff*pow(v, 5);
    Asymm = 0;
    Modes[59] = Symm + Asymm;
    Modes[45] = -std::conj(Symm - Asymm);
    // (ell, m) = (8, +/- 0)
    Symm = 0;
    Asymm = 0;
    Modes[68] = Symm + Asymm;
    // (ell, m) = (8, +/- 1)
    Symm = 0;
    Asymm = 0;
    Modes[69] = Symm + Asymm;
    Modes[67] = std::conj(Symm - Asymm);
    // (ell, m) = (8, +/- 2)
    Symm = hHat_8_2_6*rhOverM_coeff*pow(v, 6);
    Asymm = 0;
    Modes[70] = Symm + Asymm;
    Modes[66] = std::conj(Symm - Asymm);
    // (ell, m) = (8, +/- 3)
    Symm = 0;
    Asymm = 0;
    Modes[71] = Symm + Asymm;
    Modes[65] = std::conj(Symm - Asymm);
    // (ell, m) = (8, +/- 4)
    Symm = hHat_8_4_6*rhOverM_coeff*pow(v, 6);
    Asymm = 0;
    Modes[72] = Symm + Asymm;
    Modes[64] = std::conj(Symm - Asymm);
    // (ell, m) = (8, +/- 5)
    Symm = 0;
    Asymm = 0;
    Modes[73] = Symm + Asymm;
    Modes[63] = std::conj(Symm - Asymm);
    // (ell, m) = (8, +/- 6)
    Symm = hHat_8_6_6*rhOverM_coeff*pow(v, 6);
    Asymm = 0;
    Modes[74] = Symm + Asymm;
    Modes[62] = std::conj(Symm - Asymm);
    // (ell, m) = (8, +/- 7)
    Symm = 0;
    Asymm = 0;
    Modes[75] = Symm + Asymm;
    Modes[61] = std::conj(Symm - Asymm);
    // (ell, m) = (8, +/- 8)
    Symm = hHat_8_8_6*rhOverM_coeff*pow(v, 6);
    Asymm = 0;
    Modes[76] = Symm + Asymm;
    Modes[60] = std::conj(Symm - Asymm);

    return Modes;
  }

}; // class WaveformModes_3p5PN : public WaveformModes_Base
