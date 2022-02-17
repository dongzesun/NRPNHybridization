# File produced automatically by PNCodeGen.ipynb
from scipy.integrate import solve_ivp
import numpy as np
from numpy import conjugate, dot, cross, exp, log, sqrt, pi
from numpy import euler_gamma as EulerGamma
from numba import jit, float64, int32, boolean
from numba.experimental import jitclass
from scipy.interpolate import InterpolatedUnivariateSpline as Spline
import quaternionic

@jit
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

@jit
def FrameFromAngularVelocityIntegrand(rfrak, Omega):
    rfrakMag = np.sqrt(rfrak[0] * rfrak[0] + rfrak[1] * rfrak[1] + rfrak[2] * rfrak[2])
    OmegaMag = np.sqrt(Omega[0] * Omega[0] + Omega[1] * Omega[1] + Omega[2] * Omega[2])
    # If the matrix is really close to the identity, return
    if rfrakMag < 1e-12:
        return np.array([Omega[0] / 2.0, Omega[1] / 2.0, Omega[2] / 2.0])
    # If the matrix is really close to singular, it's equivalent to the identity, so return
    if abs(np.sin(rfrakMag)) < 1e-12:
        return np.array([Omega[0] / 2.0, Omega[1] / 2.0, Omega[2] / 2.0])
    OmegaOver2 = np.array([Omega[0] / 2.0, Omega[1] / 2.0, Omega[2] / 2.0])
    rfrakHat = np.array([rfrak[0] / rfrakMag, rfrak[1] / rfrakMag, rfrak[2] / rfrakMag])
    return ((OmegaOver2 - rfrakHat * np.dot(rfrakHat, OmegaOver2)) * (rfrakMag / np.tan(rfrakMag))\
        + rfrakHat * np.dot(rfrakHat, OmegaOver2) + np.cross(OmegaOver2, rfrak))

spec=[("xHat", float64),("yHat", float64),("zHat", float64),("M1", float64),("M2", float64),("v", float64),("S_chi1", float64),("S_chi2", float64),("rfrak_chi1_x", float64),("rfrak_chi1_y", float64),("rfrak_chi2_x", float64),("rfrak_chi2_y", float64),("rfrak_frame_x", float64),("rfrak_frame_y", float64),("rfrak_frame_z", float64),("M", float64),("delta", float64),("nu", float64),("R", float64),("nHat", float64),("ellHat", float64),("R_S1", float64),("R_S2", float64),("chiVec1", float64),("chiVec2", float64),("chi1chi1", float64),("chi1chi2", float64),("chi2chi2", float64),("chi1_n", float64),("chi1_ell", float64),("chi2_n", float64),("chi2_ell", float64),("S_ell", float64),("S_n", float64),("Sigma_ell", float64),("Sigma_n", float64),("chi_s_ell", float64),("chi_a_ell", float64),("Fcal_coeff", float64),("Fcal_0", float64),("E_0", float64),("xHat_i", float64),("yHat_i", float64),("zHat_i", float64),("M1_i", float64),("M2_i", float64),("v_i", float64),("S_chi1_i", float64),("S_chi2_i", float64),("rfrak_chi1_x_i", float64),("rfrak_chi1_y_i", float64),("rfrak_chi2_x_i", float64),("rfrak_chi2_y_i", float64),("rfrak_frame_x_i", float64),("rfrak_frame_y_i", float64),("rfrak_frame_z_i", float64),("t_PN1",float64),("t_PN2",float64),("StepsPerOrbit",int32),("ForwardInTime",boolean),("tol",float64),("MinStep",float64),("EvolveSpin1",boolean),("EvolveSpin2",boolean),("terminal1",boolean),("terminal2",boolean),("NotForward",boolean)]
@jitclass(spec)
class TaylorTn_0PN_Q : 
    def __init__(self, xHat_i, yHat_i, zHat_i, M1_i, M2_i, v_i, S_chi1_i, S_chi2_i, rfrak_chi1_x_i, rfrak_chi1_y_i, rfrak_chi2_x_i, rfrak_chi2_y_i, rfrak_frame_x_i, rfrak_frame_y_i, rfrak_frame_z_i, t_PN1=False, t_PN2=False, StepsPerOrbit=32, ForwardInTime=True, tol=1e-12, MinStep=1e-7) :
        self.t_PN1,self.t_PN2,self.StepsPerOrbit,self.ForwardInTime,self.tol,self.MinStep=t_PN1, t_PN2, StepsPerOrbit, ForwardInTime, tol, MinStep
        terminal1=True
        terminal2=True
        NotForward=True
        self.xHat_i=xHat_i
        self.yHat_i=yHat_i
        self.zHat_i=zHat_i
        self.M1_i=M1_i
        self.M2_i=M2_i
        self.v_i=v_i
        self.S_chi1_i=S_chi1_i
        self.S_chi2_i=S_chi2_i
        self.rfrak_chi1_x_i=rfrak_chi1_x_i
        self.rfrak_chi1_y_i=rfrak_chi1_y_i
        self.rfrak_chi2_x_i=rfrak_chi2_x_i
        self.rfrak_chi2_y_i=rfrak_chi2_y_i
        self.rfrak_frame_x_i=rfrak_frame_x_i
        self.rfrak_frame_y_i=rfrak_frame_y_i
        self.rfrak_frame_z_i=rfrak_frame_z_i
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
        chi1chi1=dot(chiVec1.vector,chiVec1.vector)
        chi1chi2=dot(chiVec1.vector,chiVec2.vector)
        chi2chi2=dot(chiVec2.vector,chiVec2.vector)
        chi1_n=dot(chiVec1.vector,nHat.vector)
        chi1_ell=dot(chiVec1.vector,ellHat.vector)
        chi2_n=dot(chiVec2.vector,nHat.vector)
        chi2_ell=dot(chiVec2.vector,ellHat.vector)
        S_ell=M1**2*chi1_ell + M2**2*chi2_ell
        S_n=M1**2*chi1_n + M2**2*chi2_n
        Sigma_ell=M*(-M1*chi1_ell + M2*chi2_ell)
        Sigma_n=M*(-M1*chi1_n + M2*chi2_n)
        chi_s_ell=chi1_ell/2 + chi2_ell/2
        chi_a_ell=chi1_ell/2 - chi2_ell/2
        Fcal_coeff=32*nu**2*v**10/5
        Fcal_0=1
        E_0=1
        self.xHat=xHat
        self.yHat=yHat
        self.zHat=zHat
        self.M1=M1
        self.M2=M2
        self.v=v
        self.S_chi1=S_chi1
        self.S_chi2=S_chi2
        self.rfrak_chi1_x=rfrak_chi1_x
        self.rfrak_chi1_y=rfrak_chi1_y
        self.rfrak_chi2_x=rfrak_chi2_x
        self.rfrak_chi2_y=rfrak_chi2_y
        self.rfrak_frame_x=rfrak_frame_x
        self.rfrak_frame_y=rfrak_frame_y
        self.rfrak_frame_z=rfrak_frame_z
        self.M=M
        self.delta=delta
        self.nu=nu
        self.R=R
        self.nHat=nHat
        self.ellHat=ellHat
        self.R_S1=R_S1
        self.R_S2=R_S2
        self.chiVec1=chiVec1
        self.chiVec2=chiVec2
        self.chi1chi1=chi1chi1
        self.chi1chi2=chi1chi2
        self.chi2chi2=chi2chi2
        self.chi1_n=chi1_n
        self.chi1_ell=chi1_ell
        self.chi2_n=chi2_n
        self.chi2_ell=chi2_ell
        self.S_ell=S_ell
        self.S_n=S_n
        self.Sigma_ell=Sigma_ell
        self.Sigma_n=Sigma_n
        self.chi_s_ell=chi_s_ell
        self.chi_a_ell=chi_a_ell
        self.Fcal_coeff=Fcal_coeff
        self.Fcal_0=Fcal_0
        self.E_0=E_0
        Phi=0.0
        EvolveSpin1=(S_chi1*conjugate(S_chi1)).norm>1e-12
        EvolveSpin2=(S_chi2*conjugate(S_chi2)).norm>1e-12
        self.v, self.rfrak_chi1_x, self.rfrak_chi1_y, self.rfrak_chi2_x, self.rfrak_chi2_y=v, rfrak_chi1_x, rfrak_chi1_y, rfrak_chi2_x, rfrak_chi2_y
        self.rfrak_frame_x, self.rfrak_frame_y, self.rfrak_frame_z, self.Phi=rfrak_frame_x, rfrak_frame_y, rfrak_frame_z, Phi
        self.EvolveSpin1=EvolveSpin1
        self.EvolveSpin2=EvolveSpin2
        self.terminal1,self.terminal2,self.NotForward=terminal1,terminal2,NotForward

    def __Recalculate__(self,t, y):
        v = y[0]
        rfrak_chi1_x = y[1]
        rfrak_chi1_y = y[2]
        rfrak_chi2_x = y[3]
        rfrak_chi2_y = y[4]
        rfrak_frame_x = y[5]
        rfrak_frame_y = y[6]
        rfrak_frame_z = y[7]
        Phi = y[8]
        self.v, self.rfrak_chi1_x, self.rfrak_chi1_y, self.rfrak_chi2_x, self.rfrak_chi2_y=v, rfrak_chi1_x, rfrak_chi1_y, rfrak_chi2_x, rfrak_chi2_y
        self.rfrak_frame_x, self.rfrak_frame_y, self.rfrak_frame_z, self.Phi=rfrak_frame_x, rfrak_frame_y, rfrak_frame_z, Phi
        xHat=self.xHat
        yHat=self.yHat
        zHat=self.zHat
        M1=self.M1
        M2=self.M2
        v=self.v
        S_chi1=self.S_chi1
        S_chi2=self.S_chi2
        rfrak_chi1_x=self.rfrak_chi1_x
        rfrak_chi1_y=self.rfrak_chi1_y
        rfrak_chi2_x=self.rfrak_chi2_x
        rfrak_chi2_y=self.rfrak_chi2_y
        rfrak_frame_x=self.rfrak_frame_x
        rfrak_frame_y=self.rfrak_frame_y
        rfrak_frame_z=self.rfrak_frame_z
        M=self.M
        delta=self.delta
        nu=self.nu
        R=self.R
        nHat=self.nHat
        ellHat=self.ellHat
        R_S1=self.R_S1
        R_S2=self.R_S2
        chiVec1=self.chiVec1
        chiVec2=self.chiVec2
        chi1chi1=self.chi1chi1
        chi1chi2=self.chi1chi2
        chi2chi2=self.chi2chi2
        chi1_n=self.chi1_n
        chi1_ell=self.chi1_ell
        chi2_n=self.chi2_n
        chi2_ell=self.chi2_ell
        S_ell=self.S_ell
        S_n=self.S_n
        Sigma_ell=self.Sigma_ell
        Sigma_n=self.Sigma_n
        chi_s_ell=self.chi_s_ell
        chi_a_ell=self.chi_a_ell
        Fcal_coeff=self.Fcal_coeff
        Fcal_0=self.Fcal_0
        E_0=self.E_0
        R = exp(rfrak_frame_x*xHat + rfrak_frame_y*yHat + rfrak_frame_z*zHat)
        nHat = R*xHat*conjugate(R)
        ellHat = R*zHat*conjugate(R)
        R_S1 = exp(rfrak_chi1_x*xHat + rfrak_chi1_y*yHat)
        R_S2 = exp(rfrak_chi2_x*xHat + rfrak_chi2_y*yHat)
        chiVec1 = S_chi1*R_S1*zHat*conjugate(R_S1)*conjugate(S_chi1)
        chiVec2 = S_chi2*R_S2*zHat*conjugate(R_S2)*conjugate(S_chi2)
        chi1_n = dot(chiVec1.vector,nHat.vector)
        chi1_ell = dot(chiVec1.vector,ellHat.vector)
        chi2_n = dot(chiVec2.vector,nHat.vector)
        chi2_ell = dot(chiVec2.vector,ellHat.vector)
        S_ell = M1**2*chi1_ell + M2**2*chi2_ell
        S_n = M1**2*chi1_n + M2**2*chi2_n
        Sigma_ell = M*(-M1*chi1_ell + M2*chi2_ell)
        Sigma_n = M*(-M1*chi1_n + M2*chi2_n)
        chi_s_ell = chi1_ell/2 + chi2_ell/2
        chi_a_ell = chi1_ell/2 - chi2_ell/2
        Fcal_coeff = 32*nu**2*v**10/5
        self.R=R
        self.nHat=nHat
        self.ellHat=ellHat
        self.R_S1=R_S1
        self.R_S2=R_S2
        self.chiVec1=chiVec1
        self.chiVec2=chiVec2
        self.chi1_n=chi1_n
        self.chi1_ell=chi1_ell
        self.chi2_n=chi2_n
        self.chi2_ell=chi2_ell
        self.S_ell=S_ell
        self.S_n=S_n
        self.Sigma_ell=Sigma_ell
        self.Sigma_n=Sigma_n
        self.chi_s_ell=chi_s_ell
        self.chi_a_ell=chi_a_ell
        self.Fcal_coeff=Fcal_coeff

    def OmegaVec_chiVec_1(self):
        xHat=self.xHat
        yHat=self.yHat
        zHat=self.zHat
        M1=self.M1
        M2=self.M2
        v=self.v
        S_chi1=self.S_chi1
        S_chi2=self.S_chi2
        rfrak_chi1_x=self.rfrak_chi1_x
        rfrak_chi1_y=self.rfrak_chi1_y
        rfrak_chi2_x=self.rfrak_chi2_x
        rfrak_chi2_y=self.rfrak_chi2_y
        rfrak_frame_x=self.rfrak_frame_x
        rfrak_frame_y=self.rfrak_frame_y
        rfrak_frame_z=self.rfrak_frame_z
        M=self.M
        delta=self.delta
        nu=self.nu
        R=self.R
        nHat=self.nHat
        ellHat=self.ellHat
        R_S1=self.R_S1
        R_S2=self.R_S2
        chiVec1=self.chiVec1
        chiVec2=self.chiVec2
        chi1_n=self.chi1_n
        chi2_n=self.chi2_n
        S_n=self.S_n
        Sigma_n=self.Sigma_n
        Omega1_coeff = v**5/M
        return Omega1_coeff*ellHat*(-0.75*delta + 0.5*nu + 0.75)

    def OmegaVec_chiVec_2(self):
        xHat=self.xHat
        yHat=self.yHat
        zHat=self.zHat
        M1=self.M1
        M2=self.M2
        v=self.v
        S_chi1=self.S_chi1
        S_chi2=self.S_chi2
        rfrak_chi1_x=self.rfrak_chi1_x
        rfrak_chi1_y=self.rfrak_chi1_y
        rfrak_chi2_x=self.rfrak_chi2_x
        rfrak_chi2_y=self.rfrak_chi2_y
        rfrak_frame_x=self.rfrak_frame_x
        rfrak_frame_y=self.rfrak_frame_y
        rfrak_frame_z=self.rfrak_frame_z
        M=self.M
        delta=self.delta
        nu=self.nu
        R=self.R
        nHat=self.nHat
        ellHat=self.ellHat
        R_S1=self.R_S1
        R_S2=self.R_S2
        chiVec1=self.chiVec1
        chiVec2=self.chiVec2
        chi1_n=self.chi1_n
        chi2_n=self.chi2_n
        S_n=self.S_n
        Sigma_n=self.Sigma_n
        Omega2_coeff = v**5/M
        return Omega2_coeff*ellHat*(0.75*delta + 0.5*nu + 0.75)

    def OmegaVec(self):
        xHat=self.xHat
        yHat=self.yHat
        zHat=self.zHat
        M1=self.M1
        M2=self.M2
        v=self.v
        S_chi1=self.S_chi1
        S_chi2=self.S_chi2
        rfrak_chi1_x=self.rfrak_chi1_x
        rfrak_chi1_y=self.rfrak_chi1_y
        rfrak_chi2_x=self.rfrak_chi2_x
        rfrak_chi2_y=self.rfrak_chi2_y
        rfrak_frame_x=self.rfrak_frame_x
        rfrak_frame_y=self.rfrak_frame_y
        rfrak_frame_z=self.rfrak_frame_z
        M=self.M
        delta=self.delta
        nu=self.nu
        R=self.R
        nHat=self.nHat
        ellHat=self.ellHat
        R_S1=self.R_S1
        R_S2=self.R_S2
        chiVec1=self.chiVec1
        chiVec2=self.chiVec2
        chi1_n=self.chi1_n
        chi2_n=self.chi2_n
        S_n=self.S_n
        Sigma_n=self.Sigma_n
        gamma_PN_0 = 1.00000000000000
        a_ell_0 = 7.0*S_n + 3.0*Sigma_n*delta
        return ellHat*v**3/M + a_ell_0*gamma_PN_0*nHat*v**6/M**3

    def OrbitalAngularMomentum(self):
        xHat=self.xHat
        yHat=self.yHat
        zHat=self.zHat
        M1=self.M1
        M2=self.M2
        v=self.v
        S_chi1=self.S_chi1
        S_chi2=self.S_chi2
        rfrak_chi1_x=self.rfrak_chi1_x
        rfrak_chi1_y=self.rfrak_chi1_y
        rfrak_chi2_x=self.rfrak_chi2_x
        rfrak_chi2_y=self.rfrak_chi2_y
        rfrak_frame_x=self.rfrak_frame_x
        rfrak_frame_y=self.rfrak_frame_y
        rfrak_frame_z=self.rfrak_frame_z
        M=self.M
        delta=self.delta
        nu=self.nu
        R=self.R
        nHat=self.nHat
        ellHat=self.ellHat
        R_S1=self.R_S1
        R_S2=self.R_S2
        chiVec1=self.chiVec1
        chiVec2=self.chiVec2
        chi1_n=self.chi1_n
        chi2_n=self.chi2_n
        S_n=self.S_n
        Sigma_n=self.Sigma_n
        L_0 = ellHat
        L_coeff = M**2*nu/v
        return L_0*L_coeff

               

    #Evolve PN
    def __TaylorT1__(self,t, y):
        self.__Recalculate__(t, y)
        xHat=self.xHat
        yHat=self.yHat
        zHat=self.zHat
        M1=self.M1
        M2=self.M2
        v=self.v
        S_chi1=self.S_chi1
        S_chi2=self.S_chi2
        rfrak_chi1_x=self.rfrak_chi1_x
        rfrak_chi1_y=self.rfrak_chi1_y
        rfrak_chi2_x=self.rfrak_chi2_x
        rfrak_chi2_y=self.rfrak_chi2_y
        rfrak_frame_x=self.rfrak_frame_x
        rfrak_frame_y=self.rfrak_frame_y
        rfrak_frame_z=self.rfrak_frame_z
        M=self.M
        delta=self.delta
        nu=self.nu
        R=self.R
        nHat=self.nHat
        ellHat=self.ellHat
        R_S1=self.R_S1
        R_S2=self.R_S2
        chiVec1=self.chiVec1
        chiVec2=self.chiVec2
        chi1chi1=self.chi1chi1
        chi1chi2=self.chi1chi2
        chi2chi2=self.chi2chi2
        chi1_n=self.chi1_n
        chi1_ell=self.chi1_ell
        chi2_n=self.chi2_n
        chi2_ell=self.chi2_ell
        S_ell=self.S_ell
        S_n=self.S_n
        Sigma_ell=self.Sigma_ell
        Sigma_n=self.Sigma_n
        chi_s_ell=self.chi_s_ell
        chi_a_ell=self.chi_a_ell
        Fcal_coeff=self.Fcal_coeff
        Fcal_0=self.Fcal_0
        E_0=self.E_0
        v, rfrak_chi1_x, rfrak_chi1_y, rfrak_chi2_x, rfrak_chi2_y=self.v, self.rfrak_chi1_x, self.rfrak_chi1_y, self.rfrak_chi2_x, self.rfrak_chi2_y
        rfrak_frame_x, rfrak_frame_y, rfrak_frame_z, Phi=self.rfrak_frame_x, self.rfrak_frame_y, self.rfrak_frame_z, self.Phi
        dvdt_T1=0
        if v>=1.0 and self.NotForward:
            print("Beyond domain of PN validity, this is a good way to terminate.")
            self.terminal1=False
        Flux = Fcal_0*Fcal_coeff
        dEdv = -E_0*M*nu*v
        Absorption = 0
        dvdt_T1 = (-Absorption - Flux)/dEdv
        if dvdt_T1<1.0e-12 and self.NotForward:
            print("v is decreasing, which is not an uncommon way to stop.")
            self.terminal2=False
        return self.CommonRHS(dvdt_T1, y)

    def CommonRHS(self,dvdt, y):
        dydt=np.zeros(9)
        rfrak_frame=np.zeros(3)
        rfrak_frame[0] = y[5]
        rfrak_frame[1] = y[6]
        rfrak_frame[2] = y[7]
        rfrakdot_frame = FrameFromAngularVelocityIntegrand(rfrak_frame, self.OmegaVec().vector)
        dydt[0] = dvdt
        if(self.EvolveSpin1):
            dydt[1], dydt[2]= FrameFromAngularVelocity_2D_Integrand(y[1], y[2],\
                (self.S_chi1.inverse*self.OmegaVec_chiVec_1()*self.S_chi1).vector)
        else:
            dydt[1] = 0.0
            dydt[2] = 0.0
        if(self.EvolveSpin2):
            dydt[3], dydt[4] = FrameFromAngularVelocity_2D_Integrand(y[3], y[4],\
                (self.S_chi2.inverse*self.OmegaVec_chiVec_2()*self.S_chi2).vector)
        else:
            dydt[3] = 0.0
            dydt[4] = 0.0
        dydt[5] = rfrakdot_frame[0]
        dydt[6] = rfrakdot_frame[1]
        dydt[7] = rfrakdot_frame[2]
        dydt[8] = self.v*self.v*self.v/self.M

        return dydt

    # Prepare for integration
    def Evolution(self): 
        v_i,rfrak_chi1_x_i,rfrak_chi1_y_i,rfrak_chi2_x_i,rfrak_chi2_y_i,rfrak_frame_x_i,rfrak_frame_y_i,rfrak_frame_z_i\
            =self.v_i,self.rfrak_chi1_x_i,self.rfrak_chi1_y_i,self.rfrak_chi2_x_i,self.rfrak_chi2_y_i,self.rfrak_frame_x_i,self.rfrak_frame_y_i,self.rfrak_frame_z_i
        def terminate(t,y):
            return 1.0*self.terminal1*self.terminal2
        terminate.terminal=True
        TMerger=5.0/(256.0*self.nu*v_i**8)
        TEnd=TMerger
        if self.t_PN2:
            TEnd=self.t_PN2
        time=[0.0]
        while time[-1]<TEnd and 2*self.M*(256*self.nu*(TMerger-time[-1])/5)**(3/8)/self.StepsPerOrbit>self.MinStep:
            time.append(time[-1]+2*self.M*(256*self.nu*(TMerger-time[-1])/5)**(3/8)/self.StepsPerOrbit)
        time=np.delete(time, -1)
        
        # Integrate
        yy=solve_ivp(self.__TaylorT1__, [time[0],time[-1]], [v_i,rfrak_chi1_x_i,\
            rfrak_chi1_y_i,rfrak_chi2_x_i,rfrak_chi2_y_i,rfrak_frame_x_i,\
            rfrak_frame_y_i,rfrak_frame_z_i,0.0], method='DOP853',\
            t_eval=time, dense_output=True, events=terminate, rtol=self.tol, atol=self.tol)           
        if self.ForwardInTime:
            self.NotForward=False
            time=[0.0]
            TStart=-3*TMerger
            if self.t_PN1:
                TStart=self.t_PN1
            while time[-1]>TStart:
                time.append(time[-1]-2*self.M*(256*self.nu*(TMerger-time[-1])/5)**(3/8)/self.StepsPerOrbit)
            yyForward=solve_ivp(self.__TaylorT1__, [time[0],time[-1]], [v_i,rfrak_chi1_x_i,\
                rfrak_chi1_y_i,rfrak_chi2_x_i,rfrak_chi2_y_i,rfrak_frame_x_i,\
                rfrak_frame_y_i,rfrak_frame_z_i,0.0], method='DOP853',\
                t_eval=time, dense_output=True, rtol=self.tol, atol=self.tol)
            yy.t=np.append(yyForward.t[1:][::-1],yy.t)
            data=np.empty((9,len(yy.t)))
            for i in range(9):
                data[i]=np.append(yyForward.y[i][1:][::-1],yy.y[i])
            yy.y=data
            
        return yy

spec=[("xHat", float64),("yHat", float64),("zHat", float64),("M1", float64),("M2", float64),("v", float64),("S_chi1", float64),("S_chi2", float64),("rfrak_chi1_x", float64),("rfrak_chi1_y", float64),("rfrak_chi2_x", float64),("rfrak_chi2_y", float64),("rfrak_frame_x", float64),("rfrak_frame_y", float64),("rfrak_frame_z", float64),("M", float64),("delta", float64),("nu", float64),("R", float64),("nHat", float64),("ellHat", float64),("R_S1", float64),("R_S2", float64),("chiVec1", float64),("chiVec2", float64),("chi1chi1", float64),("chi1chi2", float64),("chi2chi2", float64),("chi1_n", float64),("chi1_ell", float64),("chi2_n", float64),("chi2_ell", float64),("S_ell", float64),("S_n", float64),("Sigma_ell", float64),("Sigma_n", float64),("chi_s_ell", float64),("chi_a_ell", float64),("Fcal_coeff", float64),("Fcal_0", float64),("E_0", float64),("xHat_i", float64),("yHat_i", float64),("zHat_i", float64),("M1_i", float64),("M2_i", float64),("v_i", float64),("S_chi1_i", float64),("S_chi2_i", float64),("rfrak_chi1_x_i", float64),("rfrak_chi1_y_i", float64),("rfrak_chi2_x_i", float64),("rfrak_chi2_y_i", float64),("rfrak_frame_x_i", float64),("rfrak_frame_y_i", float64),("rfrak_frame_z_i", float64),("t_PN1",float64),("t_PN2",float64),("StepsPerOrbit",int32),("ForwardInTime",boolean),("tol",float64),("MinStep",float64),("EvolveSpin1",boolean),("EvolveSpin2",boolean),("terminal1",boolean),("terminal2",boolean),("NotForward",boolean)]
@jitclass(spec)
class TaylorTn_0p50PN_Q : 
    def __init__(self, xHat_i, yHat_i, zHat_i, M1_i, M2_i, v_i, S_chi1_i, S_chi2_i, rfrak_chi1_x_i, rfrak_chi1_y_i, rfrak_chi2_x_i, rfrak_chi2_y_i, rfrak_frame_x_i, rfrak_frame_y_i, rfrak_frame_z_i, t_PN1=False, t_PN2=False, StepsPerOrbit=32, ForwardInTime=True, tol=1e-12, MinStep=1e-7) :
        self.t_PN1,self.t_PN2,self.StepsPerOrbit,self.ForwardInTime,self.tol,self.MinStep=t_PN1, t_PN2, StepsPerOrbit, ForwardInTime, tol, MinStep
        terminal1=True
        terminal2=True
        NotForward=True
        self.xHat_i=xHat_i
        self.yHat_i=yHat_i
        self.zHat_i=zHat_i
        self.M1_i=M1_i
        self.M2_i=M2_i
        self.v_i=v_i
        self.S_chi1_i=S_chi1_i
        self.S_chi2_i=S_chi2_i
        self.rfrak_chi1_x_i=rfrak_chi1_x_i
        self.rfrak_chi1_y_i=rfrak_chi1_y_i
        self.rfrak_chi2_x_i=rfrak_chi2_x_i
        self.rfrak_chi2_y_i=rfrak_chi2_y_i
        self.rfrak_frame_x_i=rfrak_frame_x_i
        self.rfrak_frame_y_i=rfrak_frame_y_i
        self.rfrak_frame_z_i=rfrak_frame_z_i
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
        chi1chi1=dot(chiVec1.vector,chiVec1.vector)
        chi1chi2=dot(chiVec1.vector,chiVec2.vector)
        chi2chi2=dot(chiVec2.vector,chiVec2.vector)
        chi1_n=dot(chiVec1.vector,nHat.vector)
        chi1_ell=dot(chiVec1.vector,ellHat.vector)
        chi2_n=dot(chiVec2.vector,nHat.vector)
        chi2_ell=dot(chiVec2.vector,ellHat.vector)
        S_ell=M1**2*chi1_ell + M2**2*chi2_ell
        S_n=M1**2*chi1_n + M2**2*chi2_n
        Sigma_ell=M*(-M1*chi1_ell + M2*chi2_ell)
        Sigma_n=M*(-M1*chi1_n + M2*chi2_n)
        chi_s_ell=chi1_ell/2 + chi2_ell/2
        chi_a_ell=chi1_ell/2 - chi2_ell/2
        Fcal_coeff=32*nu**2*v**10/5
        Fcal_0=1
        E_0=1
        self.xHat=xHat
        self.yHat=yHat
        self.zHat=zHat
        self.M1=M1
        self.M2=M2
        self.v=v
        self.S_chi1=S_chi1
        self.S_chi2=S_chi2
        self.rfrak_chi1_x=rfrak_chi1_x
        self.rfrak_chi1_y=rfrak_chi1_y
        self.rfrak_chi2_x=rfrak_chi2_x
        self.rfrak_chi2_y=rfrak_chi2_y
        self.rfrak_frame_x=rfrak_frame_x
        self.rfrak_frame_y=rfrak_frame_y
        self.rfrak_frame_z=rfrak_frame_z
        self.M=M
        self.delta=delta
        self.nu=nu
        self.R=R
        self.nHat=nHat
        self.ellHat=ellHat
        self.R_S1=R_S1
        self.R_S2=R_S2
        self.chiVec1=chiVec1
        self.chiVec2=chiVec2
        self.chi1chi1=chi1chi1
        self.chi1chi2=chi1chi2
        self.chi2chi2=chi2chi2
        self.chi1_n=chi1_n
        self.chi1_ell=chi1_ell
        self.chi2_n=chi2_n
        self.chi2_ell=chi2_ell
        self.S_ell=S_ell
        self.S_n=S_n
        self.Sigma_ell=Sigma_ell
        self.Sigma_n=Sigma_n
        self.chi_s_ell=chi_s_ell
        self.chi_a_ell=chi_a_ell
        self.Fcal_coeff=Fcal_coeff
        self.Fcal_0=Fcal_0
        self.E_0=E_0
        Phi=0.0
        EvolveSpin1=(S_chi1*conjugate(S_chi1)).norm>1e-12
        EvolveSpin2=(S_chi2*conjugate(S_chi2)).norm>1e-12
        self.v, self.rfrak_chi1_x, self.rfrak_chi1_y, self.rfrak_chi2_x, self.rfrak_chi2_y=v, rfrak_chi1_x, rfrak_chi1_y, rfrak_chi2_x, rfrak_chi2_y
        self.rfrak_frame_x, self.rfrak_frame_y, self.rfrak_frame_z, self.Phi=rfrak_frame_x, rfrak_frame_y, rfrak_frame_z, Phi
        self.EvolveSpin1=EvolveSpin1
        self.EvolveSpin2=EvolveSpin2
        self.terminal1,self.terminal2,self.NotForward=terminal1,terminal2,NotForward

    def __Recalculate__(self,t, y):
        v = y[0]
        rfrak_chi1_x = y[1]
        rfrak_chi1_y = y[2]
        rfrak_chi2_x = y[3]
        rfrak_chi2_y = y[4]
        rfrak_frame_x = y[5]
        rfrak_frame_y = y[6]
        rfrak_frame_z = y[7]
        Phi = y[8]
        self.v, self.rfrak_chi1_x, self.rfrak_chi1_y, self.rfrak_chi2_x, self.rfrak_chi2_y=v, rfrak_chi1_x, rfrak_chi1_y, rfrak_chi2_x, rfrak_chi2_y
        self.rfrak_frame_x, self.rfrak_frame_y, self.rfrak_frame_z, self.Phi=rfrak_frame_x, rfrak_frame_y, rfrak_frame_z, Phi
        xHat=self.xHat
        yHat=self.yHat
        zHat=self.zHat
        M1=self.M1
        M2=self.M2
        v=self.v
        S_chi1=self.S_chi1
        S_chi2=self.S_chi2
        rfrak_chi1_x=self.rfrak_chi1_x
        rfrak_chi1_y=self.rfrak_chi1_y
        rfrak_chi2_x=self.rfrak_chi2_x
        rfrak_chi2_y=self.rfrak_chi2_y
        rfrak_frame_x=self.rfrak_frame_x
        rfrak_frame_y=self.rfrak_frame_y
        rfrak_frame_z=self.rfrak_frame_z
        M=self.M
        delta=self.delta
        nu=self.nu
        R=self.R
        nHat=self.nHat
        ellHat=self.ellHat
        R_S1=self.R_S1
        R_S2=self.R_S2
        chiVec1=self.chiVec1
        chiVec2=self.chiVec2
        chi1chi1=self.chi1chi1
        chi1chi2=self.chi1chi2
        chi2chi2=self.chi2chi2
        chi1_n=self.chi1_n
        chi1_ell=self.chi1_ell
        chi2_n=self.chi2_n
        chi2_ell=self.chi2_ell
        S_ell=self.S_ell
        S_n=self.S_n
        Sigma_ell=self.Sigma_ell
        Sigma_n=self.Sigma_n
        chi_s_ell=self.chi_s_ell
        chi_a_ell=self.chi_a_ell
        Fcal_coeff=self.Fcal_coeff
        Fcal_0=self.Fcal_0
        E_0=self.E_0
        R = exp(rfrak_frame_x*xHat + rfrak_frame_y*yHat + rfrak_frame_z*zHat)
        nHat = R*xHat*conjugate(R)
        ellHat = R*zHat*conjugate(R)
        R_S1 = exp(rfrak_chi1_x*xHat + rfrak_chi1_y*yHat)
        R_S2 = exp(rfrak_chi2_x*xHat + rfrak_chi2_y*yHat)
        chiVec1 = S_chi1*R_S1*zHat*conjugate(R_S1)*conjugate(S_chi1)
        chiVec2 = S_chi2*R_S2*zHat*conjugate(R_S2)*conjugate(S_chi2)
        chi1_n = dot(chiVec1.vector,nHat.vector)
        chi1_ell = dot(chiVec1.vector,ellHat.vector)
        chi2_n = dot(chiVec2.vector,nHat.vector)
        chi2_ell = dot(chiVec2.vector,ellHat.vector)
        S_ell = M1**2*chi1_ell + M2**2*chi2_ell
        S_n = M1**2*chi1_n + M2**2*chi2_n
        Sigma_ell = M*(-M1*chi1_ell + M2*chi2_ell)
        Sigma_n = M*(-M1*chi1_n + M2*chi2_n)
        chi_s_ell = chi1_ell/2 + chi2_ell/2
        chi_a_ell = chi1_ell/2 - chi2_ell/2
        Fcal_coeff = 32*nu**2*v**10/5
        self.R=R
        self.nHat=nHat
        self.ellHat=ellHat
        self.R_S1=R_S1
        self.R_S2=R_S2
        self.chiVec1=chiVec1
        self.chiVec2=chiVec2
        self.chi1_n=chi1_n
        self.chi1_ell=chi1_ell
        self.chi2_n=chi2_n
        self.chi2_ell=chi2_ell
        self.S_ell=S_ell
        self.S_n=S_n
        self.Sigma_ell=Sigma_ell
        self.Sigma_n=Sigma_n
        self.chi_s_ell=chi_s_ell
        self.chi_a_ell=chi_a_ell
        self.Fcal_coeff=Fcal_coeff

    def OmegaVec_chiVec_1(self):
        xHat=self.xHat
        yHat=self.yHat
        zHat=self.zHat
        M1=self.M1
        M2=self.M2
        v=self.v
        S_chi1=self.S_chi1
        S_chi2=self.S_chi2
        rfrak_chi1_x=self.rfrak_chi1_x
        rfrak_chi1_y=self.rfrak_chi1_y
        rfrak_chi2_x=self.rfrak_chi2_x
        rfrak_chi2_y=self.rfrak_chi2_y
        rfrak_frame_x=self.rfrak_frame_x
        rfrak_frame_y=self.rfrak_frame_y
        rfrak_frame_z=self.rfrak_frame_z
        M=self.M
        delta=self.delta
        nu=self.nu
        R=self.R
        nHat=self.nHat
        ellHat=self.ellHat
        R_S1=self.R_S1
        R_S2=self.R_S2
        chiVec1=self.chiVec1
        chiVec2=self.chiVec2
        chi1_n=self.chi1_n
        chi2_n=self.chi2_n
        S_n=self.S_n
        Sigma_n=self.Sigma_n
        Omega1_coeff = v**5/M
        return Omega1_coeff*(ellHat*(-0.75*delta + 0.5*nu + 0.75) + nHat*v*(3.0*chi1_n*nu + 3.0*M2**2*chi2_n/M**2) - M2**2*chiVec2*v/M**2)

    def OmegaVec_chiVec_2(self):
        xHat=self.xHat
        yHat=self.yHat
        zHat=self.zHat
        M1=self.M1
        M2=self.M2
        v=self.v
        S_chi1=self.S_chi1
        S_chi2=self.S_chi2
        rfrak_chi1_x=self.rfrak_chi1_x
        rfrak_chi1_y=self.rfrak_chi1_y
        rfrak_chi2_x=self.rfrak_chi2_x
        rfrak_chi2_y=self.rfrak_chi2_y
        rfrak_frame_x=self.rfrak_frame_x
        rfrak_frame_y=self.rfrak_frame_y
        rfrak_frame_z=self.rfrak_frame_z
        M=self.M
        delta=self.delta
        nu=self.nu
        R=self.R
        nHat=self.nHat
        ellHat=self.ellHat
        R_S1=self.R_S1
        R_S2=self.R_S2
        chiVec1=self.chiVec1
        chiVec2=self.chiVec2
        chi1_n=self.chi1_n
        chi2_n=self.chi2_n
        S_n=self.S_n
        Sigma_n=self.Sigma_n
        Omega2_coeff = v**5/M
        return Omega2_coeff*(ellHat*(0.75*delta + 0.5*nu + 0.75) + nHat*v*(3.0*chi2_n*nu + 3.0*M1**2*chi1_n/M**2) - M1**2*chiVec1*v/M**2)

    def OmegaVec(self):
        xHat=self.xHat
        yHat=self.yHat
        zHat=self.zHat
        M1=self.M1
        M2=self.M2
        v=self.v
        S_chi1=self.S_chi1
        S_chi2=self.S_chi2
        rfrak_chi1_x=self.rfrak_chi1_x
        rfrak_chi1_y=self.rfrak_chi1_y
        rfrak_chi2_x=self.rfrak_chi2_x
        rfrak_chi2_y=self.rfrak_chi2_y
        rfrak_frame_x=self.rfrak_frame_x
        rfrak_frame_y=self.rfrak_frame_y
        rfrak_frame_z=self.rfrak_frame_z
        M=self.M
        delta=self.delta
        nu=self.nu
        R=self.R
        nHat=self.nHat
        ellHat=self.ellHat
        R_S1=self.R_S1
        R_S2=self.R_S2
        chiVec1=self.chiVec1
        chiVec2=self.chiVec2
        chi1_n=self.chi1_n
        chi2_n=self.chi2_n
        S_n=self.S_n
        Sigma_n=self.Sigma_n
        gamma_PN_0 = 1.00000000000000
        a_ell_0 = 7.0*S_n + 3.0*Sigma_n*delta
        return ellHat*v**3/M + a_ell_0*gamma_PN_0*nHat*v**6/M**3

    def OrbitalAngularMomentum(self):
        xHat=self.xHat
        yHat=self.yHat
        zHat=self.zHat
        M1=self.M1
        M2=self.M2
        v=self.v
        S_chi1=self.S_chi1
        S_chi2=self.S_chi2
        rfrak_chi1_x=self.rfrak_chi1_x
        rfrak_chi1_y=self.rfrak_chi1_y
        rfrak_chi2_x=self.rfrak_chi2_x
        rfrak_chi2_y=self.rfrak_chi2_y
        rfrak_frame_x=self.rfrak_frame_x
        rfrak_frame_y=self.rfrak_frame_y
        rfrak_frame_z=self.rfrak_frame_z
        M=self.M
        delta=self.delta
        nu=self.nu
        R=self.R
        nHat=self.nHat
        ellHat=self.ellHat
        R_S1=self.R_S1
        R_S2=self.R_S2
        chiVec1=self.chiVec1
        chiVec2=self.chiVec2
        chi1_n=self.chi1_n
        chi2_n=self.chi2_n
        S_n=self.S_n
        Sigma_n=self.Sigma_n
        L_0 = ellHat
        L_coeff = M**2*nu/v
        return L_0*L_coeff

               

    #Evolve PN
    def __TaylorT1__(self,t, y):
        self.__Recalculate__(t, y)
        xHat=self.xHat
        yHat=self.yHat
        zHat=self.zHat
        M1=self.M1
        M2=self.M2
        v=self.v
        S_chi1=self.S_chi1
        S_chi2=self.S_chi2
        rfrak_chi1_x=self.rfrak_chi1_x
        rfrak_chi1_y=self.rfrak_chi1_y
        rfrak_chi2_x=self.rfrak_chi2_x
        rfrak_chi2_y=self.rfrak_chi2_y
        rfrak_frame_x=self.rfrak_frame_x
        rfrak_frame_y=self.rfrak_frame_y
        rfrak_frame_z=self.rfrak_frame_z
        M=self.M
        delta=self.delta
        nu=self.nu
        R=self.R
        nHat=self.nHat
        ellHat=self.ellHat
        R_S1=self.R_S1
        R_S2=self.R_S2
        chiVec1=self.chiVec1
        chiVec2=self.chiVec2
        chi1chi1=self.chi1chi1
        chi1chi2=self.chi1chi2
        chi2chi2=self.chi2chi2
        chi1_n=self.chi1_n
        chi1_ell=self.chi1_ell
        chi2_n=self.chi2_n
        chi2_ell=self.chi2_ell
        S_ell=self.S_ell
        S_n=self.S_n
        Sigma_ell=self.Sigma_ell
        Sigma_n=self.Sigma_n
        chi_s_ell=self.chi_s_ell
        chi_a_ell=self.chi_a_ell
        Fcal_coeff=self.Fcal_coeff
        Fcal_0=self.Fcal_0
        E_0=self.E_0
        v, rfrak_chi1_x, rfrak_chi1_y, rfrak_chi2_x, rfrak_chi2_y=self.v, self.rfrak_chi1_x, self.rfrak_chi1_y, self.rfrak_chi2_x, self.rfrak_chi2_y
        rfrak_frame_x, rfrak_frame_y, rfrak_frame_z, Phi=self.rfrak_frame_x, self.rfrak_frame_y, self.rfrak_frame_z, self.Phi
        dvdt_T1=0
        if v>=1.0 and self.NotForward:
            print("Beyond domain of PN validity, this is a good way to terminate.")
            self.terminal1=False
        Flux = Fcal_0*Fcal_coeff
        dEdv = -E_0*M*nu*v
        Absorption = 0
        dvdt_T1 = (-Absorption - Flux)/dEdv
        if dvdt_T1<1.0e-12 and self.NotForward:
            print("v is decreasing, which is not an uncommon way to stop.")
            self.terminal2=False
        return self.CommonRHS(dvdt_T1, y)

    def CommonRHS(self,dvdt, y):
        dydt=np.zeros(9)
        rfrak_frame=np.zeros(3)
        rfrak_frame[0] = y[5]
        rfrak_frame[1] = y[6]
        rfrak_frame[2] = y[7]
        rfrakdot_frame = FrameFromAngularVelocityIntegrand(rfrak_frame, self.OmegaVec().vector)
        dydt[0] = dvdt
        if(self.EvolveSpin1):
            dydt[1], dydt[2]= FrameFromAngularVelocity_2D_Integrand(y[1], y[2],\
                (self.S_chi1.inverse*self.OmegaVec_chiVec_1()*self.S_chi1).vector)
        else:
            dydt[1] = 0.0
            dydt[2] = 0.0
        if(self.EvolveSpin2):
            dydt[3], dydt[4] = FrameFromAngularVelocity_2D_Integrand(y[3], y[4],\
                (self.S_chi2.inverse*self.OmegaVec_chiVec_2()*self.S_chi2).vector)
        else:
            dydt[3] = 0.0
            dydt[4] = 0.0
        dydt[5] = rfrakdot_frame[0]
        dydt[6] = rfrakdot_frame[1]
        dydt[7] = rfrakdot_frame[2]
        dydt[8] = self.v*self.v*self.v/self.M

        return dydt

    # Prepare for integration
    def Evolution(self): 
        v_i,rfrak_chi1_x_i,rfrak_chi1_y_i,rfrak_chi2_x_i,rfrak_chi2_y_i,rfrak_frame_x_i,rfrak_frame_y_i,rfrak_frame_z_i\
            =self.v_i,self.rfrak_chi1_x_i,self.rfrak_chi1_y_i,self.rfrak_chi2_x_i,self.rfrak_chi2_y_i,self.rfrak_frame_x_i,self.rfrak_frame_y_i,self.rfrak_frame_z_i
        def terminate(t,y):
            return 1.0*self.terminal1*self.terminal2
        terminate.terminal=True
        TMerger=5.0/(256.0*self.nu*v_i**8)
        TEnd=TMerger
        if self.t_PN2:
            TEnd=self.t_PN2
        time=[0.0]
        while time[-1]<TEnd and 2*self.M*(256*self.nu*(TMerger-time[-1])/5)**(3/8)/self.StepsPerOrbit>self.MinStep:
            time.append(time[-1]+2*self.M*(256*self.nu*(TMerger-time[-1])/5)**(3/8)/self.StepsPerOrbit)
        time=np.delete(time, -1)
        
        # Integrate
        yy=solve_ivp(self.__TaylorT1__, [time[0],time[-1]], [v_i,rfrak_chi1_x_i,\
            rfrak_chi1_y_i,rfrak_chi2_x_i,rfrak_chi2_y_i,rfrak_frame_x_i,\
            rfrak_frame_y_i,rfrak_frame_z_i,0.0], method='DOP853',\
            t_eval=time, dense_output=True, events=terminate, rtol=self.tol, atol=self.tol)           
        if self.ForwardInTime:
            self.NotForward=False
            time=[0.0]
            TStart=-3*TMerger
            if self.t_PN1:
                TStart=self.t_PN1
            while time[-1]>TStart:
                time.append(time[-1]-2*self.M*(256*self.nu*(TMerger-time[-1])/5)**(3/8)/self.StepsPerOrbit)
            yyForward=solve_ivp(self.__TaylorT1__, [time[0],time[-1]], [v_i,rfrak_chi1_x_i,\
                rfrak_chi1_y_i,rfrak_chi2_x_i,rfrak_chi2_y_i,rfrak_frame_x_i,\
                rfrak_frame_y_i,rfrak_frame_z_i,0.0], method='DOP853',\
                t_eval=time, dense_output=True, rtol=self.tol, atol=self.tol)
            yy.t=np.append(yyForward.t[1:][::-1],yy.t)
            data=np.empty((9,len(yy.t)))
            for i in range(9):
                data[i]=np.append(yyForward.y[i][1:][::-1],yy.y[i])
            yy.y=data
            
        return yy

spec=[("xHat", float64),("yHat", float64),("zHat", float64),("M1", float64),("M2", float64),("v", float64),("S_chi1", float64),("S_chi2", float64),("rfrak_chi1_x", float64),("rfrak_chi1_y", float64),("rfrak_chi2_x", float64),("rfrak_chi2_y", float64),("rfrak_frame_x", float64),("rfrak_frame_y", float64),("rfrak_frame_z", float64),("M", float64),("delta", float64),("nu", float64),("R", float64),("nHat", float64),("ellHat", float64),("R_S1", float64),("R_S2", float64),("chiVec1", float64),("chiVec2", float64),("chi1chi1", float64),("chi1chi2", float64),("chi2chi2", float64),("chi1_n", float64),("chi1_ell", float64),("chi2_n", float64),("chi2_ell", float64),("S_ell", float64),("S_n", float64),("Sigma_ell", float64),("Sigma_n", float64),("chi_s_ell", float64),("chi_a_ell", float64),("Fcal_coeff", float64),("Fcal_0", float64),("Fcal_2", float64),("E_0", float64),("E_2", float64),("xHat_i", float64),("yHat_i", float64),("zHat_i", float64),("M1_i", float64),("M2_i", float64),("v_i", float64),("S_chi1_i", float64),("S_chi2_i", float64),("rfrak_chi1_x_i", float64),("rfrak_chi1_y_i", float64),("rfrak_chi2_x_i", float64),("rfrak_chi2_y_i", float64),("rfrak_frame_x_i", float64),("rfrak_frame_y_i", float64),("rfrak_frame_z_i", float64),("t_PN1",float64),("t_PN2",float64),("StepsPerOrbit",int32),("ForwardInTime",boolean),("tol",float64),("MinStep",float64),("EvolveSpin1",boolean),("EvolveSpin2",boolean),("terminal1",boolean),("terminal2",boolean),("NotForward",boolean)]
@jitclass(spec)
class TaylorTn_1p0PN_Q : 
    def __init__(self, xHat_i, yHat_i, zHat_i, M1_i, M2_i, v_i, S_chi1_i, S_chi2_i, rfrak_chi1_x_i, rfrak_chi1_y_i, rfrak_chi2_x_i, rfrak_chi2_y_i, rfrak_frame_x_i, rfrak_frame_y_i, rfrak_frame_z_i, t_PN1=False, t_PN2=False, StepsPerOrbit=32, ForwardInTime=True, tol=1e-12, MinStep=1e-7) :
        self.t_PN1,self.t_PN2,self.StepsPerOrbit,self.ForwardInTime,self.tol,self.MinStep=t_PN1, t_PN2, StepsPerOrbit, ForwardInTime, tol, MinStep
        terminal1=True
        terminal2=True
        NotForward=True
        self.xHat_i=xHat_i
        self.yHat_i=yHat_i
        self.zHat_i=zHat_i
        self.M1_i=M1_i
        self.M2_i=M2_i
        self.v_i=v_i
        self.S_chi1_i=S_chi1_i
        self.S_chi2_i=S_chi2_i
        self.rfrak_chi1_x_i=rfrak_chi1_x_i
        self.rfrak_chi1_y_i=rfrak_chi1_y_i
        self.rfrak_chi2_x_i=rfrak_chi2_x_i
        self.rfrak_chi2_y_i=rfrak_chi2_y_i
        self.rfrak_frame_x_i=rfrak_frame_x_i
        self.rfrak_frame_y_i=rfrak_frame_y_i
        self.rfrak_frame_z_i=rfrak_frame_z_i
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
        chi1chi1=dot(chiVec1.vector,chiVec1.vector)
        chi1chi2=dot(chiVec1.vector,chiVec2.vector)
        chi2chi2=dot(chiVec2.vector,chiVec2.vector)
        chi1_n=dot(chiVec1.vector,nHat.vector)
        chi1_ell=dot(chiVec1.vector,ellHat.vector)
        chi2_n=dot(chiVec2.vector,nHat.vector)
        chi2_ell=dot(chiVec2.vector,ellHat.vector)
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
        self.xHat=xHat
        self.yHat=yHat
        self.zHat=zHat
        self.M1=M1
        self.M2=M2
        self.v=v
        self.S_chi1=S_chi1
        self.S_chi2=S_chi2
        self.rfrak_chi1_x=rfrak_chi1_x
        self.rfrak_chi1_y=rfrak_chi1_y
        self.rfrak_chi2_x=rfrak_chi2_x
        self.rfrak_chi2_y=rfrak_chi2_y
        self.rfrak_frame_x=rfrak_frame_x
        self.rfrak_frame_y=rfrak_frame_y
        self.rfrak_frame_z=rfrak_frame_z
        self.M=M
        self.delta=delta
        self.nu=nu
        self.R=R
        self.nHat=nHat
        self.ellHat=ellHat
        self.R_S1=R_S1
        self.R_S2=R_S2
        self.chiVec1=chiVec1
        self.chiVec2=chiVec2
        self.chi1chi1=chi1chi1
        self.chi1chi2=chi1chi2
        self.chi2chi2=chi2chi2
        self.chi1_n=chi1_n
        self.chi1_ell=chi1_ell
        self.chi2_n=chi2_n
        self.chi2_ell=chi2_ell
        self.S_ell=S_ell
        self.S_n=S_n
        self.Sigma_ell=Sigma_ell
        self.Sigma_n=Sigma_n
        self.chi_s_ell=chi_s_ell
        self.chi_a_ell=chi_a_ell
        self.Fcal_coeff=Fcal_coeff
        self.Fcal_0=Fcal_0
        self.Fcal_2=Fcal_2
        self.E_0=E_0
        self.E_2=E_2
        Phi=0.0
        EvolveSpin1=(S_chi1*conjugate(S_chi1)).norm>1e-12
        EvolveSpin2=(S_chi2*conjugate(S_chi2)).norm>1e-12
        self.v, self.rfrak_chi1_x, self.rfrak_chi1_y, self.rfrak_chi2_x, self.rfrak_chi2_y=v, rfrak_chi1_x, rfrak_chi1_y, rfrak_chi2_x, rfrak_chi2_y
        self.rfrak_frame_x, self.rfrak_frame_y, self.rfrak_frame_z, self.Phi=rfrak_frame_x, rfrak_frame_y, rfrak_frame_z, Phi
        self.EvolveSpin1=EvolveSpin1
        self.EvolveSpin2=EvolveSpin2
        self.terminal1,self.terminal2,self.NotForward=terminal1,terminal2,NotForward

    def __Recalculate__(self,t, y):
        v = y[0]
        rfrak_chi1_x = y[1]
        rfrak_chi1_y = y[2]
        rfrak_chi2_x = y[3]
        rfrak_chi2_y = y[4]
        rfrak_frame_x = y[5]
        rfrak_frame_y = y[6]
        rfrak_frame_z = y[7]
        Phi = y[8]
        self.v, self.rfrak_chi1_x, self.rfrak_chi1_y, self.rfrak_chi2_x, self.rfrak_chi2_y=v, rfrak_chi1_x, rfrak_chi1_y, rfrak_chi2_x, rfrak_chi2_y
        self.rfrak_frame_x, self.rfrak_frame_y, self.rfrak_frame_z, self.Phi=rfrak_frame_x, rfrak_frame_y, rfrak_frame_z, Phi
        xHat=self.xHat
        yHat=self.yHat
        zHat=self.zHat
        M1=self.M1
        M2=self.M2
        v=self.v
        S_chi1=self.S_chi1
        S_chi2=self.S_chi2
        rfrak_chi1_x=self.rfrak_chi1_x
        rfrak_chi1_y=self.rfrak_chi1_y
        rfrak_chi2_x=self.rfrak_chi2_x
        rfrak_chi2_y=self.rfrak_chi2_y
        rfrak_frame_x=self.rfrak_frame_x
        rfrak_frame_y=self.rfrak_frame_y
        rfrak_frame_z=self.rfrak_frame_z
        M=self.M
        delta=self.delta
        nu=self.nu
        R=self.R
        nHat=self.nHat
        ellHat=self.ellHat
        R_S1=self.R_S1
        R_S2=self.R_S2
        chiVec1=self.chiVec1
        chiVec2=self.chiVec2
        chi1chi1=self.chi1chi1
        chi1chi2=self.chi1chi2
        chi2chi2=self.chi2chi2
        chi1_n=self.chi1_n
        chi1_ell=self.chi1_ell
        chi2_n=self.chi2_n
        chi2_ell=self.chi2_ell
        S_ell=self.S_ell
        S_n=self.S_n
        Sigma_ell=self.Sigma_ell
        Sigma_n=self.Sigma_n
        chi_s_ell=self.chi_s_ell
        chi_a_ell=self.chi_a_ell
        Fcal_coeff=self.Fcal_coeff
        Fcal_0=self.Fcal_0
        Fcal_2=self.Fcal_2
        E_0=self.E_0
        E_2=self.E_2
        R = exp(rfrak_frame_x*xHat + rfrak_frame_y*yHat + rfrak_frame_z*zHat)
        nHat = R*xHat*conjugate(R)
        ellHat = R*zHat*conjugate(R)
        R_S1 = exp(rfrak_chi1_x*xHat + rfrak_chi1_y*yHat)
        R_S2 = exp(rfrak_chi2_x*xHat + rfrak_chi2_y*yHat)
        chiVec1 = S_chi1*R_S1*zHat*conjugate(R_S1)*conjugate(S_chi1)
        chiVec2 = S_chi2*R_S2*zHat*conjugate(R_S2)*conjugate(S_chi2)
        chi1_n = dot(chiVec1.vector,nHat.vector)
        chi1_ell = dot(chiVec1.vector,ellHat.vector)
        chi2_n = dot(chiVec2.vector,nHat.vector)
        chi2_ell = dot(chiVec2.vector,ellHat.vector)
        S_ell = M1**2*chi1_ell + M2**2*chi2_ell
        S_n = M1**2*chi1_n + M2**2*chi2_n
        Sigma_ell = M*(-M1*chi1_ell + M2*chi2_ell)
        Sigma_n = M*(-M1*chi1_n + M2*chi2_n)
        chi_s_ell = chi1_ell/2 + chi2_ell/2
        chi_a_ell = chi1_ell/2 - chi2_ell/2
        Fcal_coeff = 32*nu**2*v**10/5
        self.R=R
        self.nHat=nHat
        self.ellHat=ellHat
        self.R_S1=R_S1
        self.R_S2=R_S2
        self.chiVec1=chiVec1
        self.chiVec2=chiVec2
        self.chi1_n=chi1_n
        self.chi1_ell=chi1_ell
        self.chi2_n=chi2_n
        self.chi2_ell=chi2_ell
        self.S_ell=S_ell
        self.S_n=S_n
        self.Sigma_ell=Sigma_ell
        self.Sigma_n=Sigma_n
        self.chi_s_ell=chi_s_ell
        self.chi_a_ell=chi_a_ell
        self.Fcal_coeff=Fcal_coeff

    def OmegaVec_chiVec_1(self):
        xHat=self.xHat
        yHat=self.yHat
        zHat=self.zHat
        M1=self.M1
        M2=self.M2
        v=self.v
        S_chi1=self.S_chi1
        S_chi2=self.S_chi2
        rfrak_chi1_x=self.rfrak_chi1_x
        rfrak_chi1_y=self.rfrak_chi1_y
        rfrak_chi2_x=self.rfrak_chi2_x
        rfrak_chi2_y=self.rfrak_chi2_y
        rfrak_frame_x=self.rfrak_frame_x
        rfrak_frame_y=self.rfrak_frame_y
        rfrak_frame_z=self.rfrak_frame_z
        M=self.M
        delta=self.delta
        nu=self.nu
        R=self.R
        nHat=self.nHat
        ellHat=self.ellHat
        R_S1=self.R_S1
        R_S2=self.R_S2
        chiVec1=self.chiVec1
        chiVec2=self.chiVec2
        chi1_n=self.chi1_n
        chi2_n=self.chi2_n
        S_n=self.S_n
        Sigma_n=self.Sigma_n
        Omega1_coeff = v**5/M
        return Omega1_coeff*(ellHat*(-0.75*delta + 0.5*nu + v**2*(delta*(0.625*nu - 0.5625) + nu*(1.25 - 0.0416666666666667*nu) + 0.5625) + 0.75) + nHat*v*(3.0*chi1_n*nu + 3.0*M2**2*chi2_n/M**2) - M2**2*chiVec2*v/M**2)

    def OmegaVec_chiVec_2(self):
        xHat=self.xHat
        yHat=self.yHat
        zHat=self.zHat
        M1=self.M1
        M2=self.M2
        v=self.v
        S_chi1=self.S_chi1
        S_chi2=self.S_chi2
        rfrak_chi1_x=self.rfrak_chi1_x
        rfrak_chi1_y=self.rfrak_chi1_y
        rfrak_chi2_x=self.rfrak_chi2_x
        rfrak_chi2_y=self.rfrak_chi2_y
        rfrak_frame_x=self.rfrak_frame_x
        rfrak_frame_y=self.rfrak_frame_y
        rfrak_frame_z=self.rfrak_frame_z
        M=self.M
        delta=self.delta
        nu=self.nu
        R=self.R
        nHat=self.nHat
        ellHat=self.ellHat
        R_S1=self.R_S1
        R_S2=self.R_S2
        chiVec1=self.chiVec1
        chiVec2=self.chiVec2
        chi1_n=self.chi1_n
        chi2_n=self.chi2_n
        S_n=self.S_n
        Sigma_n=self.Sigma_n
        Omega2_coeff = v**5/M
        return Omega2_coeff*(ellHat*(0.75*delta + 0.5*nu + v**2*(delta*(0.5625 - 0.625*nu) + nu*(1.25 - 0.0416666666666667*nu) + 0.5625) + 0.75) + nHat*v*(3.0*chi2_n*nu + 3.0*M1**2*chi1_n/M**2) - M1**2*chiVec1*v/M**2)

    def OmegaVec(self):
        xHat=self.xHat
        yHat=self.yHat
        zHat=self.zHat
        M1=self.M1
        M2=self.M2
        v=self.v
        S_chi1=self.S_chi1
        S_chi2=self.S_chi2
        rfrak_chi1_x=self.rfrak_chi1_x
        rfrak_chi1_y=self.rfrak_chi1_y
        rfrak_chi2_x=self.rfrak_chi2_x
        rfrak_chi2_y=self.rfrak_chi2_y
        rfrak_frame_x=self.rfrak_frame_x
        rfrak_frame_y=self.rfrak_frame_y
        rfrak_frame_z=self.rfrak_frame_z
        M=self.M
        delta=self.delta
        nu=self.nu
        R=self.R
        nHat=self.nHat
        ellHat=self.ellHat
        R_S1=self.R_S1
        R_S2=self.R_S2
        chiVec1=self.chiVec1
        chiVec2=self.chiVec2
        chi1_n=self.chi1_n
        chi2_n=self.chi2_n
        S_n=self.S_n
        Sigma_n=self.Sigma_n
        gamma_PN_2 = 1.0 - 0.333333333333333*nu
        gamma_PN_0 = 1.00000000000000
        a_ell_0 = 7.0*S_n + 3.0*Sigma_n*delta
        a_ell_2 = S_n*(-9.66666666666667*nu - 10.0) + Sigma_n*delta*(-4.5*nu - 6.0)
        return ellHat*v**3/M + nHat*v**6*(a_ell_0 + a_ell_2*v**2)*(gamma_PN_0 + gamma_PN_2*v**2)/M**3

    def OrbitalAngularMomentum(self):
        xHat=self.xHat
        yHat=self.yHat
        zHat=self.zHat
        M1=self.M1
        M2=self.M2
        v=self.v
        S_chi1=self.S_chi1
        S_chi2=self.S_chi2
        rfrak_chi1_x=self.rfrak_chi1_x
        rfrak_chi1_y=self.rfrak_chi1_y
        rfrak_chi2_x=self.rfrak_chi2_x
        rfrak_chi2_y=self.rfrak_chi2_y
        rfrak_frame_x=self.rfrak_frame_x
        rfrak_frame_y=self.rfrak_frame_y
        rfrak_frame_z=self.rfrak_frame_z
        M=self.M
        delta=self.delta
        nu=self.nu
        R=self.R
        nHat=self.nHat
        ellHat=self.ellHat
        R_S1=self.R_S1
        R_S2=self.R_S2
        chiVec1=self.chiVec1
        chiVec2=self.chiVec2
        chi1_n=self.chi1_n
        chi2_n=self.chi2_n
        S_n=self.S_n
        Sigma_n=self.Sigma_n
        L_coeff = M**2*nu/v
        L_0 = ellHat
        L_2 = ellHat*(0.166666666666667*nu + 1.5)
        return L_coeff*(L_0 + L_2*v**2)

               

    #Evolve PN
    def __TaylorT1__(self,t, y):
        self.__Recalculate__(t, y)
        xHat=self.xHat
        yHat=self.yHat
        zHat=self.zHat
        M1=self.M1
        M2=self.M2
        v=self.v
        S_chi1=self.S_chi1
        S_chi2=self.S_chi2
        rfrak_chi1_x=self.rfrak_chi1_x
        rfrak_chi1_y=self.rfrak_chi1_y
        rfrak_chi2_x=self.rfrak_chi2_x
        rfrak_chi2_y=self.rfrak_chi2_y
        rfrak_frame_x=self.rfrak_frame_x
        rfrak_frame_y=self.rfrak_frame_y
        rfrak_frame_z=self.rfrak_frame_z
        M=self.M
        delta=self.delta
        nu=self.nu
        R=self.R
        nHat=self.nHat
        ellHat=self.ellHat
        R_S1=self.R_S1
        R_S2=self.R_S2
        chiVec1=self.chiVec1
        chiVec2=self.chiVec2
        chi1chi1=self.chi1chi1
        chi1chi2=self.chi1chi2
        chi2chi2=self.chi2chi2
        chi1_n=self.chi1_n
        chi1_ell=self.chi1_ell
        chi2_n=self.chi2_n
        chi2_ell=self.chi2_ell
        S_ell=self.S_ell
        S_n=self.S_n
        Sigma_ell=self.Sigma_ell
        Sigma_n=self.Sigma_n
        chi_s_ell=self.chi_s_ell
        chi_a_ell=self.chi_a_ell
        Fcal_coeff=self.Fcal_coeff
        Fcal_0=self.Fcal_0
        Fcal_2=self.Fcal_2
        E_0=self.E_0
        E_2=self.E_2
        v, rfrak_chi1_x, rfrak_chi1_y, rfrak_chi2_x, rfrak_chi2_y=self.v, self.rfrak_chi1_x, self.rfrak_chi1_y, self.rfrak_chi2_x, self.rfrak_chi2_y
        rfrak_frame_x, rfrak_frame_y, rfrak_frame_z, Phi=self.rfrak_frame_x, self.rfrak_frame_y, self.rfrak_frame_z, self.Phi
        dvdt_T1=0
        if v>=1.0 and self.NotForward:
            print("Beyond domain of PN validity, this is a good way to terminate.")
            self.terminal1=False
        Flux = Fcal_coeff*(Fcal_0 + Fcal_2*v**2)
        dEdv = -M*nu*v*(E_0 + 2.0*E_2*v**2)
        Absorption = 0
        dvdt_T1 = (-Absorption - Flux)/dEdv
        if dvdt_T1<1.0e-12 and self.NotForward:
            print("v is decreasing, which is not an uncommon way to stop.")
            self.terminal2=False
        return self.CommonRHS(dvdt_T1, y)

    def CommonRHS(self,dvdt, y):
        dydt=np.zeros(9)
        rfrak_frame=np.zeros(3)
        rfrak_frame[0] = y[5]
        rfrak_frame[1] = y[6]
        rfrak_frame[2] = y[7]
        rfrakdot_frame = FrameFromAngularVelocityIntegrand(rfrak_frame, self.OmegaVec().vector)
        dydt[0] = dvdt
        if(self.EvolveSpin1):
            dydt[1], dydt[2]= FrameFromAngularVelocity_2D_Integrand(y[1], y[2],\
                (self.S_chi1.inverse*self.OmegaVec_chiVec_1()*self.S_chi1).vector)
        else:
            dydt[1] = 0.0
            dydt[2] = 0.0
        if(self.EvolveSpin2):
            dydt[3], dydt[4] = FrameFromAngularVelocity_2D_Integrand(y[3], y[4],\
                (self.S_chi2.inverse*self.OmegaVec_chiVec_2()*self.S_chi2).vector)
        else:
            dydt[3] = 0.0
            dydt[4] = 0.0
        dydt[5] = rfrakdot_frame[0]
        dydt[6] = rfrakdot_frame[1]
        dydt[7] = rfrakdot_frame[2]
        dydt[8] = self.v*self.v*self.v/self.M

        return dydt

    # Prepare for integration
    def Evolution(self): 
        v_i,rfrak_chi1_x_i,rfrak_chi1_y_i,rfrak_chi2_x_i,rfrak_chi2_y_i,rfrak_frame_x_i,rfrak_frame_y_i,rfrak_frame_z_i\
            =self.v_i,self.rfrak_chi1_x_i,self.rfrak_chi1_y_i,self.rfrak_chi2_x_i,self.rfrak_chi2_y_i,self.rfrak_frame_x_i,self.rfrak_frame_y_i,self.rfrak_frame_z_i
        def terminate(t,y):
            return 1.0*self.terminal1*self.terminal2
        terminate.terminal=True
        TMerger=5.0/(256.0*self.nu*v_i**8)
        TEnd=TMerger
        if self.t_PN2:
            TEnd=self.t_PN2
        time=[0.0]
        while time[-1]<TEnd and 2*self.M*(256*self.nu*(TMerger-time[-1])/5)**(3/8)/self.StepsPerOrbit>self.MinStep:
            time.append(time[-1]+2*self.M*(256*self.nu*(TMerger-time[-1])/5)**(3/8)/self.StepsPerOrbit)
        time=np.delete(time, -1)
        
        # Integrate
        yy=solve_ivp(self.__TaylorT1__, [time[0],time[-1]], [v_i,rfrak_chi1_x_i,\
            rfrak_chi1_y_i,rfrak_chi2_x_i,rfrak_chi2_y_i,rfrak_frame_x_i,\
            rfrak_frame_y_i,rfrak_frame_z_i,0.0], method='DOP853',\
            t_eval=time, dense_output=True, events=terminate, rtol=self.tol, atol=self.tol)           
        if self.ForwardInTime:
            self.NotForward=False
            time=[0.0]
            TStart=-3*TMerger
            if self.t_PN1:
                TStart=self.t_PN1
            while time[-1]>TStart:
                time.append(time[-1]-2*self.M*(256*self.nu*(TMerger-time[-1])/5)**(3/8)/self.StepsPerOrbit)
            yyForward=solve_ivp(self.__TaylorT1__, [time[0],time[-1]], [v_i,rfrak_chi1_x_i,\
                rfrak_chi1_y_i,rfrak_chi2_x_i,rfrak_chi2_y_i,rfrak_frame_x_i,\
                rfrak_frame_y_i,rfrak_frame_z_i,0.0], method='DOP853',\
                t_eval=time, dense_output=True, rtol=self.tol, atol=self.tol)
            yy.t=np.append(yyForward.t[1:][::-1],yy.t)
            data=np.empty((9,len(yy.t)))
            for i in range(9):
                data[i]=np.append(yyForward.y[i][1:][::-1],yy.y[i])
            yy.y=data
            
        return yy

spec=[("xHat", float64),("yHat", float64),("zHat", float64),("M1", float64),("M2", float64),("v", float64),("S_chi1", float64),("S_chi2", float64),("rfrak_chi1_x", float64),("rfrak_chi1_y", float64),("rfrak_chi2_x", float64),("rfrak_chi2_y", float64),("rfrak_frame_x", float64),("rfrak_frame_y", float64),("rfrak_frame_z", float64),("M", float64),("delta", float64),("nu", float64),("R", float64),("nHat", float64),("lambdaHat", float64),("ellHat", float64),("R_S1", float64),("R_S2", float64),("chiVec1", float64),("chiVec2", float64),("chi1chi1", float64),("chi1chi2", float64),("chi2chi2", float64),("chi1_n", float64),("chi1_lambda", float64),("chi1_ell", float64),("chi2_n", float64),("chi2_lambda", float64),("chi2_ell", float64),("S_ell", float64),("S_n", float64),("S_lambda", float64),("Sigma_ell", float64),("Sigma_n", float64),("Sigma_lambda", float64),("chi_s_ell", float64),("chi_a_ell", float64),("Fcal_coeff", float64),("Fcal_0", float64),("Fcal_2", float64),("Fcal_3", float64),("Fcal_SO_3", float64),("E_0", float64),("E_2", float64),("E_SO_3", float64),("xHat_i", float64),("yHat_i", float64),("zHat_i", float64),("M1_i", float64),("M2_i", float64),("v_i", float64),("S_chi1_i", float64),("S_chi2_i", float64),("rfrak_chi1_x_i", float64),("rfrak_chi1_y_i", float64),("rfrak_chi2_x_i", float64),("rfrak_chi2_y_i", float64),("rfrak_frame_x_i", float64),("rfrak_frame_y_i", float64),("rfrak_frame_z_i", float64),("t_PN1",float64),("t_PN2",float64),("StepsPerOrbit",int32),("ForwardInTime",boolean),("tol",float64),("MinStep",float64),("EvolveSpin1",boolean),("EvolveSpin2",boolean),("terminal1",boolean),("terminal2",boolean),("NotForward",boolean)]
@jitclass(spec)
class TaylorTn_1p5PN_Q : 
    def __init__(self, xHat_i, yHat_i, zHat_i, M1_i, M2_i, v_i, S_chi1_i, S_chi2_i, rfrak_chi1_x_i, rfrak_chi1_y_i, rfrak_chi2_x_i, rfrak_chi2_y_i, rfrak_frame_x_i, rfrak_frame_y_i, rfrak_frame_z_i, t_PN1=False, t_PN2=False, StepsPerOrbit=32, ForwardInTime=True, tol=1e-12, MinStep=1e-7) :
        self.t_PN1,self.t_PN2,self.StepsPerOrbit,self.ForwardInTime,self.tol,self.MinStep=t_PN1, t_PN2, StepsPerOrbit, ForwardInTime, tol, MinStep
        terminal1=True
        terminal2=True
        NotForward=True
        self.xHat_i=xHat_i
        self.yHat_i=yHat_i
        self.zHat_i=zHat_i
        self.M1_i=M1_i
        self.M2_i=M2_i
        self.v_i=v_i
        self.S_chi1_i=S_chi1_i
        self.S_chi2_i=S_chi2_i
        self.rfrak_chi1_x_i=rfrak_chi1_x_i
        self.rfrak_chi1_y_i=rfrak_chi1_y_i
        self.rfrak_chi2_x_i=rfrak_chi2_x_i
        self.rfrak_chi2_y_i=rfrak_chi2_y_i
        self.rfrak_frame_x_i=rfrak_frame_x_i
        self.rfrak_frame_y_i=rfrak_frame_y_i
        self.rfrak_frame_z_i=rfrak_frame_z_i
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
        chi1chi1=dot(chiVec1.vector,chiVec1.vector)
        chi1chi2=dot(chiVec1.vector,chiVec2.vector)
        chi2chi2=dot(chiVec2.vector,chiVec2.vector)
        chi1_n=dot(chiVec1.vector,nHat.vector)
        chi1_lambda=dot(chiVec1.vector,lambdaHat.vector)
        chi1_ell=dot(chiVec1.vector,ellHat.vector)
        chi2_n=dot(chiVec2.vector,nHat.vector)
        chi2_lambda=dot(chiVec2.vector,lambdaHat.vector)
        chi2_ell=dot(chiVec2.vector,ellHat.vector)
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
        self.xHat=xHat
        self.yHat=yHat
        self.zHat=zHat
        self.M1=M1
        self.M2=M2
        self.v=v
        self.S_chi1=S_chi1
        self.S_chi2=S_chi2
        self.rfrak_chi1_x=rfrak_chi1_x
        self.rfrak_chi1_y=rfrak_chi1_y
        self.rfrak_chi2_x=rfrak_chi2_x
        self.rfrak_chi2_y=rfrak_chi2_y
        self.rfrak_frame_x=rfrak_frame_x
        self.rfrak_frame_y=rfrak_frame_y
        self.rfrak_frame_z=rfrak_frame_z
        self.M=M
        self.delta=delta
        self.nu=nu
        self.R=R
        self.nHat=nHat
        self.lambdaHat=lambdaHat
        self.ellHat=ellHat
        self.R_S1=R_S1
        self.R_S2=R_S2
        self.chiVec1=chiVec1
        self.chiVec2=chiVec2
        self.chi1chi1=chi1chi1
        self.chi1chi2=chi1chi2
        self.chi2chi2=chi2chi2
        self.chi1_n=chi1_n
        self.chi1_lambda=chi1_lambda
        self.chi1_ell=chi1_ell
        self.chi2_n=chi2_n
        self.chi2_lambda=chi2_lambda
        self.chi2_ell=chi2_ell
        self.S_ell=S_ell
        self.S_n=S_n
        self.S_lambda=S_lambda
        self.Sigma_ell=Sigma_ell
        self.Sigma_n=Sigma_n
        self.Sigma_lambda=Sigma_lambda
        self.chi_s_ell=chi_s_ell
        self.chi_a_ell=chi_a_ell
        self.Fcal_coeff=Fcal_coeff
        self.Fcal_0=Fcal_0
        self.Fcal_2=Fcal_2
        self.Fcal_3=Fcal_3
        self.Fcal_SO_3=Fcal_SO_3
        self.E_0=E_0
        self.E_2=E_2
        self.E_SO_3=E_SO_3
        Phi=0.0
        EvolveSpin1=(S_chi1*conjugate(S_chi1)).norm>1e-12
        EvolveSpin2=(S_chi2*conjugate(S_chi2)).norm>1e-12
        self.v, self.rfrak_chi1_x, self.rfrak_chi1_y, self.rfrak_chi2_x, self.rfrak_chi2_y=v, rfrak_chi1_x, rfrak_chi1_y, rfrak_chi2_x, rfrak_chi2_y
        self.rfrak_frame_x, self.rfrak_frame_y, self.rfrak_frame_z, self.Phi=rfrak_frame_x, rfrak_frame_y, rfrak_frame_z, Phi
        self.EvolveSpin1=EvolveSpin1
        self.EvolveSpin2=EvolveSpin2
        self.terminal1,self.terminal2,self.NotForward=terminal1,terminal2,NotForward

    def __Recalculate__(self,t, y):
        v = y[0]
        rfrak_chi1_x = y[1]
        rfrak_chi1_y = y[2]
        rfrak_chi2_x = y[3]
        rfrak_chi2_y = y[4]
        rfrak_frame_x = y[5]
        rfrak_frame_y = y[6]
        rfrak_frame_z = y[7]
        Phi = y[8]
        self.v, self.rfrak_chi1_x, self.rfrak_chi1_y, self.rfrak_chi2_x, self.rfrak_chi2_y=v, rfrak_chi1_x, rfrak_chi1_y, rfrak_chi2_x, rfrak_chi2_y
        self.rfrak_frame_x, self.rfrak_frame_y, self.rfrak_frame_z, self.Phi=rfrak_frame_x, rfrak_frame_y, rfrak_frame_z, Phi
        xHat=self.xHat
        yHat=self.yHat
        zHat=self.zHat
        M1=self.M1
        M2=self.M2
        v=self.v
        S_chi1=self.S_chi1
        S_chi2=self.S_chi2
        rfrak_chi1_x=self.rfrak_chi1_x
        rfrak_chi1_y=self.rfrak_chi1_y
        rfrak_chi2_x=self.rfrak_chi2_x
        rfrak_chi2_y=self.rfrak_chi2_y
        rfrak_frame_x=self.rfrak_frame_x
        rfrak_frame_y=self.rfrak_frame_y
        rfrak_frame_z=self.rfrak_frame_z
        M=self.M
        delta=self.delta
        nu=self.nu
        R=self.R
        nHat=self.nHat
        lambdaHat=self.lambdaHat
        ellHat=self.ellHat
        R_S1=self.R_S1
        R_S2=self.R_S2
        chiVec1=self.chiVec1
        chiVec2=self.chiVec2
        chi1chi1=self.chi1chi1
        chi1chi2=self.chi1chi2
        chi2chi2=self.chi2chi2
        chi1_n=self.chi1_n
        chi1_lambda=self.chi1_lambda
        chi1_ell=self.chi1_ell
        chi2_n=self.chi2_n
        chi2_lambda=self.chi2_lambda
        chi2_ell=self.chi2_ell
        S_ell=self.S_ell
        S_n=self.S_n
        S_lambda=self.S_lambda
        Sigma_ell=self.Sigma_ell
        Sigma_n=self.Sigma_n
        Sigma_lambda=self.Sigma_lambda
        chi_s_ell=self.chi_s_ell
        chi_a_ell=self.chi_a_ell
        Fcal_coeff=self.Fcal_coeff
        Fcal_0=self.Fcal_0
        Fcal_2=self.Fcal_2
        Fcal_3=self.Fcal_3
        Fcal_SO_3=self.Fcal_SO_3
        E_0=self.E_0
        E_2=self.E_2
        E_SO_3=self.E_SO_3
        R = exp(rfrak_frame_x*xHat + rfrak_frame_y*yHat + rfrak_frame_z*zHat)
        nHat = R*xHat*conjugate(R)
        lambdaHat = R*yHat*conjugate(R)
        ellHat = R*zHat*conjugate(R)
        R_S1 = exp(rfrak_chi1_x*xHat + rfrak_chi1_y*yHat)
        R_S2 = exp(rfrak_chi2_x*xHat + rfrak_chi2_y*yHat)
        chiVec1 = S_chi1*R_S1*zHat*conjugate(R_S1)*conjugate(S_chi1)
        chiVec2 = S_chi2*R_S2*zHat*conjugate(R_S2)*conjugate(S_chi2)
        chi1_n = dot(chiVec1.vector,nHat.vector)
        chi1_lambda = dot(chiVec1.vector,lambdaHat.vector)
        chi1_ell = dot(chiVec1.vector,ellHat.vector)
        chi2_n = dot(chiVec2.vector,nHat.vector)
        chi2_lambda = dot(chiVec2.vector,lambdaHat.vector)
        chi2_ell = dot(chiVec2.vector,ellHat.vector)
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
        self.R=R
        self.nHat=nHat
        self.lambdaHat=lambdaHat
        self.ellHat=ellHat
        self.R_S1=R_S1
        self.R_S2=R_S2
        self.chiVec1=chiVec1
        self.chiVec2=chiVec2
        self.chi1_n=chi1_n
        self.chi1_lambda=chi1_lambda
        self.chi1_ell=chi1_ell
        self.chi2_n=chi2_n
        self.chi2_lambda=chi2_lambda
        self.chi2_ell=chi2_ell
        self.S_ell=S_ell
        self.S_n=S_n
        self.S_lambda=S_lambda
        self.Sigma_ell=Sigma_ell
        self.Sigma_n=Sigma_n
        self.Sigma_lambda=Sigma_lambda
        self.chi_s_ell=chi_s_ell
        self.chi_a_ell=chi_a_ell
        self.Fcal_coeff=Fcal_coeff
        self.Fcal_SO_3=Fcal_SO_3
        self.E_SO_3=E_SO_3

    def OmegaVec_chiVec_1(self):
        xHat=self.xHat
        yHat=self.yHat
        zHat=self.zHat
        M1=self.M1
        M2=self.M2
        v=self.v
        S_chi1=self.S_chi1
        S_chi2=self.S_chi2
        rfrak_chi1_x=self.rfrak_chi1_x
        rfrak_chi1_y=self.rfrak_chi1_y
        rfrak_chi2_x=self.rfrak_chi2_x
        rfrak_chi2_y=self.rfrak_chi2_y
        rfrak_frame_x=self.rfrak_frame_x
        rfrak_frame_y=self.rfrak_frame_y
        rfrak_frame_z=self.rfrak_frame_z
        M=self.M
        delta=self.delta
        nu=self.nu
        R=self.R
        nHat=self.nHat
        lambdaHat=self.lambdaHat
        ellHat=self.ellHat
        R_S1=self.R_S1
        R_S2=self.R_S2
        chiVec1=self.chiVec1
        chiVec2=self.chiVec2
        chi1_n=self.chi1_n
        chi1_lambda=self.chi1_lambda
        chi1_ell=self.chi1_ell
        chi2_n=self.chi2_n
        chi2_lambda=self.chi2_lambda
        chi2_ell=self.chi2_ell
        S_ell=self.S_ell
        S_n=self.S_n
        S_lambda=self.S_lambda
        Sigma_ell=self.Sigma_ell
        Sigma_n=self.Sigma_n
        Sigma_lambda=self.Sigma_lambda
        Omega1_coeff = v**5/M
        return Omega1_coeff*(ellHat*(-0.75*delta + 0.5*nu + v**2*(delta*(0.625*nu - 0.5625) + nu*(1.25 - 0.0416666666666667*nu) + 0.5625) + 0.75) + nHat*v*(3.0*chi1_n*nu + 3.0*M2**2*chi2_n/M**2) - M2**2*chiVec2*v/M**2)

    def OmegaVec_chiVec_2(self):
        xHat=self.xHat
        yHat=self.yHat
        zHat=self.zHat
        M1=self.M1
        M2=self.M2
        v=self.v
        S_chi1=self.S_chi1
        S_chi2=self.S_chi2
        rfrak_chi1_x=self.rfrak_chi1_x
        rfrak_chi1_y=self.rfrak_chi1_y
        rfrak_chi2_x=self.rfrak_chi2_x
        rfrak_chi2_y=self.rfrak_chi2_y
        rfrak_frame_x=self.rfrak_frame_x
        rfrak_frame_y=self.rfrak_frame_y
        rfrak_frame_z=self.rfrak_frame_z
        M=self.M
        delta=self.delta
        nu=self.nu
        R=self.R
        nHat=self.nHat
        lambdaHat=self.lambdaHat
        ellHat=self.ellHat
        R_S1=self.R_S1
        R_S2=self.R_S2
        chiVec1=self.chiVec1
        chiVec2=self.chiVec2
        chi1_n=self.chi1_n
        chi1_lambda=self.chi1_lambda
        chi1_ell=self.chi1_ell
        chi2_n=self.chi2_n
        chi2_lambda=self.chi2_lambda
        chi2_ell=self.chi2_ell
        S_ell=self.S_ell
        S_n=self.S_n
        S_lambda=self.S_lambda
        Sigma_ell=self.Sigma_ell
        Sigma_n=self.Sigma_n
        Sigma_lambda=self.Sigma_lambda
        Omega2_coeff = v**5/M
        return Omega2_coeff*(ellHat*(0.75*delta + 0.5*nu + v**2*(delta*(0.5625 - 0.625*nu) + nu*(1.25 - 0.0416666666666667*nu) + 0.5625) + 0.75) + nHat*v*(3.0*chi2_n*nu + 3.0*M1**2*chi1_n/M**2) - M1**2*chiVec1*v/M**2)

    def OmegaVec(self):
        xHat=self.xHat
        yHat=self.yHat
        zHat=self.zHat
        M1=self.M1
        M2=self.M2
        v=self.v
        S_chi1=self.S_chi1
        S_chi2=self.S_chi2
        rfrak_chi1_x=self.rfrak_chi1_x
        rfrak_chi1_y=self.rfrak_chi1_y
        rfrak_chi2_x=self.rfrak_chi2_x
        rfrak_chi2_y=self.rfrak_chi2_y
        rfrak_frame_x=self.rfrak_frame_x
        rfrak_frame_y=self.rfrak_frame_y
        rfrak_frame_z=self.rfrak_frame_z
        M=self.M
        delta=self.delta
        nu=self.nu
        R=self.R
        nHat=self.nHat
        lambdaHat=self.lambdaHat
        ellHat=self.ellHat
        R_S1=self.R_S1
        R_S2=self.R_S2
        chiVec1=self.chiVec1
        chiVec2=self.chiVec2
        chi1_n=self.chi1_n
        chi1_lambda=self.chi1_lambda
        chi1_ell=self.chi1_ell
        chi2_n=self.chi2_n
        chi2_lambda=self.chi2_lambda
        chi2_ell=self.chi2_ell
        S_ell=self.S_ell
        S_n=self.S_n
        S_lambda=self.S_lambda
        Sigma_ell=self.Sigma_ell
        Sigma_n=self.Sigma_n
        Sigma_lambda=self.Sigma_lambda
        gamma_PN_0 = 1.00000000000000
        a_ell_2 = S_n*(-9.66666666666667*nu - 10.0) + Sigma_n*delta*(-4.5*nu - 6.0)
        gamma_PN_3 = (1.66666666666667*S_ell + Sigma_ell*delta)/M**2
        a_ell_0 = 7.0*S_n + 3.0*Sigma_n*delta
        gamma_PN_2 = 1.0 - 0.333333333333333*nu
        return ellHat*v**3/M + nHat*v**6*(a_ell_0 + a_ell_2*v**2)*(gamma_PN_0 + v**2*(gamma_PN_2 + gamma_PN_3*v))/M**3

    def OrbitalAngularMomentum(self):
        xHat=self.xHat
        yHat=self.yHat
        zHat=self.zHat
        M1=self.M1
        M2=self.M2
        v=self.v
        S_chi1=self.S_chi1
        S_chi2=self.S_chi2
        rfrak_chi1_x=self.rfrak_chi1_x
        rfrak_chi1_y=self.rfrak_chi1_y
        rfrak_chi2_x=self.rfrak_chi2_x
        rfrak_chi2_y=self.rfrak_chi2_y
        rfrak_frame_x=self.rfrak_frame_x
        rfrak_frame_y=self.rfrak_frame_y
        rfrak_frame_z=self.rfrak_frame_z
        M=self.M
        delta=self.delta
        nu=self.nu
        R=self.R
        nHat=self.nHat
        lambdaHat=self.lambdaHat
        ellHat=self.ellHat
        R_S1=self.R_S1
        R_S2=self.R_S2
        chiVec1=self.chiVec1
        chiVec2=self.chiVec2
        chi1_n=self.chi1_n
        chi1_lambda=self.chi1_lambda
        chi1_ell=self.chi1_ell
        chi2_n=self.chi2_n
        chi2_lambda=self.chi2_lambda
        chi2_ell=self.chi2_ell
        S_ell=self.S_ell
        S_n=self.S_n
        S_lambda=self.S_lambda
        Sigma_ell=self.Sigma_ell
        Sigma_n=self.Sigma_n
        Sigma_lambda=self.Sigma_lambda
        L_SO_3 = 0.166666666666667*ellHat*(-35.0*S_ell - 15.0*Sigma_ell*delta)/M**2 + lambdaHat*(-3.0*S_lambda - Sigma_lambda*delta)/M**2 + 0.5*nHat*(S_n + Sigma_n*delta)/M**2
        L_2 = ellHat*(0.166666666666667*nu + 1.5)
        L_0 = ellHat
        L_coeff = M**2*nu/v
        return L_coeff*(L_0 + v**2*(L_2 + L_SO_3*v))

               

    #Evolve PN
    def __TaylorT1__(self,t, y):
        self.__Recalculate__(t, y)
        xHat=self.xHat
        yHat=self.yHat
        zHat=self.zHat
        M1=self.M1
        M2=self.M2
        v=self.v
        S_chi1=self.S_chi1
        S_chi2=self.S_chi2
        rfrak_chi1_x=self.rfrak_chi1_x
        rfrak_chi1_y=self.rfrak_chi1_y
        rfrak_chi2_x=self.rfrak_chi2_x
        rfrak_chi2_y=self.rfrak_chi2_y
        rfrak_frame_x=self.rfrak_frame_x
        rfrak_frame_y=self.rfrak_frame_y
        rfrak_frame_z=self.rfrak_frame_z
        M=self.M
        delta=self.delta
        nu=self.nu
        R=self.R
        nHat=self.nHat
        lambdaHat=self.lambdaHat
        ellHat=self.ellHat
        R_S1=self.R_S1
        R_S2=self.R_S2
        chiVec1=self.chiVec1
        chiVec2=self.chiVec2
        chi1chi1=self.chi1chi1
        chi1chi2=self.chi1chi2
        chi2chi2=self.chi2chi2
        chi1_n=self.chi1_n
        chi1_lambda=self.chi1_lambda
        chi1_ell=self.chi1_ell
        chi2_n=self.chi2_n
        chi2_lambda=self.chi2_lambda
        chi2_ell=self.chi2_ell
        S_ell=self.S_ell
        S_n=self.S_n
        S_lambda=self.S_lambda
        Sigma_ell=self.Sigma_ell
        Sigma_n=self.Sigma_n
        Sigma_lambda=self.Sigma_lambda
        chi_s_ell=self.chi_s_ell
        chi_a_ell=self.chi_a_ell
        Fcal_coeff=self.Fcal_coeff
        Fcal_0=self.Fcal_0
        Fcal_2=self.Fcal_2
        Fcal_3=self.Fcal_3
        Fcal_SO_3=self.Fcal_SO_3
        E_0=self.E_0
        E_2=self.E_2
        E_SO_3=self.E_SO_3
        v, rfrak_chi1_x, rfrak_chi1_y, rfrak_chi2_x, rfrak_chi2_y=self.v, self.rfrak_chi1_x, self.rfrak_chi1_y, self.rfrak_chi2_x, self.rfrak_chi2_y
        rfrak_frame_x, rfrak_frame_y, rfrak_frame_z, Phi=self.rfrak_frame_x, self.rfrak_frame_y, self.rfrak_frame_z, self.Phi
        dvdt_T1=0
        if v>=1.0 and self.NotForward:
            print("Beyond domain of PN validity, this is a good way to terminate.")
            self.terminal1=False
        Flux = Fcal_coeff*(Fcal_0 + v**2*(Fcal_2 + v*(Fcal_3 + Fcal_SO_3)))
        dEdv = -0.5*M*nu*v*(2.0*E_0 + v**2*(4.0*E_2 + 5.0*E_SO_3*v))
        Absorption = 0
        dvdt_T1 = (-Absorption - Flux)/dEdv
        if dvdt_T1<1.0e-12 and self.NotForward:
            print("v is decreasing, which is not an uncommon way to stop.")
            self.terminal2=False
        return self.CommonRHS(dvdt_T1, y)

    def CommonRHS(self,dvdt, y):
        dydt=np.zeros(9)
        rfrak_frame=np.zeros(3)
        rfrak_frame[0] = y[5]
        rfrak_frame[1] = y[6]
        rfrak_frame[2] = y[7]
        rfrakdot_frame = FrameFromAngularVelocityIntegrand(rfrak_frame, self.OmegaVec().vector)
        dydt[0] = dvdt
        if(self.EvolveSpin1):
            dydt[1], dydt[2]= FrameFromAngularVelocity_2D_Integrand(y[1], y[2],\
                (self.S_chi1.inverse*self.OmegaVec_chiVec_1()*self.S_chi1).vector)
        else:
            dydt[1] = 0.0
            dydt[2] = 0.0
        if(self.EvolveSpin2):
            dydt[3], dydt[4] = FrameFromAngularVelocity_2D_Integrand(y[3], y[4],\
                (self.S_chi2.inverse*self.OmegaVec_chiVec_2()*self.S_chi2).vector)
        else:
            dydt[3] = 0.0
            dydt[4] = 0.0
        dydt[5] = rfrakdot_frame[0]
        dydt[6] = rfrakdot_frame[1]
        dydt[7] = rfrakdot_frame[2]
        dydt[8] = self.v*self.v*self.v/self.M

        return dydt

    # Prepare for integration
    def Evolution(self): 
        v_i,rfrak_chi1_x_i,rfrak_chi1_y_i,rfrak_chi2_x_i,rfrak_chi2_y_i,rfrak_frame_x_i,rfrak_frame_y_i,rfrak_frame_z_i\
            =self.v_i,self.rfrak_chi1_x_i,self.rfrak_chi1_y_i,self.rfrak_chi2_x_i,self.rfrak_chi2_y_i,self.rfrak_frame_x_i,self.rfrak_frame_y_i,self.rfrak_frame_z_i
        def terminate(t,y):
            return 1.0*self.terminal1*self.terminal2
        terminate.terminal=True
        TMerger=5.0/(256.0*self.nu*v_i**8)
        TEnd=TMerger
        if self.t_PN2:
            TEnd=self.t_PN2
        time=[0.0]
        while time[-1]<TEnd and 2*self.M*(256*self.nu*(TMerger-time[-1])/5)**(3/8)/self.StepsPerOrbit>self.MinStep:
            time.append(time[-1]+2*self.M*(256*self.nu*(TMerger-time[-1])/5)**(3/8)/self.StepsPerOrbit)
        time=np.delete(time, -1)
        
        # Integrate
        yy=solve_ivp(self.__TaylorT1__, [time[0],time[-1]], [v_i,rfrak_chi1_x_i,\
            rfrak_chi1_y_i,rfrak_chi2_x_i,rfrak_chi2_y_i,rfrak_frame_x_i,\
            rfrak_frame_y_i,rfrak_frame_z_i,0.0], method='DOP853',\
            t_eval=time, dense_output=True, events=terminate, rtol=self.tol, atol=self.tol)           
        if self.ForwardInTime:
            self.NotForward=False
            time=[0.0]
            TStart=-3*TMerger
            if self.t_PN1:
                TStart=self.t_PN1
            while time[-1]>TStart:
                time.append(time[-1]-2*self.M*(256*self.nu*(TMerger-time[-1])/5)**(3/8)/self.StepsPerOrbit)
            yyForward=solve_ivp(self.__TaylorT1__, [time[0],time[-1]], [v_i,rfrak_chi1_x_i,\
                rfrak_chi1_y_i,rfrak_chi2_x_i,rfrak_chi2_y_i,rfrak_frame_x_i,\
                rfrak_frame_y_i,rfrak_frame_z_i,0.0], method='DOP853',\
                t_eval=time, dense_output=True, rtol=self.tol, atol=self.tol)
            yy.t=np.append(yyForward.t[1:][::-1],yy.t)
            data=np.empty((9,len(yy.t)))
            for i in range(9):
                data[i]=np.append(yyForward.y[i][1:][::-1],yy.y[i])
            yy.y=data
            
        return yy

spec=[("xHat", float64),("yHat", float64),("zHat", float64),("M1", float64),("M2", float64),("v", float64),("S_chi1", float64),("S_chi2", float64),("rfrak_chi1_x", float64),("rfrak_chi1_y", float64),("rfrak_chi2_x", float64),("rfrak_chi2_y", float64),("rfrak_frame_x", float64),("rfrak_frame_y", float64),("rfrak_frame_z", float64),("M", float64),("delta", float64),("nu", float64),("R", float64),("nHat", float64),("lambdaHat", float64),("ellHat", float64),("R_S1", float64),("R_S2", float64),("chiVec1", float64),("chiVec2", float64),("chi1chi1", float64),("chi1chi2", float64),("chi2chi2", float64),("chi1_n", float64),("chi1_lambda", float64),("chi1_ell", float64),("chi2_n", float64),("chi2_lambda", float64),("chi2_ell", float64),("S_ell", float64),("S_n", float64),("S_lambda", float64),("Sigma_ell", float64),("Sigma_n", float64),("Sigma_lambda", float64),("chi_s_ell", float64),("chi_a_ell", float64),("Fcal_coeff", float64),("Fcal_0", float64),("Fcal_2", float64),("Fcal_3", float64),("Fcal_4", float64),("Fcal_SQ_4", float64),("Fcal_SO_3", float64),("E_0", float64),("E_2", float64),("E_4", float64),("E_SQ_4", float64),("E_SO_3", float64),("xHat_i", float64),("yHat_i", float64),("zHat_i", float64),("M1_i", float64),("M2_i", float64),("v_i", float64),("S_chi1_i", float64),("S_chi2_i", float64),("rfrak_chi1_x_i", float64),("rfrak_chi1_y_i", float64),("rfrak_chi2_x_i", float64),("rfrak_chi2_y_i", float64),("rfrak_frame_x_i", float64),("rfrak_frame_y_i", float64),("rfrak_frame_z_i", float64),("t_PN1",float64),("t_PN2",float64),("StepsPerOrbit",int32),("ForwardInTime",boolean),("tol",float64),("MinStep",float64),("EvolveSpin1",boolean),("EvolveSpin2",boolean),("terminal1",boolean),("terminal2",boolean),("NotForward",boolean)]
@jitclass(spec)
class TaylorTn_2p0PN_Q : 
    def __init__(self, xHat_i, yHat_i, zHat_i, M1_i, M2_i, v_i, S_chi1_i, S_chi2_i, rfrak_chi1_x_i, rfrak_chi1_y_i, rfrak_chi2_x_i, rfrak_chi2_y_i, rfrak_frame_x_i, rfrak_frame_y_i, rfrak_frame_z_i, t_PN1=False, t_PN2=False, StepsPerOrbit=32, ForwardInTime=True, tol=1e-12, MinStep=1e-7) :
        self.t_PN1,self.t_PN2,self.StepsPerOrbit,self.ForwardInTime,self.tol,self.MinStep=t_PN1, t_PN2, StepsPerOrbit, ForwardInTime, tol, MinStep
        terminal1=True
        terminal2=True
        NotForward=True
        self.xHat_i=xHat_i
        self.yHat_i=yHat_i
        self.zHat_i=zHat_i
        self.M1_i=M1_i
        self.M2_i=M2_i
        self.v_i=v_i
        self.S_chi1_i=S_chi1_i
        self.S_chi2_i=S_chi2_i
        self.rfrak_chi1_x_i=rfrak_chi1_x_i
        self.rfrak_chi1_y_i=rfrak_chi1_y_i
        self.rfrak_chi2_x_i=rfrak_chi2_x_i
        self.rfrak_chi2_y_i=rfrak_chi2_y_i
        self.rfrak_frame_x_i=rfrak_frame_x_i
        self.rfrak_frame_y_i=rfrak_frame_y_i
        self.rfrak_frame_z_i=rfrak_frame_z_i
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
        chi1chi1=dot(chiVec1.vector,chiVec1.vector)
        chi1chi2=dot(chiVec1.vector,chiVec2.vector)
        chi2chi2=dot(chiVec2.vector,chiVec2.vector)
        chi1_n=dot(chiVec1.vector,nHat.vector)
        chi1_lambda=dot(chiVec1.vector,lambdaHat.vector)
        chi1_ell=dot(chiVec1.vector,ellHat.vector)
        chi2_n=dot(chiVec2.vector,nHat.vector)
        chi2_lambda=dot(chiVec2.vector,lambdaHat.vector)
        chi2_ell=dot(chiVec2.vector,ellHat.vector)
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
        self.xHat=xHat
        self.yHat=yHat
        self.zHat=zHat
        self.M1=M1
        self.M2=M2
        self.v=v
        self.S_chi1=S_chi1
        self.S_chi2=S_chi2
        self.rfrak_chi1_x=rfrak_chi1_x
        self.rfrak_chi1_y=rfrak_chi1_y
        self.rfrak_chi2_x=rfrak_chi2_x
        self.rfrak_chi2_y=rfrak_chi2_y
        self.rfrak_frame_x=rfrak_frame_x
        self.rfrak_frame_y=rfrak_frame_y
        self.rfrak_frame_z=rfrak_frame_z
        self.M=M
        self.delta=delta
        self.nu=nu
        self.R=R
        self.nHat=nHat
        self.lambdaHat=lambdaHat
        self.ellHat=ellHat
        self.R_S1=R_S1
        self.R_S2=R_S2
        self.chiVec1=chiVec1
        self.chiVec2=chiVec2
        self.chi1chi1=chi1chi1
        self.chi1chi2=chi1chi2
        self.chi2chi2=chi2chi2
        self.chi1_n=chi1_n
        self.chi1_lambda=chi1_lambda
        self.chi1_ell=chi1_ell
        self.chi2_n=chi2_n
        self.chi2_lambda=chi2_lambda
        self.chi2_ell=chi2_ell
        self.S_ell=S_ell
        self.S_n=S_n
        self.S_lambda=S_lambda
        self.Sigma_ell=Sigma_ell
        self.Sigma_n=Sigma_n
        self.Sigma_lambda=Sigma_lambda
        self.chi_s_ell=chi_s_ell
        self.chi_a_ell=chi_a_ell
        self.Fcal_coeff=Fcal_coeff
        self.Fcal_0=Fcal_0
        self.Fcal_2=Fcal_2
        self.Fcal_3=Fcal_3
        self.Fcal_4=Fcal_4
        self.Fcal_SQ_4=Fcal_SQ_4
        self.Fcal_SO_3=Fcal_SO_3
        self.E_0=E_0
        self.E_2=E_2
        self.E_4=E_4
        self.E_SQ_4=E_SQ_4
        self.E_SO_3=E_SO_3
        Phi=0.0
        EvolveSpin1=(S_chi1*conjugate(S_chi1)).norm>1e-12
        EvolveSpin2=(S_chi2*conjugate(S_chi2)).norm>1e-12
        self.v, self.rfrak_chi1_x, self.rfrak_chi1_y, self.rfrak_chi2_x, self.rfrak_chi2_y=v, rfrak_chi1_x, rfrak_chi1_y, rfrak_chi2_x, rfrak_chi2_y
        self.rfrak_frame_x, self.rfrak_frame_y, self.rfrak_frame_z, self.Phi=rfrak_frame_x, rfrak_frame_y, rfrak_frame_z, Phi
        self.EvolveSpin1=EvolveSpin1
        self.EvolveSpin2=EvolveSpin2
        self.terminal1,self.terminal2,self.NotForward=terminal1,terminal2,NotForward

    def __Recalculate__(self,t, y):
        v = y[0]
        rfrak_chi1_x = y[1]
        rfrak_chi1_y = y[2]
        rfrak_chi2_x = y[3]
        rfrak_chi2_y = y[4]
        rfrak_frame_x = y[5]
        rfrak_frame_y = y[6]
        rfrak_frame_z = y[7]
        Phi = y[8]
        self.v, self.rfrak_chi1_x, self.rfrak_chi1_y, self.rfrak_chi2_x, self.rfrak_chi2_y=v, rfrak_chi1_x, rfrak_chi1_y, rfrak_chi2_x, rfrak_chi2_y
        self.rfrak_frame_x, self.rfrak_frame_y, self.rfrak_frame_z, self.Phi=rfrak_frame_x, rfrak_frame_y, rfrak_frame_z, Phi
        xHat=self.xHat
        yHat=self.yHat
        zHat=self.zHat
        M1=self.M1
        M2=self.M2
        v=self.v
        S_chi1=self.S_chi1
        S_chi2=self.S_chi2
        rfrak_chi1_x=self.rfrak_chi1_x
        rfrak_chi1_y=self.rfrak_chi1_y
        rfrak_chi2_x=self.rfrak_chi2_x
        rfrak_chi2_y=self.rfrak_chi2_y
        rfrak_frame_x=self.rfrak_frame_x
        rfrak_frame_y=self.rfrak_frame_y
        rfrak_frame_z=self.rfrak_frame_z
        M=self.M
        delta=self.delta
        nu=self.nu
        R=self.R
        nHat=self.nHat
        lambdaHat=self.lambdaHat
        ellHat=self.ellHat
        R_S1=self.R_S1
        R_S2=self.R_S2
        chiVec1=self.chiVec1
        chiVec2=self.chiVec2
        chi1chi1=self.chi1chi1
        chi1chi2=self.chi1chi2
        chi2chi2=self.chi2chi2
        chi1_n=self.chi1_n
        chi1_lambda=self.chi1_lambda
        chi1_ell=self.chi1_ell
        chi2_n=self.chi2_n
        chi2_lambda=self.chi2_lambda
        chi2_ell=self.chi2_ell
        S_ell=self.S_ell
        S_n=self.S_n
        S_lambda=self.S_lambda
        Sigma_ell=self.Sigma_ell
        Sigma_n=self.Sigma_n
        Sigma_lambda=self.Sigma_lambda
        chi_s_ell=self.chi_s_ell
        chi_a_ell=self.chi_a_ell
        Fcal_coeff=self.Fcal_coeff
        Fcal_0=self.Fcal_0
        Fcal_2=self.Fcal_2
        Fcal_3=self.Fcal_3
        Fcal_4=self.Fcal_4
        Fcal_SQ_4=self.Fcal_SQ_4
        Fcal_SO_3=self.Fcal_SO_3
        E_0=self.E_0
        E_2=self.E_2
        E_4=self.E_4
        E_SQ_4=self.E_SQ_4
        E_SO_3=self.E_SO_3
        R = exp(rfrak_frame_x*xHat + rfrak_frame_y*yHat + rfrak_frame_z*zHat)
        nHat = R*xHat*conjugate(R)
        lambdaHat = R*yHat*conjugate(R)
        ellHat = R*zHat*conjugate(R)
        R_S1 = exp(rfrak_chi1_x*xHat + rfrak_chi1_y*yHat)
        R_S2 = exp(rfrak_chi2_x*xHat + rfrak_chi2_y*yHat)
        chiVec1 = S_chi1*R_S1*zHat*conjugate(R_S1)*conjugate(S_chi1)
        chiVec2 = S_chi2*R_S2*zHat*conjugate(R_S2)*conjugate(S_chi2)
        chi1_n = dot(chiVec1.vector,nHat.vector)
        chi1_lambda = dot(chiVec1.vector,lambdaHat.vector)
        chi1_ell = dot(chiVec1.vector,ellHat.vector)
        chi2_n = dot(chiVec2.vector,nHat.vector)
        chi2_lambda = dot(chiVec2.vector,lambdaHat.vector)
        chi2_ell = dot(chiVec2.vector,ellHat.vector)
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
        self.R=R
        self.nHat=nHat
        self.lambdaHat=lambdaHat
        self.ellHat=ellHat
        self.R_S1=R_S1
        self.R_S2=R_S2
        self.chiVec1=chiVec1
        self.chiVec2=chiVec2
        self.chi1_n=chi1_n
        self.chi1_lambda=chi1_lambda
        self.chi1_ell=chi1_ell
        self.chi2_n=chi2_n
        self.chi2_lambda=chi2_lambda
        self.chi2_ell=chi2_ell
        self.S_ell=S_ell
        self.S_n=S_n
        self.S_lambda=S_lambda
        self.Sigma_ell=Sigma_ell
        self.Sigma_n=Sigma_n
        self.Sigma_lambda=Sigma_lambda
        self.chi_s_ell=chi_s_ell
        self.chi_a_ell=chi_a_ell
        self.Fcal_coeff=Fcal_coeff
        self.Fcal_SQ_4=Fcal_SQ_4
        self.Fcal_SO_3=Fcal_SO_3
        self.E_SQ_4=E_SQ_4
        self.E_SO_3=E_SO_3

    def OmegaVec_chiVec_1(self):
        xHat=self.xHat
        yHat=self.yHat
        zHat=self.zHat
        M1=self.M1
        M2=self.M2
        v=self.v
        S_chi1=self.S_chi1
        S_chi2=self.S_chi2
        rfrak_chi1_x=self.rfrak_chi1_x
        rfrak_chi1_y=self.rfrak_chi1_y
        rfrak_chi2_x=self.rfrak_chi2_x
        rfrak_chi2_y=self.rfrak_chi2_y
        rfrak_frame_x=self.rfrak_frame_x
        rfrak_frame_y=self.rfrak_frame_y
        rfrak_frame_z=self.rfrak_frame_z
        M=self.M
        delta=self.delta
        nu=self.nu
        R=self.R
        nHat=self.nHat
        lambdaHat=self.lambdaHat
        ellHat=self.ellHat
        R_S1=self.R_S1
        R_S2=self.R_S2
        chiVec1=self.chiVec1
        chiVec2=self.chiVec2
        chi1_n=self.chi1_n
        chi1_lambda=self.chi1_lambda
        chi1_ell=self.chi1_ell
        chi2_n=self.chi2_n
        chi2_lambda=self.chi2_lambda
        chi2_ell=self.chi2_ell
        S_ell=self.S_ell
        S_n=self.S_n
        S_lambda=self.S_lambda
        Sigma_ell=self.Sigma_ell
        Sigma_n=self.Sigma_n
        Sigma_lambda=self.Sigma_lambda
        Omega1_coeff = v**5/M
        return Omega1_coeff*(ellHat*(-0.75*delta + 0.5*nu + v**2*(delta*(0.625*nu - 0.5625) + nu*(1.25 - 0.0416666666666667*nu) + v**2*(delta*(nu*(4.875 - 0.15625*nu) - 0.84375) + nu*(nu*(-0.0208333333333333*nu - 3.28125) + 0.1875) + 0.84375) + 0.5625) + 0.75) + nHat*v*(3.0*chi1_n*nu + 3.0*M2**2*chi2_n/M**2) - M2**2*chiVec2*v/M**2)

    def OmegaVec_chiVec_2(self):
        xHat=self.xHat
        yHat=self.yHat
        zHat=self.zHat
        M1=self.M1
        M2=self.M2
        v=self.v
        S_chi1=self.S_chi1
        S_chi2=self.S_chi2
        rfrak_chi1_x=self.rfrak_chi1_x
        rfrak_chi1_y=self.rfrak_chi1_y
        rfrak_chi2_x=self.rfrak_chi2_x
        rfrak_chi2_y=self.rfrak_chi2_y
        rfrak_frame_x=self.rfrak_frame_x
        rfrak_frame_y=self.rfrak_frame_y
        rfrak_frame_z=self.rfrak_frame_z
        M=self.M
        delta=self.delta
        nu=self.nu
        R=self.R
        nHat=self.nHat
        lambdaHat=self.lambdaHat
        ellHat=self.ellHat
        R_S1=self.R_S1
        R_S2=self.R_S2
        chiVec1=self.chiVec1
        chiVec2=self.chiVec2
        chi1_n=self.chi1_n
        chi1_lambda=self.chi1_lambda
        chi1_ell=self.chi1_ell
        chi2_n=self.chi2_n
        chi2_lambda=self.chi2_lambda
        chi2_ell=self.chi2_ell
        S_ell=self.S_ell
        S_n=self.S_n
        S_lambda=self.S_lambda
        Sigma_ell=self.Sigma_ell
        Sigma_n=self.Sigma_n
        Sigma_lambda=self.Sigma_lambda
        Omega2_coeff = v**5/M
        return Omega2_coeff*(ellHat*(0.75*delta + 0.5*nu + v**2*(delta*(0.5625 - 0.625*nu) + nu*(1.25 - 0.0416666666666667*nu) + v**2*(delta*(nu*(0.15625*nu - 4.875) + 0.84375) + nu*(nu*(-0.0208333333333333*nu - 3.28125) + 0.1875) + 0.84375) + 0.5625) + 0.75) + nHat*v*(3.0*chi2_n*nu + 3.0*M1**2*chi1_n/M**2) - M1**2*chiVec1*v/M**2)

    def OmegaVec(self):
        xHat=self.xHat
        yHat=self.yHat
        zHat=self.zHat
        M1=self.M1
        M2=self.M2
        v=self.v
        S_chi1=self.S_chi1
        S_chi2=self.S_chi2
        rfrak_chi1_x=self.rfrak_chi1_x
        rfrak_chi1_y=self.rfrak_chi1_y
        rfrak_chi2_x=self.rfrak_chi2_x
        rfrak_chi2_y=self.rfrak_chi2_y
        rfrak_frame_x=self.rfrak_frame_x
        rfrak_frame_y=self.rfrak_frame_y
        rfrak_frame_z=self.rfrak_frame_z
        M=self.M
        delta=self.delta
        nu=self.nu
        R=self.R
        nHat=self.nHat
        lambdaHat=self.lambdaHat
        ellHat=self.ellHat
        R_S1=self.R_S1
        R_S2=self.R_S2
        chiVec1=self.chiVec1
        chiVec2=self.chiVec2
        chi1_n=self.chi1_n
        chi1_lambda=self.chi1_lambda
        chi1_ell=self.chi1_ell
        chi2_n=self.chi2_n
        chi2_lambda=self.chi2_lambda
        chi2_ell=self.chi2_ell
        S_ell=self.S_ell
        S_n=self.S_n
        S_lambda=self.S_lambda
        Sigma_ell=self.Sigma_ell
        Sigma_n=self.Sigma_n
        Sigma_lambda=self.Sigma_lambda
        gamma_PN_4 = 1.0 - 5.41666666666667*nu
        gamma_PN_2 = 1.0 - 0.333333333333333*nu
        a_ell_4 = S_n*(5.77777777777778*nu**2 + 14.75*nu + 1.5) + Sigma_n*delta*(2.83333333333333*nu**2 + 9.125*nu + 1.5)
        gamma_PN_0 = 1.00000000000000
        gamma_PN_3 = (1.66666666666667*S_ell + Sigma_ell*delta)/M**2
        a_ell_0 = 7.0*S_n + 3.0*Sigma_n*delta
        a_ell_2 = S_n*(-9.66666666666667*nu - 10.0) + Sigma_n*delta*(-4.5*nu - 6.0)
        return ellHat*v**3/M + nHat*v**6*(a_ell_0 + v**2*(a_ell_2 + a_ell_4*v**2))*(gamma_PN_0 + v**2*(gamma_PN_2 + v*(gamma_PN_3 + gamma_PN_4*v)))/M**3

    def OrbitalAngularMomentum(self):
        xHat=self.xHat
        yHat=self.yHat
        zHat=self.zHat
        M1=self.M1
        M2=self.M2
        v=self.v
        S_chi1=self.S_chi1
        S_chi2=self.S_chi2
        rfrak_chi1_x=self.rfrak_chi1_x
        rfrak_chi1_y=self.rfrak_chi1_y
        rfrak_chi2_x=self.rfrak_chi2_x
        rfrak_chi2_y=self.rfrak_chi2_y
        rfrak_frame_x=self.rfrak_frame_x
        rfrak_frame_y=self.rfrak_frame_y
        rfrak_frame_z=self.rfrak_frame_z
        M=self.M
        delta=self.delta
        nu=self.nu
        R=self.R
        nHat=self.nHat
        lambdaHat=self.lambdaHat
        ellHat=self.ellHat
        R_S1=self.R_S1
        R_S2=self.R_S2
        chiVec1=self.chiVec1
        chiVec2=self.chiVec2
        chi1_n=self.chi1_n
        chi1_lambda=self.chi1_lambda
        chi1_ell=self.chi1_ell
        chi2_n=self.chi2_n
        chi2_lambda=self.chi2_lambda
        chi2_ell=self.chi2_ell
        S_ell=self.S_ell
        S_n=self.S_n
        S_lambda=self.S_lambda
        Sigma_ell=self.Sigma_ell
        Sigma_n=self.Sigma_n
        Sigma_lambda=self.Sigma_lambda
        L_SO_3 = 0.166666666666667*ellHat*(-35.0*S_ell - 15.0*Sigma_ell*delta)/M**2 + lambdaHat*(-3.0*S_lambda - Sigma_lambda*delta)/M**2 + 0.5*nHat*(S_n + Sigma_n*delta)/M**2
        L_2 = ellHat*(0.166666666666667*nu + 1.5)
        L_4 = ellHat*(0.0416666666666667*nu**2 - 2.375*nu + 3.375)
        L_0 = ellHat
        L_coeff = M**2*nu/v
        return L_coeff*(L_0 + v**2*(L_2 + v*(L_4*v + L_SO_3)))

               

    #Evolve PN
    def __TaylorT1__(self,t, y):
        self.__Recalculate__(t, y)
        xHat=self.xHat
        yHat=self.yHat
        zHat=self.zHat
        M1=self.M1
        M2=self.M2
        v=self.v
        S_chi1=self.S_chi1
        S_chi2=self.S_chi2
        rfrak_chi1_x=self.rfrak_chi1_x
        rfrak_chi1_y=self.rfrak_chi1_y
        rfrak_chi2_x=self.rfrak_chi2_x
        rfrak_chi2_y=self.rfrak_chi2_y
        rfrak_frame_x=self.rfrak_frame_x
        rfrak_frame_y=self.rfrak_frame_y
        rfrak_frame_z=self.rfrak_frame_z
        M=self.M
        delta=self.delta
        nu=self.nu
        R=self.R
        nHat=self.nHat
        lambdaHat=self.lambdaHat
        ellHat=self.ellHat
        R_S1=self.R_S1
        R_S2=self.R_S2
        chiVec1=self.chiVec1
        chiVec2=self.chiVec2
        chi1chi1=self.chi1chi1
        chi1chi2=self.chi1chi2
        chi2chi2=self.chi2chi2
        chi1_n=self.chi1_n
        chi1_lambda=self.chi1_lambda
        chi1_ell=self.chi1_ell
        chi2_n=self.chi2_n
        chi2_lambda=self.chi2_lambda
        chi2_ell=self.chi2_ell
        S_ell=self.S_ell
        S_n=self.S_n
        S_lambda=self.S_lambda
        Sigma_ell=self.Sigma_ell
        Sigma_n=self.Sigma_n
        Sigma_lambda=self.Sigma_lambda
        chi_s_ell=self.chi_s_ell
        chi_a_ell=self.chi_a_ell
        Fcal_coeff=self.Fcal_coeff
        Fcal_0=self.Fcal_0
        Fcal_2=self.Fcal_2
        Fcal_3=self.Fcal_3
        Fcal_4=self.Fcal_4
        Fcal_SQ_4=self.Fcal_SQ_4
        Fcal_SO_3=self.Fcal_SO_3
        E_0=self.E_0
        E_2=self.E_2
        E_4=self.E_4
        E_SQ_4=self.E_SQ_4
        E_SO_3=self.E_SO_3
        v, rfrak_chi1_x, rfrak_chi1_y, rfrak_chi2_x, rfrak_chi2_y=self.v, self.rfrak_chi1_x, self.rfrak_chi1_y, self.rfrak_chi2_x, self.rfrak_chi2_y
        rfrak_frame_x, rfrak_frame_y, rfrak_frame_z, Phi=self.rfrak_frame_x, self.rfrak_frame_y, self.rfrak_frame_z, self.Phi
        dvdt_T1=0
        if v>=1.0 and self.NotForward:
            print("Beyond domain of PN validity, this is a good way to terminate.")
            self.terminal1=False
        Flux = Fcal_coeff*(Fcal_0 + v**2*(Fcal_2 + v*(Fcal_3 + Fcal_SO_3 + v*(Fcal_4 + Fcal_SQ_4))))
        dEdv = -0.5*M*nu*v*(2.0*E_0 + v**2*(4.0*E_2 + v*(5.0*E_SO_3 + 6.0*v*(E_4 + E_SQ_4))))
        Absorption = 0
        dvdt_T1 = (-Absorption - Flux)/dEdv
        if dvdt_T1<1.0e-12 and self.NotForward:
            print("v is decreasing, which is not an uncommon way to stop.")
            self.terminal2=False
        return self.CommonRHS(dvdt_T1, y)

    def CommonRHS(self,dvdt, y):
        dydt=np.zeros(9)
        rfrak_frame=np.zeros(3)
        rfrak_frame[0] = y[5]
        rfrak_frame[1] = y[6]
        rfrak_frame[2] = y[7]
        rfrakdot_frame = FrameFromAngularVelocityIntegrand(rfrak_frame, self.OmegaVec().vector)
        dydt[0] = dvdt
        if(self.EvolveSpin1):
            dydt[1], dydt[2]= FrameFromAngularVelocity_2D_Integrand(y[1], y[2],\
                (self.S_chi1.inverse*self.OmegaVec_chiVec_1()*self.S_chi1).vector)
        else:
            dydt[1] = 0.0
            dydt[2] = 0.0
        if(self.EvolveSpin2):
            dydt[3], dydt[4] = FrameFromAngularVelocity_2D_Integrand(y[3], y[4],\
                (self.S_chi2.inverse*self.OmegaVec_chiVec_2()*self.S_chi2).vector)
        else:
            dydt[3] = 0.0
            dydt[4] = 0.0
        dydt[5] = rfrakdot_frame[0]
        dydt[6] = rfrakdot_frame[1]
        dydt[7] = rfrakdot_frame[2]
        dydt[8] = self.v*self.v*self.v/self.M

        return dydt

    # Prepare for integration
    def Evolution(self): 
        v_i,rfrak_chi1_x_i,rfrak_chi1_y_i,rfrak_chi2_x_i,rfrak_chi2_y_i,rfrak_frame_x_i,rfrak_frame_y_i,rfrak_frame_z_i\
            =self.v_i,self.rfrak_chi1_x_i,self.rfrak_chi1_y_i,self.rfrak_chi2_x_i,self.rfrak_chi2_y_i,self.rfrak_frame_x_i,self.rfrak_frame_y_i,self.rfrak_frame_z_i
        def terminate(t,y):
            return 1.0*self.terminal1*self.terminal2
        terminate.terminal=True
        TMerger=5.0/(256.0*self.nu*v_i**8)
        TEnd=TMerger
        if self.t_PN2:
            TEnd=self.t_PN2
        time=[0.0]
        while time[-1]<TEnd and 2*self.M*(256*self.nu*(TMerger-time[-1])/5)**(3/8)/self.StepsPerOrbit>self.MinStep:
            time.append(time[-1]+2*self.M*(256*self.nu*(TMerger-time[-1])/5)**(3/8)/self.StepsPerOrbit)
        time=np.delete(time, -1)
        
        # Integrate
        yy=solve_ivp(self.__TaylorT1__, [time[0],time[-1]], [v_i,rfrak_chi1_x_i,\
            rfrak_chi1_y_i,rfrak_chi2_x_i,rfrak_chi2_y_i,rfrak_frame_x_i,\
            rfrak_frame_y_i,rfrak_frame_z_i,0.0], method='DOP853',\
            t_eval=time, dense_output=True, events=terminate, rtol=self.tol, atol=self.tol)           
        if self.ForwardInTime:
            self.NotForward=False
            time=[0.0]
            TStart=-3*TMerger
            if self.t_PN1:
                TStart=self.t_PN1
            while time[-1]>TStart:
                time.append(time[-1]-2*self.M*(256*self.nu*(TMerger-time[-1])/5)**(3/8)/self.StepsPerOrbit)
            yyForward=solve_ivp(self.__TaylorT1__, [time[0],time[-1]], [v_i,rfrak_chi1_x_i,\
                rfrak_chi1_y_i,rfrak_chi2_x_i,rfrak_chi2_y_i,rfrak_frame_x_i,\
                rfrak_frame_y_i,rfrak_frame_z_i,0.0], method='DOP853',\
                t_eval=time, dense_output=True, rtol=self.tol, atol=self.tol)
            yy.t=np.append(yyForward.t[1:][::-1],yy.t)
            data=np.empty((9,len(yy.t)))
            for i in range(9):
                data[i]=np.append(yyForward.y[i][1:][::-1],yy.y[i])
            yy.y=data
            
        return yy

spec=[("xHat", float64),("yHat", float64),("zHat", float64),("M1", float64),("M2", float64),("v", float64),("S_chi1", float64),("S_chi2", float64),("rfrak_chi1_x", float64),("rfrak_chi1_y", float64),("rfrak_chi2_x", float64),("rfrak_chi2_y", float64),("rfrak_frame_x", float64),("rfrak_frame_y", float64),("rfrak_frame_z", float64),("M", float64),("delta", float64),("nu", float64),("R", float64),("nHat", float64),("lambdaHat", float64),("ellHat", float64),("R_S1", float64),("R_S2", float64),("chiVec1", float64),("chiVec2", float64),("chi1chi1", float64),("chi1chi2", float64),("chi2chi2", float64),("chi1_n", float64),("chi1_lambda", float64),("chi1_ell", float64),("chi2_n", float64),("chi2_lambda", float64),("chi2_ell", float64),("S_ell", float64),("S_n", float64),("S_lambda", float64),("Sigma_ell", float64),("Sigma_n", float64),("Sigma_lambda", float64),("chi_s_ell", float64),("chi_a_ell", float64),("Fcal_coeff", float64),("Fcal_0", float64),("Fcal_2", float64),("Fcal_3", float64),("Fcal_4", float64),("Fcal_5", float64),("Fcal_SQ_4", float64),("Fcal_SO_3", float64),("Fcal_SO_5", float64),("E_0", float64),("E_2", float64),("E_4", float64),("E_SQ_4", float64),("E_SO_3", float64),("E_SO_5", float64),("xHat_i", float64),("yHat_i", float64),("zHat_i", float64),("M1_i", float64),("M2_i", float64),("v_i", float64),("S_chi1_i", float64),("S_chi2_i", float64),("rfrak_chi1_x_i", float64),("rfrak_chi1_y_i", float64),("rfrak_chi2_x_i", float64),("rfrak_chi2_y_i", float64),("rfrak_frame_x_i", float64),("rfrak_frame_y_i", float64),("rfrak_frame_z_i", float64),("t_PN1",float64),("t_PN2",float64),("StepsPerOrbit",int32),("ForwardInTime",boolean),("tol",float64),("MinStep",float64),("EvolveSpin1",boolean),("EvolveSpin2",boolean),("terminal1",boolean),("terminal2",boolean),("NotForward",boolean)]
@jitclass(spec)
class TaylorTn_2p5PN_Q : 
    def __init__(self, xHat_i, yHat_i, zHat_i, M1_i, M2_i, v_i, S_chi1_i, S_chi2_i, rfrak_chi1_x_i, rfrak_chi1_y_i, rfrak_chi2_x_i, rfrak_chi2_y_i, rfrak_frame_x_i, rfrak_frame_y_i, rfrak_frame_z_i, t_PN1=False, t_PN2=False, StepsPerOrbit=32, ForwardInTime=True, tol=1e-12, MinStep=1e-7) :
        self.t_PN1,self.t_PN2,self.StepsPerOrbit,self.ForwardInTime,self.tol,self.MinStep=t_PN1, t_PN2, StepsPerOrbit, ForwardInTime, tol, MinStep
        terminal1=True
        terminal2=True
        NotForward=True
        self.xHat_i=xHat_i
        self.yHat_i=yHat_i
        self.zHat_i=zHat_i
        self.M1_i=M1_i
        self.M2_i=M2_i
        self.v_i=v_i
        self.S_chi1_i=S_chi1_i
        self.S_chi2_i=S_chi2_i
        self.rfrak_chi1_x_i=rfrak_chi1_x_i
        self.rfrak_chi1_y_i=rfrak_chi1_y_i
        self.rfrak_chi2_x_i=rfrak_chi2_x_i
        self.rfrak_chi2_y_i=rfrak_chi2_y_i
        self.rfrak_frame_x_i=rfrak_frame_x_i
        self.rfrak_frame_y_i=rfrak_frame_y_i
        self.rfrak_frame_z_i=rfrak_frame_z_i
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
        chi1chi1=dot(chiVec1.vector,chiVec1.vector)
        chi1chi2=dot(chiVec1.vector,chiVec2.vector)
        chi2chi2=dot(chiVec2.vector,chiVec2.vector)
        chi1_n=dot(chiVec1.vector,nHat.vector)
        chi1_lambda=dot(chiVec1.vector,lambdaHat.vector)
        chi1_ell=dot(chiVec1.vector,ellHat.vector)
        chi2_n=dot(chiVec2.vector,nHat.vector)
        chi2_lambda=dot(chiVec2.vector,lambdaHat.vector)
        chi2_ell=dot(chiVec2.vector,ellHat.vector)
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
        self.xHat=xHat
        self.yHat=yHat
        self.zHat=zHat
        self.M1=M1
        self.M2=M2
        self.v=v
        self.S_chi1=S_chi1
        self.S_chi2=S_chi2
        self.rfrak_chi1_x=rfrak_chi1_x
        self.rfrak_chi1_y=rfrak_chi1_y
        self.rfrak_chi2_x=rfrak_chi2_x
        self.rfrak_chi2_y=rfrak_chi2_y
        self.rfrak_frame_x=rfrak_frame_x
        self.rfrak_frame_y=rfrak_frame_y
        self.rfrak_frame_z=rfrak_frame_z
        self.M=M
        self.delta=delta
        self.nu=nu
        self.R=R
        self.nHat=nHat
        self.lambdaHat=lambdaHat
        self.ellHat=ellHat
        self.R_S1=R_S1
        self.R_S2=R_S2
        self.chiVec1=chiVec1
        self.chiVec2=chiVec2
        self.chi1chi1=chi1chi1
        self.chi1chi2=chi1chi2
        self.chi2chi2=chi2chi2
        self.chi1_n=chi1_n
        self.chi1_lambda=chi1_lambda
        self.chi1_ell=chi1_ell
        self.chi2_n=chi2_n
        self.chi2_lambda=chi2_lambda
        self.chi2_ell=chi2_ell
        self.S_ell=S_ell
        self.S_n=S_n
        self.S_lambda=S_lambda
        self.Sigma_ell=Sigma_ell
        self.Sigma_n=Sigma_n
        self.Sigma_lambda=Sigma_lambda
        self.chi_s_ell=chi_s_ell
        self.chi_a_ell=chi_a_ell
        self.Fcal_coeff=Fcal_coeff
        self.Fcal_0=Fcal_0
        self.Fcal_2=Fcal_2
        self.Fcal_3=Fcal_3
        self.Fcal_4=Fcal_4
        self.Fcal_5=Fcal_5
        self.Fcal_SQ_4=Fcal_SQ_4
        self.Fcal_SO_3=Fcal_SO_3
        self.Fcal_SO_5=Fcal_SO_5
        self.E_0=E_0
        self.E_2=E_2
        self.E_4=E_4
        self.E_SQ_4=E_SQ_4
        self.E_SO_3=E_SO_3
        self.E_SO_5=E_SO_5
        Phi=0.0
        EvolveSpin1=(S_chi1*conjugate(S_chi1)).norm>1e-12
        EvolveSpin2=(S_chi2*conjugate(S_chi2)).norm>1e-12
        self.v, self.rfrak_chi1_x, self.rfrak_chi1_y, self.rfrak_chi2_x, self.rfrak_chi2_y=v, rfrak_chi1_x, rfrak_chi1_y, rfrak_chi2_x, rfrak_chi2_y
        self.rfrak_frame_x, self.rfrak_frame_y, self.rfrak_frame_z, self.Phi=rfrak_frame_x, rfrak_frame_y, rfrak_frame_z, Phi
        self.EvolveSpin1=EvolveSpin1
        self.EvolveSpin2=EvolveSpin2
        self.terminal1,self.terminal2,self.NotForward=terminal1,terminal2,NotForward

    def __Recalculate__(self,t, y):
        v = y[0]
        rfrak_chi1_x = y[1]
        rfrak_chi1_y = y[2]
        rfrak_chi2_x = y[3]
        rfrak_chi2_y = y[4]
        rfrak_frame_x = y[5]
        rfrak_frame_y = y[6]
        rfrak_frame_z = y[7]
        Phi = y[8]
        self.v, self.rfrak_chi1_x, self.rfrak_chi1_y, self.rfrak_chi2_x, self.rfrak_chi2_y=v, rfrak_chi1_x, rfrak_chi1_y, rfrak_chi2_x, rfrak_chi2_y
        self.rfrak_frame_x, self.rfrak_frame_y, self.rfrak_frame_z, self.Phi=rfrak_frame_x, rfrak_frame_y, rfrak_frame_z, Phi
        xHat=self.xHat
        yHat=self.yHat
        zHat=self.zHat
        M1=self.M1
        M2=self.M2
        v=self.v
        S_chi1=self.S_chi1
        S_chi2=self.S_chi2
        rfrak_chi1_x=self.rfrak_chi1_x
        rfrak_chi1_y=self.rfrak_chi1_y
        rfrak_chi2_x=self.rfrak_chi2_x
        rfrak_chi2_y=self.rfrak_chi2_y
        rfrak_frame_x=self.rfrak_frame_x
        rfrak_frame_y=self.rfrak_frame_y
        rfrak_frame_z=self.rfrak_frame_z
        M=self.M
        delta=self.delta
        nu=self.nu
        R=self.R
        nHat=self.nHat
        lambdaHat=self.lambdaHat
        ellHat=self.ellHat
        R_S1=self.R_S1
        R_S2=self.R_S2
        chiVec1=self.chiVec1
        chiVec2=self.chiVec2
        chi1chi1=self.chi1chi1
        chi1chi2=self.chi1chi2
        chi2chi2=self.chi2chi2
        chi1_n=self.chi1_n
        chi1_lambda=self.chi1_lambda
        chi1_ell=self.chi1_ell
        chi2_n=self.chi2_n
        chi2_lambda=self.chi2_lambda
        chi2_ell=self.chi2_ell
        S_ell=self.S_ell
        S_n=self.S_n
        S_lambda=self.S_lambda
        Sigma_ell=self.Sigma_ell
        Sigma_n=self.Sigma_n
        Sigma_lambda=self.Sigma_lambda
        chi_s_ell=self.chi_s_ell
        chi_a_ell=self.chi_a_ell
        Fcal_coeff=self.Fcal_coeff
        Fcal_0=self.Fcal_0
        Fcal_2=self.Fcal_2
        Fcal_3=self.Fcal_3
        Fcal_4=self.Fcal_4
        Fcal_5=self.Fcal_5
        Fcal_SQ_4=self.Fcal_SQ_4
        Fcal_SO_3=self.Fcal_SO_3
        Fcal_SO_5=self.Fcal_SO_5
        E_0=self.E_0
        E_2=self.E_2
        E_4=self.E_4
        E_SQ_4=self.E_SQ_4
        E_SO_3=self.E_SO_3
        E_SO_5=self.E_SO_5
        R = exp(rfrak_frame_x*xHat + rfrak_frame_y*yHat + rfrak_frame_z*zHat)
        nHat = R*xHat*conjugate(R)
        lambdaHat = R*yHat*conjugate(R)
        ellHat = R*zHat*conjugate(R)
        R_S1 = exp(rfrak_chi1_x*xHat + rfrak_chi1_y*yHat)
        R_S2 = exp(rfrak_chi2_x*xHat + rfrak_chi2_y*yHat)
        chiVec1 = S_chi1*R_S1*zHat*conjugate(R_S1)*conjugate(S_chi1)
        chiVec2 = S_chi2*R_S2*zHat*conjugate(R_S2)*conjugate(S_chi2)
        chi1_n = dot(chiVec1.vector,nHat.vector)
        chi1_lambda = dot(chiVec1.vector,lambdaHat.vector)
        chi1_ell = dot(chiVec1.vector,ellHat.vector)
        chi2_n = dot(chiVec2.vector,nHat.vector)
        chi2_lambda = dot(chiVec2.vector,lambdaHat.vector)
        chi2_ell = dot(chiVec2.vector,ellHat.vector)
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
        self.R=R
        self.nHat=nHat
        self.lambdaHat=lambdaHat
        self.ellHat=ellHat
        self.R_S1=R_S1
        self.R_S2=R_S2
        self.chiVec1=chiVec1
        self.chiVec2=chiVec2
        self.chi1_n=chi1_n
        self.chi1_lambda=chi1_lambda
        self.chi1_ell=chi1_ell
        self.chi2_n=chi2_n
        self.chi2_lambda=chi2_lambda
        self.chi2_ell=chi2_ell
        self.S_ell=S_ell
        self.S_n=S_n
        self.S_lambda=S_lambda
        self.Sigma_ell=Sigma_ell
        self.Sigma_n=Sigma_n
        self.Sigma_lambda=Sigma_lambda
        self.chi_s_ell=chi_s_ell
        self.chi_a_ell=chi_a_ell
        self.Fcal_coeff=Fcal_coeff
        self.Fcal_SQ_4=Fcal_SQ_4
        self.Fcal_SO_3=Fcal_SO_3
        self.Fcal_SO_5=Fcal_SO_5
        self.E_SQ_4=E_SQ_4
        self.E_SO_3=E_SO_3
        self.E_SO_5=E_SO_5

    def OmegaVec_chiVec_1(self):
        xHat=self.xHat
        yHat=self.yHat
        zHat=self.zHat
        M1=self.M1
        M2=self.M2
        v=self.v
        S_chi1=self.S_chi1
        S_chi2=self.S_chi2
        rfrak_chi1_x=self.rfrak_chi1_x
        rfrak_chi1_y=self.rfrak_chi1_y
        rfrak_chi2_x=self.rfrak_chi2_x
        rfrak_chi2_y=self.rfrak_chi2_y
        rfrak_frame_x=self.rfrak_frame_x
        rfrak_frame_y=self.rfrak_frame_y
        rfrak_frame_z=self.rfrak_frame_z
        M=self.M
        delta=self.delta
        nu=self.nu
        R=self.R
        nHat=self.nHat
        lambdaHat=self.lambdaHat
        ellHat=self.ellHat
        R_S1=self.R_S1
        R_S2=self.R_S2
        chiVec1=self.chiVec1
        chiVec2=self.chiVec2
        chi1_n=self.chi1_n
        chi1_lambda=self.chi1_lambda
        chi1_ell=self.chi1_ell
        chi2_n=self.chi2_n
        chi2_lambda=self.chi2_lambda
        chi2_ell=self.chi2_ell
        S_ell=self.S_ell
        S_n=self.S_n
        S_lambda=self.S_lambda
        Sigma_ell=self.Sigma_ell
        Sigma_n=self.Sigma_n
        Sigma_lambda=self.Sigma_lambda
        Omega1_coeff = v**5/M
        return Omega1_coeff*(ellHat*(-0.75*delta + 0.5*nu + v**2*(delta*(0.625*nu - 0.5625) + nu*(1.25 - 0.0416666666666667*nu) + v**2*(delta*(nu*(4.875 - 0.15625*nu) - 0.84375) + nu*(nu*(-0.0208333333333333*nu - 3.28125) + 0.1875) + 0.84375) + 0.5625) + 0.75) + nHat*v*(3.0*chi1_n*nu + 3.0*M2**2*chi2_n/M**2) - M2**2*chiVec2*v/M**2)

    def OmegaVec_chiVec_2(self):
        xHat=self.xHat
        yHat=self.yHat
        zHat=self.zHat
        M1=self.M1
        M2=self.M2
        v=self.v
        S_chi1=self.S_chi1
        S_chi2=self.S_chi2
        rfrak_chi1_x=self.rfrak_chi1_x
        rfrak_chi1_y=self.rfrak_chi1_y
        rfrak_chi2_x=self.rfrak_chi2_x
        rfrak_chi2_y=self.rfrak_chi2_y
        rfrak_frame_x=self.rfrak_frame_x
        rfrak_frame_y=self.rfrak_frame_y
        rfrak_frame_z=self.rfrak_frame_z
        M=self.M
        delta=self.delta
        nu=self.nu
        R=self.R
        nHat=self.nHat
        lambdaHat=self.lambdaHat
        ellHat=self.ellHat
        R_S1=self.R_S1
        R_S2=self.R_S2
        chiVec1=self.chiVec1
        chiVec2=self.chiVec2
        chi1_n=self.chi1_n
        chi1_lambda=self.chi1_lambda
        chi1_ell=self.chi1_ell
        chi2_n=self.chi2_n
        chi2_lambda=self.chi2_lambda
        chi2_ell=self.chi2_ell
        S_ell=self.S_ell
        S_n=self.S_n
        S_lambda=self.S_lambda
        Sigma_ell=self.Sigma_ell
        Sigma_n=self.Sigma_n
        Sigma_lambda=self.Sigma_lambda
        Omega2_coeff = v**5/M
        return Omega2_coeff*(ellHat*(0.75*delta + 0.5*nu + v**2*(delta*(0.5625 - 0.625*nu) + nu*(1.25 - 0.0416666666666667*nu) + v**2*(delta*(nu*(0.15625*nu - 4.875) + 0.84375) + nu*(nu*(-0.0208333333333333*nu - 3.28125) + 0.1875) + 0.84375) + 0.5625) + 0.75) + nHat*v*(3.0*chi2_n*nu + 3.0*M1**2*chi1_n/M**2) - M1**2*chiVec1*v/M**2)

    def OmegaVec(self):
        xHat=self.xHat
        yHat=self.yHat
        zHat=self.zHat
        M1=self.M1
        M2=self.M2
        v=self.v
        S_chi1=self.S_chi1
        S_chi2=self.S_chi2
        rfrak_chi1_x=self.rfrak_chi1_x
        rfrak_chi1_y=self.rfrak_chi1_y
        rfrak_chi2_x=self.rfrak_chi2_x
        rfrak_chi2_y=self.rfrak_chi2_y
        rfrak_frame_x=self.rfrak_frame_x
        rfrak_frame_y=self.rfrak_frame_y
        rfrak_frame_z=self.rfrak_frame_z
        M=self.M
        delta=self.delta
        nu=self.nu
        R=self.R
        nHat=self.nHat
        lambdaHat=self.lambdaHat
        ellHat=self.ellHat
        R_S1=self.R_S1
        R_S2=self.R_S2
        chiVec1=self.chiVec1
        chiVec2=self.chiVec2
        chi1_n=self.chi1_n
        chi1_lambda=self.chi1_lambda
        chi1_ell=self.chi1_ell
        chi2_n=self.chi2_n
        chi2_lambda=self.chi2_lambda
        chi2_ell=self.chi2_ell
        S_ell=self.S_ell
        S_n=self.S_n
        S_lambda=self.S_lambda
        Sigma_ell=self.Sigma_ell
        Sigma_n=self.Sigma_n
        Sigma_lambda=self.Sigma_lambda
        gamma_PN_4 = 1.0 - 5.41666666666667*nu
        gamma_PN_5 = (S_ell*(0.888888888888889*nu + 3.33333333333333) + 2.0*Sigma_ell*delta)/M**2
        gamma_PN_2 = 1.0 - 0.333333333333333*nu
        a_ell_4 = S_n*(5.77777777777778*nu**2 + 14.75*nu + 1.5) + Sigma_n*delta*(2.83333333333333*nu**2 + 9.125*nu + 1.5)
        gamma_PN_0 = 1.00000000000000
        gamma_PN_3 = (1.66666666666667*S_ell + Sigma_ell*delta)/M**2
        a_ell_0 = 7.0*S_n + 3.0*Sigma_n*delta
        a_ell_2 = S_n*(-9.66666666666667*nu - 10.0) + Sigma_n*delta*(-4.5*nu - 6.0)
        return ellHat*v**3/M + nHat*v**6*(a_ell_0 + v**2*(a_ell_2 + a_ell_4*v**2))*(gamma_PN_0 + v**2*(gamma_PN_2 + v*(gamma_PN_3 + v*(gamma_PN_4 + gamma_PN_5*v))))/M**3

    def OrbitalAngularMomentum(self):
        xHat=self.xHat
        yHat=self.yHat
        zHat=self.zHat
        M1=self.M1
        M2=self.M2
        v=self.v
        S_chi1=self.S_chi1
        S_chi2=self.S_chi2
        rfrak_chi1_x=self.rfrak_chi1_x
        rfrak_chi1_y=self.rfrak_chi1_y
        rfrak_chi2_x=self.rfrak_chi2_x
        rfrak_chi2_y=self.rfrak_chi2_y
        rfrak_frame_x=self.rfrak_frame_x
        rfrak_frame_y=self.rfrak_frame_y
        rfrak_frame_z=self.rfrak_frame_z
        M=self.M
        delta=self.delta
        nu=self.nu
        R=self.R
        nHat=self.nHat
        lambdaHat=self.lambdaHat
        ellHat=self.ellHat
        R_S1=self.R_S1
        R_S2=self.R_S2
        chiVec1=self.chiVec1
        chiVec2=self.chiVec2
        chi1_n=self.chi1_n
        chi1_lambda=self.chi1_lambda
        chi1_ell=self.chi1_ell
        chi2_n=self.chi2_n
        chi2_lambda=self.chi2_lambda
        chi2_ell=self.chi2_ell
        S_ell=self.S_ell
        S_n=self.S_n
        S_lambda=self.S_lambda
        Sigma_ell=self.Sigma_ell
        Sigma_n=self.Sigma_n
        Sigma_lambda=self.Sigma_lambda
        L_SO_3 = 0.166666666666667*ellHat*(-35.0*S_ell - 15.0*Sigma_ell*delta)/M**2 + lambdaHat*(-3.0*S_lambda - Sigma_lambda*delta)/M**2 + 0.5*nHat*(S_n + Sigma_n*delta)/M**2
        L_SO_5 = 0.0138888888888889*ellHat*(427.0*S_ell*nu - 693.0*S_ell + 210.0*Sigma_ell*delta*nu - 189.0*Sigma_ell*delta)/M**2 + 0.166666666666667*lambdaHat*(18.0*S_lambda*nu - 21.0*S_lambda + 8.0*Sigma_lambda*delta*nu - 3.0*Sigma_lambda*delta)/M**2 + 0.0416666666666667*nHat*(-19.0*S_n*nu + 33.0*S_n - 10.0*Sigma_n*delta*nu + 33.0*Sigma_n*delta)/M**2
        L_2 = ellHat*(0.166666666666667*nu + 1.5)
        L_4 = ellHat*(0.0416666666666667*nu**2 - 2.375*nu + 3.375)
        L_0 = ellHat
        L_coeff = M**2*nu/v
        return L_coeff*(L_0 + v**2*(L_2 + v*(L_SO_3 + v*(L_4 + L_SO_5*v))))

               

    #Evolve PN
    def __TaylorT1__(self,t, y):
        self.__Recalculate__(t, y)
        xHat=self.xHat
        yHat=self.yHat
        zHat=self.zHat
        M1=self.M1
        M2=self.M2
        v=self.v
        S_chi1=self.S_chi1
        S_chi2=self.S_chi2
        rfrak_chi1_x=self.rfrak_chi1_x
        rfrak_chi1_y=self.rfrak_chi1_y
        rfrak_chi2_x=self.rfrak_chi2_x
        rfrak_chi2_y=self.rfrak_chi2_y
        rfrak_frame_x=self.rfrak_frame_x
        rfrak_frame_y=self.rfrak_frame_y
        rfrak_frame_z=self.rfrak_frame_z
        M=self.M
        delta=self.delta
        nu=self.nu
        R=self.R
        nHat=self.nHat
        lambdaHat=self.lambdaHat
        ellHat=self.ellHat
        R_S1=self.R_S1
        R_S2=self.R_S2
        chiVec1=self.chiVec1
        chiVec2=self.chiVec2
        chi1chi1=self.chi1chi1
        chi1chi2=self.chi1chi2
        chi2chi2=self.chi2chi2
        chi1_n=self.chi1_n
        chi1_lambda=self.chi1_lambda
        chi1_ell=self.chi1_ell
        chi2_n=self.chi2_n
        chi2_lambda=self.chi2_lambda
        chi2_ell=self.chi2_ell
        S_ell=self.S_ell
        S_n=self.S_n
        S_lambda=self.S_lambda
        Sigma_ell=self.Sigma_ell
        Sigma_n=self.Sigma_n
        Sigma_lambda=self.Sigma_lambda
        chi_s_ell=self.chi_s_ell
        chi_a_ell=self.chi_a_ell
        Fcal_coeff=self.Fcal_coeff
        Fcal_0=self.Fcal_0
        Fcal_2=self.Fcal_2
        Fcal_3=self.Fcal_3
        Fcal_4=self.Fcal_4
        Fcal_5=self.Fcal_5
        Fcal_SQ_4=self.Fcal_SQ_4
        Fcal_SO_3=self.Fcal_SO_3
        Fcal_SO_5=self.Fcal_SO_5
        E_0=self.E_0
        E_2=self.E_2
        E_4=self.E_4
        E_SQ_4=self.E_SQ_4
        E_SO_3=self.E_SO_3
        E_SO_5=self.E_SO_5
        v, rfrak_chi1_x, rfrak_chi1_y, rfrak_chi2_x, rfrak_chi2_y=self.v, self.rfrak_chi1_x, self.rfrak_chi1_y, self.rfrak_chi2_x, self.rfrak_chi2_y
        rfrak_frame_x, rfrak_frame_y, rfrak_frame_z, Phi=self.rfrak_frame_x, self.rfrak_frame_y, self.rfrak_frame_z, self.Phi
        dvdt_T1=0
        if v>=1.0 and self.NotForward:
            print("Beyond domain of PN validity, this is a good way to terminate.")
            self.terminal1=False
        Flux = Fcal_coeff*(Fcal_0 + v**2*(Fcal_2 + v*(Fcal_3 + Fcal_SO_3 + v*(Fcal_4 + Fcal_SQ_4 + v*(Fcal_5 + Fcal_SO_5)))))
        dEdv = -0.5*M*nu*v*(2.0*E_0 + v**2*(4.0*E_2 + v*(5.0*E_SO_3 + v*(6.0*E_4 + 7.0*E_SO_5*v + 6.0*E_SQ_4))))
        Absorption = 0
        dvdt_T1 = (-Absorption - Flux)/dEdv
        if dvdt_T1<1.0e-12 and self.NotForward:
            print("v is decreasing, which is not an uncommon way to stop.")
            self.terminal2=False
        return self.CommonRHS(dvdt_T1, y)

    def CommonRHS(self,dvdt, y):
        dydt=np.zeros(9)
        rfrak_frame=np.zeros(3)
        rfrak_frame[0] = y[5]
        rfrak_frame[1] = y[6]
        rfrak_frame[2] = y[7]
        rfrakdot_frame = FrameFromAngularVelocityIntegrand(rfrak_frame, self.OmegaVec().vector)
        dydt[0] = dvdt
        if(self.EvolveSpin1):
            dydt[1], dydt[2]= FrameFromAngularVelocity_2D_Integrand(y[1], y[2],\
                (self.S_chi1.inverse*self.OmegaVec_chiVec_1()*self.S_chi1).vector)
        else:
            dydt[1] = 0.0
            dydt[2] = 0.0
        if(self.EvolveSpin2):
            dydt[3], dydt[4] = FrameFromAngularVelocity_2D_Integrand(y[3], y[4],\
                (self.S_chi2.inverse*self.OmegaVec_chiVec_2()*self.S_chi2).vector)
        else:
            dydt[3] = 0.0
            dydt[4] = 0.0
        dydt[5] = rfrakdot_frame[0]
        dydt[6] = rfrakdot_frame[1]
        dydt[7] = rfrakdot_frame[2]
        dydt[8] = self.v*self.v*self.v/self.M

        return dydt

    # Prepare for integration
    def Evolution(self): 
        v_i,rfrak_chi1_x_i,rfrak_chi1_y_i,rfrak_chi2_x_i,rfrak_chi2_y_i,rfrak_frame_x_i,rfrak_frame_y_i,rfrak_frame_z_i\
            =self.v_i,self.rfrak_chi1_x_i,self.rfrak_chi1_y_i,self.rfrak_chi2_x_i,self.rfrak_chi2_y_i,self.rfrak_frame_x_i,self.rfrak_frame_y_i,self.rfrak_frame_z_i
        def terminate(t,y):
            return 1.0*self.terminal1*self.terminal2
        terminate.terminal=True
        TMerger=5.0/(256.0*self.nu*v_i**8)
        TEnd=TMerger
        if self.t_PN2:
            TEnd=self.t_PN2
        time=[0.0]
        while time[-1]<TEnd and 2*self.M*(256*self.nu*(TMerger-time[-1])/5)**(3/8)/self.StepsPerOrbit>self.MinStep:
            time.append(time[-1]+2*self.M*(256*self.nu*(TMerger-time[-1])/5)**(3/8)/self.StepsPerOrbit)
        time=np.delete(time, -1)
        
        # Integrate
        yy=solve_ivp(self.__TaylorT1__, [time[0],time[-1]], [v_i,rfrak_chi1_x_i,\
            rfrak_chi1_y_i,rfrak_chi2_x_i,rfrak_chi2_y_i,rfrak_frame_x_i,\
            rfrak_frame_y_i,rfrak_frame_z_i,0.0], method='DOP853',\
            t_eval=time, dense_output=True, events=terminate, rtol=self.tol, atol=self.tol)           
        if self.ForwardInTime:
            self.NotForward=False
            time=[0.0]
            TStart=-3*TMerger
            if self.t_PN1:
                TStart=self.t_PN1
            while time[-1]>TStart:
                time.append(time[-1]-2*self.M*(256*self.nu*(TMerger-time[-1])/5)**(3/8)/self.StepsPerOrbit)
            yyForward=solve_ivp(self.__TaylorT1__, [time[0],time[-1]], [v_i,rfrak_chi1_x_i,\
                rfrak_chi1_y_i,rfrak_chi2_x_i,rfrak_chi2_y_i,rfrak_frame_x_i,\
                rfrak_frame_y_i,rfrak_frame_z_i,0.0], method='DOP853',\
                t_eval=time, dense_output=True, rtol=self.tol, atol=self.tol)
            yy.t=np.append(yyForward.t[1:][::-1],yy.t)
            data=np.empty((9,len(yy.t)))
            for i in range(9):
                data[i]=np.append(yyForward.y[i][1:][::-1],yy.y[i])
            yy.y=data
            
        return yy

spec=[("xHat", float64),("yHat", float64),("zHat", float64),("M1", float64),("M2", float64),("v", float64),("S_chi1", float64),("S_chi2", float64),("rfrak_chi1_x", float64),("rfrak_chi1_y", float64),("rfrak_chi2_x", float64),("rfrak_chi2_y", float64),("rfrak_frame_x", float64),("rfrak_frame_y", float64),("rfrak_frame_z", float64),("M", float64),("delta", float64),("nu", float64),("R", float64),("nHat", float64),("lambdaHat", float64),("ellHat", float64),("R_S1", float64),("R_S2", float64),("chiVec1", float64),("chiVec2", float64),("chi1chi1", float64),("chi1chi2", float64),("chi2chi2", float64),("chi1_n", float64),("chi1_lambda", float64),("chi1_ell", float64),("chi2_n", float64),("chi2_lambda", float64),("chi2_ell", float64),("S_ell", float64),("S_n", float64),("S_lambda", float64),("Sigma_ell", float64),("Sigma_n", float64),("Sigma_lambda", float64),("chi_s_ell", float64),("chi_a_ell", float64),("logv", float64),("Fcal_coeff", float64),("Fcal_0", float64),("Fcal_2", float64),("Fcal_3", float64),("Fcal_4", float64),("Fcal_5", float64),("Fcal_6", float64),("Fcal_lnv_6", float64),("Fcal_SQ_4", float64),("Fcal_SO_3", float64),("Fcal_SO_5", float64),("Fcal_SO_6", float64),("E_0", float64),("E_2", float64),("E_4", float64),("E_6", float64),("E_SQ_4", float64),("E_SO_3", float64),("E_SO_5", float64),("xHat_i", float64),("yHat_i", float64),("zHat_i", float64),("M1_i", float64),("M2_i", float64),("v_i", float64),("S_chi1_i", float64),("S_chi2_i", float64),("rfrak_chi1_x_i", float64),("rfrak_chi1_y_i", float64),("rfrak_chi2_x_i", float64),("rfrak_chi2_y_i", float64),("rfrak_frame_x_i", float64),("rfrak_frame_y_i", float64),("rfrak_frame_z_i", float64),("t_PN1",float64),("t_PN2",float64),("StepsPerOrbit",int32),("ForwardInTime",boolean),("tol",float64),("MinStep",float64),("EvolveSpin1",boolean),("EvolveSpin2",boolean),("terminal1",boolean),("terminal2",boolean),("NotForward",boolean)]
@jitclass(spec)
class TaylorTn_3p0PN_Q : 
    def __init__(self, xHat_i, yHat_i, zHat_i, M1_i, M2_i, v_i, S_chi1_i, S_chi2_i, rfrak_chi1_x_i, rfrak_chi1_y_i, rfrak_chi2_x_i, rfrak_chi2_y_i, rfrak_frame_x_i, rfrak_frame_y_i, rfrak_frame_z_i, t_PN1=False, t_PN2=False, StepsPerOrbit=32, ForwardInTime=True, tol=1e-12, MinStep=1e-7) :
        self.t_PN1,self.t_PN2,self.StepsPerOrbit,self.ForwardInTime,self.tol,self.MinStep=t_PN1, t_PN2, StepsPerOrbit, ForwardInTime, tol, MinStep
        terminal1=True
        terminal2=True
        NotForward=True
        self.xHat_i=xHat_i
        self.yHat_i=yHat_i
        self.zHat_i=zHat_i
        self.M1_i=M1_i
        self.M2_i=M2_i
        self.v_i=v_i
        self.S_chi1_i=S_chi1_i
        self.S_chi2_i=S_chi2_i
        self.rfrak_chi1_x_i=rfrak_chi1_x_i
        self.rfrak_chi1_y_i=rfrak_chi1_y_i
        self.rfrak_chi2_x_i=rfrak_chi2_x_i
        self.rfrak_chi2_y_i=rfrak_chi2_y_i
        self.rfrak_frame_x_i=rfrak_frame_x_i
        self.rfrak_frame_y_i=rfrak_frame_y_i
        self.rfrak_frame_z_i=rfrak_frame_z_i
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
        chi1chi1=dot(chiVec1.vector,chiVec1.vector)
        chi1chi2=dot(chiVec1.vector,chiVec2.vector)
        chi2chi2=dot(chiVec2.vector,chiVec2.vector)
        chi1_n=dot(chiVec1.vector,nHat.vector)
        chi1_lambda=dot(chiVec1.vector,lambdaHat.vector)
        chi1_ell=dot(chiVec1.vector,ellHat.vector)
        chi2_n=dot(chiVec2.vector,nHat.vector)
        chi2_lambda=dot(chiVec2.vector,lambdaHat.vector)
        chi2_ell=dot(chiVec2.vector,ellHat.vector)
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
        self.xHat=xHat
        self.yHat=yHat
        self.zHat=zHat
        self.M1=M1
        self.M2=M2
        self.v=v
        self.S_chi1=S_chi1
        self.S_chi2=S_chi2
        self.rfrak_chi1_x=rfrak_chi1_x
        self.rfrak_chi1_y=rfrak_chi1_y
        self.rfrak_chi2_x=rfrak_chi2_x
        self.rfrak_chi2_y=rfrak_chi2_y
        self.rfrak_frame_x=rfrak_frame_x
        self.rfrak_frame_y=rfrak_frame_y
        self.rfrak_frame_z=rfrak_frame_z
        self.M=M
        self.delta=delta
        self.nu=nu
        self.R=R
        self.nHat=nHat
        self.lambdaHat=lambdaHat
        self.ellHat=ellHat
        self.R_S1=R_S1
        self.R_S2=R_S2
        self.chiVec1=chiVec1
        self.chiVec2=chiVec2
        self.chi1chi1=chi1chi1
        self.chi1chi2=chi1chi2
        self.chi2chi2=chi2chi2
        self.chi1_n=chi1_n
        self.chi1_lambda=chi1_lambda
        self.chi1_ell=chi1_ell
        self.chi2_n=chi2_n
        self.chi2_lambda=chi2_lambda
        self.chi2_ell=chi2_ell
        self.S_ell=S_ell
        self.S_n=S_n
        self.S_lambda=S_lambda
        self.Sigma_ell=Sigma_ell
        self.Sigma_n=Sigma_n
        self.Sigma_lambda=Sigma_lambda
        self.chi_s_ell=chi_s_ell
        self.chi_a_ell=chi_a_ell
        self.logv=logv
        self.Fcal_coeff=Fcal_coeff
        self.Fcal_0=Fcal_0
        self.Fcal_2=Fcal_2
        self.Fcal_3=Fcal_3
        self.Fcal_4=Fcal_4
        self.Fcal_5=Fcal_5
        self.Fcal_6=Fcal_6
        self.Fcal_lnv_6=Fcal_lnv_6
        self.Fcal_SQ_4=Fcal_SQ_4
        self.Fcal_SO_3=Fcal_SO_3
        self.Fcal_SO_5=Fcal_SO_5
        self.Fcal_SO_6=Fcal_SO_6
        self.E_0=E_0
        self.E_2=E_2
        self.E_4=E_4
        self.E_6=E_6
        self.E_SQ_4=E_SQ_4
        self.E_SO_3=E_SO_3
        self.E_SO_5=E_SO_5
        Phi=0.0
        EvolveSpin1=(S_chi1*conjugate(S_chi1)).norm>1e-12
        EvolveSpin2=(S_chi2*conjugate(S_chi2)).norm>1e-12
        self.v, self.rfrak_chi1_x, self.rfrak_chi1_y, self.rfrak_chi2_x, self.rfrak_chi2_y=v, rfrak_chi1_x, rfrak_chi1_y, rfrak_chi2_x, rfrak_chi2_y
        self.rfrak_frame_x, self.rfrak_frame_y, self.rfrak_frame_z, self.Phi=rfrak_frame_x, rfrak_frame_y, rfrak_frame_z, Phi
        self.EvolveSpin1=EvolveSpin1
        self.EvolveSpin2=EvolveSpin2
        self.terminal1,self.terminal2,self.NotForward=terminal1,terminal2,NotForward

    def __Recalculate__(self,t, y):
        v = y[0]
        rfrak_chi1_x = y[1]
        rfrak_chi1_y = y[2]
        rfrak_chi2_x = y[3]
        rfrak_chi2_y = y[4]
        rfrak_frame_x = y[5]
        rfrak_frame_y = y[6]
        rfrak_frame_z = y[7]
        Phi = y[8]
        self.v, self.rfrak_chi1_x, self.rfrak_chi1_y, self.rfrak_chi2_x, self.rfrak_chi2_y=v, rfrak_chi1_x, rfrak_chi1_y, rfrak_chi2_x, rfrak_chi2_y
        self.rfrak_frame_x, self.rfrak_frame_y, self.rfrak_frame_z, self.Phi=rfrak_frame_x, rfrak_frame_y, rfrak_frame_z, Phi
        xHat=self.xHat
        yHat=self.yHat
        zHat=self.zHat
        M1=self.M1
        M2=self.M2
        v=self.v
        S_chi1=self.S_chi1
        S_chi2=self.S_chi2
        rfrak_chi1_x=self.rfrak_chi1_x
        rfrak_chi1_y=self.rfrak_chi1_y
        rfrak_chi2_x=self.rfrak_chi2_x
        rfrak_chi2_y=self.rfrak_chi2_y
        rfrak_frame_x=self.rfrak_frame_x
        rfrak_frame_y=self.rfrak_frame_y
        rfrak_frame_z=self.rfrak_frame_z
        M=self.M
        delta=self.delta
        nu=self.nu
        R=self.R
        nHat=self.nHat
        lambdaHat=self.lambdaHat
        ellHat=self.ellHat
        R_S1=self.R_S1
        R_S2=self.R_S2
        chiVec1=self.chiVec1
        chiVec2=self.chiVec2
        chi1chi1=self.chi1chi1
        chi1chi2=self.chi1chi2
        chi2chi2=self.chi2chi2
        chi1_n=self.chi1_n
        chi1_lambda=self.chi1_lambda
        chi1_ell=self.chi1_ell
        chi2_n=self.chi2_n
        chi2_lambda=self.chi2_lambda
        chi2_ell=self.chi2_ell
        S_ell=self.S_ell
        S_n=self.S_n
        S_lambda=self.S_lambda
        Sigma_ell=self.Sigma_ell
        Sigma_n=self.Sigma_n
        Sigma_lambda=self.Sigma_lambda
        chi_s_ell=self.chi_s_ell
        chi_a_ell=self.chi_a_ell
        logv=self.logv
        Fcal_coeff=self.Fcal_coeff
        Fcal_0=self.Fcal_0
        Fcal_2=self.Fcal_2
        Fcal_3=self.Fcal_3
        Fcal_4=self.Fcal_4
        Fcal_5=self.Fcal_5
        Fcal_6=self.Fcal_6
        Fcal_lnv_6=self.Fcal_lnv_6
        Fcal_SQ_4=self.Fcal_SQ_4
        Fcal_SO_3=self.Fcal_SO_3
        Fcal_SO_5=self.Fcal_SO_5
        Fcal_SO_6=self.Fcal_SO_6
        E_0=self.E_0
        E_2=self.E_2
        E_4=self.E_4
        E_6=self.E_6
        E_SQ_4=self.E_SQ_4
        E_SO_3=self.E_SO_3
        E_SO_5=self.E_SO_5
        R = exp(rfrak_frame_x*xHat + rfrak_frame_y*yHat + rfrak_frame_z*zHat)
        nHat = R*xHat*conjugate(R)
        lambdaHat = R*yHat*conjugate(R)
        ellHat = R*zHat*conjugate(R)
        R_S1 = exp(rfrak_chi1_x*xHat + rfrak_chi1_y*yHat)
        R_S2 = exp(rfrak_chi2_x*xHat + rfrak_chi2_y*yHat)
        chiVec1 = S_chi1*R_S1*zHat*conjugate(R_S1)*conjugate(S_chi1)
        chiVec2 = S_chi2*R_S2*zHat*conjugate(R_S2)*conjugate(S_chi2)
        chi1_n = dot(chiVec1.vector,nHat.vector)
        chi1_lambda = dot(chiVec1.vector,lambdaHat.vector)
        chi1_ell = dot(chiVec1.vector,ellHat.vector)
        chi2_n = dot(chiVec2.vector,nHat.vector)
        chi2_lambda = dot(chiVec2.vector,lambdaHat.vector)
        chi2_ell = dot(chiVec2.vector,ellHat.vector)
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
        self.R=R
        self.nHat=nHat
        self.lambdaHat=lambdaHat
        self.ellHat=ellHat
        self.R_S1=R_S1
        self.R_S2=R_S2
        self.chiVec1=chiVec1
        self.chiVec2=chiVec2
        self.chi1_n=chi1_n
        self.chi1_lambda=chi1_lambda
        self.chi1_ell=chi1_ell
        self.chi2_n=chi2_n
        self.chi2_lambda=chi2_lambda
        self.chi2_ell=chi2_ell
        self.S_ell=S_ell
        self.S_n=S_n
        self.S_lambda=S_lambda
        self.Sigma_ell=Sigma_ell
        self.Sigma_n=Sigma_n
        self.Sigma_lambda=Sigma_lambda
        self.chi_s_ell=chi_s_ell
        self.chi_a_ell=chi_a_ell
        self.logv=logv
        self.Fcal_coeff=Fcal_coeff
        self.Fcal_SQ_4=Fcal_SQ_4
        self.Fcal_SO_3=Fcal_SO_3
        self.Fcal_SO_5=Fcal_SO_5
        self.Fcal_SO_6=Fcal_SO_6
        self.E_SQ_4=E_SQ_4
        self.E_SO_3=E_SO_3
        self.E_SO_5=E_SO_5

    def OmegaVec_chiVec_1(self):
        xHat=self.xHat
        yHat=self.yHat
        zHat=self.zHat
        M1=self.M1
        M2=self.M2
        v=self.v
        S_chi1=self.S_chi1
        S_chi2=self.S_chi2
        rfrak_chi1_x=self.rfrak_chi1_x
        rfrak_chi1_y=self.rfrak_chi1_y
        rfrak_chi2_x=self.rfrak_chi2_x
        rfrak_chi2_y=self.rfrak_chi2_y
        rfrak_frame_x=self.rfrak_frame_x
        rfrak_frame_y=self.rfrak_frame_y
        rfrak_frame_z=self.rfrak_frame_z
        M=self.M
        delta=self.delta
        nu=self.nu
        R=self.R
        nHat=self.nHat
        lambdaHat=self.lambdaHat
        ellHat=self.ellHat
        R_S1=self.R_S1
        R_S2=self.R_S2
        chiVec1=self.chiVec1
        chiVec2=self.chiVec2
        chi1_n=self.chi1_n
        chi1_lambda=self.chi1_lambda
        chi1_ell=self.chi1_ell
        chi2_n=self.chi2_n
        chi2_lambda=self.chi2_lambda
        chi2_ell=self.chi2_ell
        S_ell=self.S_ell
        S_n=self.S_n
        S_lambda=self.S_lambda
        Sigma_ell=self.Sigma_ell
        Sigma_n=self.Sigma_n
        Sigma_lambda=self.Sigma_lambda
        Omega1_coeff = v**5/M
        return Omega1_coeff*(ellHat*(-0.75*delta + 0.5*nu + v**2*(delta*(0.625*nu - 0.5625) + nu*(1.25 - 0.0416666666666667*nu) + v**2*(delta*(nu*(4.875 - 0.15625*nu) - 0.84375) + nu*(nu*(-0.0208333333333333*nu - 3.28125) + 0.1875) + 0.84375) + 0.5625) + 0.75) + nHat*v*(3.0*chi1_n*nu + 3.0*M2**2*chi2_n/M**2) - M2**2*chiVec2*v/M**2)

    def OmegaVec_chiVec_2(self):
        xHat=self.xHat
        yHat=self.yHat
        zHat=self.zHat
        M1=self.M1
        M2=self.M2
        v=self.v
        S_chi1=self.S_chi1
        S_chi2=self.S_chi2
        rfrak_chi1_x=self.rfrak_chi1_x
        rfrak_chi1_y=self.rfrak_chi1_y
        rfrak_chi2_x=self.rfrak_chi2_x
        rfrak_chi2_y=self.rfrak_chi2_y
        rfrak_frame_x=self.rfrak_frame_x
        rfrak_frame_y=self.rfrak_frame_y
        rfrak_frame_z=self.rfrak_frame_z
        M=self.M
        delta=self.delta
        nu=self.nu
        R=self.R
        nHat=self.nHat
        lambdaHat=self.lambdaHat
        ellHat=self.ellHat
        R_S1=self.R_S1
        R_S2=self.R_S2
        chiVec1=self.chiVec1
        chiVec2=self.chiVec2
        chi1_n=self.chi1_n
        chi1_lambda=self.chi1_lambda
        chi1_ell=self.chi1_ell
        chi2_n=self.chi2_n
        chi2_lambda=self.chi2_lambda
        chi2_ell=self.chi2_ell
        S_ell=self.S_ell
        S_n=self.S_n
        S_lambda=self.S_lambda
        Sigma_ell=self.Sigma_ell
        Sigma_n=self.Sigma_n
        Sigma_lambda=self.Sigma_lambda
        Omega2_coeff = v**5/M
        return Omega2_coeff*(ellHat*(0.75*delta + 0.5*nu + v**2*(delta*(0.5625 - 0.625*nu) + nu*(1.25 - 0.0416666666666667*nu) + v**2*(delta*(nu*(0.15625*nu - 4.875) + 0.84375) + nu*(nu*(-0.0208333333333333*nu - 3.28125) + 0.1875) + 0.84375) + 0.5625) + 0.75) + nHat*v*(3.0*chi2_n*nu + 3.0*M1**2*chi1_n/M**2) - M1**2*chiVec1*v/M**2)

    def OmegaVec(self):
        xHat=self.xHat
        yHat=self.yHat
        zHat=self.zHat
        M1=self.M1
        M2=self.M2
        v=self.v
        S_chi1=self.S_chi1
        S_chi2=self.S_chi2
        rfrak_chi1_x=self.rfrak_chi1_x
        rfrak_chi1_y=self.rfrak_chi1_y
        rfrak_chi2_x=self.rfrak_chi2_x
        rfrak_chi2_y=self.rfrak_chi2_y
        rfrak_frame_x=self.rfrak_frame_x
        rfrak_frame_y=self.rfrak_frame_y
        rfrak_frame_z=self.rfrak_frame_z
        M=self.M
        delta=self.delta
        nu=self.nu
        R=self.R
        nHat=self.nHat
        lambdaHat=self.lambdaHat
        ellHat=self.ellHat
        R_S1=self.R_S1
        R_S2=self.R_S2
        chiVec1=self.chiVec1
        chiVec2=self.chiVec2
        chi1_n=self.chi1_n
        chi1_lambda=self.chi1_lambda
        chi1_ell=self.chi1_ell
        chi2_n=self.chi2_n
        chi2_lambda=self.chi2_lambda
        chi2_ell=self.chi2_ell
        S_ell=self.S_ell
        S_n=self.S_n
        S_lambda=self.S_lambda
        Sigma_ell=self.Sigma_ell
        Sigma_n=self.Sigma_n
        Sigma_lambda=self.Sigma_lambda
        gamma_PN_6 = 0.0123456790123457*nu**3 + 6.36111111111111*nu**2 - 2.98177812235564*nu + 1.0
        gamma_PN_4 = 1.0 - 5.41666666666667*nu
        gamma_PN_5 = (S_ell*(0.888888888888889*nu + 3.33333333333333) + 2.0*Sigma_ell*delta)/M**2
        gamma_PN_2 = 1.0 - 0.333333333333333*nu
        a_ell_4 = S_n*(5.77777777777778*nu**2 + 14.75*nu + 1.5) + Sigma_n*delta*(2.83333333333333*nu**2 + 9.125*nu + 1.5)
        gamma_PN_0 = 1.00000000000000
        gamma_PN_3 = (1.66666666666667*S_ell + Sigma_ell*delta)/M**2
        a_ell_0 = 7.0*S_n + 3.0*Sigma_n*delta
        a_ell_2 = S_n*(-9.66666666666667*nu - 10.0) + Sigma_n*delta*(-4.5*nu - 6.0)
        return ellHat*v**3/M + nHat*v**6*(a_ell_0 + v**2*(a_ell_2 + a_ell_4*v**2))*(gamma_PN_0 + v**2*(gamma_PN_2 + v*(gamma_PN_3 + v*(gamma_PN_4 + v*(gamma_PN_5 + gamma_PN_6*v)))))/M**3

    def OrbitalAngularMomentum(self):
        xHat=self.xHat
        yHat=self.yHat
        zHat=self.zHat
        M1=self.M1
        M2=self.M2
        v=self.v
        S_chi1=self.S_chi1
        S_chi2=self.S_chi2
        rfrak_chi1_x=self.rfrak_chi1_x
        rfrak_chi1_y=self.rfrak_chi1_y
        rfrak_chi2_x=self.rfrak_chi2_x
        rfrak_chi2_y=self.rfrak_chi2_y
        rfrak_frame_x=self.rfrak_frame_x
        rfrak_frame_y=self.rfrak_frame_y
        rfrak_frame_z=self.rfrak_frame_z
        M=self.M
        delta=self.delta
        nu=self.nu
        R=self.R
        nHat=self.nHat
        lambdaHat=self.lambdaHat
        ellHat=self.ellHat
        R_S1=self.R_S1
        R_S2=self.R_S2
        chiVec1=self.chiVec1
        chiVec2=self.chiVec2
        chi1_n=self.chi1_n
        chi1_lambda=self.chi1_lambda
        chi1_ell=self.chi1_ell
        chi2_n=self.chi2_n
        chi2_lambda=self.chi2_lambda
        chi2_ell=self.chi2_ell
        S_ell=self.S_ell
        S_n=self.S_n
        S_lambda=self.S_lambda
        Sigma_ell=self.Sigma_ell
        Sigma_n=self.Sigma_n
        Sigma_lambda=self.Sigma_lambda
        L_6 = ellHat*(0.00540123456790123*nu**3 + 1.29166666666667*nu**2 - 30.9797035925835*nu + 8.4375)
        L_SO_3 = 0.166666666666667*ellHat*(-35.0*S_ell - 15.0*Sigma_ell*delta)/M**2 + lambdaHat*(-3.0*S_lambda - Sigma_lambda*delta)/M**2 + 0.5*nHat*(S_n + Sigma_n*delta)/M**2
        L_SO_5 = 0.0138888888888889*ellHat*(427.0*S_ell*nu - 693.0*S_ell + 210.0*Sigma_ell*delta*nu - 189.0*Sigma_ell*delta)/M**2 + 0.166666666666667*lambdaHat*(18.0*S_lambda*nu - 21.0*S_lambda + 8.0*Sigma_lambda*delta*nu - 3.0*Sigma_lambda*delta)/M**2 + 0.0416666666666667*nHat*(-19.0*S_n*nu + 33.0*S_n - 10.0*Sigma_n*delta*nu + 33.0*Sigma_n*delta)/M**2
        L_2 = ellHat*(0.166666666666667*nu + 1.5)
        L_4 = ellHat*(0.0416666666666667*nu**2 - 2.375*nu + 3.375)
        L_0 = ellHat
        L_coeff = M**2*nu/v
        return L_coeff*(L_0 + v**2*(L_2 + v*(L_SO_3 + v*(L_4 + v*(L_6*v + L_SO_5)))))

               

    #Evolve PN
    def __TaylorT1__(self,t, y):
        self.__Recalculate__(t, y)
        xHat=self.xHat
        yHat=self.yHat
        zHat=self.zHat
        M1=self.M1
        M2=self.M2
        v=self.v
        S_chi1=self.S_chi1
        S_chi2=self.S_chi2
        rfrak_chi1_x=self.rfrak_chi1_x
        rfrak_chi1_y=self.rfrak_chi1_y
        rfrak_chi2_x=self.rfrak_chi2_x
        rfrak_chi2_y=self.rfrak_chi2_y
        rfrak_frame_x=self.rfrak_frame_x
        rfrak_frame_y=self.rfrak_frame_y
        rfrak_frame_z=self.rfrak_frame_z
        M=self.M
        delta=self.delta
        nu=self.nu
        R=self.R
        nHat=self.nHat
        lambdaHat=self.lambdaHat
        ellHat=self.ellHat
        R_S1=self.R_S1
        R_S2=self.R_S2
        chiVec1=self.chiVec1
        chiVec2=self.chiVec2
        chi1chi1=self.chi1chi1
        chi1chi2=self.chi1chi2
        chi2chi2=self.chi2chi2
        chi1_n=self.chi1_n
        chi1_lambda=self.chi1_lambda
        chi1_ell=self.chi1_ell
        chi2_n=self.chi2_n
        chi2_lambda=self.chi2_lambda
        chi2_ell=self.chi2_ell
        S_ell=self.S_ell
        S_n=self.S_n
        S_lambda=self.S_lambda
        Sigma_ell=self.Sigma_ell
        Sigma_n=self.Sigma_n
        Sigma_lambda=self.Sigma_lambda
        chi_s_ell=self.chi_s_ell
        chi_a_ell=self.chi_a_ell
        logv=self.logv
        Fcal_coeff=self.Fcal_coeff
        Fcal_0=self.Fcal_0
        Fcal_2=self.Fcal_2
        Fcal_3=self.Fcal_3
        Fcal_4=self.Fcal_4
        Fcal_5=self.Fcal_5
        Fcal_6=self.Fcal_6
        Fcal_lnv_6=self.Fcal_lnv_6
        Fcal_SQ_4=self.Fcal_SQ_4
        Fcal_SO_3=self.Fcal_SO_3
        Fcal_SO_5=self.Fcal_SO_5
        Fcal_SO_6=self.Fcal_SO_6
        E_0=self.E_0
        E_2=self.E_2
        E_4=self.E_4
        E_6=self.E_6
        E_SQ_4=self.E_SQ_4
        E_SO_3=self.E_SO_3
        E_SO_5=self.E_SO_5
        v, rfrak_chi1_x, rfrak_chi1_y, rfrak_chi2_x, rfrak_chi2_y=self.v, self.rfrak_chi1_x, self.rfrak_chi1_y, self.rfrak_chi2_x, self.rfrak_chi2_y
        rfrak_frame_x, rfrak_frame_y, rfrak_frame_z, Phi=self.rfrak_frame_x, self.rfrak_frame_y, self.rfrak_frame_z, self.Phi
        dvdt_T1=0
        if v>=1.0 and self.NotForward:
            print("Beyond domain of PN validity, this is a good way to terminate.")
            self.terminal1=False
        Flux = Fcal_coeff*(Fcal_0 + v**2*(Fcal_2 + v*(Fcal_3 + Fcal_SO_3 + v*(Fcal_4 + Fcal_SQ_4 + v*(Fcal_5 + Fcal_SO_5 + v*(Fcal_6 + Fcal_SO_6 + Fcal_lnv_6*logv))))))
        dEdv = -0.5*M*nu*v*(2.0*E_0 + v**2*(4.0*E_2 + v*(5.0*E_SO_3 + v*(6.0*E_4 + 6.0*E_SQ_4 + v*(8.0*E_6*v + 7.0*E_SO_5)))))
        Absorption = 0
        dvdt_T1 = (-Absorption - Flux)/dEdv
        if dvdt_T1<1.0e-12 and self.NotForward:
            print("v is decreasing, which is not an uncommon way to stop.")
            self.terminal2=False
        return self.CommonRHS(dvdt_T1, y)

    def CommonRHS(self,dvdt, y):
        dydt=np.zeros(9)
        rfrak_frame=np.zeros(3)
        rfrak_frame[0] = y[5]
        rfrak_frame[1] = y[6]
        rfrak_frame[2] = y[7]
        rfrakdot_frame = FrameFromAngularVelocityIntegrand(rfrak_frame, self.OmegaVec().vector)
        dydt[0] = dvdt
        if(self.EvolveSpin1):
            dydt[1], dydt[2]= FrameFromAngularVelocity_2D_Integrand(y[1], y[2],\
                (self.S_chi1.inverse*self.OmegaVec_chiVec_1()*self.S_chi1).vector)
        else:
            dydt[1] = 0.0
            dydt[2] = 0.0
        if(self.EvolveSpin2):
            dydt[3], dydt[4] = FrameFromAngularVelocity_2D_Integrand(y[3], y[4],\
                (self.S_chi2.inverse*self.OmegaVec_chiVec_2()*self.S_chi2).vector)
        else:
            dydt[3] = 0.0
            dydt[4] = 0.0
        dydt[5] = rfrakdot_frame[0]
        dydt[6] = rfrakdot_frame[1]
        dydt[7] = rfrakdot_frame[2]
        dydt[8] = self.v*self.v*self.v/self.M

        return dydt

    # Prepare for integration
    def Evolution(self): 
        v_i,rfrak_chi1_x_i,rfrak_chi1_y_i,rfrak_chi2_x_i,rfrak_chi2_y_i,rfrak_frame_x_i,rfrak_frame_y_i,rfrak_frame_z_i\
            =self.v_i,self.rfrak_chi1_x_i,self.rfrak_chi1_y_i,self.rfrak_chi2_x_i,self.rfrak_chi2_y_i,self.rfrak_frame_x_i,self.rfrak_frame_y_i,self.rfrak_frame_z_i
        def terminate(t,y):
            return 1.0*self.terminal1*self.terminal2
        terminate.terminal=True
        TMerger=5.0/(256.0*self.nu*v_i**8)
        TEnd=TMerger
        if self.t_PN2:
            TEnd=self.t_PN2
        time=[0.0]
        while time[-1]<TEnd and 2*self.M*(256*self.nu*(TMerger-time[-1])/5)**(3/8)/self.StepsPerOrbit>self.MinStep:
            time.append(time[-1]+2*self.M*(256*self.nu*(TMerger-time[-1])/5)**(3/8)/self.StepsPerOrbit)
        time=np.delete(time, -1)
        
        # Integrate
        yy=solve_ivp(self.__TaylorT1__, [time[0],time[-1]], [v_i,rfrak_chi1_x_i,\
            rfrak_chi1_y_i,rfrak_chi2_x_i,rfrak_chi2_y_i,rfrak_frame_x_i,\
            rfrak_frame_y_i,rfrak_frame_z_i,0.0], method='DOP853',\
            t_eval=time, dense_output=True, events=terminate, rtol=self.tol, atol=self.tol)           
        if self.ForwardInTime:
            self.NotForward=False
            time=[0.0]
            TStart=-3*TMerger
            if self.t_PN1:
                TStart=self.t_PN1
            while time[-1]>TStart:
                time.append(time[-1]-2*self.M*(256*self.nu*(TMerger-time[-1])/5)**(3/8)/self.StepsPerOrbit)
            yyForward=solve_ivp(self.__TaylorT1__, [time[0],time[-1]], [v_i,rfrak_chi1_x_i,\
                rfrak_chi1_y_i,rfrak_chi2_x_i,rfrak_chi2_y_i,rfrak_frame_x_i,\
                rfrak_frame_y_i,rfrak_frame_z_i,0.0], method='DOP853',\
                t_eval=time, dense_output=True, rtol=self.tol, atol=self.tol)
            yy.t=np.append(yyForward.t[1:][::-1],yy.t)
            data=np.empty((9,len(yy.t)))
            for i in range(9):
                data[i]=np.append(yyForward.y[i][1:][::-1],yy.y[i])
            yy.y=data
            
        return yy

spec=[("xHat", float64),("yHat", float64),("zHat", float64),("M1", float64),("M2", float64),("v", float64),("S_chi1", float64),("S_chi2", float64),("rfrak_chi1_x", float64),("rfrak_chi1_y", float64),("rfrak_chi2_x", float64),("rfrak_chi2_y", float64),("rfrak_frame_x", float64),("rfrak_frame_y", float64),("rfrak_frame_z", float64),("M", float64),("delta", float64),("nu", float64),("R", float64),("nHat", float64),("lambdaHat", float64),("ellHat", float64),("R_S1", float64),("R_S2", float64),("chiVec1", float64),("chiVec2", float64),("chi1chi1", float64),("chi1chi2", float64),("chi2chi2", float64),("chi1_n", float64),("chi1_lambda", float64),("chi1_ell", float64),("chi2_n", float64),("chi2_lambda", float64),("chi2_ell", float64),("S_ell", float64),("S_n", float64),("S_lambda", float64),("Sigma_ell", float64),("Sigma_n", float64),("Sigma_lambda", float64),("chi_s_ell", float64),("chi_a_ell", float64),("logv", float64),("Fcal_coeff", float64),("Fcal_0", float64),("Fcal_2", float64),("Fcal_3", float64),("Fcal_4", float64),("Fcal_5", float64),("Fcal_6", float64),("Fcal_lnv_6", float64),("Fcal_7", float64),("Fcal_SQ_4", float64),("Fcal_SO_3", float64),("Fcal_SO_5", float64),("Fcal_SO_6", float64),("Fcal_SO_7", float64),("E_0", float64),("E_2", float64),("E_4", float64),("E_6", float64),("E_SQ_4", float64),("E_SO_3", float64),("E_SO_5", float64),("E_SO_7", float64),("xHat_i", float64),("yHat_i", float64),("zHat_i", float64),("M1_i", float64),("M2_i", float64),("v_i", float64),("S_chi1_i", float64),("S_chi2_i", float64),("rfrak_chi1_x_i", float64),("rfrak_chi1_y_i", float64),("rfrak_chi2_x_i", float64),("rfrak_chi2_y_i", float64),("rfrak_frame_x_i", float64),("rfrak_frame_y_i", float64),("rfrak_frame_z_i", float64),("t_PN1",float64),("t_PN2",float64),("StepsPerOrbit",int32),("ForwardInTime",boolean),("tol",float64),("MinStep",float64),("EvolveSpin1",boolean),("EvolveSpin2",boolean),("terminal1",boolean),("terminal2",boolean),("NotForward",boolean)]
@jitclass(spec)
class TaylorTn_3p5PN_Q : 
    def __init__(self, xHat_i, yHat_i, zHat_i, M1_i, M2_i, v_i, S_chi1_i, S_chi2_i, rfrak_chi1_x_i, rfrak_chi1_y_i, rfrak_chi2_x_i, rfrak_chi2_y_i, rfrak_frame_x_i, rfrak_frame_y_i, rfrak_frame_z_i, t_PN1=False, t_PN2=False, StepsPerOrbit=32, ForwardInTime=True, tol=1e-12, MinStep=1e-7) :
        self.t_PN1,self.t_PN2,self.StepsPerOrbit,self.ForwardInTime,self.tol,self.MinStep=t_PN1, t_PN2, StepsPerOrbit, ForwardInTime, tol, MinStep
        terminal1=True
        terminal2=True
        NotForward=True
        self.xHat_i=xHat_i
        self.yHat_i=yHat_i
        self.zHat_i=zHat_i
        self.M1_i=M1_i
        self.M2_i=M2_i
        self.v_i=v_i
        self.S_chi1_i=S_chi1_i
        self.S_chi2_i=S_chi2_i
        self.rfrak_chi1_x_i=rfrak_chi1_x_i
        self.rfrak_chi1_y_i=rfrak_chi1_y_i
        self.rfrak_chi2_x_i=rfrak_chi2_x_i
        self.rfrak_chi2_y_i=rfrak_chi2_y_i
        self.rfrak_frame_x_i=rfrak_frame_x_i
        self.rfrak_frame_y_i=rfrak_frame_y_i
        self.rfrak_frame_z_i=rfrak_frame_z_i
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
        chi1chi1=dot(chiVec1.vector,chiVec1.vector)
        chi1chi2=dot(chiVec1.vector,chiVec2.vector)
        chi2chi2=dot(chiVec2.vector,chiVec2.vector)
        chi1_n=dot(chiVec1.vector,nHat.vector)
        chi1_lambda=dot(chiVec1.vector,lambdaHat.vector)
        chi1_ell=dot(chiVec1.vector,ellHat.vector)
        chi2_n=dot(chiVec2.vector,nHat.vector)
        chi2_lambda=dot(chiVec2.vector,lambdaHat.vector)
        chi2_ell=dot(chiVec2.vector,ellHat.vector)
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
        self.xHat=xHat
        self.yHat=yHat
        self.zHat=zHat
        self.M1=M1
        self.M2=M2
        self.v=v
        self.S_chi1=S_chi1
        self.S_chi2=S_chi2
        self.rfrak_chi1_x=rfrak_chi1_x
        self.rfrak_chi1_y=rfrak_chi1_y
        self.rfrak_chi2_x=rfrak_chi2_x
        self.rfrak_chi2_y=rfrak_chi2_y
        self.rfrak_frame_x=rfrak_frame_x
        self.rfrak_frame_y=rfrak_frame_y
        self.rfrak_frame_z=rfrak_frame_z
        self.M=M
        self.delta=delta
        self.nu=nu
        self.R=R
        self.nHat=nHat
        self.lambdaHat=lambdaHat
        self.ellHat=ellHat
        self.R_S1=R_S1
        self.R_S2=R_S2
        self.chiVec1=chiVec1
        self.chiVec2=chiVec2
        self.chi1chi1=chi1chi1
        self.chi1chi2=chi1chi2
        self.chi2chi2=chi2chi2
        self.chi1_n=chi1_n
        self.chi1_lambda=chi1_lambda
        self.chi1_ell=chi1_ell
        self.chi2_n=chi2_n
        self.chi2_lambda=chi2_lambda
        self.chi2_ell=chi2_ell
        self.S_ell=S_ell
        self.S_n=S_n
        self.S_lambda=S_lambda
        self.Sigma_ell=Sigma_ell
        self.Sigma_n=Sigma_n
        self.Sigma_lambda=Sigma_lambda
        self.chi_s_ell=chi_s_ell
        self.chi_a_ell=chi_a_ell
        self.logv=logv
        self.Fcal_coeff=Fcal_coeff
        self.Fcal_0=Fcal_0
        self.Fcal_2=Fcal_2
        self.Fcal_3=Fcal_3
        self.Fcal_4=Fcal_4
        self.Fcal_5=Fcal_5
        self.Fcal_6=Fcal_6
        self.Fcal_lnv_6=Fcal_lnv_6
        self.Fcal_7=Fcal_7
        self.Fcal_SQ_4=Fcal_SQ_4
        self.Fcal_SO_3=Fcal_SO_3
        self.Fcal_SO_5=Fcal_SO_5
        self.Fcal_SO_6=Fcal_SO_6
        self.Fcal_SO_7=Fcal_SO_7
        self.E_0=E_0
        self.E_2=E_2
        self.E_4=E_4
        self.E_6=E_6
        self.E_SQ_4=E_SQ_4
        self.E_SO_3=E_SO_3
        self.E_SO_5=E_SO_5
        self.E_SO_7=E_SO_7
        Phi=0.0
        EvolveSpin1=(S_chi1*conjugate(S_chi1)).norm>1e-12
        EvolveSpin2=(S_chi2*conjugate(S_chi2)).norm>1e-12
        self.v, self.rfrak_chi1_x, self.rfrak_chi1_y, self.rfrak_chi2_x, self.rfrak_chi2_y=v, rfrak_chi1_x, rfrak_chi1_y, rfrak_chi2_x, rfrak_chi2_y
        self.rfrak_frame_x, self.rfrak_frame_y, self.rfrak_frame_z, self.Phi=rfrak_frame_x, rfrak_frame_y, rfrak_frame_z, Phi
        self.EvolveSpin1=EvolveSpin1
        self.EvolveSpin2=EvolveSpin2
        self.terminal1,self.terminal2,self.NotForward=terminal1,terminal2,NotForward

    def __Recalculate__(self,t, y):
        v = y[0]
        rfrak_chi1_x = y[1]
        rfrak_chi1_y = y[2]
        rfrak_chi2_x = y[3]
        rfrak_chi2_y = y[4]
        rfrak_frame_x = y[5]
        rfrak_frame_y = y[6]
        rfrak_frame_z = y[7]
        Phi = y[8]
        self.v, self.rfrak_chi1_x, self.rfrak_chi1_y, self.rfrak_chi2_x, self.rfrak_chi2_y=v, rfrak_chi1_x, rfrak_chi1_y, rfrak_chi2_x, rfrak_chi2_y
        self.rfrak_frame_x, self.rfrak_frame_y, self.rfrak_frame_z, self.Phi=rfrak_frame_x, rfrak_frame_y, rfrak_frame_z, Phi
        xHat=self.xHat
        yHat=self.yHat
        zHat=self.zHat
        M1=self.M1
        M2=self.M2
        v=self.v
        S_chi1=self.S_chi1
        S_chi2=self.S_chi2
        rfrak_chi1_x=self.rfrak_chi1_x
        rfrak_chi1_y=self.rfrak_chi1_y
        rfrak_chi2_x=self.rfrak_chi2_x
        rfrak_chi2_y=self.rfrak_chi2_y
        rfrak_frame_x=self.rfrak_frame_x
        rfrak_frame_y=self.rfrak_frame_y
        rfrak_frame_z=self.rfrak_frame_z
        M=self.M
        delta=self.delta
        nu=self.nu
        R=self.R
        nHat=self.nHat
        lambdaHat=self.lambdaHat
        ellHat=self.ellHat
        R_S1=self.R_S1
        R_S2=self.R_S2
        chiVec1=self.chiVec1
        chiVec2=self.chiVec2
        chi1chi1=self.chi1chi1
        chi1chi2=self.chi1chi2
        chi2chi2=self.chi2chi2
        chi1_n=self.chi1_n
        chi1_lambda=self.chi1_lambda
        chi1_ell=self.chi1_ell
        chi2_n=self.chi2_n
        chi2_lambda=self.chi2_lambda
        chi2_ell=self.chi2_ell
        S_ell=self.S_ell
        S_n=self.S_n
        S_lambda=self.S_lambda
        Sigma_ell=self.Sigma_ell
        Sigma_n=self.Sigma_n
        Sigma_lambda=self.Sigma_lambda
        chi_s_ell=self.chi_s_ell
        chi_a_ell=self.chi_a_ell
        logv=self.logv
        Fcal_coeff=self.Fcal_coeff
        Fcal_0=self.Fcal_0
        Fcal_2=self.Fcal_2
        Fcal_3=self.Fcal_3
        Fcal_4=self.Fcal_4
        Fcal_5=self.Fcal_5
        Fcal_6=self.Fcal_6
        Fcal_lnv_6=self.Fcal_lnv_6
        Fcal_7=self.Fcal_7
        Fcal_SQ_4=self.Fcal_SQ_4
        Fcal_SO_3=self.Fcal_SO_3
        Fcal_SO_5=self.Fcal_SO_5
        Fcal_SO_6=self.Fcal_SO_6
        Fcal_SO_7=self.Fcal_SO_7
        E_0=self.E_0
        E_2=self.E_2
        E_4=self.E_4
        E_6=self.E_6
        E_SQ_4=self.E_SQ_4
        E_SO_3=self.E_SO_3
        E_SO_5=self.E_SO_5
        E_SO_7=self.E_SO_7
        R = exp(rfrak_frame_x*xHat + rfrak_frame_y*yHat + rfrak_frame_z*zHat)
        nHat = R*xHat*conjugate(R)
        lambdaHat = R*yHat*conjugate(R)
        ellHat = R*zHat*conjugate(R)
        R_S1 = exp(rfrak_chi1_x*xHat + rfrak_chi1_y*yHat)
        R_S2 = exp(rfrak_chi2_x*xHat + rfrak_chi2_y*yHat)
        chiVec1 = S_chi1*R_S1*zHat*conjugate(R_S1)*conjugate(S_chi1)
        chiVec2 = S_chi2*R_S2*zHat*conjugate(R_S2)*conjugate(S_chi2)
        chi1_n = dot(chiVec1.vector,nHat.vector)
        chi1_lambda = dot(chiVec1.vector,lambdaHat.vector)
        chi1_ell = dot(chiVec1.vector,ellHat.vector)
        chi2_n = dot(chiVec2.vector,nHat.vector)
        chi2_lambda = dot(chiVec2.vector,lambdaHat.vector)
        chi2_ell = dot(chiVec2.vector,ellHat.vector)
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
        self.R=R
        self.nHat=nHat
        self.lambdaHat=lambdaHat
        self.ellHat=ellHat
        self.R_S1=R_S1
        self.R_S2=R_S2
        self.chiVec1=chiVec1
        self.chiVec2=chiVec2
        self.chi1_n=chi1_n
        self.chi1_lambda=chi1_lambda
        self.chi1_ell=chi1_ell
        self.chi2_n=chi2_n
        self.chi2_lambda=chi2_lambda
        self.chi2_ell=chi2_ell
        self.S_ell=S_ell
        self.S_n=S_n
        self.S_lambda=S_lambda
        self.Sigma_ell=Sigma_ell
        self.Sigma_n=Sigma_n
        self.Sigma_lambda=Sigma_lambda
        self.chi_s_ell=chi_s_ell
        self.chi_a_ell=chi_a_ell
        self.logv=logv
        self.Fcal_coeff=Fcal_coeff
        self.Fcal_SQ_4=Fcal_SQ_4
        self.Fcal_SO_3=Fcal_SO_3
        self.Fcal_SO_5=Fcal_SO_5
        self.Fcal_SO_6=Fcal_SO_6
        self.Fcal_SO_7=Fcal_SO_7
        self.E_SQ_4=E_SQ_4
        self.E_SO_3=E_SO_3
        self.E_SO_5=E_SO_5
        self.E_SO_7=E_SO_7

    def OmegaVec_chiVec_1(self):
        xHat=self.xHat
        yHat=self.yHat
        zHat=self.zHat
        M1=self.M1
        M2=self.M2
        v=self.v
        S_chi1=self.S_chi1
        S_chi2=self.S_chi2
        rfrak_chi1_x=self.rfrak_chi1_x
        rfrak_chi1_y=self.rfrak_chi1_y
        rfrak_chi2_x=self.rfrak_chi2_x
        rfrak_chi2_y=self.rfrak_chi2_y
        rfrak_frame_x=self.rfrak_frame_x
        rfrak_frame_y=self.rfrak_frame_y
        rfrak_frame_z=self.rfrak_frame_z
        M=self.M
        delta=self.delta
        nu=self.nu
        R=self.R
        nHat=self.nHat
        lambdaHat=self.lambdaHat
        ellHat=self.ellHat
        R_S1=self.R_S1
        R_S2=self.R_S2
        chiVec1=self.chiVec1
        chiVec2=self.chiVec2
        chi1_n=self.chi1_n
        chi1_lambda=self.chi1_lambda
        chi1_ell=self.chi1_ell
        chi2_n=self.chi2_n
        chi2_lambda=self.chi2_lambda
        chi2_ell=self.chi2_ell
        S_ell=self.S_ell
        S_n=self.S_n
        S_lambda=self.S_lambda
        Sigma_ell=self.Sigma_ell
        Sigma_n=self.Sigma_n
        Sigma_lambda=self.Sigma_lambda
        Omega1_coeff = v**5/M
        return Omega1_coeff*(ellHat*(-0.75*delta + 0.5*nu + v**2*(delta*(0.625*nu - 0.5625) + nu*(1.25 - 0.0416666666666667*nu) + v**2*(delta*(nu*(4.875 - 0.15625*nu) - 0.84375) + nu*(nu*(-0.0208333333333333*nu - 3.28125) + 0.1875) + 0.84375) + 0.5625) + 0.75) + nHat*v*(3.0*chi1_n*nu + 3.0*M2**2*chi2_n/M**2) - M2**2*chiVec2*v/M**2)

    def OmegaVec_chiVec_2(self):
        xHat=self.xHat
        yHat=self.yHat
        zHat=self.zHat
        M1=self.M1
        M2=self.M2
        v=self.v
        S_chi1=self.S_chi1
        S_chi2=self.S_chi2
        rfrak_chi1_x=self.rfrak_chi1_x
        rfrak_chi1_y=self.rfrak_chi1_y
        rfrak_chi2_x=self.rfrak_chi2_x
        rfrak_chi2_y=self.rfrak_chi2_y
        rfrak_frame_x=self.rfrak_frame_x
        rfrak_frame_y=self.rfrak_frame_y
        rfrak_frame_z=self.rfrak_frame_z
        M=self.M
        delta=self.delta
        nu=self.nu
        R=self.R
        nHat=self.nHat
        lambdaHat=self.lambdaHat
        ellHat=self.ellHat
        R_S1=self.R_S1
        R_S2=self.R_S2
        chiVec1=self.chiVec1
        chiVec2=self.chiVec2
        chi1_n=self.chi1_n
        chi1_lambda=self.chi1_lambda
        chi1_ell=self.chi1_ell
        chi2_n=self.chi2_n
        chi2_lambda=self.chi2_lambda
        chi2_ell=self.chi2_ell
        S_ell=self.S_ell
        S_n=self.S_n
        S_lambda=self.S_lambda
        Sigma_ell=self.Sigma_ell
        Sigma_n=self.Sigma_n
        Sigma_lambda=self.Sigma_lambda
        Omega2_coeff = v**5/M
        return Omega2_coeff*(ellHat*(0.75*delta + 0.5*nu + v**2*(delta*(0.5625 - 0.625*nu) + nu*(1.25 - 0.0416666666666667*nu) + v**2*(delta*(nu*(0.15625*nu - 4.875) + 0.84375) + nu*(nu*(-0.0208333333333333*nu - 3.28125) + 0.1875) + 0.84375) + 0.5625) + 0.75) + nHat*v*(3.0*chi2_n*nu + 3.0*M1**2*chi1_n/M**2) - M1**2*chiVec1*v/M**2)

    def OmegaVec(self):
        xHat=self.xHat
        yHat=self.yHat
        zHat=self.zHat
        M1=self.M1
        M2=self.M2
        v=self.v
        S_chi1=self.S_chi1
        S_chi2=self.S_chi2
        rfrak_chi1_x=self.rfrak_chi1_x
        rfrak_chi1_y=self.rfrak_chi1_y
        rfrak_chi2_x=self.rfrak_chi2_x
        rfrak_chi2_y=self.rfrak_chi2_y
        rfrak_frame_x=self.rfrak_frame_x
        rfrak_frame_y=self.rfrak_frame_y
        rfrak_frame_z=self.rfrak_frame_z
        M=self.M
        delta=self.delta
        nu=self.nu
        R=self.R
        nHat=self.nHat
        lambdaHat=self.lambdaHat
        ellHat=self.ellHat
        R_S1=self.R_S1
        R_S2=self.R_S2
        chiVec1=self.chiVec1
        chiVec2=self.chiVec2
        chi1_n=self.chi1_n
        chi1_lambda=self.chi1_lambda
        chi1_ell=self.chi1_ell
        chi2_n=self.chi2_n
        chi2_lambda=self.chi2_lambda
        chi2_ell=self.chi2_ell
        S_ell=self.S_ell
        S_n=self.S_n
        S_lambda=self.S_lambda
        Sigma_ell=self.Sigma_ell
        Sigma_n=self.Sigma_n
        Sigma_lambda=self.Sigma_lambda
        gamma_PN_6 = 0.0123456790123457*nu**3 + 6.36111111111111*nu**2 - 2.98177812235564*nu + 1.0
        gamma_PN_4 = 1.0 - 5.41666666666667*nu
        gamma_PN_5 = (S_ell*(0.888888888888889*nu + 3.33333333333333) + 2.0*Sigma_ell*delta)/M**2
        gamma_PN_2 = 1.0 - 0.333333333333333*nu
        a_ell_4 = S_n*(5.77777777777778*nu**2 + 14.75*nu + 1.5) + Sigma_n*delta*(2.83333333333333*nu**2 + 9.125*nu + 1.5)
        gamma_PN_0 = 1.00000000000000
        gamma_PN_7 = (S_ell*(-6.0*nu**2 - 10.5833333333333*nu + 5.0) - 2.66666666666667*Sigma_ell*delta*nu**2 + Sigma_ell*delta*(3.0 - 10.1666666666667*nu))/M**2
        gamma_PN_3 = (1.66666666666667*S_ell + Sigma_ell*delta)/M**2
        a_ell_0 = 7.0*S_n + 3.0*Sigma_n*delta
        a_ell_2 = S_n*(-9.66666666666667*nu - 10.0) + Sigma_n*delta*(-4.5*nu - 6.0)
        return ellHat*v**3/M + nHat*v**6*(a_ell_0 + v**2*(a_ell_2 + a_ell_4*v**2))*(gamma_PN_0 + v**2*(gamma_PN_2 + v*(gamma_PN_3 + v*(gamma_PN_4 + v*(gamma_PN_5 + v*(gamma_PN_6 + gamma_PN_7*v))))))/M**3

    def OrbitalAngularMomentum(self):
        xHat=self.xHat
        yHat=self.yHat
        zHat=self.zHat
        M1=self.M1
        M2=self.M2
        v=self.v
        S_chi1=self.S_chi1
        S_chi2=self.S_chi2
        rfrak_chi1_x=self.rfrak_chi1_x
        rfrak_chi1_y=self.rfrak_chi1_y
        rfrak_chi2_x=self.rfrak_chi2_x
        rfrak_chi2_y=self.rfrak_chi2_y
        rfrak_frame_x=self.rfrak_frame_x
        rfrak_frame_y=self.rfrak_frame_y
        rfrak_frame_z=self.rfrak_frame_z
        M=self.M
        delta=self.delta
        nu=self.nu
        R=self.R
        nHat=self.nHat
        lambdaHat=self.lambdaHat
        ellHat=self.ellHat
        R_S1=self.R_S1
        R_S2=self.R_S2
        chiVec1=self.chiVec1
        chiVec2=self.chiVec2
        chi1_n=self.chi1_n
        chi1_lambda=self.chi1_lambda
        chi1_ell=self.chi1_ell
        chi2_n=self.chi2_n
        chi2_lambda=self.chi2_lambda
        chi2_ell=self.chi2_ell
        S_ell=self.S_ell
        S_n=self.S_n
        S_lambda=self.S_lambda
        Sigma_ell=self.Sigma_ell
        Sigma_n=self.Sigma_n
        Sigma_lambda=self.Sigma_lambda
        L_6 = ellHat*(0.00540123456790123*nu**3 + 1.29166666666667*nu**2 - 30.9797035925835*nu + 8.4375)
        L_SO_7 = 0.0625*ellHat*(-29.0*S_ell*nu**2 + 1101.0*S_ell*nu - 405.0*S_ell - 15.0*Sigma_ell*delta*nu**2 + 468.0*Sigma_ell*delta*nu - 81.0*Sigma_ell*delta)/M**2 + 0.0416666666666667*lambdaHat*(-32.0*S_lambda*nu**2 + 2.0*S_lambda*nu - 174.0*S_lambda - 16.0*Sigma_lambda*delta*nu**2 - 79.0*Sigma_lambda*delta*nu - 12.0*Sigma_lambda*delta)/M**2 + 0.0208333333333333*nHat*(11.0*S_n*nu**2 - 1331.0*S_n*nu + 183.0*S_n + 5.0*Sigma_n*delta*nu**2 - 734.0*Sigma_n*delta*nu + 183.0*Sigma_n*delta)/M**2
        L_SO_5 = 0.0138888888888889*ellHat*(427.0*S_ell*nu - 693.0*S_ell + 210.0*Sigma_ell*delta*nu - 189.0*Sigma_ell*delta)/M**2 + 0.166666666666667*lambdaHat*(18.0*S_lambda*nu - 21.0*S_lambda + 8.0*Sigma_lambda*delta*nu - 3.0*Sigma_lambda*delta)/M**2 + 0.0416666666666667*nHat*(-19.0*S_n*nu + 33.0*S_n - 10.0*Sigma_n*delta*nu + 33.0*Sigma_n*delta)/M**2
        L_SO_3 = 0.166666666666667*ellHat*(-35.0*S_ell - 15.0*Sigma_ell*delta)/M**2 + lambdaHat*(-3.0*S_lambda - Sigma_lambda*delta)/M**2 + 0.5*nHat*(S_n + Sigma_n*delta)/M**2
        L_2 = ellHat*(0.166666666666667*nu + 1.5)
        L_4 = ellHat*(0.0416666666666667*nu**2 - 2.375*nu + 3.375)
        L_0 = ellHat
        L_coeff = M**2*nu/v
        return L_coeff*(L_0 + v**2*(L_2 + v*(L_SO_3 + v*(L_4 + v*(L_SO_5 + v*(L_6 + L_SO_7*v))))))

               

    #Evolve PN
    def __TaylorT1__(self,t, y):
        self.__Recalculate__(t, y)
        xHat=self.xHat
        yHat=self.yHat
        zHat=self.zHat
        M1=self.M1
        M2=self.M2
        v=self.v
        S_chi1=self.S_chi1
        S_chi2=self.S_chi2
        rfrak_chi1_x=self.rfrak_chi1_x
        rfrak_chi1_y=self.rfrak_chi1_y
        rfrak_chi2_x=self.rfrak_chi2_x
        rfrak_chi2_y=self.rfrak_chi2_y
        rfrak_frame_x=self.rfrak_frame_x
        rfrak_frame_y=self.rfrak_frame_y
        rfrak_frame_z=self.rfrak_frame_z
        M=self.M
        delta=self.delta
        nu=self.nu
        R=self.R
        nHat=self.nHat
        lambdaHat=self.lambdaHat
        ellHat=self.ellHat
        R_S1=self.R_S1
        R_S2=self.R_S2
        chiVec1=self.chiVec1
        chiVec2=self.chiVec2
        chi1chi1=self.chi1chi1
        chi1chi2=self.chi1chi2
        chi2chi2=self.chi2chi2
        chi1_n=self.chi1_n
        chi1_lambda=self.chi1_lambda
        chi1_ell=self.chi1_ell
        chi2_n=self.chi2_n
        chi2_lambda=self.chi2_lambda
        chi2_ell=self.chi2_ell
        S_ell=self.S_ell
        S_n=self.S_n
        S_lambda=self.S_lambda
        Sigma_ell=self.Sigma_ell
        Sigma_n=self.Sigma_n
        Sigma_lambda=self.Sigma_lambda
        chi_s_ell=self.chi_s_ell
        chi_a_ell=self.chi_a_ell
        logv=self.logv
        Fcal_coeff=self.Fcal_coeff
        Fcal_0=self.Fcal_0
        Fcal_2=self.Fcal_2
        Fcal_3=self.Fcal_3
        Fcal_4=self.Fcal_4
        Fcal_5=self.Fcal_5
        Fcal_6=self.Fcal_6
        Fcal_lnv_6=self.Fcal_lnv_6
        Fcal_7=self.Fcal_7
        Fcal_SQ_4=self.Fcal_SQ_4
        Fcal_SO_3=self.Fcal_SO_3
        Fcal_SO_5=self.Fcal_SO_5
        Fcal_SO_6=self.Fcal_SO_6
        Fcal_SO_7=self.Fcal_SO_7
        E_0=self.E_0
        E_2=self.E_2
        E_4=self.E_4
        E_6=self.E_6
        E_SQ_4=self.E_SQ_4
        E_SO_3=self.E_SO_3
        E_SO_5=self.E_SO_5
        E_SO_7=self.E_SO_7
        v, rfrak_chi1_x, rfrak_chi1_y, rfrak_chi2_x, rfrak_chi2_y=self.v, self.rfrak_chi1_x, self.rfrak_chi1_y, self.rfrak_chi2_x, self.rfrak_chi2_y
        rfrak_frame_x, rfrak_frame_y, rfrak_frame_z, Phi=self.rfrak_frame_x, self.rfrak_frame_y, self.rfrak_frame_z, self.Phi
        dvdt_T1=0
        if v>=1.0 and self.NotForward:
            print("Beyond domain of PN validity, this is a good way to terminate.")
            self.terminal1=False
        Flux = Fcal_coeff*(Fcal_0 + v**2*(Fcal_2 + v*(Fcal_3 + Fcal_SO_3 + v*(Fcal_4 + Fcal_SQ_4 + v*(Fcal_5 + Fcal_SO_5 + v*(Fcal_6 + Fcal_SO_6 + Fcal_lnv_6*logv + v*(Fcal_7 + Fcal_SO_7)))))))
        dEdv = -0.5*M*nu*v*(2.0*E_0 + v**2*(4.0*E_2 + v*(5.0*E_SO_3 + v*(6.0*E_4 + 6.0*E_SQ_4 + v*(7.0*E_SO_5 + v*(8.0*E_6 + 9.0*E_SO_7*v))))))
        Absorption = 0
        dvdt_T1 = (-Absorption - Flux)/dEdv
        if dvdt_T1<1.0e-12 and self.NotForward:
            print("v is decreasing, which is not an uncommon way to stop.")
            self.terminal2=False
        return self.CommonRHS(dvdt_T1, y)

    def CommonRHS(self,dvdt, y):
        dydt=np.zeros(9)
        rfrak_frame=np.zeros(3)
        rfrak_frame[0] = y[5]
        rfrak_frame[1] = y[6]
        rfrak_frame[2] = y[7]
        rfrakdot_frame = FrameFromAngularVelocityIntegrand(rfrak_frame, self.OmegaVec().vector)
        dydt[0] = dvdt
        if(self.EvolveSpin1):
            dydt[1], dydt[2]= FrameFromAngularVelocity_2D_Integrand(y[1], y[2],\
                (self.S_chi1.inverse*self.OmegaVec_chiVec_1()*self.S_chi1).vector)
        else:
            dydt[1] = 0.0
            dydt[2] = 0.0
        if(self.EvolveSpin2):
            dydt[3], dydt[4] = FrameFromAngularVelocity_2D_Integrand(y[3], y[4],\
                (self.S_chi2.inverse*self.OmegaVec_chiVec_2()*self.S_chi2).vector)
        else:
            dydt[3] = 0.0
            dydt[4] = 0.0
        dydt[5] = rfrakdot_frame[0]
        dydt[6] = rfrakdot_frame[1]
        dydt[7] = rfrakdot_frame[2]
        dydt[8] = self.v*self.v*self.v/self.M

        return dydt

    # Prepare for integration
    def Evolution(self): 
        v_i,rfrak_chi1_x_i,rfrak_chi1_y_i,rfrak_chi2_x_i,rfrak_chi2_y_i,rfrak_frame_x_i,rfrak_frame_y_i,rfrak_frame_z_i\
            =self.v_i,self.rfrak_chi1_x_i,self.rfrak_chi1_y_i,self.rfrak_chi2_x_i,self.rfrak_chi2_y_i,self.rfrak_frame_x_i,self.rfrak_frame_y_i,self.rfrak_frame_z_i
        def terminate(t,y):
            return 1.0*self.terminal1*self.terminal2
        terminate.terminal=True
        TMerger=5.0/(256.0*self.nu*v_i**8)
        TEnd=TMerger
        if self.t_PN2:
            TEnd=self.t_PN2
        time=[0.0]
        while time[-1]<TEnd and 2*self.M*(256*self.nu*(TMerger-time[-1])/5)**(3/8)/self.StepsPerOrbit>self.MinStep:
            time.append(time[-1]+2*self.M*(256*self.nu*(TMerger-time[-1])/5)**(3/8)/self.StepsPerOrbit)
        time=np.delete(time, -1)
        
        # Integrate
        yy=solve_ivp(self.__TaylorT1__, [time[0],time[-1]], [v_i,rfrak_chi1_x_i,\
            rfrak_chi1_y_i,rfrak_chi2_x_i,rfrak_chi2_y_i,rfrak_frame_x_i,\
            rfrak_frame_y_i,rfrak_frame_z_i,0.0], method='DOP853',\
            t_eval=time, dense_output=True, events=terminate, rtol=self.tol, atol=self.tol)           
        if self.ForwardInTime:
            self.NotForward=False
            time=[0.0]
            TStart=-3*TMerger
            if self.t_PN1:
                TStart=self.t_PN1
            while time[-1]>TStart:
                time.append(time[-1]-2*self.M*(256*self.nu*(TMerger-time[-1])/5)**(3/8)/self.StepsPerOrbit)
            yyForward=solve_ivp(self.__TaylorT1__, [time[0],time[-1]], [v_i,rfrak_chi1_x_i,\
                rfrak_chi1_y_i,rfrak_chi2_x_i,rfrak_chi2_y_i,rfrak_frame_x_i,\
                rfrak_frame_y_i,rfrak_frame_z_i,0.0], method='DOP853',\
                t_eval=time, dense_output=True, rtol=self.tol, atol=self.tol)
            yy.t=np.append(yyForward.t[1:][::-1],yy.t)
            data=np.empty((9,len(yy.t)))
            for i in range(9):
                data[i]=np.append(yyForward.y[i][1:][::-1],yy.y[i])
            yy.y=data
            
        return yy

spec=[("xHat", float64),("yHat", float64),("zHat", float64),("M1", float64),("M2", float64),("v", float64),("S_chi1", float64),("S_chi2", float64),("rfrak_chi1_x", float64),("rfrak_chi1_y", float64),("rfrak_chi2_x", float64),("rfrak_chi2_y", float64),("rfrak_frame_x", float64),("rfrak_frame_y", float64),("rfrak_frame_z", float64),("M", float64),("delta", float64),("nu", float64),("R", float64),("nHat", float64),("lambdaHat", float64),("ellHat", float64),("R_S1", float64),("R_S2", float64),("chiVec1", float64),("chiVec2", float64),("chi1chi1", float64),("chi1chi2", float64),("chi2chi2", float64),("chi1_n", float64),("chi1_lambda", float64),("chi1_ell", float64),("chi2_n", float64),("chi2_lambda", float64),("chi2_ell", float64),("S_ell", float64),("S_n", float64),("S_lambda", float64),("Sigma_ell", float64),("Sigma_n", float64),("Sigma_lambda", float64),("chi_s_ell", float64),("chi_a_ell", float64),("logv", float64),("Fcal_coeff", float64),("Fcal_0", float64),("Fcal_2", float64),("Fcal_3", float64),("Fcal_4", float64),("Fcal_5", float64),("Fcal_6", float64),("Fcal_lnv_6", float64),("Fcal_7", float64),("Fcal_8", float64),("Fcal_lnv_8", float64),("Fcal_SQ_4", float64),("Fcal_SO_3", float64),("Fcal_SO_5", float64),("Fcal_SO_6", float64),("Fcal_SO_7", float64),("Fcal_SO_8", float64),("E_0", float64),("E_2", float64),("E_4", float64),("E_6", float64),("E_8", float64),("E_lnv_8", float64),("E_SQ_4", float64),("E_SO_3", float64),("E_SO_5", float64),("E_SO_7", float64),("xHat_i", float64),("yHat_i", float64),("zHat_i", float64),("M1_i", float64),("M2_i", float64),("v_i", float64),("S_chi1_i", float64),("S_chi2_i", float64),("rfrak_chi1_x_i", float64),("rfrak_chi1_y_i", float64),("rfrak_chi2_x_i", float64),("rfrak_chi2_y_i", float64),("rfrak_frame_x_i", float64),("rfrak_frame_y_i", float64),("rfrak_frame_z_i", float64),("t_PN1",float64),("t_PN2",float64),("StepsPerOrbit",int32),("ForwardInTime",boolean),("tol",float64),("MinStep",float64),("EvolveSpin1",boolean),("EvolveSpin2",boolean),("terminal1",boolean),("terminal2",boolean),("NotForward",boolean)]
@jitclass(spec)
class TaylorTn_4p0PN_Q : 
    def __init__(self, xHat_i, yHat_i, zHat_i, M1_i, M2_i, v_i, S_chi1_i, S_chi2_i, rfrak_chi1_x_i, rfrak_chi1_y_i, rfrak_chi2_x_i, rfrak_chi2_y_i, rfrak_frame_x_i, rfrak_frame_y_i, rfrak_frame_z_i, t_PN1=False, t_PN2=False, StepsPerOrbit=32, ForwardInTime=True, tol=1e-12, MinStep=1e-7) :
        self.t_PN1,self.t_PN2,self.StepsPerOrbit,self.ForwardInTime,self.tol,self.MinStep=t_PN1, t_PN2, StepsPerOrbit, ForwardInTime, tol, MinStep
        terminal1=True
        terminal2=True
        NotForward=True
        self.xHat_i=xHat_i
        self.yHat_i=yHat_i
        self.zHat_i=zHat_i
        self.M1_i=M1_i
        self.M2_i=M2_i
        self.v_i=v_i
        self.S_chi1_i=S_chi1_i
        self.S_chi2_i=S_chi2_i
        self.rfrak_chi1_x_i=rfrak_chi1_x_i
        self.rfrak_chi1_y_i=rfrak_chi1_y_i
        self.rfrak_chi2_x_i=rfrak_chi2_x_i
        self.rfrak_chi2_y_i=rfrak_chi2_y_i
        self.rfrak_frame_x_i=rfrak_frame_x_i
        self.rfrak_frame_y_i=rfrak_frame_y_i
        self.rfrak_frame_z_i=rfrak_frame_z_i
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
        chi1chi1=dot(chiVec1.vector,chiVec1.vector)
        chi1chi2=dot(chiVec1.vector,chiVec2.vector)
        chi2chi2=dot(chiVec2.vector,chiVec2.vector)
        chi1_n=dot(chiVec1.vector,nHat.vector)
        chi1_lambda=dot(chiVec1.vector,lambdaHat.vector)
        chi1_ell=dot(chiVec1.vector,ellHat.vector)
        chi2_n=dot(chiVec2.vector,nHat.vector)
        chi2_lambda=dot(chiVec2.vector,lambdaHat.vector)
        chi2_ell=dot(chiVec2.vector,ellHat.vector)
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
        self.xHat=xHat
        self.yHat=yHat
        self.zHat=zHat
        self.M1=M1
        self.M2=M2
        self.v=v
        self.S_chi1=S_chi1
        self.S_chi2=S_chi2
        self.rfrak_chi1_x=rfrak_chi1_x
        self.rfrak_chi1_y=rfrak_chi1_y
        self.rfrak_chi2_x=rfrak_chi2_x
        self.rfrak_chi2_y=rfrak_chi2_y
        self.rfrak_frame_x=rfrak_frame_x
        self.rfrak_frame_y=rfrak_frame_y
        self.rfrak_frame_z=rfrak_frame_z
        self.M=M
        self.delta=delta
        self.nu=nu
        self.R=R
        self.nHat=nHat
        self.lambdaHat=lambdaHat
        self.ellHat=ellHat
        self.R_S1=R_S1
        self.R_S2=R_S2
        self.chiVec1=chiVec1
        self.chiVec2=chiVec2
        self.chi1chi1=chi1chi1
        self.chi1chi2=chi1chi2
        self.chi2chi2=chi2chi2
        self.chi1_n=chi1_n
        self.chi1_lambda=chi1_lambda
        self.chi1_ell=chi1_ell
        self.chi2_n=chi2_n
        self.chi2_lambda=chi2_lambda
        self.chi2_ell=chi2_ell
        self.S_ell=S_ell
        self.S_n=S_n
        self.S_lambda=S_lambda
        self.Sigma_ell=Sigma_ell
        self.Sigma_n=Sigma_n
        self.Sigma_lambda=Sigma_lambda
        self.chi_s_ell=chi_s_ell
        self.chi_a_ell=chi_a_ell
        self.logv=logv
        self.Fcal_coeff=Fcal_coeff
        self.Fcal_0=Fcal_0
        self.Fcal_2=Fcal_2
        self.Fcal_3=Fcal_3
        self.Fcal_4=Fcal_4
        self.Fcal_5=Fcal_5
        self.Fcal_6=Fcal_6
        self.Fcal_lnv_6=Fcal_lnv_6
        self.Fcal_7=Fcal_7
        self.Fcal_8=Fcal_8
        self.Fcal_lnv_8=Fcal_lnv_8
        self.Fcal_SQ_4=Fcal_SQ_4
        self.Fcal_SO_3=Fcal_SO_3
        self.Fcal_SO_5=Fcal_SO_5
        self.Fcal_SO_6=Fcal_SO_6
        self.Fcal_SO_7=Fcal_SO_7
        self.Fcal_SO_8=Fcal_SO_8
        self.E_0=E_0
        self.E_2=E_2
        self.E_4=E_4
        self.E_6=E_6
        self.E_8=E_8
        self.E_lnv_8=E_lnv_8
        self.E_SQ_4=E_SQ_4
        self.E_SO_3=E_SO_3
        self.E_SO_5=E_SO_5
        self.E_SO_7=E_SO_7
        Phi=0.0
        EvolveSpin1=(S_chi1*conjugate(S_chi1)).norm>1e-12
        EvolveSpin2=(S_chi2*conjugate(S_chi2)).norm>1e-12
        self.v, self.rfrak_chi1_x, self.rfrak_chi1_y, self.rfrak_chi2_x, self.rfrak_chi2_y=v, rfrak_chi1_x, rfrak_chi1_y, rfrak_chi2_x, rfrak_chi2_y
        self.rfrak_frame_x, self.rfrak_frame_y, self.rfrak_frame_z, self.Phi=rfrak_frame_x, rfrak_frame_y, rfrak_frame_z, Phi
        self.EvolveSpin1=EvolveSpin1
        self.EvolveSpin2=EvolveSpin2
        self.terminal1,self.terminal2,self.NotForward=terminal1,terminal2,NotForward

    def __Recalculate__(self,t, y):
        v = y[0]
        rfrak_chi1_x = y[1]
        rfrak_chi1_y = y[2]
        rfrak_chi2_x = y[3]
        rfrak_chi2_y = y[4]
        rfrak_frame_x = y[5]
        rfrak_frame_y = y[6]
        rfrak_frame_z = y[7]
        Phi = y[8]
        self.v, self.rfrak_chi1_x, self.rfrak_chi1_y, self.rfrak_chi2_x, self.rfrak_chi2_y=v, rfrak_chi1_x, rfrak_chi1_y, rfrak_chi2_x, rfrak_chi2_y
        self.rfrak_frame_x, self.rfrak_frame_y, self.rfrak_frame_z, self.Phi=rfrak_frame_x, rfrak_frame_y, rfrak_frame_z, Phi
        xHat=self.xHat
        yHat=self.yHat
        zHat=self.zHat
        M1=self.M1
        M2=self.M2
        v=self.v
        S_chi1=self.S_chi1
        S_chi2=self.S_chi2
        rfrak_chi1_x=self.rfrak_chi1_x
        rfrak_chi1_y=self.rfrak_chi1_y
        rfrak_chi2_x=self.rfrak_chi2_x
        rfrak_chi2_y=self.rfrak_chi2_y
        rfrak_frame_x=self.rfrak_frame_x
        rfrak_frame_y=self.rfrak_frame_y
        rfrak_frame_z=self.rfrak_frame_z
        M=self.M
        delta=self.delta
        nu=self.nu
        R=self.R
        nHat=self.nHat
        lambdaHat=self.lambdaHat
        ellHat=self.ellHat
        R_S1=self.R_S1
        R_S2=self.R_S2
        chiVec1=self.chiVec1
        chiVec2=self.chiVec2
        chi1chi1=self.chi1chi1
        chi1chi2=self.chi1chi2
        chi2chi2=self.chi2chi2
        chi1_n=self.chi1_n
        chi1_lambda=self.chi1_lambda
        chi1_ell=self.chi1_ell
        chi2_n=self.chi2_n
        chi2_lambda=self.chi2_lambda
        chi2_ell=self.chi2_ell
        S_ell=self.S_ell
        S_n=self.S_n
        S_lambda=self.S_lambda
        Sigma_ell=self.Sigma_ell
        Sigma_n=self.Sigma_n
        Sigma_lambda=self.Sigma_lambda
        chi_s_ell=self.chi_s_ell
        chi_a_ell=self.chi_a_ell
        logv=self.logv
        Fcal_coeff=self.Fcal_coeff
        Fcal_0=self.Fcal_0
        Fcal_2=self.Fcal_2
        Fcal_3=self.Fcal_3
        Fcal_4=self.Fcal_4
        Fcal_5=self.Fcal_5
        Fcal_6=self.Fcal_6
        Fcal_lnv_6=self.Fcal_lnv_6
        Fcal_7=self.Fcal_7
        Fcal_8=self.Fcal_8
        Fcal_lnv_8=self.Fcal_lnv_8
        Fcal_SQ_4=self.Fcal_SQ_4
        Fcal_SO_3=self.Fcal_SO_3
        Fcal_SO_5=self.Fcal_SO_5
        Fcal_SO_6=self.Fcal_SO_6
        Fcal_SO_7=self.Fcal_SO_7
        Fcal_SO_8=self.Fcal_SO_8
        E_0=self.E_0
        E_2=self.E_2
        E_4=self.E_4
        E_6=self.E_6
        E_8=self.E_8
        E_lnv_8=self.E_lnv_8
        E_SQ_4=self.E_SQ_4
        E_SO_3=self.E_SO_3
        E_SO_5=self.E_SO_5
        E_SO_7=self.E_SO_7
        R = exp(rfrak_frame_x*xHat + rfrak_frame_y*yHat + rfrak_frame_z*zHat)
        nHat = R*xHat*conjugate(R)
        lambdaHat = R*yHat*conjugate(R)
        ellHat = R*zHat*conjugate(R)
        R_S1 = exp(rfrak_chi1_x*xHat + rfrak_chi1_y*yHat)
        R_S2 = exp(rfrak_chi2_x*xHat + rfrak_chi2_y*yHat)
        chiVec1 = S_chi1*R_S1*zHat*conjugate(R_S1)*conjugate(S_chi1)
        chiVec2 = S_chi2*R_S2*zHat*conjugate(R_S2)*conjugate(S_chi2)
        chi1_n = dot(chiVec1.vector,nHat.vector)
        chi1_lambda = dot(chiVec1.vector,lambdaHat.vector)
        chi1_ell = dot(chiVec1.vector,ellHat.vector)
        chi2_n = dot(chiVec2.vector,nHat.vector)
        chi2_lambda = dot(chiVec2.vector,lambdaHat.vector)
        chi2_ell = dot(chiVec2.vector,ellHat.vector)
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
        self.R=R
        self.nHat=nHat
        self.lambdaHat=lambdaHat
        self.ellHat=ellHat
        self.R_S1=R_S1
        self.R_S2=R_S2
        self.chiVec1=chiVec1
        self.chiVec2=chiVec2
        self.chi1_n=chi1_n
        self.chi1_lambda=chi1_lambda
        self.chi1_ell=chi1_ell
        self.chi2_n=chi2_n
        self.chi2_lambda=chi2_lambda
        self.chi2_ell=chi2_ell
        self.S_ell=S_ell
        self.S_n=S_n
        self.S_lambda=S_lambda
        self.Sigma_ell=Sigma_ell
        self.Sigma_n=Sigma_n
        self.Sigma_lambda=Sigma_lambda
        self.chi_s_ell=chi_s_ell
        self.chi_a_ell=chi_a_ell
        self.logv=logv
        self.Fcal_coeff=Fcal_coeff
        self.Fcal_SQ_4=Fcal_SQ_4
        self.Fcal_SO_3=Fcal_SO_3
        self.Fcal_SO_5=Fcal_SO_5
        self.Fcal_SO_6=Fcal_SO_6
        self.Fcal_SO_7=Fcal_SO_7
        self.Fcal_SO_8=Fcal_SO_8
        self.E_SQ_4=E_SQ_4
        self.E_SO_3=E_SO_3
        self.E_SO_5=E_SO_5
        self.E_SO_7=E_SO_7

    def OmegaVec_chiVec_1(self):
        xHat=self.xHat
        yHat=self.yHat
        zHat=self.zHat
        M1=self.M1
        M2=self.M2
        v=self.v
        S_chi1=self.S_chi1
        S_chi2=self.S_chi2
        rfrak_chi1_x=self.rfrak_chi1_x
        rfrak_chi1_y=self.rfrak_chi1_y
        rfrak_chi2_x=self.rfrak_chi2_x
        rfrak_chi2_y=self.rfrak_chi2_y
        rfrak_frame_x=self.rfrak_frame_x
        rfrak_frame_y=self.rfrak_frame_y
        rfrak_frame_z=self.rfrak_frame_z
        M=self.M
        delta=self.delta
        nu=self.nu
        R=self.R
        nHat=self.nHat
        lambdaHat=self.lambdaHat
        ellHat=self.ellHat
        R_S1=self.R_S1
        R_S2=self.R_S2
        chiVec1=self.chiVec1
        chiVec2=self.chiVec2
        chi1_n=self.chi1_n
        chi1_lambda=self.chi1_lambda
        chi1_ell=self.chi1_ell
        chi2_n=self.chi2_n
        chi2_lambda=self.chi2_lambda
        chi2_ell=self.chi2_ell
        S_ell=self.S_ell
        S_n=self.S_n
        S_lambda=self.S_lambda
        Sigma_ell=self.Sigma_ell
        Sigma_n=self.Sigma_n
        Sigma_lambda=self.Sigma_lambda
        Omega1_coeff = v**5/M
        return Omega1_coeff*(ellHat*(-0.75*delta + 0.5*nu + v**2*(delta*(0.625*nu - 0.5625) + nu*(1.25 - 0.0416666666666667*nu) + v**2*(delta*(nu*(4.875 - 0.15625*nu) - 0.84375) + nu*(nu*(-0.0208333333333333*nu - 3.28125) + 0.1875) + 0.84375) + 0.5625) + 0.75) + nHat*v*(3.0*chi1_n*nu + 3.0*M2**2*chi2_n/M**2) - M2**2*chiVec2*v/M**2)

    def OmegaVec_chiVec_2(self):
        xHat=self.xHat
        yHat=self.yHat
        zHat=self.zHat
        M1=self.M1
        M2=self.M2
        v=self.v
        S_chi1=self.S_chi1
        S_chi2=self.S_chi2
        rfrak_chi1_x=self.rfrak_chi1_x
        rfrak_chi1_y=self.rfrak_chi1_y
        rfrak_chi2_x=self.rfrak_chi2_x
        rfrak_chi2_y=self.rfrak_chi2_y
        rfrak_frame_x=self.rfrak_frame_x
        rfrak_frame_y=self.rfrak_frame_y
        rfrak_frame_z=self.rfrak_frame_z
        M=self.M
        delta=self.delta
        nu=self.nu
        R=self.R
        nHat=self.nHat
        lambdaHat=self.lambdaHat
        ellHat=self.ellHat
        R_S1=self.R_S1
        R_S2=self.R_S2
        chiVec1=self.chiVec1
        chiVec2=self.chiVec2
        chi1_n=self.chi1_n
        chi1_lambda=self.chi1_lambda
        chi1_ell=self.chi1_ell
        chi2_n=self.chi2_n
        chi2_lambda=self.chi2_lambda
        chi2_ell=self.chi2_ell
        S_ell=self.S_ell
        S_n=self.S_n
        S_lambda=self.S_lambda
        Sigma_ell=self.Sigma_ell
        Sigma_n=self.Sigma_n
        Sigma_lambda=self.Sigma_lambda
        Omega2_coeff = v**5/M
        return Omega2_coeff*(ellHat*(0.75*delta + 0.5*nu + v**2*(delta*(0.5625 - 0.625*nu) + nu*(1.25 - 0.0416666666666667*nu) + v**2*(delta*(nu*(0.15625*nu - 4.875) + 0.84375) + nu*(nu*(-0.0208333333333333*nu - 3.28125) + 0.1875) + 0.84375) + 0.5625) + 0.75) + nHat*v*(3.0*chi2_n*nu + 3.0*M1**2*chi1_n/M**2) - M1**2*chiVec1*v/M**2)

    def OmegaVec(self):
        xHat=self.xHat
        yHat=self.yHat
        zHat=self.zHat
        M1=self.M1
        M2=self.M2
        v=self.v
        S_chi1=self.S_chi1
        S_chi2=self.S_chi2
        rfrak_chi1_x=self.rfrak_chi1_x
        rfrak_chi1_y=self.rfrak_chi1_y
        rfrak_chi2_x=self.rfrak_chi2_x
        rfrak_chi2_y=self.rfrak_chi2_y
        rfrak_frame_x=self.rfrak_frame_x
        rfrak_frame_y=self.rfrak_frame_y
        rfrak_frame_z=self.rfrak_frame_z
        M=self.M
        delta=self.delta
        nu=self.nu
        R=self.R
        nHat=self.nHat
        lambdaHat=self.lambdaHat
        ellHat=self.ellHat
        R_S1=self.R_S1
        R_S2=self.R_S2
        chiVec1=self.chiVec1
        chiVec2=self.chiVec2
        chi1_n=self.chi1_n
        chi1_lambda=self.chi1_lambda
        chi1_ell=self.chi1_ell
        chi2_n=self.chi2_n
        chi2_lambda=self.chi2_lambda
        chi2_ell=self.chi2_ell
        S_ell=self.S_ell
        S_n=self.S_n
        S_lambda=self.S_lambda
        Sigma_ell=self.Sigma_ell
        Sigma_n=self.Sigma_n
        Sigma_lambda=self.Sigma_lambda
        gamma_PN_6 = 0.0123456790123457*nu**3 + 6.36111111111111*nu**2 - 2.98177812235564*nu + 1.0
        gamma_PN_4 = 1.0 - 5.41666666666667*nu
        gamma_PN_5 = (S_ell*(0.888888888888889*nu + 3.33333333333333) + 2.0*Sigma_ell*delta)/M**2
        gamma_PN_2 = 1.0 - 0.333333333333333*nu
        a_ell_4 = S_n*(5.77777777777778*nu**2 + 14.75*nu + 1.5) + Sigma_n*delta*(2.83333333333333*nu**2 + 9.125*nu + 1.5)
        gamma_PN_0 = 1.00000000000000
        gamma_PN_7 = (S_ell*(-6.0*nu**2 - 10.5833333333333*nu + 5.0) - 2.66666666666667*Sigma_ell*delta*nu**2 + Sigma_ell*delta*(3.0 - 10.1666666666667*nu))/M**2
        gamma_PN_3 = (1.66666666666667*S_ell + Sigma_ell*delta)/M**2
        a_ell_0 = 7.0*S_n + 3.0*Sigma_n*delta
        a_ell_2 = S_n*(-9.66666666666667*nu - 10.0) + Sigma_n*delta*(-4.5*nu - 6.0)
        return ellHat*v**3/M + nHat*v**6*(a_ell_0 + v**2*(a_ell_2 + a_ell_4*v**2))*(gamma_PN_0 + v**2*(gamma_PN_2 + v*(gamma_PN_3 + v*(gamma_PN_4 + v*(gamma_PN_5 + v*(gamma_PN_6 + gamma_PN_7*v))))))/M**3

    def OrbitalAngularMomentum(self):
        xHat=self.xHat
        yHat=self.yHat
        zHat=self.zHat
        M1=self.M1
        M2=self.M2
        v=self.v
        S_chi1=self.S_chi1
        S_chi2=self.S_chi2
        rfrak_chi1_x=self.rfrak_chi1_x
        rfrak_chi1_y=self.rfrak_chi1_y
        rfrak_chi2_x=self.rfrak_chi2_x
        rfrak_chi2_y=self.rfrak_chi2_y
        rfrak_frame_x=self.rfrak_frame_x
        rfrak_frame_y=self.rfrak_frame_y
        rfrak_frame_z=self.rfrak_frame_z
        M=self.M
        delta=self.delta
        nu=self.nu
        R=self.R
        nHat=self.nHat
        lambdaHat=self.lambdaHat
        ellHat=self.ellHat
        R_S1=self.R_S1
        R_S2=self.R_S2
        chiVec1=self.chiVec1
        chiVec2=self.chiVec2
        chi1_n=self.chi1_n
        chi1_lambda=self.chi1_lambda
        chi1_ell=self.chi1_ell
        chi2_n=self.chi2_n
        chi2_lambda=self.chi2_lambda
        chi2_ell=self.chi2_ell
        S_ell=self.S_ell
        S_n=self.S_n
        S_lambda=self.S_lambda
        Sigma_ell=self.Sigma_ell
        Sigma_n=self.Sigma_n
        Sigma_lambda=self.Sigma_lambda
        L_lnv_8 = -42.6666666666667*ellHat*nu
        L_6 = ellHat*(0.00540123456790123*nu**3 + 1.29166666666667*nu**2 - 30.9797035925835*nu + 8.4375)
        L_SO_7 = 0.0625*ellHat*(-29.0*S_ell*nu**2 + 1101.0*S_ell*nu - 405.0*S_ell - 15.0*Sigma_ell*delta*nu**2 + 468.0*Sigma_ell*delta*nu - 81.0*Sigma_ell*delta)/M**2 + 0.0416666666666667*lambdaHat*(-32.0*S_lambda*nu**2 + 2.0*S_lambda*nu - 174.0*S_lambda - 16.0*Sigma_lambda*delta*nu**2 - 79.0*Sigma_lambda*delta*nu - 12.0*Sigma_lambda*delta)/M**2 + 0.0208333333333333*nHat*(11.0*S_n*nu**2 - 1331.0*S_n*nu + 183.0*S_n + 5.0*Sigma_n*delta*nu**2 - 734.0*Sigma_n*delta*nu + 183.0*Sigma_n*delta)/M**2
        L_SO_5 = 0.0138888888888889*ellHat*(427.0*S_ell*nu - 693.0*S_ell + 210.0*Sigma_ell*delta*nu - 189.0*Sigma_ell*delta)/M**2 + 0.166666666666667*lambdaHat*(18.0*S_lambda*nu - 21.0*S_lambda + 8.0*Sigma_lambda*delta*nu - 3.0*Sigma_lambda*delta)/M**2 + 0.0416666666666667*nHat*(-19.0*S_n*nu + 33.0*S_n - 10.0*Sigma_n*delta*nu + 33.0*Sigma_n*delta)/M**2
        L_SO_3 = 0.166666666666667*ellHat*(-35.0*S_ell - 15.0*Sigma_ell*delta)/M**2 + lambdaHat*(-3.0*S_lambda - Sigma_lambda*delta)/M**2 + 0.5*nHat*(S_n + Sigma_n*delta)/M**2
        L_2 = ellHat*(0.166666666666667*nu + 1.5)
        L_4 = ellHat*(0.0416666666666667*nu**2 - 2.375*nu + 3.375)
        L_0 = ellHat
        L_8 = ellHat*(-0.714285714285714*nu*(0.0024755658436214*nu**3 + 0.174189814814815*nu**2 - 90.1327990262052*nu + 153.88379682994) + 1.82857142857143*nu + 22.1484375)
        L_coeff = M**2*nu/v
        return L_coeff*(L_0 + v**2*(L_2 + v*(L_SO_3 + v*(L_4 + v*(L_SO_5 + v*(L_6 + v*(L_SO_7 + v*(L_8 + L_lnv_8*math.log(v)))))))))

               

    #Evolve PN
    def __TaylorT1__(self,t, y):
        self.__Recalculate__(t, y)
        xHat=self.xHat
        yHat=self.yHat
        zHat=self.zHat
        M1=self.M1
        M2=self.M2
        v=self.v
        S_chi1=self.S_chi1
        S_chi2=self.S_chi2
        rfrak_chi1_x=self.rfrak_chi1_x
        rfrak_chi1_y=self.rfrak_chi1_y
        rfrak_chi2_x=self.rfrak_chi2_x
        rfrak_chi2_y=self.rfrak_chi2_y
        rfrak_frame_x=self.rfrak_frame_x
        rfrak_frame_y=self.rfrak_frame_y
        rfrak_frame_z=self.rfrak_frame_z
        M=self.M
        delta=self.delta
        nu=self.nu
        R=self.R
        nHat=self.nHat
        lambdaHat=self.lambdaHat
        ellHat=self.ellHat
        R_S1=self.R_S1
        R_S2=self.R_S2
        chiVec1=self.chiVec1
        chiVec2=self.chiVec2
        chi1chi1=self.chi1chi1
        chi1chi2=self.chi1chi2
        chi2chi2=self.chi2chi2
        chi1_n=self.chi1_n
        chi1_lambda=self.chi1_lambda
        chi1_ell=self.chi1_ell
        chi2_n=self.chi2_n
        chi2_lambda=self.chi2_lambda
        chi2_ell=self.chi2_ell
        S_ell=self.S_ell
        S_n=self.S_n
        S_lambda=self.S_lambda
        Sigma_ell=self.Sigma_ell
        Sigma_n=self.Sigma_n
        Sigma_lambda=self.Sigma_lambda
        chi_s_ell=self.chi_s_ell
        chi_a_ell=self.chi_a_ell
        logv=self.logv
        Fcal_coeff=self.Fcal_coeff
        Fcal_0=self.Fcal_0
        Fcal_2=self.Fcal_2
        Fcal_3=self.Fcal_3
        Fcal_4=self.Fcal_4
        Fcal_5=self.Fcal_5
        Fcal_6=self.Fcal_6
        Fcal_lnv_6=self.Fcal_lnv_6
        Fcal_7=self.Fcal_7
        Fcal_8=self.Fcal_8
        Fcal_lnv_8=self.Fcal_lnv_8
        Fcal_SQ_4=self.Fcal_SQ_4
        Fcal_SO_3=self.Fcal_SO_3
        Fcal_SO_5=self.Fcal_SO_5
        Fcal_SO_6=self.Fcal_SO_6
        Fcal_SO_7=self.Fcal_SO_7
        Fcal_SO_8=self.Fcal_SO_8
        E_0=self.E_0
        E_2=self.E_2
        E_4=self.E_4
        E_6=self.E_6
        E_8=self.E_8
        E_lnv_8=self.E_lnv_8
        E_SQ_4=self.E_SQ_4
        E_SO_3=self.E_SO_3
        E_SO_5=self.E_SO_5
        E_SO_7=self.E_SO_7
        v, rfrak_chi1_x, rfrak_chi1_y, rfrak_chi2_x, rfrak_chi2_y=self.v, self.rfrak_chi1_x, self.rfrak_chi1_y, self.rfrak_chi2_x, self.rfrak_chi2_y
        rfrak_frame_x, rfrak_frame_y, rfrak_frame_z, Phi=self.rfrak_frame_x, self.rfrak_frame_y, self.rfrak_frame_z, self.Phi
        dvdt_T1=0
        if v>=1.0 and self.NotForward:
            print("Beyond domain of PN validity, this is a good way to terminate.")
            self.terminal1=False
        Flux = Fcal_coeff*(Fcal_0 + v**2*(Fcal_2 + v*(Fcal_3 + Fcal_SO_3 + v*(Fcal_4 + Fcal_SQ_4 + v*(Fcal_5 + Fcal_SO_5 + v*(Fcal_6 + Fcal_SO_6 + Fcal_lnv_6*logv + v*(Fcal_7 + Fcal_SO_7 + v*(Fcal_8 + Fcal_SO_8 + Fcal_lnv_8*logv))))))))
        dEdv = -0.5*M*nu*v*(2.0*E_0 + v**2*(4.0*E_2 + v*(5.0*E_SO_3 + v*(6.0*E_4 + 6.0*E_SQ_4 + v*(7.0*E_SO_5 + v*(8.0*E_6 + v*(9.0*E_SO_7 + v*(10.0*E_8 + E_lnv_8*(10.0*logv + 1.0)))))))))
        Absorption = 0
        dvdt_T1 = (-Absorption - Flux)/dEdv
        if dvdt_T1<1.0e-12 and self.NotForward:
            print("v is decreasing, which is not an uncommon way to stop.")
            self.terminal2=False
        return self.CommonRHS(dvdt_T1, y)

    def CommonRHS(self,dvdt, y):
        dydt=np.zeros(9)
        rfrak_frame=np.zeros(3)
        rfrak_frame[0] = y[5]
        rfrak_frame[1] = y[6]
        rfrak_frame[2] = y[7]
        rfrakdot_frame = FrameFromAngularVelocityIntegrand(rfrak_frame, self.OmegaVec().vector)
        dydt[0] = dvdt
        if(self.EvolveSpin1):
            dydt[1], dydt[2]= FrameFromAngularVelocity_2D_Integrand(y[1], y[2],\
                (self.S_chi1.inverse*self.OmegaVec_chiVec_1()*self.S_chi1).vector)
        else:
            dydt[1] = 0.0
            dydt[2] = 0.0
        if(self.EvolveSpin2):
            dydt[3], dydt[4] = FrameFromAngularVelocity_2D_Integrand(y[3], y[4],\
                (self.S_chi2.inverse*self.OmegaVec_chiVec_2()*self.S_chi2).vector)
        else:
            dydt[3] = 0.0
            dydt[4] = 0.0
        dydt[5] = rfrakdot_frame[0]
        dydt[6] = rfrakdot_frame[1]
        dydt[7] = rfrakdot_frame[2]
        dydt[8] = self.v*self.v*self.v/self.M

        return dydt

    # Prepare for integration
    def Evolution(self): 
        v_i,rfrak_chi1_x_i,rfrak_chi1_y_i,rfrak_chi2_x_i,rfrak_chi2_y_i,rfrak_frame_x_i,rfrak_frame_y_i,rfrak_frame_z_i\
            =self.v_i,self.rfrak_chi1_x_i,self.rfrak_chi1_y_i,self.rfrak_chi2_x_i,self.rfrak_chi2_y_i,self.rfrak_frame_x_i,self.rfrak_frame_y_i,self.rfrak_frame_z_i
        def terminate(t,y):
            return 1.0*self.terminal1*self.terminal2
        terminate.terminal=True
        TMerger=5.0/(256.0*self.nu*v_i**8)
        TEnd=TMerger
        if self.t_PN2:
            TEnd=self.t_PN2
        time=[0.0]
        while time[-1]<TEnd and 2*self.M*(256*self.nu*(TMerger-time[-1])/5)**(3/8)/self.StepsPerOrbit>self.MinStep:
            time.append(time[-1]+2*self.M*(256*self.nu*(TMerger-time[-1])/5)**(3/8)/self.StepsPerOrbit)
        time=np.delete(time, -1)
        
        # Integrate
        yy=solve_ivp(self.__TaylorT1__, [time[0],time[-1]], [v_i,rfrak_chi1_x_i,\
            rfrak_chi1_y_i,rfrak_chi2_x_i,rfrak_chi2_y_i,rfrak_frame_x_i,\
            rfrak_frame_y_i,rfrak_frame_z_i,0.0], method='DOP853',\
            t_eval=time, dense_output=True, events=terminate, rtol=self.tol, atol=self.tol)           
        if self.ForwardInTime:
            self.NotForward=False
            time=[0.0]
            TStart=-3*TMerger
            if self.t_PN1:
                TStart=self.t_PN1
            while time[-1]>TStart:
                time.append(time[-1]-2*self.M*(256*self.nu*(TMerger-time[-1])/5)**(3/8)/self.StepsPerOrbit)
            yyForward=solve_ivp(self.__TaylorT1__, [time[0],time[-1]], [v_i,rfrak_chi1_x_i,\
                rfrak_chi1_y_i,rfrak_chi2_x_i,rfrak_chi2_y_i,rfrak_frame_x_i,\
                rfrak_frame_y_i,rfrak_frame_z_i,0.0], method='DOP853',\
                t_eval=time, dense_output=True, rtol=self.tol, atol=self.tol)
            yy.t=np.append(yyForward.t[1:][::-1],yy.t)
            data=np.empty((9,len(yy.t)))
            for i in range(9):
                data[i]=np.append(yyForward.y[i][1:][::-1],yy.y[i])
            yy.y=data
            
        return yy

spec=[("xHat", float64),("yHat", float64),("zHat", float64),("M1", float64),("M2", float64),("v", float64),("S_chi1", float64),("S_chi2", float64),("rfrak_chi1_x", float64),("rfrak_chi1_y", float64),("rfrak_chi2_x", float64),("rfrak_chi2_y", float64),("rfrak_frame_x", float64),("rfrak_frame_y", float64),("rfrak_frame_z", float64),("M", float64),("delta", float64),("nu", float64),("R", float64),("nHat", float64),("lambdaHat", float64),("ellHat", float64),("R_S1", float64),("R_S2", float64),("chiVec1", float64),("chiVec2", float64),("chi1chi1", float64),("chi1chi2", float64),("chi2chi2", float64),("chi1_n", float64),("chi1_lambda", float64),("chi1_ell", float64),("chi2_n", float64),("chi2_lambda", float64),("chi2_ell", float64),("S_ell", float64),("S_n", float64),("S_lambda", float64),("Sigma_ell", float64),("Sigma_n", float64),("Sigma_lambda", float64),("chi_s_ell", float64),("chi_a_ell", float64),("logv", float64),("Fcal_coeff", float64),("Fcal_0", float64),("Fcal_2", float64),("Fcal_3", float64),("Fcal_4", float64),("Fcal_5", float64),("Fcal_6", float64),("Fcal_lnv_6", float64),("Fcal_7", float64),("Fcal_8", float64),("Fcal_lnv_8", float64),("Fcal_9", float64),("Fcal_lnv_9", float64),("Fcal_SQ_4", float64),("Fcal_SO_3", float64),("Fcal_SO_5", float64),("Fcal_SO_6", float64),("Fcal_SO_7", float64),("Fcal_SO_8", float64),("E_0", float64),("E_2", float64),("E_4", float64),("E_6", float64),("E_8", float64),("E_lnv_8", float64),("E_SQ_4", float64),("E_SO_3", float64),("E_SO_5", float64),("E_SO_7", float64),("xHat_i", float64),("yHat_i", float64),("zHat_i", float64),("M1_i", float64),("M2_i", float64),("v_i", float64),("S_chi1_i", float64),("S_chi2_i", float64),("rfrak_chi1_x_i", float64),("rfrak_chi1_y_i", float64),("rfrak_chi2_x_i", float64),("rfrak_chi2_y_i", float64),("rfrak_frame_x_i", float64),("rfrak_frame_y_i", float64),("rfrak_frame_z_i", float64),("t_PN1",float64),("t_PN2",float64),("StepsPerOrbit",int32),("ForwardInTime",boolean),("tol",float64),("MinStep",float64),("EvolveSpin1",boolean),("EvolveSpin2",boolean),("terminal1",boolean),("terminal2",boolean),("NotForward",boolean)]
@jitclass(spec)
class TaylorTn_4p5PN_Q : 
    def __init__(self, xHat_i, yHat_i, zHat_i, M1_i, M2_i, v_i, S_chi1_i, S_chi2_i, rfrak_chi1_x_i, rfrak_chi1_y_i, rfrak_chi2_x_i, rfrak_chi2_y_i, rfrak_frame_x_i, rfrak_frame_y_i, rfrak_frame_z_i, t_PN1=False, t_PN2=False, StepsPerOrbit=32, ForwardInTime=True, tol=1e-12, MinStep=1e-7) :
        self.t_PN1,self.t_PN2,self.StepsPerOrbit,self.ForwardInTime,self.tol,self.MinStep=t_PN1, t_PN2, StepsPerOrbit, ForwardInTime, tol, MinStep
        terminal1=True
        terminal2=True
        NotForward=True
        self.xHat_i=xHat_i
        self.yHat_i=yHat_i
        self.zHat_i=zHat_i
        self.M1_i=M1_i
        self.M2_i=M2_i
        self.v_i=v_i
        self.S_chi1_i=S_chi1_i
        self.S_chi2_i=S_chi2_i
        self.rfrak_chi1_x_i=rfrak_chi1_x_i
        self.rfrak_chi1_y_i=rfrak_chi1_y_i
        self.rfrak_chi2_x_i=rfrak_chi2_x_i
        self.rfrak_chi2_y_i=rfrak_chi2_y_i
        self.rfrak_frame_x_i=rfrak_frame_x_i
        self.rfrak_frame_y_i=rfrak_frame_y_i
        self.rfrak_frame_z_i=rfrak_frame_z_i
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
        chi1chi1=dot(chiVec1.vector,chiVec1.vector)
        chi1chi2=dot(chiVec1.vector,chiVec2.vector)
        chi2chi2=dot(chiVec2.vector,chiVec2.vector)
        chi1_n=dot(chiVec1.vector,nHat.vector)
        chi1_lambda=dot(chiVec1.vector,lambdaHat.vector)
        chi1_ell=dot(chiVec1.vector,ellHat.vector)
        chi2_n=dot(chiVec2.vector,nHat.vector)
        chi2_lambda=dot(chiVec2.vector,lambdaHat.vector)
        chi2_ell=dot(chiVec2.vector,ellHat.vector)
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
        self.xHat=xHat
        self.yHat=yHat
        self.zHat=zHat
        self.M1=M1
        self.M2=M2
        self.v=v
        self.S_chi1=S_chi1
        self.S_chi2=S_chi2
        self.rfrak_chi1_x=rfrak_chi1_x
        self.rfrak_chi1_y=rfrak_chi1_y
        self.rfrak_chi2_x=rfrak_chi2_x
        self.rfrak_chi2_y=rfrak_chi2_y
        self.rfrak_frame_x=rfrak_frame_x
        self.rfrak_frame_y=rfrak_frame_y
        self.rfrak_frame_z=rfrak_frame_z
        self.M=M
        self.delta=delta
        self.nu=nu
        self.R=R
        self.nHat=nHat
        self.lambdaHat=lambdaHat
        self.ellHat=ellHat
        self.R_S1=R_S1
        self.R_S2=R_S2
        self.chiVec1=chiVec1
        self.chiVec2=chiVec2
        self.chi1chi1=chi1chi1
        self.chi1chi2=chi1chi2
        self.chi2chi2=chi2chi2
        self.chi1_n=chi1_n
        self.chi1_lambda=chi1_lambda
        self.chi1_ell=chi1_ell
        self.chi2_n=chi2_n
        self.chi2_lambda=chi2_lambda
        self.chi2_ell=chi2_ell
        self.S_ell=S_ell
        self.S_n=S_n
        self.S_lambda=S_lambda
        self.Sigma_ell=Sigma_ell
        self.Sigma_n=Sigma_n
        self.Sigma_lambda=Sigma_lambda
        self.chi_s_ell=chi_s_ell
        self.chi_a_ell=chi_a_ell
        self.logv=logv
        self.Fcal_coeff=Fcal_coeff
        self.Fcal_0=Fcal_0
        self.Fcal_2=Fcal_2
        self.Fcal_3=Fcal_3
        self.Fcal_4=Fcal_4
        self.Fcal_5=Fcal_5
        self.Fcal_6=Fcal_6
        self.Fcal_lnv_6=Fcal_lnv_6
        self.Fcal_7=Fcal_7
        self.Fcal_8=Fcal_8
        self.Fcal_lnv_8=Fcal_lnv_8
        self.Fcal_9=Fcal_9
        self.Fcal_lnv_9=Fcal_lnv_9
        self.Fcal_SQ_4=Fcal_SQ_4
        self.Fcal_SO_3=Fcal_SO_3
        self.Fcal_SO_5=Fcal_SO_5
        self.Fcal_SO_6=Fcal_SO_6
        self.Fcal_SO_7=Fcal_SO_7
        self.Fcal_SO_8=Fcal_SO_8
        self.E_0=E_0
        self.E_2=E_2
        self.E_4=E_4
        self.E_6=E_6
        self.E_8=E_8
        self.E_lnv_8=E_lnv_8
        self.E_SQ_4=E_SQ_4
        self.E_SO_3=E_SO_3
        self.E_SO_5=E_SO_5
        self.E_SO_7=E_SO_7
        Phi=0.0
        EvolveSpin1=(S_chi1*conjugate(S_chi1)).norm>1e-12
        EvolveSpin2=(S_chi2*conjugate(S_chi2)).norm>1e-12
        self.v, self.rfrak_chi1_x, self.rfrak_chi1_y, self.rfrak_chi2_x, self.rfrak_chi2_y=v, rfrak_chi1_x, rfrak_chi1_y, rfrak_chi2_x, rfrak_chi2_y
        self.rfrak_frame_x, self.rfrak_frame_y, self.rfrak_frame_z, self.Phi=rfrak_frame_x, rfrak_frame_y, rfrak_frame_z, Phi
        self.EvolveSpin1=EvolveSpin1
        self.EvolveSpin2=EvolveSpin2
        self.terminal1,self.terminal2,self.NotForward=terminal1,terminal2,NotForward

    def __Recalculate__(self,t, y):
        v = y[0]
        rfrak_chi1_x = y[1]
        rfrak_chi1_y = y[2]
        rfrak_chi2_x = y[3]
        rfrak_chi2_y = y[4]
        rfrak_frame_x = y[5]
        rfrak_frame_y = y[6]
        rfrak_frame_z = y[7]
        Phi = y[8]
        self.v, self.rfrak_chi1_x, self.rfrak_chi1_y, self.rfrak_chi2_x, self.rfrak_chi2_y=v, rfrak_chi1_x, rfrak_chi1_y, rfrak_chi2_x, rfrak_chi2_y
        self.rfrak_frame_x, self.rfrak_frame_y, self.rfrak_frame_z, self.Phi=rfrak_frame_x, rfrak_frame_y, rfrak_frame_z, Phi
        xHat=self.xHat
        yHat=self.yHat
        zHat=self.zHat
        M1=self.M1
        M2=self.M2
        v=self.v
        S_chi1=self.S_chi1
        S_chi2=self.S_chi2
        rfrak_chi1_x=self.rfrak_chi1_x
        rfrak_chi1_y=self.rfrak_chi1_y
        rfrak_chi2_x=self.rfrak_chi2_x
        rfrak_chi2_y=self.rfrak_chi2_y
        rfrak_frame_x=self.rfrak_frame_x
        rfrak_frame_y=self.rfrak_frame_y
        rfrak_frame_z=self.rfrak_frame_z
        M=self.M
        delta=self.delta
        nu=self.nu
        R=self.R
        nHat=self.nHat
        lambdaHat=self.lambdaHat
        ellHat=self.ellHat
        R_S1=self.R_S1
        R_S2=self.R_S2
        chiVec1=self.chiVec1
        chiVec2=self.chiVec2
        chi1chi1=self.chi1chi1
        chi1chi2=self.chi1chi2
        chi2chi2=self.chi2chi2
        chi1_n=self.chi1_n
        chi1_lambda=self.chi1_lambda
        chi1_ell=self.chi1_ell
        chi2_n=self.chi2_n
        chi2_lambda=self.chi2_lambda
        chi2_ell=self.chi2_ell
        S_ell=self.S_ell
        S_n=self.S_n
        S_lambda=self.S_lambda
        Sigma_ell=self.Sigma_ell
        Sigma_n=self.Sigma_n
        Sigma_lambda=self.Sigma_lambda
        chi_s_ell=self.chi_s_ell
        chi_a_ell=self.chi_a_ell
        logv=self.logv
        Fcal_coeff=self.Fcal_coeff
        Fcal_0=self.Fcal_0
        Fcal_2=self.Fcal_2
        Fcal_3=self.Fcal_3
        Fcal_4=self.Fcal_4
        Fcal_5=self.Fcal_5
        Fcal_6=self.Fcal_6
        Fcal_lnv_6=self.Fcal_lnv_6
        Fcal_7=self.Fcal_7
        Fcal_8=self.Fcal_8
        Fcal_lnv_8=self.Fcal_lnv_8
        Fcal_9=self.Fcal_9
        Fcal_lnv_9=self.Fcal_lnv_9
        Fcal_SQ_4=self.Fcal_SQ_4
        Fcal_SO_3=self.Fcal_SO_3
        Fcal_SO_5=self.Fcal_SO_5
        Fcal_SO_6=self.Fcal_SO_6
        Fcal_SO_7=self.Fcal_SO_7
        Fcal_SO_8=self.Fcal_SO_8
        E_0=self.E_0
        E_2=self.E_2
        E_4=self.E_4
        E_6=self.E_6
        E_8=self.E_8
        E_lnv_8=self.E_lnv_8
        E_SQ_4=self.E_SQ_4
        E_SO_3=self.E_SO_3
        E_SO_5=self.E_SO_5
        E_SO_7=self.E_SO_7
        R = exp(rfrak_frame_x*xHat + rfrak_frame_y*yHat + rfrak_frame_z*zHat)
        nHat = R*xHat*conjugate(R)
        lambdaHat = R*yHat*conjugate(R)
        ellHat = R*zHat*conjugate(R)
        R_S1 = exp(rfrak_chi1_x*xHat + rfrak_chi1_y*yHat)
        R_S2 = exp(rfrak_chi2_x*xHat + rfrak_chi2_y*yHat)
        chiVec1 = S_chi1*R_S1*zHat*conjugate(R_S1)*conjugate(S_chi1)
        chiVec2 = S_chi2*R_S2*zHat*conjugate(R_S2)*conjugate(S_chi2)
        chi1_n = dot(chiVec1.vector,nHat.vector)
        chi1_lambda = dot(chiVec1.vector,lambdaHat.vector)
        chi1_ell = dot(chiVec1.vector,ellHat.vector)
        chi2_n = dot(chiVec2.vector,nHat.vector)
        chi2_lambda = dot(chiVec2.vector,lambdaHat.vector)
        chi2_ell = dot(chiVec2.vector,ellHat.vector)
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
        self.R=R
        self.nHat=nHat
        self.lambdaHat=lambdaHat
        self.ellHat=ellHat
        self.R_S1=R_S1
        self.R_S2=R_S2
        self.chiVec1=chiVec1
        self.chiVec2=chiVec2
        self.chi1_n=chi1_n
        self.chi1_lambda=chi1_lambda
        self.chi1_ell=chi1_ell
        self.chi2_n=chi2_n
        self.chi2_lambda=chi2_lambda
        self.chi2_ell=chi2_ell
        self.S_ell=S_ell
        self.S_n=S_n
        self.S_lambda=S_lambda
        self.Sigma_ell=Sigma_ell
        self.Sigma_n=Sigma_n
        self.Sigma_lambda=Sigma_lambda
        self.chi_s_ell=chi_s_ell
        self.chi_a_ell=chi_a_ell
        self.logv=logv
        self.Fcal_coeff=Fcal_coeff
        self.Fcal_SQ_4=Fcal_SQ_4
        self.Fcal_SO_3=Fcal_SO_3
        self.Fcal_SO_5=Fcal_SO_5
        self.Fcal_SO_6=Fcal_SO_6
        self.Fcal_SO_7=Fcal_SO_7
        self.Fcal_SO_8=Fcal_SO_8
        self.E_SQ_4=E_SQ_4
        self.E_SO_3=E_SO_3
        self.E_SO_5=E_SO_5
        self.E_SO_7=E_SO_7

    def OmegaVec_chiVec_1(self):
        xHat=self.xHat
        yHat=self.yHat
        zHat=self.zHat
        M1=self.M1
        M2=self.M2
        v=self.v
        S_chi1=self.S_chi1
        S_chi2=self.S_chi2
        rfrak_chi1_x=self.rfrak_chi1_x
        rfrak_chi1_y=self.rfrak_chi1_y
        rfrak_chi2_x=self.rfrak_chi2_x
        rfrak_chi2_y=self.rfrak_chi2_y
        rfrak_frame_x=self.rfrak_frame_x
        rfrak_frame_y=self.rfrak_frame_y
        rfrak_frame_z=self.rfrak_frame_z
        M=self.M
        delta=self.delta
        nu=self.nu
        R=self.R
        nHat=self.nHat
        lambdaHat=self.lambdaHat
        ellHat=self.ellHat
        R_S1=self.R_S1
        R_S2=self.R_S2
        chiVec1=self.chiVec1
        chiVec2=self.chiVec2
        chi1_n=self.chi1_n
        chi1_lambda=self.chi1_lambda
        chi1_ell=self.chi1_ell
        chi2_n=self.chi2_n
        chi2_lambda=self.chi2_lambda
        chi2_ell=self.chi2_ell
        S_ell=self.S_ell
        S_n=self.S_n
        S_lambda=self.S_lambda
        Sigma_ell=self.Sigma_ell
        Sigma_n=self.Sigma_n
        Sigma_lambda=self.Sigma_lambda
        Omega1_coeff = v**5/M
        return Omega1_coeff*(ellHat*(-0.75*delta + 0.5*nu + v**2*(delta*(0.625*nu - 0.5625) + nu*(1.25 - 0.0416666666666667*nu) + v**2*(delta*(nu*(4.875 - 0.15625*nu) - 0.84375) + nu*(nu*(-0.0208333333333333*nu - 3.28125) + 0.1875) + 0.84375) + 0.5625) + 0.75) + nHat*v*(3.0*chi1_n*nu + 3.0*M2**2*chi2_n/M**2) - M2**2*chiVec2*v/M**2)

    def OmegaVec_chiVec_2(self):
        xHat=self.xHat
        yHat=self.yHat
        zHat=self.zHat
        M1=self.M1
        M2=self.M2
        v=self.v
        S_chi1=self.S_chi1
        S_chi2=self.S_chi2
        rfrak_chi1_x=self.rfrak_chi1_x
        rfrak_chi1_y=self.rfrak_chi1_y
        rfrak_chi2_x=self.rfrak_chi2_x
        rfrak_chi2_y=self.rfrak_chi2_y
        rfrak_frame_x=self.rfrak_frame_x
        rfrak_frame_y=self.rfrak_frame_y
        rfrak_frame_z=self.rfrak_frame_z
        M=self.M
        delta=self.delta
        nu=self.nu
        R=self.R
        nHat=self.nHat
        lambdaHat=self.lambdaHat
        ellHat=self.ellHat
        R_S1=self.R_S1
        R_S2=self.R_S2
        chiVec1=self.chiVec1
        chiVec2=self.chiVec2
        chi1_n=self.chi1_n
        chi1_lambda=self.chi1_lambda
        chi1_ell=self.chi1_ell
        chi2_n=self.chi2_n
        chi2_lambda=self.chi2_lambda
        chi2_ell=self.chi2_ell
        S_ell=self.S_ell
        S_n=self.S_n
        S_lambda=self.S_lambda
        Sigma_ell=self.Sigma_ell
        Sigma_n=self.Sigma_n
        Sigma_lambda=self.Sigma_lambda
        Omega2_coeff = v**5/M
        return Omega2_coeff*(ellHat*(0.75*delta + 0.5*nu + v**2*(delta*(0.5625 - 0.625*nu) + nu*(1.25 - 0.0416666666666667*nu) + v**2*(delta*(nu*(0.15625*nu - 4.875) + 0.84375) + nu*(nu*(-0.0208333333333333*nu - 3.28125) + 0.1875) + 0.84375) + 0.5625) + 0.75) + nHat*v*(3.0*chi2_n*nu + 3.0*M1**2*chi1_n/M**2) - M1**2*chiVec1*v/M**2)

    def OmegaVec(self):
        xHat=self.xHat
        yHat=self.yHat
        zHat=self.zHat
        M1=self.M1
        M2=self.M2
        v=self.v
        S_chi1=self.S_chi1
        S_chi2=self.S_chi2
        rfrak_chi1_x=self.rfrak_chi1_x
        rfrak_chi1_y=self.rfrak_chi1_y
        rfrak_chi2_x=self.rfrak_chi2_x
        rfrak_chi2_y=self.rfrak_chi2_y
        rfrak_frame_x=self.rfrak_frame_x
        rfrak_frame_y=self.rfrak_frame_y
        rfrak_frame_z=self.rfrak_frame_z
        M=self.M
        delta=self.delta
        nu=self.nu
        R=self.R
        nHat=self.nHat
        lambdaHat=self.lambdaHat
        ellHat=self.ellHat
        R_S1=self.R_S1
        R_S2=self.R_S2
        chiVec1=self.chiVec1
        chiVec2=self.chiVec2
        chi1_n=self.chi1_n
        chi1_lambda=self.chi1_lambda
        chi1_ell=self.chi1_ell
        chi2_n=self.chi2_n
        chi2_lambda=self.chi2_lambda
        chi2_ell=self.chi2_ell
        S_ell=self.S_ell
        S_n=self.S_n
        S_lambda=self.S_lambda
        Sigma_ell=self.Sigma_ell
        Sigma_n=self.Sigma_n
        Sigma_lambda=self.Sigma_lambda
        gamma_PN_6 = 0.0123456790123457*nu**3 + 6.36111111111111*nu**2 - 2.98177812235564*nu + 1.0
        gamma_PN_4 = 1.0 - 5.41666666666667*nu
        gamma_PN_5 = (S_ell*(0.888888888888889*nu + 3.33333333333333) + 2.0*Sigma_ell*delta)/M**2
        gamma_PN_2 = 1.0 - 0.333333333333333*nu
        a_ell_4 = S_n*(5.77777777777778*nu**2 + 14.75*nu + 1.5) + Sigma_n*delta*(2.83333333333333*nu**2 + 9.125*nu + 1.5)
        gamma_PN_0 = 1.00000000000000
        gamma_PN_7 = (S_ell*(-6.0*nu**2 - 10.5833333333333*nu + 5.0) - 2.66666666666667*Sigma_ell*delta*nu**2 + Sigma_ell*delta*(3.0 - 10.1666666666667*nu))/M**2
        gamma_PN_3 = (1.66666666666667*S_ell + Sigma_ell*delta)/M**2
        a_ell_0 = 7.0*S_n + 3.0*Sigma_n*delta
        a_ell_2 = S_n*(-9.66666666666667*nu - 10.0) + Sigma_n*delta*(-4.5*nu - 6.0)
        return ellHat*v**3/M + nHat*v**6*(a_ell_0 + v**2*(a_ell_2 + a_ell_4*v**2))*(gamma_PN_0 + v**2*(gamma_PN_2 + v*(gamma_PN_3 + v*(gamma_PN_4 + v*(gamma_PN_5 + v*(gamma_PN_6 + gamma_PN_7*v))))))/M**3

    def OrbitalAngularMomentum(self):
        xHat=self.xHat
        yHat=self.yHat
        zHat=self.zHat
        M1=self.M1
        M2=self.M2
        v=self.v
        S_chi1=self.S_chi1
        S_chi2=self.S_chi2
        rfrak_chi1_x=self.rfrak_chi1_x
        rfrak_chi1_y=self.rfrak_chi1_y
        rfrak_chi2_x=self.rfrak_chi2_x
        rfrak_chi2_y=self.rfrak_chi2_y
        rfrak_frame_x=self.rfrak_frame_x
        rfrak_frame_y=self.rfrak_frame_y
        rfrak_frame_z=self.rfrak_frame_z
        M=self.M
        delta=self.delta
        nu=self.nu
        R=self.R
        nHat=self.nHat
        lambdaHat=self.lambdaHat
        ellHat=self.ellHat
        R_S1=self.R_S1
        R_S2=self.R_S2
        chiVec1=self.chiVec1
        chiVec2=self.chiVec2
        chi1_n=self.chi1_n
        chi1_lambda=self.chi1_lambda
        chi1_ell=self.chi1_ell
        chi2_n=self.chi2_n
        chi2_lambda=self.chi2_lambda
        chi2_ell=self.chi2_ell
        S_ell=self.S_ell
        S_n=self.S_n
        S_lambda=self.S_lambda
        Sigma_ell=self.Sigma_ell
        Sigma_n=self.Sigma_n
        Sigma_lambda=self.Sigma_lambda
        L_lnv_8 = -42.6666666666667*ellHat*nu
        L_6 = ellHat*(0.00540123456790123*nu**3 + 1.29166666666667*nu**2 - 30.9797035925835*nu + 8.4375)
        L_SO_7 = 0.0625*ellHat*(-29.0*S_ell*nu**2 + 1101.0*S_ell*nu - 405.0*S_ell - 15.0*Sigma_ell*delta*nu**2 + 468.0*Sigma_ell*delta*nu - 81.0*Sigma_ell*delta)/M**2 + 0.0416666666666667*lambdaHat*(-32.0*S_lambda*nu**2 + 2.0*S_lambda*nu - 174.0*S_lambda - 16.0*Sigma_lambda*delta*nu**2 - 79.0*Sigma_lambda*delta*nu - 12.0*Sigma_lambda*delta)/M**2 + 0.0208333333333333*nHat*(11.0*S_n*nu**2 - 1331.0*S_n*nu + 183.0*S_n + 5.0*Sigma_n*delta*nu**2 - 734.0*Sigma_n*delta*nu + 183.0*Sigma_n*delta)/M**2
        L_SO_5 = 0.0138888888888889*ellHat*(427.0*S_ell*nu - 693.0*S_ell + 210.0*Sigma_ell*delta*nu - 189.0*Sigma_ell*delta)/M**2 + 0.166666666666667*lambdaHat*(18.0*S_lambda*nu - 21.0*S_lambda + 8.0*Sigma_lambda*delta*nu - 3.0*Sigma_lambda*delta)/M**2 + 0.0416666666666667*nHat*(-19.0*S_n*nu + 33.0*S_n - 10.0*Sigma_n*delta*nu + 33.0*Sigma_n*delta)/M**2
        L_SO_3 = 0.166666666666667*ellHat*(-35.0*S_ell - 15.0*Sigma_ell*delta)/M**2 + lambdaHat*(-3.0*S_lambda - Sigma_lambda*delta)/M**2 + 0.5*nHat*(S_n + Sigma_n*delta)/M**2
        L_2 = ellHat*(0.166666666666667*nu + 1.5)
        L_4 = ellHat*(0.0416666666666667*nu**2 - 2.375*nu + 3.375)
        L_0 = ellHat
        L_8 = ellHat*(-0.714285714285714*nu*(0.0024755658436214*nu**3 + 0.174189814814815*nu**2 - 90.1327990262052*nu + 153.88379682994) + 1.82857142857143*nu + 22.1484375)
        L_coeff = M**2*nu/v
        return L_coeff*(L_0 + v**2*(L_2 + v*(L_SO_3 + v*(L_4 + v*(L_SO_5 + v*(L_6 + v*(L_SO_7 + v*(L_8 + L_lnv_8*math.log(v)))))))))

               

    #Evolve PN
    def __TaylorT1__(self,t, y):
        self.__Recalculate__(t, y)
        xHat=self.xHat
        yHat=self.yHat
        zHat=self.zHat
        M1=self.M1
        M2=self.M2
        v=self.v
        S_chi1=self.S_chi1
        S_chi2=self.S_chi2
        rfrak_chi1_x=self.rfrak_chi1_x
        rfrak_chi1_y=self.rfrak_chi1_y
        rfrak_chi2_x=self.rfrak_chi2_x
        rfrak_chi2_y=self.rfrak_chi2_y
        rfrak_frame_x=self.rfrak_frame_x
        rfrak_frame_y=self.rfrak_frame_y
        rfrak_frame_z=self.rfrak_frame_z
        M=self.M
        delta=self.delta
        nu=self.nu
        R=self.R
        nHat=self.nHat
        lambdaHat=self.lambdaHat
        ellHat=self.ellHat
        R_S1=self.R_S1
        R_S2=self.R_S2
        chiVec1=self.chiVec1
        chiVec2=self.chiVec2
        chi1chi1=self.chi1chi1
        chi1chi2=self.chi1chi2
        chi2chi2=self.chi2chi2
        chi1_n=self.chi1_n
        chi1_lambda=self.chi1_lambda
        chi1_ell=self.chi1_ell
        chi2_n=self.chi2_n
        chi2_lambda=self.chi2_lambda
        chi2_ell=self.chi2_ell
        S_ell=self.S_ell
        S_n=self.S_n
        S_lambda=self.S_lambda
        Sigma_ell=self.Sigma_ell
        Sigma_n=self.Sigma_n
        Sigma_lambda=self.Sigma_lambda
        chi_s_ell=self.chi_s_ell
        chi_a_ell=self.chi_a_ell
        logv=self.logv
        Fcal_coeff=self.Fcal_coeff
        Fcal_0=self.Fcal_0
        Fcal_2=self.Fcal_2
        Fcal_3=self.Fcal_3
        Fcal_4=self.Fcal_4
        Fcal_5=self.Fcal_5
        Fcal_6=self.Fcal_6
        Fcal_lnv_6=self.Fcal_lnv_6
        Fcal_7=self.Fcal_7
        Fcal_8=self.Fcal_8
        Fcal_lnv_8=self.Fcal_lnv_8
        Fcal_9=self.Fcal_9
        Fcal_lnv_9=self.Fcal_lnv_9
        Fcal_SQ_4=self.Fcal_SQ_4
        Fcal_SO_3=self.Fcal_SO_3
        Fcal_SO_5=self.Fcal_SO_5
        Fcal_SO_6=self.Fcal_SO_6
        Fcal_SO_7=self.Fcal_SO_7
        Fcal_SO_8=self.Fcal_SO_8
        E_0=self.E_0
        E_2=self.E_2
        E_4=self.E_4
        E_6=self.E_6
        E_8=self.E_8
        E_lnv_8=self.E_lnv_8
        E_SQ_4=self.E_SQ_4
        E_SO_3=self.E_SO_3
        E_SO_5=self.E_SO_5
        E_SO_7=self.E_SO_7
        v, rfrak_chi1_x, rfrak_chi1_y, rfrak_chi2_x, rfrak_chi2_y=self.v, self.rfrak_chi1_x, self.rfrak_chi1_y, self.rfrak_chi2_x, self.rfrak_chi2_y
        rfrak_frame_x, rfrak_frame_y, rfrak_frame_z, Phi=self.rfrak_frame_x, self.rfrak_frame_y, self.rfrak_frame_z, self.Phi
        dvdt_T1=0
        if v>=1.0 and self.NotForward:
            print("Beyond domain of PN validity, this is a good way to terminate.")
            self.terminal1=False
        Flux = Fcal_coeff*(Fcal_0 + v**2*(Fcal_2 + v*(Fcal_3 + Fcal_SO_3 + v*(Fcal_4 + Fcal_SQ_4 + v*(Fcal_5 + Fcal_SO_5 + v*(Fcal_6 + Fcal_SO_6 + Fcal_lnv_6*logv + v*(Fcal_7 + Fcal_SO_7 + v*(Fcal_8 + Fcal_SO_8 + Fcal_lnv_8*logv + v*(Fcal_9 + Fcal_lnv_9*logv)))))))))
        dEdv = -0.5*M*nu*v*(2.0*E_0 + v**2*(4.0*E_2 + v*(5.0*E_SO_3 + v*(6.0*E_4 + 6.0*E_SQ_4 + v*(7.0*E_SO_5 + v*(8.0*E_6 + v*(9.0*E_SO_7 + v*(10.0*E_8 + E_lnv_8*(10.0*logv + 1.0)))))))))
        Absorption = 0
        dvdt_T1 = (-Absorption - Flux)/dEdv
        if dvdt_T1<1.0e-12 and self.NotForward:
            print("v is decreasing, which is not an uncommon way to stop.")
            self.terminal2=False
        return self.CommonRHS(dvdt_T1, y)

    def CommonRHS(self,dvdt, y):
        dydt=np.zeros(9)
        rfrak_frame=np.zeros(3)
        rfrak_frame[0] = y[5]
        rfrak_frame[1] = y[6]
        rfrak_frame[2] = y[7]
        rfrakdot_frame = FrameFromAngularVelocityIntegrand(rfrak_frame, self.OmegaVec().vector)
        dydt[0] = dvdt
        if(self.EvolveSpin1):
            dydt[1], dydt[2]= FrameFromAngularVelocity_2D_Integrand(y[1], y[2],\
                (self.S_chi1.inverse*self.OmegaVec_chiVec_1()*self.S_chi1).vector)
        else:
            dydt[1] = 0.0
            dydt[2] = 0.0
        if(self.EvolveSpin2):
            dydt[3], dydt[4] = FrameFromAngularVelocity_2D_Integrand(y[3], y[4],\
                (self.S_chi2.inverse*self.OmegaVec_chiVec_2()*self.S_chi2).vector)
        else:
            dydt[3] = 0.0
            dydt[4] = 0.0
        dydt[5] = rfrakdot_frame[0]
        dydt[6] = rfrakdot_frame[1]
        dydt[7] = rfrakdot_frame[2]
        dydt[8] = self.v*self.v*self.v/self.M

        return dydt

    # Prepare for integration
    def Evolution(self): 
        v_i,rfrak_chi1_x_i,rfrak_chi1_y_i,rfrak_chi2_x_i,rfrak_chi2_y_i,rfrak_frame_x_i,rfrak_frame_y_i,rfrak_frame_z_i\
            =self.v_i,self.rfrak_chi1_x_i,self.rfrak_chi1_y_i,self.rfrak_chi2_x_i,self.rfrak_chi2_y_i,self.rfrak_frame_x_i,self.rfrak_frame_y_i,self.rfrak_frame_z_i
        def terminate(t,y):
            return 1.0*self.terminal1*self.terminal2
        terminate.terminal=True
        TMerger=5.0/(256.0*self.nu*v_i**8)
        TEnd=TMerger
        if self.t_PN2:
            TEnd=self.t_PN2
        time=[0.0]
        while time[-1]<TEnd and 2*self.M*(256*self.nu*(TMerger-time[-1])/5)**(3/8)/self.StepsPerOrbit>self.MinStep:
            time.append(time[-1]+2*self.M*(256*self.nu*(TMerger-time[-1])/5)**(3/8)/self.StepsPerOrbit)
        time=np.delete(time, -1)
        
        # Integrate
        yy=solve_ivp(self.__TaylorT1__, [time[0],time[-1]], [v_i,rfrak_chi1_x_i,\
            rfrak_chi1_y_i,rfrak_chi2_x_i,rfrak_chi2_y_i,rfrak_frame_x_i,\
            rfrak_frame_y_i,rfrak_frame_z_i,0.0], method='DOP853',\
            t_eval=time, dense_output=True, events=terminate, rtol=self.tol, atol=self.tol)           
        if self.ForwardInTime:
            self.NotForward=False
            time=[0.0]
            TStart=-3*TMerger
            if self.t_PN1:
                TStart=self.t_PN1
            while time[-1]>TStart:
                time.append(time[-1]-2*self.M*(256*self.nu*(TMerger-time[-1])/5)**(3/8)/self.StepsPerOrbit)
            yyForward=solve_ivp(self.__TaylorT1__, [time[0],time[-1]], [v_i,rfrak_chi1_x_i,\
                rfrak_chi1_y_i,rfrak_chi2_x_i,rfrak_chi2_y_i,rfrak_frame_x_i,\
                rfrak_frame_y_i,rfrak_frame_z_i,0.0], method='DOP853',\
                t_eval=time, dense_output=True, rtol=self.tol, atol=self.tol)
            yy.t=np.append(yyForward.t[1:][::-1],yy.t)
            data=np.empty((9,len(yy.t)))
            for i in range(9):
                data[i]=np.append(yyForward.y[i][1:][::-1],yy.y[i])
            yy.y=data
            
        return yy

spec=[("xHat", float64),("yHat", float64),("zHat", float64),("M1", float64),("M2", float64),("v", float64),("S_chi1", float64),("S_chi2", float64),("rfrak_chi1_x", float64),("rfrak_chi1_y", float64),("rfrak_chi2_x", float64),("rfrak_chi2_y", float64),("rfrak_frame_x", float64),("rfrak_frame_y", float64),("rfrak_frame_z", float64),("M", float64),("delta", float64),("nu", float64),("R", float64),("nHat", float64),("lambdaHat", float64),("ellHat", float64),("R_S1", float64),("R_S2", float64),("chiVec1", float64),("chiVec2", float64),("chi1chi1", float64),("chi1chi2", float64),("chi2chi2", float64),("chi1_n", float64),("chi1_lambda", float64),("chi1_ell", float64),("chi2_n", float64),("chi2_lambda", float64),("chi2_ell", float64),("S_ell", float64),("S_n", float64),("S_lambda", float64),("Sigma_ell", float64),("Sigma_n", float64),("Sigma_lambda", float64),("chi_s_ell", float64),("chi_a_ell", float64),("logv", float64),("Fcal_coeff", float64),("Fcal_0", float64),("Fcal_2", float64),("Fcal_3", float64),("Fcal_4", float64),("Fcal_5", float64),("Fcal_6", float64),("Fcal_lnv_6", float64),("Fcal_7", float64),("Fcal_8", float64),("Fcal_lnv_8", float64),("Fcal_9", float64),("Fcal_lnv_9", float64),("Fcal_10", float64),("Fcal_lnv_10", float64),("Fcal_SQ_4", float64),("Fcal_SO_3", float64),("Fcal_SO_5", float64),("Fcal_SO_6", float64),("Fcal_SO_7", float64),("Fcal_SO_8", float64),("E_0", float64),("E_2", float64),("E_4", float64),("E_6", float64),("E_8", float64),("E_lnv_8", float64),("E_10", float64),("E_lnv_10", float64),("E_SQ_4", float64),("E_SO_3", float64),("E_SO_5", float64),("E_SO_7", float64),("xHat_i", float64),("yHat_i", float64),("zHat_i", float64),("M1_i", float64),("M2_i", float64),("v_i", float64),("S_chi1_i", float64),("S_chi2_i", float64),("rfrak_chi1_x_i", float64),("rfrak_chi1_y_i", float64),("rfrak_chi2_x_i", float64),("rfrak_chi2_y_i", float64),("rfrak_frame_x_i", float64),("rfrak_frame_y_i", float64),("rfrak_frame_z_i", float64),("t_PN1",float64),("t_PN2",float64),("StepsPerOrbit",int32),("ForwardInTime",boolean),("tol",float64),("MinStep",float64),("EvolveSpin1",boolean),("EvolveSpin2",boolean),("terminal1",boolean),("terminal2",boolean),("NotForward",boolean)]
@jitclass(spec)
class TaylorTn_5p0PN_Q : 
    def __init__(self, xHat_i, yHat_i, zHat_i, M1_i, M2_i, v_i, S_chi1_i, S_chi2_i, rfrak_chi1_x_i, rfrak_chi1_y_i, rfrak_chi2_x_i, rfrak_chi2_y_i, rfrak_frame_x_i, rfrak_frame_y_i, rfrak_frame_z_i, t_PN1=False, t_PN2=False, StepsPerOrbit=32, ForwardInTime=True, tol=1e-12, MinStep=1e-7) :
        self.t_PN1,self.t_PN2,self.StepsPerOrbit,self.ForwardInTime,self.tol,self.MinStep=t_PN1, t_PN2, StepsPerOrbit, ForwardInTime, tol, MinStep
        terminal1=True
        terminal2=True
        NotForward=True
        self.xHat_i=xHat_i
        self.yHat_i=yHat_i
        self.zHat_i=zHat_i
        self.M1_i=M1_i
        self.M2_i=M2_i
        self.v_i=v_i
        self.S_chi1_i=S_chi1_i
        self.S_chi2_i=S_chi2_i
        self.rfrak_chi1_x_i=rfrak_chi1_x_i
        self.rfrak_chi1_y_i=rfrak_chi1_y_i
        self.rfrak_chi2_x_i=rfrak_chi2_x_i
        self.rfrak_chi2_y_i=rfrak_chi2_y_i
        self.rfrak_frame_x_i=rfrak_frame_x_i
        self.rfrak_frame_y_i=rfrak_frame_y_i
        self.rfrak_frame_z_i=rfrak_frame_z_i
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
        chi1chi1=dot(chiVec1.vector,chiVec1.vector)
        chi1chi2=dot(chiVec1.vector,chiVec2.vector)
        chi2chi2=dot(chiVec2.vector,chiVec2.vector)
        chi1_n=dot(chiVec1.vector,nHat.vector)
        chi1_lambda=dot(chiVec1.vector,lambdaHat.vector)
        chi1_ell=dot(chiVec1.vector,ellHat.vector)
        chi2_n=dot(chiVec2.vector,nHat.vector)
        chi2_lambda=dot(chiVec2.vector,lambdaHat.vector)
        chi2_ell=dot(chiVec2.vector,ellHat.vector)
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
        self.xHat=xHat
        self.yHat=yHat
        self.zHat=zHat
        self.M1=M1
        self.M2=M2
        self.v=v
        self.S_chi1=S_chi1
        self.S_chi2=S_chi2
        self.rfrak_chi1_x=rfrak_chi1_x
        self.rfrak_chi1_y=rfrak_chi1_y
        self.rfrak_chi2_x=rfrak_chi2_x
        self.rfrak_chi2_y=rfrak_chi2_y
        self.rfrak_frame_x=rfrak_frame_x
        self.rfrak_frame_y=rfrak_frame_y
        self.rfrak_frame_z=rfrak_frame_z
        self.M=M
        self.delta=delta
        self.nu=nu
        self.R=R
        self.nHat=nHat
        self.lambdaHat=lambdaHat
        self.ellHat=ellHat
        self.R_S1=R_S1
        self.R_S2=R_S2
        self.chiVec1=chiVec1
        self.chiVec2=chiVec2
        self.chi1chi1=chi1chi1
        self.chi1chi2=chi1chi2
        self.chi2chi2=chi2chi2
        self.chi1_n=chi1_n
        self.chi1_lambda=chi1_lambda
        self.chi1_ell=chi1_ell
        self.chi2_n=chi2_n
        self.chi2_lambda=chi2_lambda
        self.chi2_ell=chi2_ell
        self.S_ell=S_ell
        self.S_n=S_n
        self.S_lambda=S_lambda
        self.Sigma_ell=Sigma_ell
        self.Sigma_n=Sigma_n
        self.Sigma_lambda=Sigma_lambda
        self.chi_s_ell=chi_s_ell
        self.chi_a_ell=chi_a_ell
        self.logv=logv
        self.Fcal_coeff=Fcal_coeff
        self.Fcal_0=Fcal_0
        self.Fcal_2=Fcal_2
        self.Fcal_3=Fcal_3
        self.Fcal_4=Fcal_4
        self.Fcal_5=Fcal_5
        self.Fcal_6=Fcal_6
        self.Fcal_lnv_6=Fcal_lnv_6
        self.Fcal_7=Fcal_7
        self.Fcal_8=Fcal_8
        self.Fcal_lnv_8=Fcal_lnv_8
        self.Fcal_9=Fcal_9
        self.Fcal_lnv_9=Fcal_lnv_9
        self.Fcal_10=Fcal_10
        self.Fcal_lnv_10=Fcal_lnv_10
        self.Fcal_SQ_4=Fcal_SQ_4
        self.Fcal_SO_3=Fcal_SO_3
        self.Fcal_SO_5=Fcal_SO_5
        self.Fcal_SO_6=Fcal_SO_6
        self.Fcal_SO_7=Fcal_SO_7
        self.Fcal_SO_8=Fcal_SO_8
        self.E_0=E_0
        self.E_2=E_2
        self.E_4=E_4
        self.E_6=E_6
        self.E_8=E_8
        self.E_lnv_8=E_lnv_8
        self.E_10=E_10
        self.E_lnv_10=E_lnv_10
        self.E_SQ_4=E_SQ_4
        self.E_SO_3=E_SO_3
        self.E_SO_5=E_SO_5
        self.E_SO_7=E_SO_7
        Phi=0.0
        EvolveSpin1=(S_chi1*conjugate(S_chi1)).norm>1e-12
        EvolveSpin2=(S_chi2*conjugate(S_chi2)).norm>1e-12
        self.v, self.rfrak_chi1_x, self.rfrak_chi1_y, self.rfrak_chi2_x, self.rfrak_chi2_y=v, rfrak_chi1_x, rfrak_chi1_y, rfrak_chi2_x, rfrak_chi2_y
        self.rfrak_frame_x, self.rfrak_frame_y, self.rfrak_frame_z, self.Phi=rfrak_frame_x, rfrak_frame_y, rfrak_frame_z, Phi
        self.EvolveSpin1=EvolveSpin1
        self.EvolveSpin2=EvolveSpin2
        self.terminal1,self.terminal2,self.NotForward=terminal1,terminal2,NotForward

    def __Recalculate__(self,t, y):
        v = y[0]
        rfrak_chi1_x = y[1]
        rfrak_chi1_y = y[2]
        rfrak_chi2_x = y[3]
        rfrak_chi2_y = y[4]
        rfrak_frame_x = y[5]
        rfrak_frame_y = y[6]
        rfrak_frame_z = y[7]
        Phi = y[8]
        self.v, self.rfrak_chi1_x, self.rfrak_chi1_y, self.rfrak_chi2_x, self.rfrak_chi2_y=v, rfrak_chi1_x, rfrak_chi1_y, rfrak_chi2_x, rfrak_chi2_y
        self.rfrak_frame_x, self.rfrak_frame_y, self.rfrak_frame_z, self.Phi=rfrak_frame_x, rfrak_frame_y, rfrak_frame_z, Phi
        xHat=self.xHat
        yHat=self.yHat
        zHat=self.zHat
        M1=self.M1
        M2=self.M2
        v=self.v
        S_chi1=self.S_chi1
        S_chi2=self.S_chi2
        rfrak_chi1_x=self.rfrak_chi1_x
        rfrak_chi1_y=self.rfrak_chi1_y
        rfrak_chi2_x=self.rfrak_chi2_x
        rfrak_chi2_y=self.rfrak_chi2_y
        rfrak_frame_x=self.rfrak_frame_x
        rfrak_frame_y=self.rfrak_frame_y
        rfrak_frame_z=self.rfrak_frame_z
        M=self.M
        delta=self.delta
        nu=self.nu
        R=self.R
        nHat=self.nHat
        lambdaHat=self.lambdaHat
        ellHat=self.ellHat
        R_S1=self.R_S1
        R_S2=self.R_S2
        chiVec1=self.chiVec1
        chiVec2=self.chiVec2
        chi1chi1=self.chi1chi1
        chi1chi2=self.chi1chi2
        chi2chi2=self.chi2chi2
        chi1_n=self.chi1_n
        chi1_lambda=self.chi1_lambda
        chi1_ell=self.chi1_ell
        chi2_n=self.chi2_n
        chi2_lambda=self.chi2_lambda
        chi2_ell=self.chi2_ell
        S_ell=self.S_ell
        S_n=self.S_n
        S_lambda=self.S_lambda
        Sigma_ell=self.Sigma_ell
        Sigma_n=self.Sigma_n
        Sigma_lambda=self.Sigma_lambda
        chi_s_ell=self.chi_s_ell
        chi_a_ell=self.chi_a_ell
        logv=self.logv
        Fcal_coeff=self.Fcal_coeff
        Fcal_0=self.Fcal_0
        Fcal_2=self.Fcal_2
        Fcal_3=self.Fcal_3
        Fcal_4=self.Fcal_4
        Fcal_5=self.Fcal_5
        Fcal_6=self.Fcal_6
        Fcal_lnv_6=self.Fcal_lnv_6
        Fcal_7=self.Fcal_7
        Fcal_8=self.Fcal_8
        Fcal_lnv_8=self.Fcal_lnv_8
        Fcal_9=self.Fcal_9
        Fcal_lnv_9=self.Fcal_lnv_9
        Fcal_10=self.Fcal_10
        Fcal_lnv_10=self.Fcal_lnv_10
        Fcal_SQ_4=self.Fcal_SQ_4
        Fcal_SO_3=self.Fcal_SO_3
        Fcal_SO_5=self.Fcal_SO_5
        Fcal_SO_6=self.Fcal_SO_6
        Fcal_SO_7=self.Fcal_SO_7
        Fcal_SO_8=self.Fcal_SO_8
        E_0=self.E_0
        E_2=self.E_2
        E_4=self.E_4
        E_6=self.E_6
        E_8=self.E_8
        E_lnv_8=self.E_lnv_8
        E_10=self.E_10
        E_lnv_10=self.E_lnv_10
        E_SQ_4=self.E_SQ_4
        E_SO_3=self.E_SO_3
        E_SO_5=self.E_SO_5
        E_SO_7=self.E_SO_7
        R = exp(rfrak_frame_x*xHat + rfrak_frame_y*yHat + rfrak_frame_z*zHat)
        nHat = R*xHat*conjugate(R)
        lambdaHat = R*yHat*conjugate(R)
        ellHat = R*zHat*conjugate(R)
        R_S1 = exp(rfrak_chi1_x*xHat + rfrak_chi1_y*yHat)
        R_S2 = exp(rfrak_chi2_x*xHat + rfrak_chi2_y*yHat)
        chiVec1 = S_chi1*R_S1*zHat*conjugate(R_S1)*conjugate(S_chi1)
        chiVec2 = S_chi2*R_S2*zHat*conjugate(R_S2)*conjugate(S_chi2)
        chi1_n = dot(chiVec1.vector,nHat.vector)
        chi1_lambda = dot(chiVec1.vector,lambdaHat.vector)
        chi1_ell = dot(chiVec1.vector,ellHat.vector)
        chi2_n = dot(chiVec2.vector,nHat.vector)
        chi2_lambda = dot(chiVec2.vector,lambdaHat.vector)
        chi2_ell = dot(chiVec2.vector,ellHat.vector)
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
        self.R=R
        self.nHat=nHat
        self.lambdaHat=lambdaHat
        self.ellHat=ellHat
        self.R_S1=R_S1
        self.R_S2=R_S2
        self.chiVec1=chiVec1
        self.chiVec2=chiVec2
        self.chi1_n=chi1_n
        self.chi1_lambda=chi1_lambda
        self.chi1_ell=chi1_ell
        self.chi2_n=chi2_n
        self.chi2_lambda=chi2_lambda
        self.chi2_ell=chi2_ell
        self.S_ell=S_ell
        self.S_n=S_n
        self.S_lambda=S_lambda
        self.Sigma_ell=Sigma_ell
        self.Sigma_n=Sigma_n
        self.Sigma_lambda=Sigma_lambda
        self.chi_s_ell=chi_s_ell
        self.chi_a_ell=chi_a_ell
        self.logv=logv
        self.Fcal_coeff=Fcal_coeff
        self.Fcal_SQ_4=Fcal_SQ_4
        self.Fcal_SO_3=Fcal_SO_3
        self.Fcal_SO_5=Fcal_SO_5
        self.Fcal_SO_6=Fcal_SO_6
        self.Fcal_SO_7=Fcal_SO_7
        self.Fcal_SO_8=Fcal_SO_8
        self.E_SQ_4=E_SQ_4
        self.E_SO_3=E_SO_3
        self.E_SO_5=E_SO_5
        self.E_SO_7=E_SO_7

    def OmegaVec_chiVec_1(self):
        xHat=self.xHat
        yHat=self.yHat
        zHat=self.zHat
        M1=self.M1
        M2=self.M2
        v=self.v
        S_chi1=self.S_chi1
        S_chi2=self.S_chi2
        rfrak_chi1_x=self.rfrak_chi1_x
        rfrak_chi1_y=self.rfrak_chi1_y
        rfrak_chi2_x=self.rfrak_chi2_x
        rfrak_chi2_y=self.rfrak_chi2_y
        rfrak_frame_x=self.rfrak_frame_x
        rfrak_frame_y=self.rfrak_frame_y
        rfrak_frame_z=self.rfrak_frame_z
        M=self.M
        delta=self.delta
        nu=self.nu
        R=self.R
        nHat=self.nHat
        lambdaHat=self.lambdaHat
        ellHat=self.ellHat
        R_S1=self.R_S1
        R_S2=self.R_S2
        chiVec1=self.chiVec1
        chiVec2=self.chiVec2
        chi1_n=self.chi1_n
        chi1_lambda=self.chi1_lambda
        chi1_ell=self.chi1_ell
        chi2_n=self.chi2_n
        chi2_lambda=self.chi2_lambda
        chi2_ell=self.chi2_ell
        S_ell=self.S_ell
        S_n=self.S_n
        S_lambda=self.S_lambda
        Sigma_ell=self.Sigma_ell
        Sigma_n=self.Sigma_n
        Sigma_lambda=self.Sigma_lambda
        Omega1_coeff = v**5/M
        return Omega1_coeff*(ellHat*(-0.75*delta + 0.5*nu + v**2*(delta*(0.625*nu - 0.5625) + nu*(1.25 - 0.0416666666666667*nu) + v**2*(delta*(nu*(4.875 - 0.15625*nu) - 0.84375) + nu*(nu*(-0.0208333333333333*nu - 3.28125) + 0.1875) + 0.84375) + 0.5625) + 0.75) + nHat*v*(3.0*chi1_n*nu + 3.0*M2**2*chi2_n/M**2) - M2**2*chiVec2*v/M**2)

    def OmegaVec_chiVec_2(self):
        xHat=self.xHat
        yHat=self.yHat
        zHat=self.zHat
        M1=self.M1
        M2=self.M2
        v=self.v
        S_chi1=self.S_chi1
        S_chi2=self.S_chi2
        rfrak_chi1_x=self.rfrak_chi1_x
        rfrak_chi1_y=self.rfrak_chi1_y
        rfrak_chi2_x=self.rfrak_chi2_x
        rfrak_chi2_y=self.rfrak_chi2_y
        rfrak_frame_x=self.rfrak_frame_x
        rfrak_frame_y=self.rfrak_frame_y
        rfrak_frame_z=self.rfrak_frame_z
        M=self.M
        delta=self.delta
        nu=self.nu
        R=self.R
        nHat=self.nHat
        lambdaHat=self.lambdaHat
        ellHat=self.ellHat
        R_S1=self.R_S1
        R_S2=self.R_S2
        chiVec1=self.chiVec1
        chiVec2=self.chiVec2
        chi1_n=self.chi1_n
        chi1_lambda=self.chi1_lambda
        chi1_ell=self.chi1_ell
        chi2_n=self.chi2_n
        chi2_lambda=self.chi2_lambda
        chi2_ell=self.chi2_ell
        S_ell=self.S_ell
        S_n=self.S_n
        S_lambda=self.S_lambda
        Sigma_ell=self.Sigma_ell
        Sigma_n=self.Sigma_n
        Sigma_lambda=self.Sigma_lambda
        Omega2_coeff = v**5/M
        return Omega2_coeff*(ellHat*(0.75*delta + 0.5*nu + v**2*(delta*(0.5625 - 0.625*nu) + nu*(1.25 - 0.0416666666666667*nu) + v**2*(delta*(nu*(0.15625*nu - 4.875) + 0.84375) + nu*(nu*(-0.0208333333333333*nu - 3.28125) + 0.1875) + 0.84375) + 0.5625) + 0.75) + nHat*v*(3.0*chi2_n*nu + 3.0*M1**2*chi1_n/M**2) - M1**2*chiVec1*v/M**2)

    def OmegaVec(self):
        xHat=self.xHat
        yHat=self.yHat
        zHat=self.zHat
        M1=self.M1
        M2=self.M2
        v=self.v
        S_chi1=self.S_chi1
        S_chi2=self.S_chi2
        rfrak_chi1_x=self.rfrak_chi1_x
        rfrak_chi1_y=self.rfrak_chi1_y
        rfrak_chi2_x=self.rfrak_chi2_x
        rfrak_chi2_y=self.rfrak_chi2_y
        rfrak_frame_x=self.rfrak_frame_x
        rfrak_frame_y=self.rfrak_frame_y
        rfrak_frame_z=self.rfrak_frame_z
        M=self.M
        delta=self.delta
        nu=self.nu
        R=self.R
        nHat=self.nHat
        lambdaHat=self.lambdaHat
        ellHat=self.ellHat
        R_S1=self.R_S1
        R_S2=self.R_S2
        chiVec1=self.chiVec1
        chiVec2=self.chiVec2
        chi1_n=self.chi1_n
        chi1_lambda=self.chi1_lambda
        chi1_ell=self.chi1_ell
        chi2_n=self.chi2_n
        chi2_lambda=self.chi2_lambda
        chi2_ell=self.chi2_ell
        S_ell=self.S_ell
        S_n=self.S_n
        S_lambda=self.S_lambda
        Sigma_ell=self.Sigma_ell
        Sigma_n=self.Sigma_n
        Sigma_lambda=self.Sigma_lambda
        gamma_PN_6 = 0.0123456790123457*nu**3 + 6.36111111111111*nu**2 - 2.98177812235564*nu + 1.0
        gamma_PN_4 = 1.0 - 5.41666666666667*nu
        gamma_PN_5 = (S_ell*(0.888888888888889*nu + 3.33333333333333) + 2.0*Sigma_ell*delta)/M**2
        gamma_PN_2 = 1.0 - 0.333333333333333*nu
        a_ell_4 = S_n*(5.77777777777778*nu**2 + 14.75*nu + 1.5) + Sigma_n*delta*(2.83333333333333*nu**2 + 9.125*nu + 1.5)
        gamma_PN_0 = 1.00000000000000
        gamma_PN_7 = (S_ell*(-6.0*nu**2 - 10.5833333333333*nu + 5.0) - 2.66666666666667*Sigma_ell*delta*nu**2 + Sigma_ell*delta*(3.0 - 10.1666666666667*nu))/M**2
        gamma_PN_3 = (1.66666666666667*S_ell + Sigma_ell*delta)/M**2
        a_ell_0 = 7.0*S_n + 3.0*Sigma_n*delta
        a_ell_2 = S_n*(-9.66666666666667*nu - 10.0) + Sigma_n*delta*(-4.5*nu - 6.0)
        return ellHat*v**3/M + nHat*v**6*(a_ell_0 + v**2*(a_ell_2 + a_ell_4*v**2))*(gamma_PN_0 + v**2*(gamma_PN_2 + v*(gamma_PN_3 + v*(gamma_PN_4 + v*(gamma_PN_5 + v*(gamma_PN_6 + gamma_PN_7*v))))))/M**3

    def OrbitalAngularMomentum(self):
        xHat=self.xHat
        yHat=self.yHat
        zHat=self.zHat
        M1=self.M1
        M2=self.M2
        v=self.v
        S_chi1=self.S_chi1
        S_chi2=self.S_chi2
        rfrak_chi1_x=self.rfrak_chi1_x
        rfrak_chi1_y=self.rfrak_chi1_y
        rfrak_chi2_x=self.rfrak_chi2_x
        rfrak_chi2_y=self.rfrak_chi2_y
        rfrak_frame_x=self.rfrak_frame_x
        rfrak_frame_y=self.rfrak_frame_y
        rfrak_frame_z=self.rfrak_frame_z
        M=self.M
        delta=self.delta
        nu=self.nu
        R=self.R
        nHat=self.nHat
        lambdaHat=self.lambdaHat
        ellHat=self.ellHat
        R_S1=self.R_S1
        R_S2=self.R_S2
        chiVec1=self.chiVec1
        chiVec2=self.chiVec2
        chi1_n=self.chi1_n
        chi1_lambda=self.chi1_lambda
        chi1_ell=self.chi1_ell
        chi2_n=self.chi2_n
        chi2_lambda=self.chi2_lambda
        chi2_ell=self.chi2_ell
        S_ell=self.S_ell
        S_n=self.S_n
        S_lambda=self.S_lambda
        Sigma_ell=self.Sigma_ell
        Sigma_n=self.Sigma_n
        Sigma_lambda=self.Sigma_lambda
        L_lnv_8 = -42.6666666666667*ellHat*nu
        L_6 = ellHat*(0.00540123456790123*nu**3 + 1.29166666666667*nu**2 - 30.9797035925835*nu + 8.4375)
        L_lnv_10 = 2.0*ellHat*nu*(87.4666666666667*nu + 95.0095238095238)
        L_SO_7 = 0.0625*ellHat*(-29.0*S_ell*nu**2 + 1101.0*S_ell*nu - 405.0*S_ell - 15.0*Sigma_ell*delta*nu**2 + 468.0*Sigma_ell*delta*nu - 81.0*Sigma_ell*delta)/M**2 + 0.0416666666666667*lambdaHat*(-32.0*S_lambda*nu**2 + 2.0*S_lambda*nu - 174.0*S_lambda - 16.0*Sigma_lambda*delta*nu**2 - 79.0*Sigma_lambda*delta*nu - 12.0*Sigma_lambda*delta)/M**2 + 0.0208333333333333*nHat*(11.0*S_n*nu**2 - 1331.0*S_n*nu + 183.0*S_n + 5.0*Sigma_n*delta*nu**2 - 734.0*Sigma_n*delta*nu + 183.0*Sigma_n*delta)/M**2
        L_SO_5 = 0.0138888888888889*ellHat*(427.0*S_ell*nu - 693.0*S_ell + 210.0*Sigma_ell*delta*nu - 189.0*Sigma_ell*delta)/M**2 + 0.166666666666667*lambdaHat*(18.0*S_lambda*nu - 21.0*S_lambda + 8.0*Sigma_lambda*delta*nu - 3.0*Sigma_lambda*delta)/M**2 + 0.0416666666666667*nHat*(-19.0*S_n*nu + 33.0*S_n - 10.0*Sigma_n*delta*nu + 33.0*Sigma_n*delta)/M**2
        L_SO_3 = 0.166666666666667*ellHat*(-35.0*S_ell - 15.0*Sigma_ell*delta)/M**2 + lambdaHat*(-3.0*S_lambda - Sigma_lambda*delta)/M**2 + 0.5*nHat*(S_n + Sigma_n*delta)/M**2
        L_2 = ellHat*(0.166666666666667*nu + 1.5)
        L_10 = ellHat*(nu*(-4.85925925925926*nu - 5.27830687830688) + 59.80078125)
        L_4 = ellHat*(0.0416666666666667*nu**2 - 2.375*nu + 3.375)
        L_0 = ellHat
        L_8 = ellHat*(-0.714285714285714*nu*(0.0024755658436214*nu**3 + 0.174189814814815*nu**2 - 90.1327990262052*nu + 153.88379682994) + 1.82857142857143*nu + 22.1484375)
        L_coeff = M**2*nu/v
        return L_coeff*(L_0 + v**2*(L_2 + v*(L_SO_3 + v*(L_4 + v*(L_SO_5 + v*(L_6 + v*(L_SO_7 + v*(L_8 + L_lnv_8*math.log(v) + v**2*(L_10 + L_lnv_10*math.log(v))))))))))

               

    #Evolve PN
    def __TaylorT1__(self,t, y):
        self.__Recalculate__(t, y)
        xHat=self.xHat
        yHat=self.yHat
        zHat=self.zHat
        M1=self.M1
        M2=self.M2
        v=self.v
        S_chi1=self.S_chi1
        S_chi2=self.S_chi2
        rfrak_chi1_x=self.rfrak_chi1_x
        rfrak_chi1_y=self.rfrak_chi1_y
        rfrak_chi2_x=self.rfrak_chi2_x
        rfrak_chi2_y=self.rfrak_chi2_y
        rfrak_frame_x=self.rfrak_frame_x
        rfrak_frame_y=self.rfrak_frame_y
        rfrak_frame_z=self.rfrak_frame_z
        M=self.M
        delta=self.delta
        nu=self.nu
        R=self.R
        nHat=self.nHat
        lambdaHat=self.lambdaHat
        ellHat=self.ellHat
        R_S1=self.R_S1
        R_S2=self.R_S2
        chiVec1=self.chiVec1
        chiVec2=self.chiVec2
        chi1chi1=self.chi1chi1
        chi1chi2=self.chi1chi2
        chi2chi2=self.chi2chi2
        chi1_n=self.chi1_n
        chi1_lambda=self.chi1_lambda
        chi1_ell=self.chi1_ell
        chi2_n=self.chi2_n
        chi2_lambda=self.chi2_lambda
        chi2_ell=self.chi2_ell
        S_ell=self.S_ell
        S_n=self.S_n
        S_lambda=self.S_lambda
        Sigma_ell=self.Sigma_ell
        Sigma_n=self.Sigma_n
        Sigma_lambda=self.Sigma_lambda
        chi_s_ell=self.chi_s_ell
        chi_a_ell=self.chi_a_ell
        logv=self.logv
        Fcal_coeff=self.Fcal_coeff
        Fcal_0=self.Fcal_0
        Fcal_2=self.Fcal_2
        Fcal_3=self.Fcal_3
        Fcal_4=self.Fcal_4
        Fcal_5=self.Fcal_5
        Fcal_6=self.Fcal_6
        Fcal_lnv_6=self.Fcal_lnv_6
        Fcal_7=self.Fcal_7
        Fcal_8=self.Fcal_8
        Fcal_lnv_8=self.Fcal_lnv_8
        Fcal_9=self.Fcal_9
        Fcal_lnv_9=self.Fcal_lnv_9
        Fcal_10=self.Fcal_10
        Fcal_lnv_10=self.Fcal_lnv_10
        Fcal_SQ_4=self.Fcal_SQ_4
        Fcal_SO_3=self.Fcal_SO_3
        Fcal_SO_5=self.Fcal_SO_5
        Fcal_SO_6=self.Fcal_SO_6
        Fcal_SO_7=self.Fcal_SO_7
        Fcal_SO_8=self.Fcal_SO_8
        E_0=self.E_0
        E_2=self.E_2
        E_4=self.E_4
        E_6=self.E_6
        E_8=self.E_8
        E_lnv_8=self.E_lnv_8
        E_10=self.E_10
        E_lnv_10=self.E_lnv_10
        E_SQ_4=self.E_SQ_4
        E_SO_3=self.E_SO_3
        E_SO_5=self.E_SO_5
        E_SO_7=self.E_SO_7
        v, rfrak_chi1_x, rfrak_chi1_y, rfrak_chi2_x, rfrak_chi2_y=self.v, self.rfrak_chi1_x, self.rfrak_chi1_y, self.rfrak_chi2_x, self.rfrak_chi2_y
        rfrak_frame_x, rfrak_frame_y, rfrak_frame_z, Phi=self.rfrak_frame_x, self.rfrak_frame_y, self.rfrak_frame_z, self.Phi
        dvdt_T1=0
        if v>=1.0 and self.NotForward:
            print("Beyond domain of PN validity, this is a good way to terminate.")
            self.terminal1=False
        Flux = Fcal_coeff*(Fcal_0 + v**2*(Fcal_2 + v*(Fcal_3 + Fcal_SO_3 + v*(Fcal_4 + Fcal_SQ_4 + v*(Fcal_5 + Fcal_SO_5 + v*(Fcal_6 + Fcal_SO_6 + Fcal_lnv_6*logv + v*(Fcal_7 + Fcal_SO_7 + v*(Fcal_8 + Fcal_SO_8 + Fcal_lnv_8*logv + v*(Fcal_9 + Fcal_lnv_9*logv + v*(Fcal_10 + Fcal_lnv_10*logv))))))))))
        dEdv = -0.5*M*nu*v*(2.0*E_0 + v**2*(4.0*E_2 + v*(5.0*E_SO_3 + v*(6.0*E_4 + 6.0*E_SQ_4 + v*(7.0*E_SO_5 + v*(8.0*E_6 + v*(9.0*E_SO_7 + v*(10.0*E_8 + E_lnv_8*(10.0*logv + 1.0) + v**2*(12.0*E_10 + E_lnv_10*(12.0*logv + 1.0))))))))))
        Absorption = 0
        dvdt_T1 = (-Absorption - Flux)/dEdv
        if dvdt_T1<1.0e-12 and self.NotForward:
            print("v is decreasing, which is not an uncommon way to stop.")
            self.terminal2=False
        return self.CommonRHS(dvdt_T1, y)

    def CommonRHS(self,dvdt, y):
        dydt=np.zeros(9)
        rfrak_frame=np.zeros(3)
        rfrak_frame[0] = y[5]
        rfrak_frame[1] = y[6]
        rfrak_frame[2] = y[7]
        rfrakdot_frame = FrameFromAngularVelocityIntegrand(rfrak_frame, self.OmegaVec().vector)
        dydt[0] = dvdt
        if(self.EvolveSpin1):
            dydt[1], dydt[2]= FrameFromAngularVelocity_2D_Integrand(y[1], y[2],\
                (self.S_chi1.inverse*self.OmegaVec_chiVec_1()*self.S_chi1).vector)
        else:
            dydt[1] = 0.0
            dydt[2] = 0.0
        if(self.EvolveSpin2):
            dydt[3], dydt[4] = FrameFromAngularVelocity_2D_Integrand(y[3], y[4],\
                (self.S_chi2.inverse*self.OmegaVec_chiVec_2()*self.S_chi2).vector)
        else:
            dydt[3] = 0.0
            dydt[4] = 0.0
        dydt[5] = rfrakdot_frame[0]
        dydt[6] = rfrakdot_frame[1]
        dydt[7] = rfrakdot_frame[2]
        dydt[8] = self.v*self.v*self.v/self.M

        return dydt

    # Prepare for integration
    def Evolution(self): 
        v_i,rfrak_chi1_x_i,rfrak_chi1_y_i,rfrak_chi2_x_i,rfrak_chi2_y_i,rfrak_frame_x_i,rfrak_frame_y_i,rfrak_frame_z_i\
            =self.v_i,self.rfrak_chi1_x_i,self.rfrak_chi1_y_i,self.rfrak_chi2_x_i,self.rfrak_chi2_y_i,self.rfrak_frame_x_i,self.rfrak_frame_y_i,self.rfrak_frame_z_i
        def terminate(t,y):
            return 1.0*self.terminal1*self.terminal2
        terminate.terminal=True
        TMerger=5.0/(256.0*self.nu*v_i**8)
        TEnd=TMerger
        if self.t_PN2:
            TEnd=self.t_PN2
        time=[0.0]
        while time[-1]<TEnd and 2*self.M*(256*self.nu*(TMerger-time[-1])/5)**(3/8)/self.StepsPerOrbit>self.MinStep:
            time.append(time[-1]+2*self.M*(256*self.nu*(TMerger-time[-1])/5)**(3/8)/self.StepsPerOrbit)
        time=np.delete(time, -1)
        
        # Integrate
        yy=solve_ivp(self.__TaylorT1__, [time[0],time[-1]], [v_i,rfrak_chi1_x_i,\
            rfrak_chi1_y_i,rfrak_chi2_x_i,rfrak_chi2_y_i,rfrak_frame_x_i,\
            rfrak_frame_y_i,rfrak_frame_z_i,0.0], method='DOP853',\
            t_eval=time, dense_output=True, events=terminate, rtol=self.tol, atol=self.tol)           
        if self.ForwardInTime:
            self.NotForward=False
            time=[0.0]
            TStart=-3*TMerger
            if self.t_PN1:
                TStart=self.t_PN1
            while time[-1]>TStart:
                time.append(time[-1]-2*self.M*(256*self.nu*(TMerger-time[-1])/5)**(3/8)/self.StepsPerOrbit)
            yyForward=solve_ivp(self.__TaylorT1__, [time[0],time[-1]], [v_i,rfrak_chi1_x_i,\
                rfrak_chi1_y_i,rfrak_chi2_x_i,rfrak_chi2_y_i,rfrak_frame_x_i,\
                rfrak_frame_y_i,rfrak_frame_z_i,0.0], method='DOP853',\
                t_eval=time, dense_output=True, rtol=self.tol, atol=self.tol)
            yy.t=np.append(yyForward.t[1:][::-1],yy.t)
            data=np.empty((9,len(yy.t)))
            for i in range(9):
                data[i]=np.append(yyForward.y[i][1:][::-1],yy.y[i])
            yy.y=data
            
        return yy

spec=[("xHat", float64),("yHat", float64),("zHat", float64),("M1", float64),("M2", float64),("v", float64),("S_chi1", float64),("S_chi2", float64),("rfrak_chi1_x", float64),("rfrak_chi1_y", float64),("rfrak_chi2_x", float64),("rfrak_chi2_y", float64),("rfrak_frame_x", float64),("rfrak_frame_y", float64),("rfrak_frame_z", float64),("M", float64),("delta", float64),("nu", float64),("R", float64),("nHat", float64),("lambdaHat", float64),("ellHat", float64),("R_S1", float64),("R_S2", float64),("chiVec1", float64),("chiVec2", float64),("chi1chi1", float64),("chi1chi2", float64),("chi2chi2", float64),("chi1_n", float64),("chi1_lambda", float64),("chi1_ell", float64),("chi2_n", float64),("chi2_lambda", float64),("chi2_ell", float64),("S_ell", float64),("S_n", float64),("S_lambda", float64),("Sigma_ell", float64),("Sigma_n", float64),("Sigma_lambda", float64),("chi_s_ell", float64),("chi_a_ell", float64),("logv", float64),("Fcal_coeff", float64),("Fcal_0", float64),("Fcal_2", float64),("Fcal_3", float64),("Fcal_4", float64),("Fcal_5", float64),("Fcal_6", float64),("Fcal_lnv_6", float64),("Fcal_7", float64),("Fcal_8", float64),("Fcal_lnv_8", float64),("Fcal_9", float64),("Fcal_lnv_9", float64),("Fcal_10", float64),("Fcal_lnv_10", float64),("Fcal_11", float64),("Fcal_lnv_11", float64),("Fcal_SQ_4", float64),("Fcal_SO_3", float64),("Fcal_SO_5", float64),("Fcal_SO_6", float64),("Fcal_SO_7", float64),("Fcal_SO_8", float64),("E_0", float64),("E_2", float64),("E_4", float64),("E_6", float64),("E_8", float64),("E_lnv_8", float64),("E_10", float64),("E_lnv_10", float64),("E_11", float64),("E_SQ_4", float64),("E_SO_3", float64),("E_SO_5", float64),("E_SO_7", float64),("xHat_i", float64),("yHat_i", float64),("zHat_i", float64),("M1_i", float64),("M2_i", float64),("v_i", float64),("S_chi1_i", float64),("S_chi2_i", float64),("rfrak_chi1_x_i", float64),("rfrak_chi1_y_i", float64),("rfrak_chi2_x_i", float64),("rfrak_chi2_y_i", float64),("rfrak_frame_x_i", float64),("rfrak_frame_y_i", float64),("rfrak_frame_z_i", float64),("t_PN1",float64),("t_PN2",float64),("StepsPerOrbit",int32),("ForwardInTime",boolean),("tol",float64),("MinStep",float64),("EvolveSpin1",boolean),("EvolveSpin2",boolean),("terminal1",boolean),("terminal2",boolean),("NotForward",boolean)]
@jitclass(spec)
class TaylorTn_5p5PN_Q : 
    def __init__(self, xHat_i, yHat_i, zHat_i, M1_i, M2_i, v_i, S_chi1_i, S_chi2_i, rfrak_chi1_x_i, rfrak_chi1_y_i, rfrak_chi2_x_i, rfrak_chi2_y_i, rfrak_frame_x_i, rfrak_frame_y_i, rfrak_frame_z_i, t_PN1=False, t_PN2=False, StepsPerOrbit=32, ForwardInTime=True, tol=1e-12, MinStep=1e-7) :
        self.t_PN1,self.t_PN2,self.StepsPerOrbit,self.ForwardInTime,self.tol,self.MinStep=t_PN1, t_PN2, StepsPerOrbit, ForwardInTime, tol, MinStep
        terminal1=True
        terminal2=True
        NotForward=True
        self.xHat_i=xHat_i
        self.yHat_i=yHat_i
        self.zHat_i=zHat_i
        self.M1_i=M1_i
        self.M2_i=M2_i
        self.v_i=v_i
        self.S_chi1_i=S_chi1_i
        self.S_chi2_i=S_chi2_i
        self.rfrak_chi1_x_i=rfrak_chi1_x_i
        self.rfrak_chi1_y_i=rfrak_chi1_y_i
        self.rfrak_chi2_x_i=rfrak_chi2_x_i
        self.rfrak_chi2_y_i=rfrak_chi2_y_i
        self.rfrak_frame_x_i=rfrak_frame_x_i
        self.rfrak_frame_y_i=rfrak_frame_y_i
        self.rfrak_frame_z_i=rfrak_frame_z_i
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
        chi1chi1=dot(chiVec1.vector,chiVec1.vector)
        chi1chi2=dot(chiVec1.vector,chiVec2.vector)
        chi2chi2=dot(chiVec2.vector,chiVec2.vector)
        chi1_n=dot(chiVec1.vector,nHat.vector)
        chi1_lambda=dot(chiVec1.vector,lambdaHat.vector)
        chi1_ell=dot(chiVec1.vector,ellHat.vector)
        chi2_n=dot(chiVec2.vector,nHat.vector)
        chi2_lambda=dot(chiVec2.vector,lambdaHat.vector)
        chi2_ell=dot(chiVec2.vector,ellHat.vector)
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
        self.xHat=xHat
        self.yHat=yHat
        self.zHat=zHat
        self.M1=M1
        self.M2=M2
        self.v=v
        self.S_chi1=S_chi1
        self.S_chi2=S_chi2
        self.rfrak_chi1_x=rfrak_chi1_x
        self.rfrak_chi1_y=rfrak_chi1_y
        self.rfrak_chi2_x=rfrak_chi2_x
        self.rfrak_chi2_y=rfrak_chi2_y
        self.rfrak_frame_x=rfrak_frame_x
        self.rfrak_frame_y=rfrak_frame_y
        self.rfrak_frame_z=rfrak_frame_z
        self.M=M
        self.delta=delta
        self.nu=nu
        self.R=R
        self.nHat=nHat
        self.lambdaHat=lambdaHat
        self.ellHat=ellHat
        self.R_S1=R_S1
        self.R_S2=R_S2
        self.chiVec1=chiVec1
        self.chiVec2=chiVec2
        self.chi1chi1=chi1chi1
        self.chi1chi2=chi1chi2
        self.chi2chi2=chi2chi2
        self.chi1_n=chi1_n
        self.chi1_lambda=chi1_lambda
        self.chi1_ell=chi1_ell
        self.chi2_n=chi2_n
        self.chi2_lambda=chi2_lambda
        self.chi2_ell=chi2_ell
        self.S_ell=S_ell
        self.S_n=S_n
        self.S_lambda=S_lambda
        self.Sigma_ell=Sigma_ell
        self.Sigma_n=Sigma_n
        self.Sigma_lambda=Sigma_lambda
        self.chi_s_ell=chi_s_ell
        self.chi_a_ell=chi_a_ell
        self.logv=logv
        self.Fcal_coeff=Fcal_coeff
        self.Fcal_0=Fcal_0
        self.Fcal_2=Fcal_2
        self.Fcal_3=Fcal_3
        self.Fcal_4=Fcal_4
        self.Fcal_5=Fcal_5
        self.Fcal_6=Fcal_6
        self.Fcal_lnv_6=Fcal_lnv_6
        self.Fcal_7=Fcal_7
        self.Fcal_8=Fcal_8
        self.Fcal_lnv_8=Fcal_lnv_8
        self.Fcal_9=Fcal_9
        self.Fcal_lnv_9=Fcal_lnv_9
        self.Fcal_10=Fcal_10
        self.Fcal_lnv_10=Fcal_lnv_10
        self.Fcal_11=Fcal_11
        self.Fcal_lnv_11=Fcal_lnv_11
        self.Fcal_SQ_4=Fcal_SQ_4
        self.Fcal_SO_3=Fcal_SO_3
        self.Fcal_SO_5=Fcal_SO_5
        self.Fcal_SO_6=Fcal_SO_6
        self.Fcal_SO_7=Fcal_SO_7
        self.Fcal_SO_8=Fcal_SO_8
        self.E_0=E_0
        self.E_2=E_2
        self.E_4=E_4
        self.E_6=E_6
        self.E_8=E_8
        self.E_lnv_8=E_lnv_8
        self.E_10=E_10
        self.E_lnv_10=E_lnv_10
        self.E_11=E_11
        self.E_SQ_4=E_SQ_4
        self.E_SO_3=E_SO_3
        self.E_SO_5=E_SO_5
        self.E_SO_7=E_SO_7
        Phi=0.0
        EvolveSpin1=(S_chi1*conjugate(S_chi1)).norm>1e-12
        EvolveSpin2=(S_chi2*conjugate(S_chi2)).norm>1e-12
        self.v, self.rfrak_chi1_x, self.rfrak_chi1_y, self.rfrak_chi2_x, self.rfrak_chi2_y=v, rfrak_chi1_x, rfrak_chi1_y, rfrak_chi2_x, rfrak_chi2_y
        self.rfrak_frame_x, self.rfrak_frame_y, self.rfrak_frame_z, self.Phi=rfrak_frame_x, rfrak_frame_y, rfrak_frame_z, Phi
        self.EvolveSpin1=EvolveSpin1
        self.EvolveSpin2=EvolveSpin2
        self.terminal1,self.terminal2,self.NotForward=terminal1,terminal2,NotForward

    def __Recalculate__(self,t, y):
        v = y[0]
        rfrak_chi1_x = y[1]
        rfrak_chi1_y = y[2]
        rfrak_chi2_x = y[3]
        rfrak_chi2_y = y[4]
        rfrak_frame_x = y[5]
        rfrak_frame_y = y[6]
        rfrak_frame_z = y[7]
        Phi = y[8]
        self.v, self.rfrak_chi1_x, self.rfrak_chi1_y, self.rfrak_chi2_x, self.rfrak_chi2_y=v, rfrak_chi1_x, rfrak_chi1_y, rfrak_chi2_x, rfrak_chi2_y
        self.rfrak_frame_x, self.rfrak_frame_y, self.rfrak_frame_z, self.Phi=rfrak_frame_x, rfrak_frame_y, rfrak_frame_z, Phi
        xHat=self.xHat
        yHat=self.yHat
        zHat=self.zHat
        M1=self.M1
        M2=self.M2
        v=self.v
        S_chi1=self.S_chi1
        S_chi2=self.S_chi2
        rfrak_chi1_x=self.rfrak_chi1_x
        rfrak_chi1_y=self.rfrak_chi1_y
        rfrak_chi2_x=self.rfrak_chi2_x
        rfrak_chi2_y=self.rfrak_chi2_y
        rfrak_frame_x=self.rfrak_frame_x
        rfrak_frame_y=self.rfrak_frame_y
        rfrak_frame_z=self.rfrak_frame_z
        M=self.M
        delta=self.delta
        nu=self.nu
        R=self.R
        nHat=self.nHat
        lambdaHat=self.lambdaHat
        ellHat=self.ellHat
        R_S1=self.R_S1
        R_S2=self.R_S2
        chiVec1=self.chiVec1
        chiVec2=self.chiVec2
        chi1chi1=self.chi1chi1
        chi1chi2=self.chi1chi2
        chi2chi2=self.chi2chi2
        chi1_n=self.chi1_n
        chi1_lambda=self.chi1_lambda
        chi1_ell=self.chi1_ell
        chi2_n=self.chi2_n
        chi2_lambda=self.chi2_lambda
        chi2_ell=self.chi2_ell
        S_ell=self.S_ell
        S_n=self.S_n
        S_lambda=self.S_lambda
        Sigma_ell=self.Sigma_ell
        Sigma_n=self.Sigma_n
        Sigma_lambda=self.Sigma_lambda
        chi_s_ell=self.chi_s_ell
        chi_a_ell=self.chi_a_ell
        logv=self.logv
        Fcal_coeff=self.Fcal_coeff
        Fcal_0=self.Fcal_0
        Fcal_2=self.Fcal_2
        Fcal_3=self.Fcal_3
        Fcal_4=self.Fcal_4
        Fcal_5=self.Fcal_5
        Fcal_6=self.Fcal_6
        Fcal_lnv_6=self.Fcal_lnv_6
        Fcal_7=self.Fcal_7
        Fcal_8=self.Fcal_8
        Fcal_lnv_8=self.Fcal_lnv_8
        Fcal_9=self.Fcal_9
        Fcal_lnv_9=self.Fcal_lnv_9
        Fcal_10=self.Fcal_10
        Fcal_lnv_10=self.Fcal_lnv_10
        Fcal_11=self.Fcal_11
        Fcal_lnv_11=self.Fcal_lnv_11
        Fcal_SQ_4=self.Fcal_SQ_4
        Fcal_SO_3=self.Fcal_SO_3
        Fcal_SO_5=self.Fcal_SO_5
        Fcal_SO_6=self.Fcal_SO_6
        Fcal_SO_7=self.Fcal_SO_7
        Fcal_SO_8=self.Fcal_SO_8
        E_0=self.E_0
        E_2=self.E_2
        E_4=self.E_4
        E_6=self.E_6
        E_8=self.E_8
        E_lnv_8=self.E_lnv_8
        E_10=self.E_10
        E_lnv_10=self.E_lnv_10
        E_11=self.E_11
        E_SQ_4=self.E_SQ_4
        E_SO_3=self.E_SO_3
        E_SO_5=self.E_SO_5
        E_SO_7=self.E_SO_7
        R = exp(rfrak_frame_x*xHat + rfrak_frame_y*yHat + rfrak_frame_z*zHat)
        nHat = R*xHat*conjugate(R)
        lambdaHat = R*yHat*conjugate(R)
        ellHat = R*zHat*conjugate(R)
        R_S1 = exp(rfrak_chi1_x*xHat + rfrak_chi1_y*yHat)
        R_S2 = exp(rfrak_chi2_x*xHat + rfrak_chi2_y*yHat)
        chiVec1 = S_chi1*R_S1*zHat*conjugate(R_S1)*conjugate(S_chi1)
        chiVec2 = S_chi2*R_S2*zHat*conjugate(R_S2)*conjugate(S_chi2)
        chi1_n = dot(chiVec1.vector,nHat.vector)
        chi1_lambda = dot(chiVec1.vector,lambdaHat.vector)
        chi1_ell = dot(chiVec1.vector,ellHat.vector)
        chi2_n = dot(chiVec2.vector,nHat.vector)
        chi2_lambda = dot(chiVec2.vector,lambdaHat.vector)
        chi2_ell = dot(chiVec2.vector,ellHat.vector)
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
        self.R=R
        self.nHat=nHat
        self.lambdaHat=lambdaHat
        self.ellHat=ellHat
        self.R_S1=R_S1
        self.R_S2=R_S2
        self.chiVec1=chiVec1
        self.chiVec2=chiVec2
        self.chi1_n=chi1_n
        self.chi1_lambda=chi1_lambda
        self.chi1_ell=chi1_ell
        self.chi2_n=chi2_n
        self.chi2_lambda=chi2_lambda
        self.chi2_ell=chi2_ell
        self.S_ell=S_ell
        self.S_n=S_n
        self.S_lambda=S_lambda
        self.Sigma_ell=Sigma_ell
        self.Sigma_n=Sigma_n
        self.Sigma_lambda=Sigma_lambda
        self.chi_s_ell=chi_s_ell
        self.chi_a_ell=chi_a_ell
        self.logv=logv
        self.Fcal_coeff=Fcal_coeff
        self.Fcal_SQ_4=Fcal_SQ_4
        self.Fcal_SO_3=Fcal_SO_3
        self.Fcal_SO_5=Fcal_SO_5
        self.Fcal_SO_6=Fcal_SO_6
        self.Fcal_SO_7=Fcal_SO_7
        self.Fcal_SO_8=Fcal_SO_8
        self.E_SQ_4=E_SQ_4
        self.E_SO_3=E_SO_3
        self.E_SO_5=E_SO_5
        self.E_SO_7=E_SO_7

    def OmegaVec_chiVec_1(self):
        xHat=self.xHat
        yHat=self.yHat
        zHat=self.zHat
        M1=self.M1
        M2=self.M2
        v=self.v
        S_chi1=self.S_chi1
        S_chi2=self.S_chi2
        rfrak_chi1_x=self.rfrak_chi1_x
        rfrak_chi1_y=self.rfrak_chi1_y
        rfrak_chi2_x=self.rfrak_chi2_x
        rfrak_chi2_y=self.rfrak_chi2_y
        rfrak_frame_x=self.rfrak_frame_x
        rfrak_frame_y=self.rfrak_frame_y
        rfrak_frame_z=self.rfrak_frame_z
        M=self.M
        delta=self.delta
        nu=self.nu
        R=self.R
        nHat=self.nHat
        lambdaHat=self.lambdaHat
        ellHat=self.ellHat
        R_S1=self.R_S1
        R_S2=self.R_S2
        chiVec1=self.chiVec1
        chiVec2=self.chiVec2
        chi1_n=self.chi1_n
        chi1_lambda=self.chi1_lambda
        chi1_ell=self.chi1_ell
        chi2_n=self.chi2_n
        chi2_lambda=self.chi2_lambda
        chi2_ell=self.chi2_ell
        S_ell=self.S_ell
        S_n=self.S_n
        S_lambda=self.S_lambda
        Sigma_ell=self.Sigma_ell
        Sigma_n=self.Sigma_n
        Sigma_lambda=self.Sigma_lambda
        Omega1_coeff = v**5/M
        return Omega1_coeff*(ellHat*(-0.75*delta + 0.5*nu + v**2*(delta*(0.625*nu - 0.5625) + nu*(1.25 - 0.0416666666666667*nu) + v**2*(delta*(nu*(4.875 - 0.15625*nu) - 0.84375) + nu*(nu*(-0.0208333333333333*nu - 3.28125) + 0.1875) + 0.84375) + 0.5625) + 0.75) + nHat*v*(3.0*chi1_n*nu + 3.0*M2**2*chi2_n/M**2) - M2**2*chiVec2*v/M**2)

    def OmegaVec_chiVec_2(self):
        xHat=self.xHat
        yHat=self.yHat
        zHat=self.zHat
        M1=self.M1
        M2=self.M2
        v=self.v
        S_chi1=self.S_chi1
        S_chi2=self.S_chi2
        rfrak_chi1_x=self.rfrak_chi1_x
        rfrak_chi1_y=self.rfrak_chi1_y
        rfrak_chi2_x=self.rfrak_chi2_x
        rfrak_chi2_y=self.rfrak_chi2_y
        rfrak_frame_x=self.rfrak_frame_x
        rfrak_frame_y=self.rfrak_frame_y
        rfrak_frame_z=self.rfrak_frame_z
        M=self.M
        delta=self.delta
        nu=self.nu
        R=self.R
        nHat=self.nHat
        lambdaHat=self.lambdaHat
        ellHat=self.ellHat
        R_S1=self.R_S1
        R_S2=self.R_S2
        chiVec1=self.chiVec1
        chiVec2=self.chiVec2
        chi1_n=self.chi1_n
        chi1_lambda=self.chi1_lambda
        chi1_ell=self.chi1_ell
        chi2_n=self.chi2_n
        chi2_lambda=self.chi2_lambda
        chi2_ell=self.chi2_ell
        S_ell=self.S_ell
        S_n=self.S_n
        S_lambda=self.S_lambda
        Sigma_ell=self.Sigma_ell
        Sigma_n=self.Sigma_n
        Sigma_lambda=self.Sigma_lambda
        Omega2_coeff = v**5/M
        return Omega2_coeff*(ellHat*(0.75*delta + 0.5*nu + v**2*(delta*(0.5625 - 0.625*nu) + nu*(1.25 - 0.0416666666666667*nu) + v**2*(delta*(nu*(0.15625*nu - 4.875) + 0.84375) + nu*(nu*(-0.0208333333333333*nu - 3.28125) + 0.1875) + 0.84375) + 0.5625) + 0.75) + nHat*v*(3.0*chi2_n*nu + 3.0*M1**2*chi1_n/M**2) - M1**2*chiVec1*v/M**2)

    def OmegaVec(self):
        xHat=self.xHat
        yHat=self.yHat
        zHat=self.zHat
        M1=self.M1
        M2=self.M2
        v=self.v
        S_chi1=self.S_chi1
        S_chi2=self.S_chi2
        rfrak_chi1_x=self.rfrak_chi1_x
        rfrak_chi1_y=self.rfrak_chi1_y
        rfrak_chi2_x=self.rfrak_chi2_x
        rfrak_chi2_y=self.rfrak_chi2_y
        rfrak_frame_x=self.rfrak_frame_x
        rfrak_frame_y=self.rfrak_frame_y
        rfrak_frame_z=self.rfrak_frame_z
        M=self.M
        delta=self.delta
        nu=self.nu
        R=self.R
        nHat=self.nHat
        lambdaHat=self.lambdaHat
        ellHat=self.ellHat
        R_S1=self.R_S1
        R_S2=self.R_S2
        chiVec1=self.chiVec1
        chiVec2=self.chiVec2
        chi1_n=self.chi1_n
        chi1_lambda=self.chi1_lambda
        chi1_ell=self.chi1_ell
        chi2_n=self.chi2_n
        chi2_lambda=self.chi2_lambda
        chi2_ell=self.chi2_ell
        S_ell=self.S_ell
        S_n=self.S_n
        S_lambda=self.S_lambda
        Sigma_ell=self.Sigma_ell
        Sigma_n=self.Sigma_n
        Sigma_lambda=self.Sigma_lambda
        gamma_PN_6 = 0.0123456790123457*nu**3 + 6.36111111111111*nu**2 - 2.98177812235564*nu + 1.0
        gamma_PN_4 = 1.0 - 5.41666666666667*nu
        gamma_PN_5 = (S_ell*(0.888888888888889*nu + 3.33333333333333) + 2.0*Sigma_ell*delta)/M**2
        gamma_PN_2 = 1.0 - 0.333333333333333*nu
        a_ell_4 = S_n*(5.77777777777778*nu**2 + 14.75*nu + 1.5) + Sigma_n*delta*(2.83333333333333*nu**2 + 9.125*nu + 1.5)
        gamma_PN_0 = 1.00000000000000
        gamma_PN_7 = (S_ell*(-6.0*nu**2 - 10.5833333333333*nu + 5.0) - 2.66666666666667*Sigma_ell*delta*nu**2 + Sigma_ell*delta*(3.0 - 10.1666666666667*nu))/M**2
        gamma_PN_3 = (1.66666666666667*S_ell + Sigma_ell*delta)/M**2
        a_ell_0 = 7.0*S_n + 3.0*Sigma_n*delta
        a_ell_2 = S_n*(-9.66666666666667*nu - 10.0) + Sigma_n*delta*(-4.5*nu - 6.0)
        return ellHat*v**3/M + nHat*v**6*(a_ell_0 + v**2*(a_ell_2 + a_ell_4*v**2))*(gamma_PN_0 + v**2*(gamma_PN_2 + v*(gamma_PN_3 + v*(gamma_PN_4 + v*(gamma_PN_5 + v*(gamma_PN_6 + gamma_PN_7*v))))))/M**3

    def OrbitalAngularMomentum(self):
        xHat=self.xHat
        yHat=self.yHat
        zHat=self.zHat
        M1=self.M1
        M2=self.M2
        v=self.v
        S_chi1=self.S_chi1
        S_chi2=self.S_chi2
        rfrak_chi1_x=self.rfrak_chi1_x
        rfrak_chi1_y=self.rfrak_chi1_y
        rfrak_chi2_x=self.rfrak_chi2_x
        rfrak_chi2_y=self.rfrak_chi2_y
        rfrak_frame_x=self.rfrak_frame_x
        rfrak_frame_y=self.rfrak_frame_y
        rfrak_frame_z=self.rfrak_frame_z
        M=self.M
        delta=self.delta
        nu=self.nu
        R=self.R
        nHat=self.nHat
        lambdaHat=self.lambdaHat
        ellHat=self.ellHat
        R_S1=self.R_S1
        R_S2=self.R_S2
        chiVec1=self.chiVec1
        chiVec2=self.chiVec2
        chi1_n=self.chi1_n
        chi1_lambda=self.chi1_lambda
        chi1_ell=self.chi1_ell
        chi2_n=self.chi2_n
        chi2_lambda=self.chi2_lambda
        chi2_ell=self.chi2_ell
        S_ell=self.S_ell
        S_n=self.S_n
        S_lambda=self.S_lambda
        Sigma_ell=self.Sigma_ell
        Sigma_n=self.Sigma_n
        Sigma_lambda=self.Sigma_lambda
        L_lnv_8 = -42.6666666666667*ellHat*nu
        L_6 = ellHat*(0.00540123456790123*nu**3 + 1.29166666666667*nu**2 - 30.9797035925835*nu + 8.4375)
        L_lnv_10 = 2.0*ellHat*nu*(87.4666666666667*nu + 95.0095238095238)
        L_SO_7 = 0.0625*ellHat*(-29.0*S_ell*nu**2 + 1101.0*S_ell*nu - 405.0*S_ell - 15.0*Sigma_ell*delta*nu**2 + 468.0*Sigma_ell*delta*nu - 81.0*Sigma_ell*delta)/M**2 + 0.0416666666666667*lambdaHat*(-32.0*S_lambda*nu**2 + 2.0*S_lambda*nu - 174.0*S_lambda - 16.0*Sigma_lambda*delta*nu**2 - 79.0*Sigma_lambda*delta*nu - 12.0*Sigma_lambda*delta)/M**2 + 0.0208333333333333*nHat*(11.0*S_n*nu**2 - 1331.0*S_n*nu + 183.0*S_n + 5.0*Sigma_n*delta*nu**2 - 734.0*Sigma_n*delta*nu + 183.0*Sigma_n*delta)/M**2
        L_SO_5 = 0.0138888888888889*ellHat*(427.0*S_ell*nu - 693.0*S_ell + 210.0*Sigma_ell*delta*nu - 189.0*Sigma_ell*delta)/M**2 + 0.166666666666667*lambdaHat*(18.0*S_lambda*nu - 21.0*S_lambda + 8.0*Sigma_lambda*delta*nu - 3.0*Sigma_lambda*delta)/M**2 + 0.0416666666666667*nHat*(-19.0*S_n*nu + 33.0*S_n - 10.0*Sigma_n*delta*nu + 33.0*Sigma_n*delta)/M**2
        L_SO_3 = 0.166666666666667*ellHat*(-35.0*S_ell - 15.0*Sigma_ell*delta)/M**2 + lambdaHat*(-3.0*S_lambda - Sigma_lambda*delta)/M**2 + 0.5*nHat*(S_n + Sigma_n*delta)/M**2
        L_2 = ellHat*(0.166666666666667*nu + 1.5)
        L_10 = ellHat*(nu*(-4.85925925925926*nu - 5.27830687830688) + 59.80078125)
        L_4 = ellHat*(0.0416666666666667*nu**2 - 2.375*nu + 3.375)
        L_0 = ellHat
        L_8 = ellHat*(-0.714285714285714*nu*(0.0024755658436214*nu**3 + 0.174189814814815*nu**2 - 90.1327990262052*nu + 153.88379682994) + 1.82857142857143*nu + 22.1484375)
        L_coeff = M**2*nu/v
        return L_coeff*(L_0 + v**2*(L_2 + v*(L_SO_3 + v*(L_4 + v*(L_SO_5 + v*(L_6 + v*(L_SO_7 + v*(L_8 + L_lnv_8*math.log(v) + v**2*(L_10 + L_lnv_10*math.log(v))))))))))

               

    #Evolve PN
    def __TaylorT1__(self,t, y):
        self.__Recalculate__(t, y)
        xHat=self.xHat
        yHat=self.yHat
        zHat=self.zHat
        M1=self.M1
        M2=self.M2
        v=self.v
        S_chi1=self.S_chi1
        S_chi2=self.S_chi2
        rfrak_chi1_x=self.rfrak_chi1_x
        rfrak_chi1_y=self.rfrak_chi1_y
        rfrak_chi2_x=self.rfrak_chi2_x
        rfrak_chi2_y=self.rfrak_chi2_y
        rfrak_frame_x=self.rfrak_frame_x
        rfrak_frame_y=self.rfrak_frame_y
        rfrak_frame_z=self.rfrak_frame_z
        M=self.M
        delta=self.delta
        nu=self.nu
        R=self.R
        nHat=self.nHat
        lambdaHat=self.lambdaHat
        ellHat=self.ellHat
        R_S1=self.R_S1
        R_S2=self.R_S2
        chiVec1=self.chiVec1
        chiVec2=self.chiVec2
        chi1chi1=self.chi1chi1
        chi1chi2=self.chi1chi2
        chi2chi2=self.chi2chi2
        chi1_n=self.chi1_n
        chi1_lambda=self.chi1_lambda
        chi1_ell=self.chi1_ell
        chi2_n=self.chi2_n
        chi2_lambda=self.chi2_lambda
        chi2_ell=self.chi2_ell
        S_ell=self.S_ell
        S_n=self.S_n
        S_lambda=self.S_lambda
        Sigma_ell=self.Sigma_ell
        Sigma_n=self.Sigma_n
        Sigma_lambda=self.Sigma_lambda
        chi_s_ell=self.chi_s_ell
        chi_a_ell=self.chi_a_ell
        logv=self.logv
        Fcal_coeff=self.Fcal_coeff
        Fcal_0=self.Fcal_0
        Fcal_2=self.Fcal_2
        Fcal_3=self.Fcal_3
        Fcal_4=self.Fcal_4
        Fcal_5=self.Fcal_5
        Fcal_6=self.Fcal_6
        Fcal_lnv_6=self.Fcal_lnv_6
        Fcal_7=self.Fcal_7
        Fcal_8=self.Fcal_8
        Fcal_lnv_8=self.Fcal_lnv_8
        Fcal_9=self.Fcal_9
        Fcal_lnv_9=self.Fcal_lnv_9
        Fcal_10=self.Fcal_10
        Fcal_lnv_10=self.Fcal_lnv_10
        Fcal_11=self.Fcal_11
        Fcal_lnv_11=self.Fcal_lnv_11
        Fcal_SQ_4=self.Fcal_SQ_4
        Fcal_SO_3=self.Fcal_SO_3
        Fcal_SO_5=self.Fcal_SO_5
        Fcal_SO_6=self.Fcal_SO_6
        Fcal_SO_7=self.Fcal_SO_7
        Fcal_SO_8=self.Fcal_SO_8
        E_0=self.E_0
        E_2=self.E_2
        E_4=self.E_4
        E_6=self.E_6
        E_8=self.E_8
        E_lnv_8=self.E_lnv_8
        E_10=self.E_10
        E_lnv_10=self.E_lnv_10
        E_11=self.E_11
        E_SQ_4=self.E_SQ_4
        E_SO_3=self.E_SO_3
        E_SO_5=self.E_SO_5
        E_SO_7=self.E_SO_7
        v, rfrak_chi1_x, rfrak_chi1_y, rfrak_chi2_x, rfrak_chi2_y=self.v, self.rfrak_chi1_x, self.rfrak_chi1_y, self.rfrak_chi2_x, self.rfrak_chi2_y
        rfrak_frame_x, rfrak_frame_y, rfrak_frame_z, Phi=self.rfrak_frame_x, self.rfrak_frame_y, self.rfrak_frame_z, self.Phi
        dvdt_T1=0
        if v>=1.0 and self.NotForward:
            print("Beyond domain of PN validity, this is a good way to terminate.")
            self.terminal1=False
        Flux = Fcal_coeff*(Fcal_0 + v**2*(Fcal_2 + v*(Fcal_3 + Fcal_SO_3 + v*(Fcal_4 + Fcal_SQ_4 + v*(Fcal_5 + Fcal_SO_5 + v*(Fcal_6 + Fcal_SO_6 + Fcal_lnv_6*logv + v*(Fcal_7 + Fcal_SO_7 + v*(Fcal_8 + Fcal_SO_8 + Fcal_lnv_8*logv + v*(Fcal_9 + Fcal_lnv_9*logv + v*(Fcal_10 + Fcal_lnv_10*logv + v*(Fcal_11 + Fcal_lnv_11*logv)))))))))))
        dEdv = -0.5*M*nu*v*(2.0*E_0 + v**2*(4.0*E_2 + v*(5.0*E_SO_3 + v*(6.0*E_4 + 6.0*E_SQ_4 + v*(7.0*E_SO_5 + v*(8.0*E_6 + v*(9.0*E_SO_7 + v*(10.0*E_8 + E_lnv_8*(10.0*logv + 1.0) + v**2*(12.0*E_10 + 13.0*E_11*v + E_lnv_10*(12.0*logv + 1.0))))))))))
        Absorption = 0
        dvdt_T1 = (-Absorption - Flux)/dEdv
        if dvdt_T1<1.0e-12 and self.NotForward:
            print("v is decreasing, which is not an uncommon way to stop.")
            self.terminal2=False
        return self.CommonRHS(dvdt_T1, y)

    def CommonRHS(self,dvdt, y):
        dydt=np.zeros(9)
        rfrak_frame=np.zeros(3)
        rfrak_frame[0] = y[5]
        rfrak_frame[1] = y[6]
        rfrak_frame[2] = y[7]
        rfrakdot_frame = FrameFromAngularVelocityIntegrand(rfrak_frame, self.OmegaVec().vector)
        dydt[0] = dvdt
        if(self.EvolveSpin1):
            dydt[1], dydt[2]= FrameFromAngularVelocity_2D_Integrand(y[1], y[2],\
                (self.S_chi1.inverse*self.OmegaVec_chiVec_1()*self.S_chi1).vector)
        else:
            dydt[1] = 0.0
            dydt[2] = 0.0
        if(self.EvolveSpin2):
            dydt[3], dydt[4] = FrameFromAngularVelocity_2D_Integrand(y[3], y[4],\
                (self.S_chi2.inverse*self.OmegaVec_chiVec_2()*self.S_chi2).vector)
        else:
            dydt[3] = 0.0
            dydt[4] = 0.0
        dydt[5] = rfrakdot_frame[0]
        dydt[6] = rfrakdot_frame[1]
        dydt[7] = rfrakdot_frame[2]
        dydt[8] = self.v*self.v*self.v/self.M

        return dydt

    # Prepare for integration
    def Evolution(self): 
        v_i,rfrak_chi1_x_i,rfrak_chi1_y_i,rfrak_chi2_x_i,rfrak_chi2_y_i,rfrak_frame_x_i,rfrak_frame_y_i,rfrak_frame_z_i\
            =self.v_i,self.rfrak_chi1_x_i,self.rfrak_chi1_y_i,self.rfrak_chi2_x_i,self.rfrak_chi2_y_i,self.rfrak_frame_x_i,self.rfrak_frame_y_i,self.rfrak_frame_z_i
        def terminate(t,y):
            return 1.0*self.terminal1*self.terminal2
        terminate.terminal=True
        TMerger=5.0/(256.0*self.nu*v_i**8)
        TEnd=TMerger
        if self.t_PN2:
            TEnd=self.t_PN2
        time=[0.0]
        while time[-1]<TEnd and 2*self.M*(256*self.nu*(TMerger-time[-1])/5)**(3/8)/self.StepsPerOrbit>self.MinStep:
            time.append(time[-1]+2*self.M*(256*self.nu*(TMerger-time[-1])/5)**(3/8)/self.StepsPerOrbit)
        time=np.delete(time, -1)
        
        # Integrate
        yy=solve_ivp(self.__TaylorT1__, [time[0],time[-1]], [v_i,rfrak_chi1_x_i,\
            rfrak_chi1_y_i,rfrak_chi2_x_i,rfrak_chi2_y_i,rfrak_frame_x_i,\
            rfrak_frame_y_i,rfrak_frame_z_i,0.0], method='DOP853',\
            t_eval=time, dense_output=True, events=terminate, rtol=self.tol, atol=self.tol)           
        if self.ForwardInTime:
            self.NotForward=False
            time=[0.0]
            TStart=-3*TMerger
            if self.t_PN1:
                TStart=self.t_PN1
            while time[-1]>TStart:
                time.append(time[-1]-2*self.M*(256*self.nu*(TMerger-time[-1])/5)**(3/8)/self.StepsPerOrbit)
            yyForward=solve_ivp(self.__TaylorT1__, [time[0],time[-1]], [v_i,rfrak_chi1_x_i,\
                rfrak_chi1_y_i,rfrak_chi2_x_i,rfrak_chi2_y_i,rfrak_frame_x_i,\
                rfrak_frame_y_i,rfrak_frame_z_i,0.0], method='DOP853',\
                t_eval=time, dense_output=True, rtol=self.tol, atol=self.tol)
            yy.t=np.append(yyForward.t[1:][::-1],yy.t)
            data=np.empty((9,len(yy.t)))
            for i in range(9):
                data[i]=np.append(yyForward.y[i][1:][::-1],yy.y[i])
            yy.y=data
            
        return yy

spec=[("xHat", float64),("yHat", float64),("zHat", float64),("M1", float64),("M2", float64),("v", float64),("S_chi1", float64),("S_chi2", float64),("rfrak_chi1_x", float64),("rfrak_chi1_y", float64),("rfrak_chi2_x", float64),("rfrak_chi2_y", float64),("rfrak_frame_x", float64),("rfrak_frame_y", float64),("rfrak_frame_z", float64),("M", float64),("delta", float64),("nu", float64),("R", float64),("nHat", float64),("lambdaHat", float64),("ellHat", float64),("R_S1", float64),("R_S2", float64),("chiVec1", float64),("chiVec2", float64),("chi1chi1", float64),("chi1chi2", float64),("chi2chi2", float64),("chi1_n", float64),("chi1_lambda", float64),("chi1_ell", float64),("chi2_n", float64),("chi2_lambda", float64),("chi2_ell", float64),("S_ell", float64),("S_n", float64),("S_lambda", float64),("Sigma_ell", float64),("Sigma_n", float64),("Sigma_lambda", float64),("chi_s_ell", float64),("chi_a_ell", float64),("logv", float64),("Fcal_coeff", float64),("Fcal_0", float64),("Fcal_2", float64),("Fcal_3", float64),("Fcal_4", float64),("Fcal_5", float64),("Fcal_6", float64),("Fcal_lnv_6", float64),("Fcal_7", float64),("Fcal_8", float64),("Fcal_lnv_8", float64),("Fcal_9", float64),("Fcal_lnv_9", float64),("Fcal_10", float64),("Fcal_lnv_10", float64),("Fcal_11", float64),("Fcal_lnv_11", float64),("Fcal_12", float64),("Fcal_lnv_12", float64),("Fcal_lnv2_12", float64),("Fcal_SQ_4", float64),("Fcal_SO_3", float64),("Fcal_SO_5", float64),("Fcal_SO_6", float64),("Fcal_SO_7", float64),("Fcal_SO_8", float64),("E_0", float64),("E_2", float64),("E_4", float64),("E_6", float64),("E_8", float64),("E_lnv_8", float64),("E_10", float64),("E_lnv_10", float64),("E_11", float64),("E_12", float64),("E_lnv_12", float64),("E_SQ_4", float64),("E_SO_3", float64),("E_SO_5", float64),("E_SO_7", float64),("xHat_i", float64),("yHat_i", float64),("zHat_i", float64),("M1_i", float64),("M2_i", float64),("v_i", float64),("S_chi1_i", float64),("S_chi2_i", float64),("rfrak_chi1_x_i", float64),("rfrak_chi1_y_i", float64),("rfrak_chi2_x_i", float64),("rfrak_chi2_y_i", float64),("rfrak_frame_x_i", float64),("rfrak_frame_y_i", float64),("rfrak_frame_z_i", float64),("t_PN1",float64),("t_PN2",float64),("StepsPerOrbit",int32),("ForwardInTime",boolean),("tol",float64),("MinStep",float64),("EvolveSpin1",boolean),("EvolveSpin2",boolean),("terminal1",boolean),("terminal2",boolean),("NotForward",boolean)]
@jitclass(spec)
class TaylorTn_6p0PN_Q : 
    def __init__(self, xHat_i, yHat_i, zHat_i, M1_i, M2_i, v_i, S_chi1_i, S_chi2_i, rfrak_chi1_x_i, rfrak_chi1_y_i, rfrak_chi2_x_i, rfrak_chi2_y_i, rfrak_frame_x_i, rfrak_frame_y_i, rfrak_frame_z_i, t_PN1=False, t_PN2=False, StepsPerOrbit=32, ForwardInTime=True, tol=1e-12, MinStep=1e-7) :
        self.t_PN1,self.t_PN2,self.StepsPerOrbit,self.ForwardInTime,self.tol,self.MinStep=t_PN1, t_PN2, StepsPerOrbit, ForwardInTime, tol, MinStep
        terminal1=True
        terminal2=True
        NotForward=True
        self.xHat_i=xHat_i
        self.yHat_i=yHat_i
        self.zHat_i=zHat_i
        self.M1_i=M1_i
        self.M2_i=M2_i
        self.v_i=v_i
        self.S_chi1_i=S_chi1_i
        self.S_chi2_i=S_chi2_i
        self.rfrak_chi1_x_i=rfrak_chi1_x_i
        self.rfrak_chi1_y_i=rfrak_chi1_y_i
        self.rfrak_chi2_x_i=rfrak_chi2_x_i
        self.rfrak_chi2_y_i=rfrak_chi2_y_i
        self.rfrak_frame_x_i=rfrak_frame_x_i
        self.rfrak_frame_y_i=rfrak_frame_y_i
        self.rfrak_frame_z_i=rfrak_frame_z_i
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
        chi1chi1=dot(chiVec1.vector,chiVec1.vector)
        chi1chi2=dot(chiVec1.vector,chiVec2.vector)
        chi2chi2=dot(chiVec2.vector,chiVec2.vector)
        chi1_n=dot(chiVec1.vector,nHat.vector)
        chi1_lambda=dot(chiVec1.vector,lambdaHat.vector)
        chi1_ell=dot(chiVec1.vector,ellHat.vector)
        chi2_n=dot(chiVec2.vector,nHat.vector)
        chi2_lambda=dot(chiVec2.vector,lambdaHat.vector)
        chi2_ell=dot(chiVec2.vector,ellHat.vector)
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
        self.xHat=xHat
        self.yHat=yHat
        self.zHat=zHat
        self.M1=M1
        self.M2=M2
        self.v=v
        self.S_chi1=S_chi1
        self.S_chi2=S_chi2
        self.rfrak_chi1_x=rfrak_chi1_x
        self.rfrak_chi1_y=rfrak_chi1_y
        self.rfrak_chi2_x=rfrak_chi2_x
        self.rfrak_chi2_y=rfrak_chi2_y
        self.rfrak_frame_x=rfrak_frame_x
        self.rfrak_frame_y=rfrak_frame_y
        self.rfrak_frame_z=rfrak_frame_z
        self.M=M
        self.delta=delta
        self.nu=nu
        self.R=R
        self.nHat=nHat
        self.lambdaHat=lambdaHat
        self.ellHat=ellHat
        self.R_S1=R_S1
        self.R_S2=R_S2
        self.chiVec1=chiVec1
        self.chiVec2=chiVec2
        self.chi1chi1=chi1chi1
        self.chi1chi2=chi1chi2
        self.chi2chi2=chi2chi2
        self.chi1_n=chi1_n
        self.chi1_lambda=chi1_lambda
        self.chi1_ell=chi1_ell
        self.chi2_n=chi2_n
        self.chi2_lambda=chi2_lambda
        self.chi2_ell=chi2_ell
        self.S_ell=S_ell
        self.S_n=S_n
        self.S_lambda=S_lambda
        self.Sigma_ell=Sigma_ell
        self.Sigma_n=Sigma_n
        self.Sigma_lambda=Sigma_lambda
        self.chi_s_ell=chi_s_ell
        self.chi_a_ell=chi_a_ell
        self.logv=logv
        self.Fcal_coeff=Fcal_coeff
        self.Fcal_0=Fcal_0
        self.Fcal_2=Fcal_2
        self.Fcal_3=Fcal_3
        self.Fcal_4=Fcal_4
        self.Fcal_5=Fcal_5
        self.Fcal_6=Fcal_6
        self.Fcal_lnv_6=Fcal_lnv_6
        self.Fcal_7=Fcal_7
        self.Fcal_8=Fcal_8
        self.Fcal_lnv_8=Fcal_lnv_8
        self.Fcal_9=Fcal_9
        self.Fcal_lnv_9=Fcal_lnv_9
        self.Fcal_10=Fcal_10
        self.Fcal_lnv_10=Fcal_lnv_10
        self.Fcal_11=Fcal_11
        self.Fcal_lnv_11=Fcal_lnv_11
        self.Fcal_12=Fcal_12
        self.Fcal_lnv_12=Fcal_lnv_12
        self.Fcal_lnv2_12=Fcal_lnv2_12
        self.Fcal_SQ_4=Fcal_SQ_4
        self.Fcal_SO_3=Fcal_SO_3
        self.Fcal_SO_5=Fcal_SO_5
        self.Fcal_SO_6=Fcal_SO_6
        self.Fcal_SO_7=Fcal_SO_7
        self.Fcal_SO_8=Fcal_SO_8
        self.E_0=E_0
        self.E_2=E_2
        self.E_4=E_4
        self.E_6=E_6
        self.E_8=E_8
        self.E_lnv_8=E_lnv_8
        self.E_10=E_10
        self.E_lnv_10=E_lnv_10
        self.E_11=E_11
        self.E_12=E_12
        self.E_lnv_12=E_lnv_12
        self.E_SQ_4=E_SQ_4
        self.E_SO_3=E_SO_3
        self.E_SO_5=E_SO_5
        self.E_SO_7=E_SO_7
        Phi=0.0
        EvolveSpin1=(S_chi1*conjugate(S_chi1)).norm>1e-12
        EvolveSpin2=(S_chi2*conjugate(S_chi2)).norm>1e-12
        self.v, self.rfrak_chi1_x, self.rfrak_chi1_y, self.rfrak_chi2_x, self.rfrak_chi2_y=v, rfrak_chi1_x, rfrak_chi1_y, rfrak_chi2_x, rfrak_chi2_y
        self.rfrak_frame_x, self.rfrak_frame_y, self.rfrak_frame_z, self.Phi=rfrak_frame_x, rfrak_frame_y, rfrak_frame_z, Phi
        self.EvolveSpin1=EvolveSpin1
        self.EvolveSpin2=EvolveSpin2
        self.terminal1,self.terminal2,self.NotForward=terminal1,terminal2,NotForward

    def __Recalculate__(self,t, y):
        v = y[0]
        rfrak_chi1_x = y[1]
        rfrak_chi1_y = y[2]
        rfrak_chi2_x = y[3]
        rfrak_chi2_y = y[4]
        rfrak_frame_x = y[5]
        rfrak_frame_y = y[6]
        rfrak_frame_z = y[7]
        Phi = y[8]
        self.v, self.rfrak_chi1_x, self.rfrak_chi1_y, self.rfrak_chi2_x, self.rfrak_chi2_y=v, rfrak_chi1_x, rfrak_chi1_y, rfrak_chi2_x, rfrak_chi2_y
        self.rfrak_frame_x, self.rfrak_frame_y, self.rfrak_frame_z, self.Phi=rfrak_frame_x, rfrak_frame_y, rfrak_frame_z, Phi
        xHat=self.xHat
        yHat=self.yHat
        zHat=self.zHat
        M1=self.M1
        M2=self.M2
        v=self.v
        S_chi1=self.S_chi1
        S_chi2=self.S_chi2
        rfrak_chi1_x=self.rfrak_chi1_x
        rfrak_chi1_y=self.rfrak_chi1_y
        rfrak_chi2_x=self.rfrak_chi2_x
        rfrak_chi2_y=self.rfrak_chi2_y
        rfrak_frame_x=self.rfrak_frame_x
        rfrak_frame_y=self.rfrak_frame_y
        rfrak_frame_z=self.rfrak_frame_z
        M=self.M
        delta=self.delta
        nu=self.nu
        R=self.R
        nHat=self.nHat
        lambdaHat=self.lambdaHat
        ellHat=self.ellHat
        R_S1=self.R_S1
        R_S2=self.R_S2
        chiVec1=self.chiVec1
        chiVec2=self.chiVec2
        chi1chi1=self.chi1chi1
        chi1chi2=self.chi1chi2
        chi2chi2=self.chi2chi2
        chi1_n=self.chi1_n
        chi1_lambda=self.chi1_lambda
        chi1_ell=self.chi1_ell
        chi2_n=self.chi2_n
        chi2_lambda=self.chi2_lambda
        chi2_ell=self.chi2_ell
        S_ell=self.S_ell
        S_n=self.S_n
        S_lambda=self.S_lambda
        Sigma_ell=self.Sigma_ell
        Sigma_n=self.Sigma_n
        Sigma_lambda=self.Sigma_lambda
        chi_s_ell=self.chi_s_ell
        chi_a_ell=self.chi_a_ell
        logv=self.logv
        Fcal_coeff=self.Fcal_coeff
        Fcal_0=self.Fcal_0
        Fcal_2=self.Fcal_2
        Fcal_3=self.Fcal_3
        Fcal_4=self.Fcal_4
        Fcal_5=self.Fcal_5
        Fcal_6=self.Fcal_6
        Fcal_lnv_6=self.Fcal_lnv_6
        Fcal_7=self.Fcal_7
        Fcal_8=self.Fcal_8
        Fcal_lnv_8=self.Fcal_lnv_8
        Fcal_9=self.Fcal_9
        Fcal_lnv_9=self.Fcal_lnv_9
        Fcal_10=self.Fcal_10
        Fcal_lnv_10=self.Fcal_lnv_10
        Fcal_11=self.Fcal_11
        Fcal_lnv_11=self.Fcal_lnv_11
        Fcal_12=self.Fcal_12
        Fcal_lnv_12=self.Fcal_lnv_12
        Fcal_lnv2_12=self.Fcal_lnv2_12
        Fcal_SQ_4=self.Fcal_SQ_4
        Fcal_SO_3=self.Fcal_SO_3
        Fcal_SO_5=self.Fcal_SO_5
        Fcal_SO_6=self.Fcal_SO_6
        Fcal_SO_7=self.Fcal_SO_7
        Fcal_SO_8=self.Fcal_SO_8
        E_0=self.E_0
        E_2=self.E_2
        E_4=self.E_4
        E_6=self.E_6
        E_8=self.E_8
        E_lnv_8=self.E_lnv_8
        E_10=self.E_10
        E_lnv_10=self.E_lnv_10
        E_11=self.E_11
        E_12=self.E_12
        E_lnv_12=self.E_lnv_12
        E_SQ_4=self.E_SQ_4
        E_SO_3=self.E_SO_3
        E_SO_5=self.E_SO_5
        E_SO_7=self.E_SO_7
        R = exp(rfrak_frame_x*xHat + rfrak_frame_y*yHat + rfrak_frame_z*zHat)
        nHat = R*xHat*conjugate(R)
        lambdaHat = R*yHat*conjugate(R)
        ellHat = R*zHat*conjugate(R)
        R_S1 = exp(rfrak_chi1_x*xHat + rfrak_chi1_y*yHat)
        R_S2 = exp(rfrak_chi2_x*xHat + rfrak_chi2_y*yHat)
        chiVec1 = S_chi1*R_S1*zHat*conjugate(R_S1)*conjugate(S_chi1)
        chiVec2 = S_chi2*R_S2*zHat*conjugate(R_S2)*conjugate(S_chi2)
        chi1_n = dot(chiVec1.vector,nHat.vector)
        chi1_lambda = dot(chiVec1.vector,lambdaHat.vector)
        chi1_ell = dot(chiVec1.vector,ellHat.vector)
        chi2_n = dot(chiVec2.vector,nHat.vector)
        chi2_lambda = dot(chiVec2.vector,lambdaHat.vector)
        chi2_ell = dot(chiVec2.vector,ellHat.vector)
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
        self.R=R
        self.nHat=nHat
        self.lambdaHat=lambdaHat
        self.ellHat=ellHat
        self.R_S1=R_S1
        self.R_S2=R_S2
        self.chiVec1=chiVec1
        self.chiVec2=chiVec2
        self.chi1_n=chi1_n
        self.chi1_lambda=chi1_lambda
        self.chi1_ell=chi1_ell
        self.chi2_n=chi2_n
        self.chi2_lambda=chi2_lambda
        self.chi2_ell=chi2_ell
        self.S_ell=S_ell
        self.S_n=S_n
        self.S_lambda=S_lambda
        self.Sigma_ell=Sigma_ell
        self.Sigma_n=Sigma_n
        self.Sigma_lambda=Sigma_lambda
        self.chi_s_ell=chi_s_ell
        self.chi_a_ell=chi_a_ell
        self.logv=logv
        self.Fcal_coeff=Fcal_coeff
        self.Fcal_SQ_4=Fcal_SQ_4
        self.Fcal_SO_3=Fcal_SO_3
        self.Fcal_SO_5=Fcal_SO_5
        self.Fcal_SO_6=Fcal_SO_6
        self.Fcal_SO_7=Fcal_SO_7
        self.Fcal_SO_8=Fcal_SO_8
        self.E_SQ_4=E_SQ_4
        self.E_SO_3=E_SO_3
        self.E_SO_5=E_SO_5
        self.E_SO_7=E_SO_7

    def OmegaVec_chiVec_1(self):
        xHat=self.xHat
        yHat=self.yHat
        zHat=self.zHat
        M1=self.M1
        M2=self.M2
        v=self.v
        S_chi1=self.S_chi1
        S_chi2=self.S_chi2
        rfrak_chi1_x=self.rfrak_chi1_x
        rfrak_chi1_y=self.rfrak_chi1_y
        rfrak_chi2_x=self.rfrak_chi2_x
        rfrak_chi2_y=self.rfrak_chi2_y
        rfrak_frame_x=self.rfrak_frame_x
        rfrak_frame_y=self.rfrak_frame_y
        rfrak_frame_z=self.rfrak_frame_z
        M=self.M
        delta=self.delta
        nu=self.nu
        R=self.R
        nHat=self.nHat
        lambdaHat=self.lambdaHat
        ellHat=self.ellHat
        R_S1=self.R_S1
        R_S2=self.R_S2
        chiVec1=self.chiVec1
        chiVec2=self.chiVec2
        chi1_n=self.chi1_n
        chi1_lambda=self.chi1_lambda
        chi1_ell=self.chi1_ell
        chi2_n=self.chi2_n
        chi2_lambda=self.chi2_lambda
        chi2_ell=self.chi2_ell
        S_ell=self.S_ell
        S_n=self.S_n
        S_lambda=self.S_lambda
        Sigma_ell=self.Sigma_ell
        Sigma_n=self.Sigma_n
        Sigma_lambda=self.Sigma_lambda
        Omega1_coeff = v**5/M
        return Omega1_coeff*(ellHat*(-0.75*delta + 0.5*nu + v**2*(delta*(0.625*nu - 0.5625) + nu*(1.25 - 0.0416666666666667*nu) + v**2*(delta*(nu*(4.875 - 0.15625*nu) - 0.84375) + nu*(nu*(-0.0208333333333333*nu - 3.28125) + 0.1875) + 0.84375) + 0.5625) + 0.75) + nHat*v*(3.0*chi1_n*nu + 3.0*M2**2*chi2_n/M**2) - M2**2*chiVec2*v/M**2)

    def OmegaVec_chiVec_2(self):
        xHat=self.xHat
        yHat=self.yHat
        zHat=self.zHat
        M1=self.M1
        M2=self.M2
        v=self.v
        S_chi1=self.S_chi1
        S_chi2=self.S_chi2
        rfrak_chi1_x=self.rfrak_chi1_x
        rfrak_chi1_y=self.rfrak_chi1_y
        rfrak_chi2_x=self.rfrak_chi2_x
        rfrak_chi2_y=self.rfrak_chi2_y
        rfrak_frame_x=self.rfrak_frame_x
        rfrak_frame_y=self.rfrak_frame_y
        rfrak_frame_z=self.rfrak_frame_z
        M=self.M
        delta=self.delta
        nu=self.nu
        R=self.R
        nHat=self.nHat
        lambdaHat=self.lambdaHat
        ellHat=self.ellHat
        R_S1=self.R_S1
        R_S2=self.R_S2
        chiVec1=self.chiVec1
        chiVec2=self.chiVec2
        chi1_n=self.chi1_n
        chi1_lambda=self.chi1_lambda
        chi1_ell=self.chi1_ell
        chi2_n=self.chi2_n
        chi2_lambda=self.chi2_lambda
        chi2_ell=self.chi2_ell
        S_ell=self.S_ell
        S_n=self.S_n
        S_lambda=self.S_lambda
        Sigma_ell=self.Sigma_ell
        Sigma_n=self.Sigma_n
        Sigma_lambda=self.Sigma_lambda
        Omega2_coeff = v**5/M
        return Omega2_coeff*(ellHat*(0.75*delta + 0.5*nu + v**2*(delta*(0.5625 - 0.625*nu) + nu*(1.25 - 0.0416666666666667*nu) + v**2*(delta*(nu*(0.15625*nu - 4.875) + 0.84375) + nu*(nu*(-0.0208333333333333*nu - 3.28125) + 0.1875) + 0.84375) + 0.5625) + 0.75) + nHat*v*(3.0*chi2_n*nu + 3.0*M1**2*chi1_n/M**2) - M1**2*chiVec1*v/M**2)

    def OmegaVec(self):
        xHat=self.xHat
        yHat=self.yHat
        zHat=self.zHat
        M1=self.M1
        M2=self.M2
        v=self.v
        S_chi1=self.S_chi1
        S_chi2=self.S_chi2
        rfrak_chi1_x=self.rfrak_chi1_x
        rfrak_chi1_y=self.rfrak_chi1_y
        rfrak_chi2_x=self.rfrak_chi2_x
        rfrak_chi2_y=self.rfrak_chi2_y
        rfrak_frame_x=self.rfrak_frame_x
        rfrak_frame_y=self.rfrak_frame_y
        rfrak_frame_z=self.rfrak_frame_z
        M=self.M
        delta=self.delta
        nu=self.nu
        R=self.R
        nHat=self.nHat
        lambdaHat=self.lambdaHat
        ellHat=self.ellHat
        R_S1=self.R_S1
        R_S2=self.R_S2
        chiVec1=self.chiVec1
        chiVec2=self.chiVec2
        chi1_n=self.chi1_n
        chi1_lambda=self.chi1_lambda
        chi1_ell=self.chi1_ell
        chi2_n=self.chi2_n
        chi2_lambda=self.chi2_lambda
        chi2_ell=self.chi2_ell
        S_ell=self.S_ell
        S_n=self.S_n
        S_lambda=self.S_lambda
        Sigma_ell=self.Sigma_ell
        Sigma_n=self.Sigma_n
        Sigma_lambda=self.Sigma_lambda
        gamma_PN_6 = 0.0123456790123457*nu**3 + 6.36111111111111*nu**2 - 2.98177812235564*nu + 1.0
        gamma_PN_4 = 1.0 - 5.41666666666667*nu
        gamma_PN_5 = (S_ell*(0.888888888888889*nu + 3.33333333333333) + 2.0*Sigma_ell*delta)/M**2
        gamma_PN_2 = 1.0 - 0.333333333333333*nu
        a_ell_4 = S_n*(5.77777777777778*nu**2 + 14.75*nu + 1.5) + Sigma_n*delta*(2.83333333333333*nu**2 + 9.125*nu + 1.5)
        gamma_PN_0 = 1.00000000000000
        gamma_PN_7 = (S_ell*(-6.0*nu**2 - 10.5833333333333*nu + 5.0) - 2.66666666666667*Sigma_ell*delta*nu**2 + Sigma_ell*delta*(3.0 - 10.1666666666667*nu))/M**2
        gamma_PN_3 = (1.66666666666667*S_ell + Sigma_ell*delta)/M**2
        a_ell_0 = 7.0*S_n + 3.0*Sigma_n*delta
        a_ell_2 = S_n*(-9.66666666666667*nu - 10.0) + Sigma_n*delta*(-4.5*nu - 6.0)
        return ellHat*v**3/M + nHat*v**6*(a_ell_0 + v**2*(a_ell_2 + a_ell_4*v**2))*(gamma_PN_0 + v**2*(gamma_PN_2 + v*(gamma_PN_3 + v*(gamma_PN_4 + v*(gamma_PN_5 + v*(gamma_PN_6 + gamma_PN_7*v))))))/M**3

    def OrbitalAngularMomentum(self):
        xHat=self.xHat
        yHat=self.yHat
        zHat=self.zHat
        M1=self.M1
        M2=self.M2
        v=self.v
        S_chi1=self.S_chi1
        S_chi2=self.S_chi2
        rfrak_chi1_x=self.rfrak_chi1_x
        rfrak_chi1_y=self.rfrak_chi1_y
        rfrak_chi2_x=self.rfrak_chi2_x
        rfrak_chi2_y=self.rfrak_chi2_y
        rfrak_frame_x=self.rfrak_frame_x
        rfrak_frame_y=self.rfrak_frame_y
        rfrak_frame_z=self.rfrak_frame_z
        M=self.M
        delta=self.delta
        nu=self.nu
        R=self.R
        nHat=self.nHat
        lambdaHat=self.lambdaHat
        ellHat=self.ellHat
        R_S1=self.R_S1
        R_S2=self.R_S2
        chiVec1=self.chiVec1
        chiVec2=self.chiVec2
        chi1_n=self.chi1_n
        chi1_lambda=self.chi1_lambda
        chi1_ell=self.chi1_ell
        chi2_n=self.chi2_n
        chi2_lambda=self.chi2_lambda
        chi2_ell=self.chi2_ell
        S_ell=self.S_ell
        S_n=self.S_n
        S_lambda=self.S_lambda
        Sigma_ell=self.Sigma_ell
        Sigma_n=self.Sigma_n
        Sigma_lambda=self.Sigma_lambda
        L_lnv_8 = -42.6666666666667*ellHat*nu
        L_6 = ellHat*(0.00540123456790123*nu**3 + 1.29166666666667*nu**2 - 30.9797035925835*nu + 8.4375)
        L_lnv_10 = 2.0*ellHat*nu*(87.4666666666667*nu + 95.0095238095238)
        L_SO_7 = 0.0625*ellHat*(-29.0*S_ell*nu**2 + 1101.0*S_ell*nu - 405.0*S_ell - 15.0*Sigma_ell*delta*nu**2 + 468.0*Sigma_ell*delta*nu - 81.0*Sigma_ell*delta)/M**2 + 0.0416666666666667*lambdaHat*(-32.0*S_lambda*nu**2 + 2.0*S_lambda*nu - 174.0*S_lambda - 16.0*Sigma_lambda*delta*nu**2 - 79.0*Sigma_lambda*delta*nu - 12.0*Sigma_lambda*delta)/M**2 + 0.0208333333333333*nHat*(11.0*S_n*nu**2 - 1331.0*S_n*nu + 183.0*S_n + 5.0*Sigma_n*delta*nu**2 - 734.0*Sigma_n*delta*nu + 183.0*Sigma_n*delta)/M**2
        L_SO_5 = 0.0138888888888889*ellHat*(427.0*S_ell*nu - 693.0*S_ell + 210.0*Sigma_ell*delta*nu - 189.0*Sigma_ell*delta)/M**2 + 0.166666666666667*lambdaHat*(18.0*S_lambda*nu - 21.0*S_lambda + 8.0*Sigma_lambda*delta*nu - 3.0*Sigma_lambda*delta)/M**2 + 0.0416666666666667*nHat*(-19.0*S_n*nu + 33.0*S_n - 10.0*Sigma_n*delta*nu + 33.0*Sigma_n*delta)/M**2
        L_SO_3 = 0.166666666666667*ellHat*(-35.0*S_ell - 15.0*Sigma_ell*delta)/M**2 + lambdaHat*(-3.0*S_lambda - Sigma_lambda*delta)/M**2 + 0.5*nHat*(S_n + Sigma_n*delta)/M**2
        L_2 = ellHat*(0.166666666666667*nu + 1.5)
        L_10 = ellHat*(nu*(-4.85925925925926*nu - 5.27830687830688) + 59.80078125)
        L_4 = ellHat*(0.0416666666666667*nu**2 - 2.375*nu + 3.375)
        L_0 = ellHat
        L_8 = ellHat*(-0.714285714285714*nu*(0.0024755658436214*nu**3 + 0.174189814814815*nu**2 - 90.1327990262052*nu + 153.88379682994) + 1.82857142857143*nu + 22.1484375)
        L_coeff = M**2*nu/v
        return L_coeff*(L_0 + v**2*(L_2 + v*(L_SO_3 + v*(L_4 + v*(L_SO_5 + v*(L_6 + v*(L_SO_7 + v*(L_8 + L_lnv_8*math.log(v) + v**2*(L_10 + L_lnv_10*math.log(v))))))))))

               

    #Evolve PN
    def __TaylorT1__(self,t, y):
        self.__Recalculate__(t, y)
        xHat=self.xHat
        yHat=self.yHat
        zHat=self.zHat
        M1=self.M1
        M2=self.M2
        v=self.v
        S_chi1=self.S_chi1
        S_chi2=self.S_chi2
        rfrak_chi1_x=self.rfrak_chi1_x
        rfrak_chi1_y=self.rfrak_chi1_y
        rfrak_chi2_x=self.rfrak_chi2_x
        rfrak_chi2_y=self.rfrak_chi2_y
        rfrak_frame_x=self.rfrak_frame_x
        rfrak_frame_y=self.rfrak_frame_y
        rfrak_frame_z=self.rfrak_frame_z
        M=self.M
        delta=self.delta
        nu=self.nu
        R=self.R
        nHat=self.nHat
        lambdaHat=self.lambdaHat
        ellHat=self.ellHat
        R_S1=self.R_S1
        R_S2=self.R_S2
        chiVec1=self.chiVec1
        chiVec2=self.chiVec2
        chi1chi1=self.chi1chi1
        chi1chi2=self.chi1chi2
        chi2chi2=self.chi2chi2
        chi1_n=self.chi1_n
        chi1_lambda=self.chi1_lambda
        chi1_ell=self.chi1_ell
        chi2_n=self.chi2_n
        chi2_lambda=self.chi2_lambda
        chi2_ell=self.chi2_ell
        S_ell=self.S_ell
        S_n=self.S_n
        S_lambda=self.S_lambda
        Sigma_ell=self.Sigma_ell
        Sigma_n=self.Sigma_n
        Sigma_lambda=self.Sigma_lambda
        chi_s_ell=self.chi_s_ell
        chi_a_ell=self.chi_a_ell
        logv=self.logv
        Fcal_coeff=self.Fcal_coeff
        Fcal_0=self.Fcal_0
        Fcal_2=self.Fcal_2
        Fcal_3=self.Fcal_3
        Fcal_4=self.Fcal_4
        Fcal_5=self.Fcal_5
        Fcal_6=self.Fcal_6
        Fcal_lnv_6=self.Fcal_lnv_6
        Fcal_7=self.Fcal_7
        Fcal_8=self.Fcal_8
        Fcal_lnv_8=self.Fcal_lnv_8
        Fcal_9=self.Fcal_9
        Fcal_lnv_9=self.Fcal_lnv_9
        Fcal_10=self.Fcal_10
        Fcal_lnv_10=self.Fcal_lnv_10
        Fcal_11=self.Fcal_11
        Fcal_lnv_11=self.Fcal_lnv_11
        Fcal_12=self.Fcal_12
        Fcal_lnv_12=self.Fcal_lnv_12
        Fcal_lnv2_12=self.Fcal_lnv2_12
        Fcal_SQ_4=self.Fcal_SQ_4
        Fcal_SO_3=self.Fcal_SO_3
        Fcal_SO_5=self.Fcal_SO_5
        Fcal_SO_6=self.Fcal_SO_6
        Fcal_SO_7=self.Fcal_SO_7
        Fcal_SO_8=self.Fcal_SO_8
        E_0=self.E_0
        E_2=self.E_2
        E_4=self.E_4
        E_6=self.E_6
        E_8=self.E_8
        E_lnv_8=self.E_lnv_8
        E_10=self.E_10
        E_lnv_10=self.E_lnv_10
        E_11=self.E_11
        E_12=self.E_12
        E_lnv_12=self.E_lnv_12
        E_SQ_4=self.E_SQ_4
        E_SO_3=self.E_SO_3
        E_SO_5=self.E_SO_5
        E_SO_7=self.E_SO_7
        v, rfrak_chi1_x, rfrak_chi1_y, rfrak_chi2_x, rfrak_chi2_y=self.v, self.rfrak_chi1_x, self.rfrak_chi1_y, self.rfrak_chi2_x, self.rfrak_chi2_y
        rfrak_frame_x, rfrak_frame_y, rfrak_frame_z, Phi=self.rfrak_frame_x, self.rfrak_frame_y, self.rfrak_frame_z, self.Phi
        dvdt_T1=0
        if v>=1.0 and self.NotForward:
            print("Beyond domain of PN validity, this is a good way to terminate.")
            self.terminal1=False
        Flux = Fcal_coeff*(Fcal_0 + v**2*(Fcal_2 + v*(Fcal_3 + Fcal_SO_3 + v*(Fcal_4 + Fcal_SQ_4 + v*(Fcal_5 + Fcal_SO_5 + v*(Fcal_6 + Fcal_SO_6 + Fcal_lnv_6*logv + v*(Fcal_7 + Fcal_SO_7 + v*(Fcal_8 + Fcal_SO_8 + Fcal_lnv_8*logv + v*(Fcal_9 + Fcal_lnv_9*logv + v*(Fcal_10 + Fcal_lnv_10*logv + v*(Fcal_11 + Fcal_lnv_11*logv + v*(Fcal_12 + Fcal_lnv2_12*logv**2 + Fcal_lnv_12*logv))))))))))))
        dEdv = -0.5*M*nu*v*(2.0*E_0 + v**2*(4.0*E_2 + v*(5.0*E_SO_3 + v*(6.0*E_4 + 6.0*E_SQ_4 + v*(7.0*E_SO_5 + v*(8.0*E_6 + v*(9.0*E_SO_7 + v*(10.0*E_8 + E_lnv_8*(10.0*logv + 1.0) + v**2*(12.0*E_10 + E_lnv_10*(12.0*logv + 1.0) + v*(13.0*E_11 + v*(14.0*E_12 + E_lnv_12*(14.0*logv + 1.0))))))))))))
        Absorption = 0
        dvdt_T1 = (-Absorption - Flux)/dEdv
        if dvdt_T1<1.0e-12 and self.NotForward:
            print("v is decreasing, which is not an uncommon way to stop.")
            self.terminal2=False
        return self.CommonRHS(dvdt_T1, y)

    def CommonRHS(self,dvdt, y):
        dydt=np.zeros(9)
        rfrak_frame=np.zeros(3)
        rfrak_frame[0] = y[5]
        rfrak_frame[1] = y[6]
        rfrak_frame[2] = y[7]
        rfrakdot_frame = FrameFromAngularVelocityIntegrand(rfrak_frame, self.OmegaVec().vector)
        dydt[0] = dvdt
        if(self.EvolveSpin1):
            dydt[1], dydt[2]= FrameFromAngularVelocity_2D_Integrand(y[1], y[2],\
                (self.S_chi1.inverse*self.OmegaVec_chiVec_1()*self.S_chi1).vector)
        else:
            dydt[1] = 0.0
            dydt[2] = 0.0
        if(self.EvolveSpin2):
            dydt[3], dydt[4] = FrameFromAngularVelocity_2D_Integrand(y[3], y[4],\
                (self.S_chi2.inverse*self.OmegaVec_chiVec_2()*self.S_chi2).vector)
        else:
            dydt[3] = 0.0
            dydt[4] = 0.0
        dydt[5] = rfrakdot_frame[0]
        dydt[6] = rfrakdot_frame[1]
        dydt[7] = rfrakdot_frame[2]
        dydt[8] = self.v*self.v*self.v/self.M

        return dydt

    # Prepare for integration
    def Evolution(self): 
        v_i,rfrak_chi1_x_i,rfrak_chi1_y_i,rfrak_chi2_x_i,rfrak_chi2_y_i,rfrak_frame_x_i,rfrak_frame_y_i,rfrak_frame_z_i\
            =self.v_i,self.rfrak_chi1_x_i,self.rfrak_chi1_y_i,self.rfrak_chi2_x_i,self.rfrak_chi2_y_i,self.rfrak_frame_x_i,self.rfrak_frame_y_i,self.rfrak_frame_z_i
        def terminate(t,y):
            return 1.0*self.terminal1*self.terminal2
        terminate.terminal=True
        TMerger=5.0/(256.0*self.nu*v_i**8)
        TEnd=TMerger
        if self.t_PN2:
            TEnd=self.t_PN2
        time=[0.0]
        while time[-1]<TEnd and 2*self.M*(256*self.nu*(TMerger-time[-1])/5)**(3/8)/self.StepsPerOrbit>self.MinStep:
            time.append(time[-1]+2*self.M*(256*self.nu*(TMerger-time[-1])/5)**(3/8)/self.StepsPerOrbit)
        time=np.delete(time, -1)
        
        # Integrate
        yy=solve_ivp(self.__TaylorT1__, [time[0],time[-1]], [v_i,rfrak_chi1_x_i,\
            rfrak_chi1_y_i,rfrak_chi2_x_i,rfrak_chi2_y_i,rfrak_frame_x_i,\
            rfrak_frame_y_i,rfrak_frame_z_i,0.0], method='DOP853',\
            t_eval=time, dense_output=True, events=terminate, rtol=self.tol, atol=self.tol)           
        if self.ForwardInTime:
            self.NotForward=False
            time=[0.0]
            TStart=-3*TMerger
            if self.t_PN1:
                TStart=self.t_PN1
            while time[-1]>TStart:
                time.append(time[-1]-2*self.M*(256*self.nu*(TMerger-time[-1])/5)**(3/8)/self.StepsPerOrbit)
            yyForward=solve_ivp(self.__TaylorT1__, [time[0],time[-1]], [v_i,rfrak_chi1_x_i,\
                rfrak_chi1_y_i,rfrak_chi2_x_i,rfrak_chi2_y_i,rfrak_frame_x_i,\
                rfrak_frame_y_i,rfrak_frame_z_i,0.0], method='DOP853',\
                t_eval=time, dense_output=True, rtol=self.tol, atol=self.tol)
            yy.t=np.append(yyForward.t[1:][::-1],yy.t)
            data=np.empty((9,len(yy.t)))
            for i in range(9):
                data[i]=np.append(yyForward.y[i][1:][::-1],yy.y[i])
            yy.y=data
            
        return yy
