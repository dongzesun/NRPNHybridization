# File produced automatically by PNCodeGen.ipynb
from scipy.integrate import solve_ivp
import numpy as np
from numpy import dot, cross, log, sqrt, pi
from numpy import euler_gamma as EulerGamma
from numba import jit, njit
from numba.typed import List, Dict
from scipy.interpolate import InterpolatedUnivariateSpline as Spline
from scipy.special import zeta
import quaternionic

qmul = njit(quaternionic.algebra.multiply)
qexp=njit(quaternionic.algebra.exp)
qconj=njit(quaternionic.algebra.conj)
qinverse=njit(quaternionic.algebra.reciprocal)

@njit
def mul(A,B):
    C=np.empty(4)
    qmul(A,B,C)
    return C
    
@njit
def exp(A):
    B=np.empty(4)
    qexp(A,B)
    return B
    
@njit
def conjugate(A):
    B=np.empty(4)
    qconj(A,B)
    return B
    
@njit
def inverse(A):
    B=np.empty(4)
    qinverse(A,B)
    return B

@njit
def FrameFromAngularVelocity_2D_Integrand(rfrak_x, rfrak_y, Omega):
    rfrakMag = np.sqrt(rfrak_x*rfrak_x+rfrak_y*rfrak_y)
    rfrakDot_x = Omega[0]/2.0
    rfrakDot_y = Omega[1]/2.0
    if np.abs(np.sin(rfrakMag)) > 1e-12 and np.abs(np.cos(rfrakMag)) > 1e-12:
        omega_v = (Omega[0]*(-rfrak_y/rfrakMag)+Omega[1]*(rfrak_x/rfrakMag))*np.tan(rfrakMag)-Omega[2]
        Omega[0] += -omega_v*np.sin(2*rfrakMag)*(-rfrak_y/rfrakMag)
        Omega[1] += -omega_v*np.sin(2*rfrakMag)*(rfrak_x/rfrakMag)
        Omega[2] +=  omega_v*np.cos(2*rfrakMag)
        dotTerm = (rfrak_x*Omega[0]+rfrak_y*Omega[1])/(rfrakMag*rfrakMag)
        cotTerm = rfrakMag/(2*np.tan(rfrakMag))
        rfrakDot_x = (Omega[0] - rfrak_x*dotTerm)*cotTerm + rfrak_x*dotTerm/2. - 0.5*Omega[2]*rfrak_y
        rfrakDot_y = (Omega[1] - rfrak_y*dotTerm)*cotTerm + rfrak_y*dotTerm/2. + 0.5*Omega[2]*rfrak_x
    return rfrakDot_x, rfrakDot_y

@njit
def FrameFromAngularVelocityIntegrand(rfrak, Omega):
    rfrakMag = np.sqrt(rfrak[0] * rfrak[0] + rfrak[1] * rfrak[1] + rfrak[2] * rfrak[2])
    OmegaMag = np.sqrt(Omega[0] * Omega[0] + Omega[1] * Omega[1] + Omega[2] * Omega[2])
    # If the matrix is really close to the identity, return
    if rfrakMag < 1e-12*OmegaMag:
        return np.array([Omega[0] / 2.0, Omega[1] / 2.0, Omega[2] / 2.0])
    # If the matrix is really close to singular, it's equivalent to the identity, so return
    if np.abs(np.sin(rfrakMag)) < 1e-12:
        return np.array([Omega[0] / 2.0, Omega[1] / 2.0, Omega[2] / 2.0])
    OmegaOver2 = np.array([Omega[0] / 2.0, Omega[1] / 2.0, Omega[2] / 2.0])
    rfrakHat = np.array([rfrak[0] / rfrakMag, rfrak[1] / rfrakMag, rfrak[2] / rfrakMag])
    return ((OmegaOver2 - rfrakHat * np.dot(rfrakHat, OmegaOver2)) * (rfrakMag / np.tan(rfrakMag))\
        + rfrakHat * np.dot(rfrakHat, OmegaOver2) + np.cross(OmegaOver2, rfrak))  

@njit
def Initialization(xHat_i, yHat_i, zHat_i, M1_i, M2_i, v_i, S_chi1_i, S_chi2_i, rfrak_chi1_x_i, rfrak_chi1_y_i, rfrak_chi2_x_i, rfrak_chi2_y_i, rfrak_frame_x_i, rfrak_frame_y_i, rfrak_frame_z_i): 
    Cons=Dict()
    Cons['xHat']=xHat_i
    Cons['yHat']=yHat_i
    Cons['zHat']=zHat_i
    Cons['M1']=np.array([M1_i])
    Cons['M2']=np.array([M2_i])
    v=np.array([v_i])
    Cons['S_chi1']=S_chi1_i
    Cons['S_chi2']=S_chi2_i
    rfrak_chi1_x=np.array([rfrak_chi1_x_i])
    rfrak_chi1_y=np.array([rfrak_chi1_y_i])
    rfrak_chi2_x=np.array([rfrak_chi2_x_i])
    rfrak_chi2_y=np.array([rfrak_chi2_y_i])
    rfrak_frame_x=np.array([rfrak_frame_x_i])
    rfrak_frame_y=np.array([rfrak_frame_y_i])
    rfrak_frame_z=np.array([rfrak_frame_z_i])
    Cons['M']=Cons['M1'] + Cons['M2']
    Cons['delta']=(Cons['M1'] - Cons['M2'])/Cons['M']
    Cons['nu']=Cons['M1']*Cons['M2']/Cons['M']**2
    R=exp(rfrak_frame_x*Cons['xHat'] + rfrak_frame_y*Cons['yHat'] + rfrak_frame_z*Cons['zHat'])
    nHat=mul(mul(R,Cons['xHat']),conjugate(R))
    ellHat=mul(mul(R,Cons['zHat']),conjugate(R))
    R_S1=exp(rfrak_chi1_x*Cons['xHat'] + rfrak_chi1_y*Cons['yHat'])
    R_S2=exp(rfrak_chi2_x*Cons['xHat'] + rfrak_chi2_y*Cons['yHat'])
    chiVec1=mul(mul(mul(Cons['S_chi1'],R_S1),Cons['zHat']),mul(conjugate(R_S1),conjugate(Cons['S_chi1'])))
    chiVec2=mul(mul(mul(Cons['S_chi2'],R_S2),Cons['zHat']),mul(conjugate(R_S2),conjugate(Cons['S_chi2'])))
    Cons['chi1chi1']=np.array([dot(chiVec1[1:],chiVec1[1:])])
    Cons['chi1chi2']=np.array([dot(chiVec1[1:],chiVec2[1:])])
    Cons['chi2chi2']=np.array([dot(chiVec2[1:],chiVec2[1:])])
    chi1_n=np.array([dot(chiVec1[1:],nHat[1:])])
    chi1_ell=np.array([dot(chiVec1[1:],ellHat[1:])])
    chi2_n=np.array([dot(chiVec2[1:],nHat[1:])])
    chi2_ell=np.array([dot(chiVec2[1:],ellHat[1:])])
    S_ell=Cons['M1']**2*chi1_ell + Cons['M2']**2*chi2_ell
    S_n=Cons['M1']**2*chi1_n + Cons['M2']**2*chi2_n
    Sigma_ell=Cons['M']*(-Cons['M1']*chi1_ell + Cons['M2']*chi2_ell)
    Sigma_n=Cons['M']*(-Cons['M1']*chi1_n + Cons['M2']*chi2_n)
    chi_s_ell=chi1_ell/2 + chi2_ell/2
    chi_a_ell=chi1_ell/2 - chi2_ell/2
    Fcal_coeff=32*Cons['nu']**2*v**10/5
    Cons['Fcal_0']=np.array([1.0])
    Cons['Fcal_2']=-35*Cons['nu']/12 - 1247/336
    Cons['E_0']=np.array([1.0])
    Cons['E_2']=-Cons['nu']/12 - 3/4
    return Cons

@njit
def Recalculate_0(Cons,y):
    Vars=Dict()
    Vars['v'] = np.array([y[0]])
    Vars['rfrak_chi1_x'] = np.array([y[1]])
    Vars['rfrak_chi1_y'] = np.array([y[2]])
    Vars['rfrak_chi2_x'] = np.array([y[3]])
    Vars['rfrak_chi2_y'] = np.array([y[4]])
    Vars['rfrak_frame_x'] = np.array([y[5]])
    Vars['rfrak_frame_y'] = np.array([y[6]])
    Vars['rfrak_frame_z'] = np.array([y[7]])
    Vars['Phi'] = np.array([y[8]])
    Vars['R'] = exp(Vars['rfrak_frame_x']*Cons['xHat'] + Vars['rfrak_frame_y']*Cons['yHat'] + Vars['rfrak_frame_z']*Cons['zHat'])
    Vars['nHat'] = mul(mul(Vars['R'],Cons['xHat']),conjugate(Vars['R']))
    Vars['ellHat'] = mul(mul(Vars['R'],Cons['zHat']),conjugate(Vars['R']))
    Vars['R_S1'] = exp(Vars['rfrak_chi1_x']*Cons['xHat'] + Vars['rfrak_chi1_y']*Cons['yHat'])
    Vars['R_S2'] = exp(Vars['rfrak_chi2_x']*Cons['xHat'] + Vars['rfrak_chi2_y']*Cons['yHat'])
    Vars['chiVec1'] = mul(mul(mul(Cons['S_chi1'],Vars['R_S1']),Cons['zHat']),mul(conjugate(Vars['R_S1']),conjugate(Cons['S_chi1'])))
    Vars['chiVec2'] = mul(mul(mul(Cons['S_chi2'],Vars['R_S2']),Cons['zHat']),mul(conjugate(Vars['R_S2']),conjugate(Cons['S_chi2'])))
    Vars['chi1_n'] = np.array([dot(Vars['chiVec1'][1:],Vars['nHat'][1:])])
    Vars['chi1_ell'] = np.array([dot(Vars['chiVec1'][1:],Vars['ellHat'][1:])])
    Vars['chi2_n'] = np.array([dot(Vars['chiVec2'][1:],Vars['nHat'][1:])])
    Vars['chi2_ell'] = np.array([dot(Vars['chiVec2'][1:],Vars['ellHat'][1:])])
    Vars['S_ell'] = Cons['M1']**2*Vars['chi1_ell'] + Cons['M2']**2*Vars['chi2_ell']
    Vars['S_n'] = Cons['M1']**2*Vars['chi1_n'] + Cons['M2']**2*Vars['chi2_n']
    Vars['Sigma_ell'] = Cons['M']*(-Cons['M1']*Vars['chi1_ell'] + Cons['M2']*Vars['chi2_ell'])
    Vars['Sigma_n'] = Cons['M']*(-Cons['M1']*Vars['chi1_n'] + Cons['M2']*Vars['chi2_n'])
    Vars['chi_s_ell'] = Vars['chi1_ell']/2 + Vars['chi2_ell']/2
    Vars['chi_a_ell'] = Vars['chi1_ell']/2 - Vars['chi2_ell']/2
    Vars['Fcal_coeff'] = 32*Cons['nu']**2*Vars['v']**10/5
    return Vars

@njit
def OmegaVec_chiVec_1_0(Cons,Vars):
    Omega1_coeff = Vars['v']**5/Cons['M']
    return Omega1_coeff*Vars['ellHat']*(-0.75*Cons['delta'] + 0.5*Cons['nu'] + 0.75)

@njit
def OmegaVec_chiVec_2_0(Cons,Vars):
    Omega2_coeff = Vars['v']**5/Cons['M']
    return Omega2_coeff*Vars['ellHat']*(0.75*Cons['delta'] + 0.5*Cons['nu'] + 0.75)

@njit
def OmegaVec_0(Cons,Vars):
    a_ell_0 = 7.0*Vars['S_n'] + 3.0*Vars['Sigma_n']*Cons['delta']
    gamma_PN_0 = 1.00000000000000
    return Vars['ellHat']*Vars['v']**3/Cons['M'] + a_ell_0*gamma_PN_0*Vars['nHat']*Vars['v']**6/Cons['M']**3


@njit
def TaylorT1_0(Cons,Vars):
    Flux = Cons['Fcal_0']*Vars['Fcal_coeff']
    dEdV = -Cons['E_0']*Cons['M']*Cons['nu']*Vars['v']
    Absorption = 0
    dvdt = (-Absorption - Flux)/dEdV
    dydt=np.zeros(9)
    rfrak_frame=[Vars['rfrak_frame_x'][0],Vars['rfrak_frame_y'][0],Vars['rfrak_frame_z'][0]]
    [dydt[5],dydt[6],dydt[7]] = FrameFromAngularVelocityIntegrand(rfrak_frame, OmegaVec_0(Cons,Vars)[1:])
    dydt[0] = dvdt[0]
    if(Cons['EvolveSpin1']):
        dydt[1], dydt[2]=FrameFromAngularVelocity_2D_Integrand(Vars['rfrak_chi1_x'][0], Vars['rfrak_chi1_y'][0],(mul(mul(inverse(Cons['S_chi1']),OmegaVec_chiVec_1_0(Cons,Vars)),Cons['S_chi1']))[1:])
    else:
        dydt[1] = 0.0
        dydt[2] = 0.0
    if(Cons['EvolveSpin2']):
        dydt[3], dydt[4]=FrameFromAngularVelocity_2D_Integrand(Vars['rfrak_chi2_x'][0], Vars['rfrak_chi2_y'][0],(mul(mul(inverse(Cons['S_chi2']),OmegaVec_chiVec_2_0(Cons,Vars)),Cons['S_chi2']))[1:])
    else:
        dydt[3] = 0.0
        dydt[4] = 0.0
    dydt[8] = (Vars['v']*Vars['v']*Vars['v']/Cons['M'])[0]
    return dydt      

@njit
def Recalculate_0p50(Cons,y):
    Vars=Dict()
    Vars['v'] = np.array([y[0]])
    Vars['rfrak_chi1_x'] = np.array([y[1]])
    Vars['rfrak_chi1_y'] = np.array([y[2]])
    Vars['rfrak_chi2_x'] = np.array([y[3]])
    Vars['rfrak_chi2_y'] = np.array([y[4]])
    Vars['rfrak_frame_x'] = np.array([y[5]])
    Vars['rfrak_frame_y'] = np.array([y[6]])
    Vars['rfrak_frame_z'] = np.array([y[7]])
    Vars['Phi'] = np.array([y[8]])
    Vars['R'] = exp(Vars['rfrak_frame_x']*Cons['xHat'] + Vars['rfrak_frame_y']*Cons['yHat'] + Vars['rfrak_frame_z']*Cons['zHat'])
    Vars['nHat'] = mul(mul(Vars['R'],Cons['xHat']),conjugate(Vars['R']))
    Vars['ellHat'] = mul(mul(Vars['R'],Cons['zHat']),conjugate(Vars['R']))
    Vars['R_S1'] = exp(Vars['rfrak_chi1_x']*Cons['xHat'] + Vars['rfrak_chi1_y']*Cons['yHat'])
    Vars['R_S2'] = exp(Vars['rfrak_chi2_x']*Cons['xHat'] + Vars['rfrak_chi2_y']*Cons['yHat'])
    Vars['chiVec1'] = mul(mul(mul(Cons['S_chi1'],Vars['R_S1']),Cons['zHat']),mul(conjugate(Vars['R_S1']),conjugate(Cons['S_chi1'])))
    Vars['chiVec2'] = mul(mul(mul(Cons['S_chi2'],Vars['R_S2']),Cons['zHat']),mul(conjugate(Vars['R_S2']),conjugate(Cons['S_chi2'])))
    Vars['chi1_n'] = np.array([dot(Vars['chiVec1'][1:],Vars['nHat'][1:])])
    Vars['chi1_ell'] = np.array([dot(Vars['chiVec1'][1:],Vars['ellHat'][1:])])
    Vars['chi2_n'] = np.array([dot(Vars['chiVec2'][1:],Vars['nHat'][1:])])
    Vars['chi2_ell'] = np.array([dot(Vars['chiVec2'][1:],Vars['ellHat'][1:])])
    Vars['S_ell'] = Cons['M1']**2*Vars['chi1_ell'] + Cons['M2']**2*Vars['chi2_ell']
    Vars['S_n'] = Cons['M1']**2*Vars['chi1_n'] + Cons['M2']**2*Vars['chi2_n']
    Vars['Sigma_ell'] = Cons['M']*(-Cons['M1']*Vars['chi1_ell'] + Cons['M2']*Vars['chi2_ell'])
    Vars['Sigma_n'] = Cons['M']*(-Cons['M1']*Vars['chi1_n'] + Cons['M2']*Vars['chi2_n'])
    Vars['chi_s_ell'] = Vars['chi1_ell']/2 + Vars['chi2_ell']/2
    Vars['chi_a_ell'] = Vars['chi1_ell']/2 - Vars['chi2_ell']/2
    Vars['Fcal_coeff'] = 32*Cons['nu']**2*Vars['v']**10/5
    return Vars

@njit
def OmegaVec_chiVec_1_0p50(Cons,Vars):
    Omega1_coeff = Vars['v']**5/Cons['M']
    return Omega1_coeff*(Vars['ellHat']*(-0.75*Cons['delta'] + 0.5*Cons['nu'] + 0.75) + Vars['nHat']*Vars['v']*(3.0*Vars['chi1_n']*Cons['nu'] + 3.0*Cons['M2']**2*Vars['chi2_n']/Cons['M']**2) - Cons['M2']**2*Vars['chiVec2']*Vars['v']/Cons['M']**2)

@njit
def OmegaVec_chiVec_2_0p50(Cons,Vars):
    Omega2_coeff = Vars['v']**5/Cons['M']
    return Omega2_coeff*(Vars['ellHat']*(0.75*Cons['delta'] + 0.5*Cons['nu'] + 0.75) + Vars['nHat']*Vars['v']*(3.0*Vars['chi2_n']*Cons['nu'] + 3.0*Cons['M1']**2*Vars['chi1_n']/Cons['M']**2) - Cons['M1']**2*Vars['chiVec1']*Vars['v']/Cons['M']**2)

@njit
def OmegaVec_0p50(Cons,Vars):
    a_ell_0 = 7.0*Vars['S_n'] + 3.0*Vars['Sigma_n']*Cons['delta']
    gamma_PN_0 = 1.00000000000000
    return Vars['ellHat']*Vars['v']**3/Cons['M'] + a_ell_0*gamma_PN_0*Vars['nHat']*Vars['v']**6/Cons['M']**3


@njit
def TaylorT1_0p50(Cons,Vars):
    Flux = Cons['Fcal_0']*Vars['Fcal_coeff']
    dEdV = -Cons['E_0']*Cons['M']*Cons['nu']*Vars['v']
    Absorption = 0
    dvdt = (-Absorption - Flux)/dEdV
    dydt=np.zeros(9)
    rfrak_frame=[Vars['rfrak_frame_x'][0],Vars['rfrak_frame_y'][0],Vars['rfrak_frame_z'][0]]
    [dydt[5],dydt[6],dydt[7]] = FrameFromAngularVelocityIntegrand(rfrak_frame, OmegaVec_0p50(Cons,Vars)[1:])
    dydt[0] = dvdt[0]
    if(Cons['EvolveSpin1']):
        dydt[1], dydt[2]=FrameFromAngularVelocity_2D_Integrand(Vars['rfrak_chi1_x'][0], Vars['rfrak_chi1_y'][0],(mul(mul(inverse(Cons['S_chi1']),OmegaVec_chiVec_1_0p50(Cons,Vars)),Cons['S_chi1']))[1:])
    else:
        dydt[1] = 0.0
        dydt[2] = 0.0
    if(Cons['EvolveSpin2']):
        dydt[3], dydt[4]=FrameFromAngularVelocity_2D_Integrand(Vars['rfrak_chi2_x'][0], Vars['rfrak_chi2_y'][0],(mul(mul(inverse(Cons['S_chi2']),OmegaVec_chiVec_2_0p50(Cons,Vars)),Cons['S_chi2']))[1:])
    else:
        dydt[3] = 0.0
        dydt[4] = 0.0
    dydt[8] = (Vars['v']*Vars['v']*Vars['v']/Cons['M'])[0]
    return dydt      

@njit
def Recalculate_1p0(Cons,y):
    Vars=Dict()
    Vars['v'] = np.array([y[0]])
    Vars['rfrak_chi1_x'] = np.array([y[1]])
    Vars['rfrak_chi1_y'] = np.array([y[2]])
    Vars['rfrak_chi2_x'] = np.array([y[3]])
    Vars['rfrak_chi2_y'] = np.array([y[4]])
    Vars['rfrak_frame_x'] = np.array([y[5]])
    Vars['rfrak_frame_y'] = np.array([y[6]])
    Vars['rfrak_frame_z'] = np.array([y[7]])
    Vars['Phi'] = np.array([y[8]])
    Vars['R'] = exp(Vars['rfrak_frame_x']*Cons['xHat'] + Vars['rfrak_frame_y']*Cons['yHat'] + Vars['rfrak_frame_z']*Cons['zHat'])
    Vars['nHat'] = mul(mul(Vars['R'],Cons['xHat']),conjugate(Vars['R']))
    Vars['ellHat'] = mul(mul(Vars['R'],Cons['zHat']),conjugate(Vars['R']))
    Vars['R_S1'] = exp(Vars['rfrak_chi1_x']*Cons['xHat'] + Vars['rfrak_chi1_y']*Cons['yHat'])
    Vars['R_S2'] = exp(Vars['rfrak_chi2_x']*Cons['xHat'] + Vars['rfrak_chi2_y']*Cons['yHat'])
    Vars['chiVec1'] = mul(mul(mul(Cons['S_chi1'],Vars['R_S1']),Cons['zHat']),mul(conjugate(Vars['R_S1']),conjugate(Cons['S_chi1'])))
    Vars['chiVec2'] = mul(mul(mul(Cons['S_chi2'],Vars['R_S2']),Cons['zHat']),mul(conjugate(Vars['R_S2']),conjugate(Cons['S_chi2'])))
    Vars['chi1_n'] = np.array([dot(Vars['chiVec1'][1:],Vars['nHat'][1:])])
    Vars['chi1_ell'] = np.array([dot(Vars['chiVec1'][1:],Vars['ellHat'][1:])])
    Vars['chi2_n'] = np.array([dot(Vars['chiVec2'][1:],Vars['nHat'][1:])])
    Vars['chi2_ell'] = np.array([dot(Vars['chiVec2'][1:],Vars['ellHat'][1:])])
    Vars['S_ell'] = Cons['M1']**2*Vars['chi1_ell'] + Cons['M2']**2*Vars['chi2_ell']
    Vars['S_n'] = Cons['M1']**2*Vars['chi1_n'] + Cons['M2']**2*Vars['chi2_n']
    Vars['Sigma_ell'] = Cons['M']*(-Cons['M1']*Vars['chi1_ell'] + Cons['M2']*Vars['chi2_ell'])
    Vars['Sigma_n'] = Cons['M']*(-Cons['M1']*Vars['chi1_n'] + Cons['M2']*Vars['chi2_n'])
    Vars['chi_s_ell'] = Vars['chi1_ell']/2 + Vars['chi2_ell']/2
    Vars['chi_a_ell'] = Vars['chi1_ell']/2 - Vars['chi2_ell']/2
    Vars['Fcal_coeff'] = 32*Cons['nu']**2*Vars['v']**10/5
    return Vars

@njit
def OmegaVec_chiVec_1_1p0(Cons,Vars):
    Omega1_coeff = Vars['v']**5/Cons['M']
    return Omega1_coeff*(Vars['ellHat']*(-0.75*Cons['delta'] + 0.5*Cons['nu'] + Vars['v']**2*(Cons['delta']*(0.625*Cons['nu'] - 0.5625) + Cons['nu']*(1.25 - 0.0416666666666667*Cons['nu']) + 0.5625) + 0.75) + Vars['nHat']*Vars['v']*(3.0*Vars['chi1_n']*Cons['nu'] + 3.0*Cons['M2']**2*Vars['chi2_n']/Cons['M']**2) - Cons['M2']**2*Vars['chiVec2']*Vars['v']/Cons['M']**2)

@njit
def OmegaVec_chiVec_2_1p0(Cons,Vars):
    Omega2_coeff = Vars['v']**5/Cons['M']
    return Omega2_coeff*(Vars['ellHat']*(0.75*Cons['delta'] + 0.5*Cons['nu'] + Vars['v']**2*(Cons['delta']*(0.5625 - 0.625*Cons['nu']) + Cons['nu']*(1.25 - 0.0416666666666667*Cons['nu']) + 0.5625) + 0.75) + Vars['nHat']*Vars['v']*(3.0*Vars['chi2_n']*Cons['nu'] + 3.0*Cons['M1']**2*Vars['chi1_n']/Cons['M']**2) - Cons['M1']**2*Vars['chiVec1']*Vars['v']/Cons['M']**2)

@njit
def OmegaVec_1p0(Cons,Vars):
    a_ell_0 = 7.0*Vars['S_n'] + 3.0*Vars['Sigma_n']*Cons['delta']
    gamma_PN_2 = 1.0 - 0.333333333333333*Cons['nu']
    gamma_PN_0 = 1.00000000000000
    a_ell_2 = Vars['S_n']*(-9.66666666666667*Cons['nu'] - 10.0) + Vars['Sigma_n']*Cons['delta']*(-4.5*Cons['nu'] - 6.0)
    return Vars['ellHat']*Vars['v']**3/Cons['M'] + Vars['nHat']*Vars['v']**6*(a_ell_0 + a_ell_2*Vars['v']**2)*(gamma_PN_0 + gamma_PN_2*Vars['v']**2)/Cons['M']**3


@njit
def TaylorT1_1p0(Cons,Vars):
    Flux = Vars['Fcal_coeff']*(Cons['Fcal_0'] + Cons['Fcal_2']*Vars['v']**2)
    dEdV = -Cons['M']*Cons['nu']*Vars['v']*(Cons['E_0'] + 2.0*Cons['E_2']*Vars['v']**2)
    Absorption = 0
    dvdt = (-Absorption - Flux)/dEdV
    dydt=np.zeros(9)
    rfrak_frame=[Vars['rfrak_frame_x'][0],Vars['rfrak_frame_y'][0],Vars['rfrak_frame_z'][0]]
    [dydt[5],dydt[6],dydt[7]] = FrameFromAngularVelocityIntegrand(rfrak_frame, OmegaVec_1p0(Cons,Vars)[1:])
    dydt[0] = dvdt[0]
    if(Cons['EvolveSpin1']):
        dydt[1], dydt[2]=FrameFromAngularVelocity_2D_Integrand(Vars['rfrak_chi1_x'][0], Vars['rfrak_chi1_y'][0],(mul(mul(inverse(Cons['S_chi1']),OmegaVec_chiVec_1_1p0(Cons,Vars)),Cons['S_chi1']))[1:])
    else:
        dydt[1] = 0.0
        dydt[2] = 0.0
    if(Cons['EvolveSpin2']):
        dydt[3], dydt[4]=FrameFromAngularVelocity_2D_Integrand(Vars['rfrak_chi2_x'][0], Vars['rfrak_chi2_y'][0],(mul(mul(inverse(Cons['S_chi2']),OmegaVec_chiVec_2_1p0(Cons,Vars)),Cons['S_chi2']))[1:])
    else:
        dydt[3] = 0.0
        dydt[4] = 0.0
    dydt[8] = (Vars['v']*Vars['v']*Vars['v']/Cons['M'])[0]
    return dydt      

@njit
def Recalculate_1p5(Cons,y):
    Vars=Dict()
    Vars['v'] = np.array([y[0]])
    Vars['rfrak_chi1_x'] = np.array([y[1]])
    Vars['rfrak_chi1_y'] = np.array([y[2]])
    Vars['rfrak_chi2_x'] = np.array([y[3]])
    Vars['rfrak_chi2_y'] = np.array([y[4]])
    Vars['rfrak_frame_x'] = np.array([y[5]])
    Vars['rfrak_frame_y'] = np.array([y[6]])
    Vars['rfrak_frame_z'] = np.array([y[7]])
    Vars['Phi'] = np.array([y[8]])
    Vars['R'] = exp(Vars['rfrak_frame_x']*Cons['xHat'] + Vars['rfrak_frame_y']*Cons['yHat'] + Vars['rfrak_frame_z']*Cons['zHat'])
    Vars['nHat'] = mul(mul(Vars['R'],Cons['xHat']),conjugate(Vars['R']))
    Vars['lambdaHat'] = mul(mul(Vars['R'],Cons['yHat']),conjugate(Vars['R']))
    Vars['ellHat'] = mul(mul(Vars['R'],Cons['zHat']),conjugate(Vars['R']))
    Vars['R_S1'] = exp(Vars['rfrak_chi1_x']*Cons['xHat'] + Vars['rfrak_chi1_y']*Cons['yHat'])
    Vars['R_S2'] = exp(Vars['rfrak_chi2_x']*Cons['xHat'] + Vars['rfrak_chi2_y']*Cons['yHat'])
    Vars['chiVec1'] = mul(mul(mul(Cons['S_chi1'],Vars['R_S1']),Cons['zHat']),mul(conjugate(Vars['R_S1']),conjugate(Cons['S_chi1'])))
    Vars['chiVec2'] = mul(mul(mul(Cons['S_chi2'],Vars['R_S2']),Cons['zHat']),mul(conjugate(Vars['R_S2']),conjugate(Cons['S_chi2'])))
    Vars['chi1_n'] = np.array([dot(Vars['chiVec1'][1:],Vars['nHat'][1:])])
    Vars['chi1_lambda'] = np.array([dot(Vars['chiVec1'][1:],Vars['lambdaHat'][1:])])
    Vars['chi1_ell'] = np.array([dot(Vars['chiVec1'][1:],Vars['ellHat'][1:])])
    Vars['chi2_n'] = np.array([dot(Vars['chiVec2'][1:],Vars['nHat'][1:])])
    Vars['chi2_lambda'] = np.array([dot(Vars['chiVec2'][1:],Vars['lambdaHat'][1:])])
    Vars['chi2_ell'] = np.array([dot(Vars['chiVec2'][1:],Vars['ellHat'][1:])])
    Vars['S_ell'] = Cons['M1']**2*Vars['chi1_ell'] + Cons['M2']**2*Vars['chi2_ell']
    Vars['S_n'] = Cons['M1']**2*Vars['chi1_n'] + Cons['M2']**2*Vars['chi2_n']
    Vars['S_lambda'] = Cons['M1']**2*Vars['chi1_lambda'] + Cons['M2']**2*Vars['chi2_lambda']
    Vars['Sigma_ell'] = Cons['M']*(-Cons['M1']*Vars['chi1_ell'] + Cons['M2']*Vars['chi2_ell'])
    Vars['Sigma_n'] = Cons['M']*(-Cons['M1']*Vars['chi1_n'] + Cons['M2']*Vars['chi2_n'])
    Vars['Sigma_lambda'] = Cons['M']*(-Cons['M1']*Vars['chi1_lambda'] + Cons['M2']*Vars['chi2_lambda'])
    Vars['chi_s_ell'] = Vars['chi1_ell']/2 + Vars['chi2_ell']/2
    Vars['chi_a_ell'] = Vars['chi1_ell']/2 - Vars['chi2_ell']/2
    Vars['Fcal_coeff'] = 32*Cons['nu']**2*Vars['v']**10/5
    Vars['Fcal_SO_3'] = (-4*Vars['S_ell'] - 5*Vars['Sigma_ell']*Cons['delta']/4)/Cons['M']**2
    Vars['E_SO_3'] = (14*Vars['S_ell']/3 + 2*Vars['Sigma_ell']*Cons['delta'])/Cons['M']**2
    return Vars

@njit
def OmegaVec_chiVec_1_1p5(Cons,Vars):
    Omega1_coeff = Vars['v']**5/Cons['M']
    return Omega1_coeff*(Vars['ellHat']*(-0.75*Cons['delta'] + 0.5*Cons['nu'] + Vars['v']**2*(Cons['delta']*(0.625*Cons['nu'] - 0.5625) + Cons['nu']*(1.25 - 0.0416666666666667*Cons['nu']) + 0.5625) + 0.75) + Vars['nHat']*Vars['v']*(3.0*Vars['chi1_n']*Cons['nu'] + 3.0*Cons['M2']**2*Vars['chi2_n']/Cons['M']**2) - Cons['M2']**2*Vars['chiVec2']*Vars['v']/Cons['M']**2)

@njit
def OmegaVec_chiVec_2_1p5(Cons,Vars):
    Omega2_coeff = Vars['v']**5/Cons['M']
    return Omega2_coeff*(Vars['ellHat']*(0.75*Cons['delta'] + 0.5*Cons['nu'] + Vars['v']**2*(Cons['delta']*(0.5625 - 0.625*Cons['nu']) + Cons['nu']*(1.25 - 0.0416666666666667*Cons['nu']) + 0.5625) + 0.75) + Vars['nHat']*Vars['v']*(3.0*Vars['chi2_n']*Cons['nu'] + 3.0*Cons['M1']**2*Vars['chi1_n']/Cons['M']**2) - Cons['M1']**2*Vars['chiVec1']*Vars['v']/Cons['M']**2)

@njit
def OmegaVec_1p5(Cons,Vars):
    gamma_PN_3 = (1.66666666666667*Vars['S_ell'] + Vars['Sigma_ell']*Cons['delta'])/Cons['M']**2
    a_ell_0 = 7.0*Vars['S_n'] + 3.0*Vars['Sigma_n']*Cons['delta']
    gamma_PN_2 = 1.0 - 0.333333333333333*Cons['nu']
    gamma_PN_0 = 1.00000000000000
    a_ell_2 = Vars['S_n']*(-9.66666666666667*Cons['nu'] - 10.0) + Vars['Sigma_n']*Cons['delta']*(-4.5*Cons['nu'] - 6.0)
    return Vars['ellHat']*Vars['v']**3/Cons['M'] + Vars['nHat']*Vars['v']**6*(a_ell_0 + a_ell_2*Vars['v']**2)*(gamma_PN_0 + Vars['v']**2*(gamma_PN_2 + gamma_PN_3*Vars['v']))/Cons['M']**3


@njit
def TaylorT1_1p5(Cons,Vars):
    Flux = Vars['Fcal_coeff']*(Cons['Fcal_0'] + Vars['v']**2*(Cons['Fcal_2'] + Vars['v']*(Cons['Fcal_3'] + Vars['Fcal_SO_3'])))
    dEdV = -0.5*Cons['M']*Cons['nu']*Vars['v']*(2.0*Cons['E_0'] + Vars['v']**2*(4.0*Cons['E_2'] + 5.0*Vars['E_SO_3']*Vars['v']))
    Absorption = 0
    dvdt = (-Absorption - Flux)/dEdV
    dydt=np.zeros(9)
    rfrak_frame=[Vars['rfrak_frame_x'][0],Vars['rfrak_frame_y'][0],Vars['rfrak_frame_z'][0]]
    [dydt[5],dydt[6],dydt[7]] = FrameFromAngularVelocityIntegrand(rfrak_frame, OmegaVec_1p5(Cons,Vars)[1:])
    dydt[0] = dvdt[0]
    if(Cons['EvolveSpin1']):
        dydt[1], dydt[2]=FrameFromAngularVelocity_2D_Integrand(Vars['rfrak_chi1_x'][0], Vars['rfrak_chi1_y'][0],(mul(mul(inverse(Cons['S_chi1']),OmegaVec_chiVec_1_1p5(Cons,Vars)),Cons['S_chi1']))[1:])
    else:
        dydt[1] = 0.0
        dydt[2] = 0.0
    if(Cons['EvolveSpin2']):
        dydt[3], dydt[4]=FrameFromAngularVelocity_2D_Integrand(Vars['rfrak_chi2_x'][0], Vars['rfrak_chi2_y'][0],(mul(mul(inverse(Cons['S_chi2']),OmegaVec_chiVec_2_1p5(Cons,Vars)),Cons['S_chi2']))[1:])
    else:
        dydt[3] = 0.0
        dydt[4] = 0.0
    dydt[8] = (Vars['v']*Vars['v']*Vars['v']/Cons['M'])[0]
    return dydt      

@njit
def Recalculate_2p0(Cons,y):
    Vars=Dict()
    Vars['v'] = np.array([y[0]])
    Vars['rfrak_chi1_x'] = np.array([y[1]])
    Vars['rfrak_chi1_y'] = np.array([y[2]])
    Vars['rfrak_chi2_x'] = np.array([y[3]])
    Vars['rfrak_chi2_y'] = np.array([y[4]])
    Vars['rfrak_frame_x'] = np.array([y[5]])
    Vars['rfrak_frame_y'] = np.array([y[6]])
    Vars['rfrak_frame_z'] = np.array([y[7]])
    Vars['Phi'] = np.array([y[8]])
    Vars['R'] = exp(Vars['rfrak_frame_x']*Cons['xHat'] + Vars['rfrak_frame_y']*Cons['yHat'] + Vars['rfrak_frame_z']*Cons['zHat'])
    Vars['nHat'] = mul(mul(Vars['R'],Cons['xHat']),conjugate(Vars['R']))
    Vars['lambdaHat'] = mul(mul(Vars['R'],Cons['yHat']),conjugate(Vars['R']))
    Vars['ellHat'] = mul(mul(Vars['R'],Cons['zHat']),conjugate(Vars['R']))
    Vars['R_S1'] = exp(Vars['rfrak_chi1_x']*Cons['xHat'] + Vars['rfrak_chi1_y']*Cons['yHat'])
    Vars['R_S2'] = exp(Vars['rfrak_chi2_x']*Cons['xHat'] + Vars['rfrak_chi2_y']*Cons['yHat'])
    Vars['chiVec1'] = mul(mul(mul(Cons['S_chi1'],Vars['R_S1']),Cons['zHat']),mul(conjugate(Vars['R_S1']),conjugate(Cons['S_chi1'])))
    Vars['chiVec2'] = mul(mul(mul(Cons['S_chi2'],Vars['R_S2']),Cons['zHat']),mul(conjugate(Vars['R_S2']),conjugate(Cons['S_chi2'])))
    Vars['chi1_n'] = np.array([dot(Vars['chiVec1'][1:],Vars['nHat'][1:])])
    Vars['chi1_lambda'] = np.array([dot(Vars['chiVec1'][1:],Vars['lambdaHat'][1:])])
    Vars['chi1_ell'] = np.array([dot(Vars['chiVec1'][1:],Vars['ellHat'][1:])])
    Vars['chi2_n'] = np.array([dot(Vars['chiVec2'][1:],Vars['nHat'][1:])])
    Vars['chi2_lambda'] = np.array([dot(Vars['chiVec2'][1:],Vars['lambdaHat'][1:])])
    Vars['chi2_ell'] = np.array([dot(Vars['chiVec2'][1:],Vars['ellHat'][1:])])
    Vars['S_ell'] = Cons['M1']**2*Vars['chi1_ell'] + Cons['M2']**2*Vars['chi2_ell']
    Vars['S_n'] = Cons['M1']**2*Vars['chi1_n'] + Cons['M2']**2*Vars['chi2_n']
    Vars['S_lambda'] = Cons['M1']**2*Vars['chi1_lambda'] + Cons['M2']**2*Vars['chi2_lambda']
    Vars['Sigma_ell'] = Cons['M']*(-Cons['M1']*Vars['chi1_ell'] + Cons['M2']*Vars['chi2_ell'])
    Vars['Sigma_n'] = Cons['M']*(-Cons['M1']*Vars['chi1_n'] + Cons['M2']*Vars['chi2_n'])
    Vars['Sigma_lambda'] = Cons['M']*(-Cons['M1']*Vars['chi1_lambda'] + Cons['M2']*Vars['chi2_lambda'])
    Vars['chi_s_ell'] = Vars['chi1_ell']/2 + Vars['chi2_ell']/2
    Vars['chi_a_ell'] = Vars['chi1_ell']/2 - Vars['chi2_ell']/2
    Vars['Fcal_coeff'] = 32*Cons['nu']**2*Vars['v']**10/5
    Vars['Fcal_SQ_4'] = Cons['chi1chi1']*(-89*Cons['delta']/192 + 89*Cons['nu']/96 - 89/192) - 103*Cons['chi1chi2']*Cons['nu']/48 + Cons['chi2chi2']*(89*Cons['delta']/192 + 89*Cons['nu']/96 - 89/192) + Vars['chi_a_ell']*(Vars['chi_a_ell']*(287/96 - 12*Cons['nu']) + 287*Vars['chi_s_ell']*Cons['delta']/48) + Vars['chi_s_ell']**2*(Cons['nu']/24 + 287/96)
    Vars['Fcal_SO_3'] = (-4*Vars['S_ell'] - 5*Vars['Sigma_ell']*Cons['delta']/4)/Cons['M']**2
    Vars['E_SQ_4'] = -3*Vars['chi_a_ell']**2/2 - 3*Vars['chi_s_ell']**2/2 - Cons['delta']*(Cons['chi2chi2']/2 + 3*Vars['chi_a_ell']*Vars['chi_s_ell']) + Cons['nu']*(Cons['chi1chi2'] + 6*Vars['chi_a_ell']**2) + (Cons['chi1chi1'] + Cons['chi2chi2'])*(Cons['delta'] - 2*Cons['nu'] + 1)/4
    Vars['E_SO_3'] = (14*Vars['S_ell']/3 + 2*Vars['Sigma_ell']*Cons['delta'])/Cons['M']**2
    return Vars

@njit
def OmegaVec_chiVec_1_2p0(Cons,Vars):
    Omega1_coeff = Vars['v']**5/Cons['M']
    return Omega1_coeff*(Vars['ellHat']*(-0.75*Cons['delta'] + 0.5*Cons['nu'] + Vars['v']**2*(Cons['delta']*(0.625*Cons['nu'] - 0.5625) + Cons['nu']*(1.25 - 0.0416666666666667*Cons['nu']) + Vars['v']**2*(Cons['delta']*(Cons['nu']*(4.875 - 0.15625*Cons['nu']) - 0.84375) + Cons['nu']*(Cons['nu']*(-0.0208333333333333*Cons['nu'] - 3.28125) + 0.1875) + 0.84375) + 0.5625) + 0.75) + Vars['nHat']*Vars['v']*(3.0*Vars['chi1_n']*Cons['nu'] + 3.0*Cons['M2']**2*Vars['chi2_n']/Cons['M']**2) - Cons['M2']**2*Vars['chiVec2']*Vars['v']/Cons['M']**2)

@njit
def OmegaVec_chiVec_2_2p0(Cons,Vars):
    Omega2_coeff = Vars['v']**5/Cons['M']
    return Omega2_coeff*(Vars['ellHat']*(0.75*Cons['delta'] + 0.5*Cons['nu'] + Vars['v']**2*(Cons['delta']*(0.5625 - 0.625*Cons['nu']) + Cons['nu']*(1.25 - 0.0416666666666667*Cons['nu']) + Vars['v']**2*(Cons['delta']*(Cons['nu']*(0.15625*Cons['nu'] - 4.875) + 0.84375) + Cons['nu']*(Cons['nu']*(-0.0208333333333333*Cons['nu'] - 3.28125) + 0.1875) + 0.84375) + 0.5625) + 0.75) + Vars['nHat']*Vars['v']*(3.0*Vars['chi2_n']*Cons['nu'] + 3.0*Cons['M1']**2*Vars['chi1_n']/Cons['M']**2) - Cons['M1']**2*Vars['chiVec1']*Vars['v']/Cons['M']**2)

@njit
def OmegaVec_2p0(Cons,Vars):
    a_ell_0 = 7.0*Vars['S_n'] + 3.0*Vars['Sigma_n']*Cons['delta']
    gamma_PN_3 = (1.66666666666667*Vars['S_ell'] + Vars['Sigma_ell']*Cons['delta'])/Cons['M']**2
    gamma_PN_4 = 1.0 - 5.41666666666667*Cons['nu']
    gamma_PN_2 = 1.0 - 0.333333333333333*Cons['nu']
    gamma_PN_0 = 1.00000000000000
    a_ell_2 = Vars['S_n']*(-9.66666666666667*Cons['nu'] - 10.0) + Vars['Sigma_n']*Cons['delta']*(-4.5*Cons['nu'] - 6.0)
    a_ell_4 = Vars['S_n']*(5.77777777777778*Cons['nu']**2 + 14.75*Cons['nu'] + 1.5) + Vars['Sigma_n']*Cons['delta']*(2.83333333333333*Cons['nu']**2 + 9.125*Cons['nu'] + 1.5)
    return Vars['ellHat']*Vars['v']**3/Cons['M'] + Vars['nHat']*Vars['v']**6*(a_ell_0 + Vars['v']**2*(a_ell_2 + a_ell_4*Vars['v']**2))*(gamma_PN_0 + Vars['v']**2*(gamma_PN_2 + Vars['v']*(gamma_PN_3 + gamma_PN_4*Vars['v'])))/Cons['M']**3


@njit
def TaylorT1_2p0(Cons,Vars):
    Flux = Vars['Fcal_coeff']*(Cons['Fcal_0'] + Vars['v']**2*(Cons['Fcal_2'] + Vars['v']*(Cons['Fcal_3'] + Vars['Fcal_SO_3'] + Vars['v']*(Cons['Fcal_4'] + Vars['Fcal_SQ_4']))))
    dEdV = -0.5*Cons['M']*Cons['nu']*Vars['v']*(2.0*Cons['E_0'] + Vars['v']**2*(4.0*Cons['E_2'] + Vars['v']*(5.0*Vars['E_SO_3'] + 6.0*Vars['v']*(Cons['E_4'] + Vars['E_SQ_4']))))
    Absorption = 0
    dvdt = (-Absorption - Flux)/dEdV
    dydt=np.zeros(9)
    rfrak_frame=[Vars['rfrak_frame_x'][0],Vars['rfrak_frame_y'][0],Vars['rfrak_frame_z'][0]]
    [dydt[5],dydt[6],dydt[7]] = FrameFromAngularVelocityIntegrand(rfrak_frame, OmegaVec_2p0(Cons,Vars)[1:])
    dydt[0] = dvdt[0]
    if(Cons['EvolveSpin1']):
        dydt[1], dydt[2]=FrameFromAngularVelocity_2D_Integrand(Vars['rfrak_chi1_x'][0], Vars['rfrak_chi1_y'][0],(mul(mul(inverse(Cons['S_chi1']),OmegaVec_chiVec_1_2p0(Cons,Vars)),Cons['S_chi1']))[1:])
    else:
        dydt[1] = 0.0
        dydt[2] = 0.0
    if(Cons['EvolveSpin2']):
        dydt[3], dydt[4]=FrameFromAngularVelocity_2D_Integrand(Vars['rfrak_chi2_x'][0], Vars['rfrak_chi2_y'][0],(mul(mul(inverse(Cons['S_chi2']),OmegaVec_chiVec_2_2p0(Cons,Vars)),Cons['S_chi2']))[1:])
    else:
        dydt[3] = 0.0
        dydt[4] = 0.0
    dydt[8] = (Vars['v']*Vars['v']*Vars['v']/Cons['M'])[0]
    return dydt      

@njit
def Recalculate_2p5(Cons,y):
    Vars=Dict()
    Vars['v'] = np.array([y[0]])
    Vars['rfrak_chi1_x'] = np.array([y[1]])
    Vars['rfrak_chi1_y'] = np.array([y[2]])
    Vars['rfrak_chi2_x'] = np.array([y[3]])
    Vars['rfrak_chi2_y'] = np.array([y[4]])
    Vars['rfrak_frame_x'] = np.array([y[5]])
    Vars['rfrak_frame_y'] = np.array([y[6]])
    Vars['rfrak_frame_z'] = np.array([y[7]])
    Vars['Phi'] = np.array([y[8]])
    Vars['R'] = exp(Vars['rfrak_frame_x']*Cons['xHat'] + Vars['rfrak_frame_y']*Cons['yHat'] + Vars['rfrak_frame_z']*Cons['zHat'])
    Vars['nHat'] = mul(mul(Vars['R'],Cons['xHat']),conjugate(Vars['R']))
    Vars['lambdaHat'] = mul(mul(Vars['R'],Cons['yHat']),conjugate(Vars['R']))
    Vars['ellHat'] = mul(mul(Vars['R'],Cons['zHat']),conjugate(Vars['R']))
    Vars['R_S1'] = exp(Vars['rfrak_chi1_x']*Cons['xHat'] + Vars['rfrak_chi1_y']*Cons['yHat'])
    Vars['R_S2'] = exp(Vars['rfrak_chi2_x']*Cons['xHat'] + Vars['rfrak_chi2_y']*Cons['yHat'])
    Vars['chiVec1'] = mul(mul(mul(Cons['S_chi1'],Vars['R_S1']),Cons['zHat']),mul(conjugate(Vars['R_S1']),conjugate(Cons['S_chi1'])))
    Vars['chiVec2'] = mul(mul(mul(Cons['S_chi2'],Vars['R_S2']),Cons['zHat']),mul(conjugate(Vars['R_S2']),conjugate(Cons['S_chi2'])))
    Vars['chi1_n'] = np.array([dot(Vars['chiVec1'][1:],Vars['nHat'][1:])])
    Vars['chi1_lambda'] = np.array([dot(Vars['chiVec1'][1:],Vars['lambdaHat'][1:])])
    Vars['chi1_ell'] = np.array([dot(Vars['chiVec1'][1:],Vars['ellHat'][1:])])
    Vars['chi2_n'] = np.array([dot(Vars['chiVec2'][1:],Vars['nHat'][1:])])
    Vars['chi2_lambda'] = np.array([dot(Vars['chiVec2'][1:],Vars['lambdaHat'][1:])])
    Vars['chi2_ell'] = np.array([dot(Vars['chiVec2'][1:],Vars['ellHat'][1:])])
    Vars['S_ell'] = Cons['M1']**2*Vars['chi1_ell'] + Cons['M2']**2*Vars['chi2_ell']
    Vars['S_n'] = Cons['M1']**2*Vars['chi1_n'] + Cons['M2']**2*Vars['chi2_n']
    Vars['S_lambda'] = Cons['M1']**2*Vars['chi1_lambda'] + Cons['M2']**2*Vars['chi2_lambda']
    Vars['Sigma_ell'] = Cons['M']*(-Cons['M1']*Vars['chi1_ell'] + Cons['M2']*Vars['chi2_ell'])
    Vars['Sigma_n'] = Cons['M']*(-Cons['M1']*Vars['chi1_n'] + Cons['M2']*Vars['chi2_n'])
    Vars['Sigma_lambda'] = Cons['M']*(-Cons['M1']*Vars['chi1_lambda'] + Cons['M2']*Vars['chi2_lambda'])
    Vars['chi_s_ell'] = Vars['chi1_ell']/2 + Vars['chi2_ell']/2
    Vars['chi_a_ell'] = Vars['chi1_ell']/2 - Vars['chi2_ell']/2
    Vars['Fcal_coeff'] = 32*Cons['nu']**2*Vars['v']**10/5
    Vars['Fcal_SQ_4'] = Cons['chi1chi1']*(-89*Cons['delta']/192 + 89*Cons['nu']/96 - 89/192) - 103*Cons['chi1chi2']*Cons['nu']/48 + Cons['chi2chi2']*(89*Cons['delta']/192 + 89*Cons['nu']/96 - 89/192) + Vars['chi_a_ell']*(Vars['chi_a_ell']*(287/96 - 12*Cons['nu']) + 287*Vars['chi_s_ell']*Cons['delta']/48) + Vars['chi_s_ell']**2*(Cons['nu']/24 + 287/96)
    Vars['Fcal_SO_3'] = (-4*Vars['S_ell'] - 5*Vars['Sigma_ell']*Cons['delta']/4)/Cons['M']**2
    Vars['Fcal_SO_5'] = (Vars['S_ell']*(272*Cons['nu']/9 - 9/2) + Vars['Sigma_ell']*Cons['delta']*(43*Cons['nu']/4 - 13/16))/Cons['M']**2
    Vars['E_SQ_4'] = -3*Vars['chi_a_ell']**2/2 - 3*Vars['chi_s_ell']**2/2 - Cons['delta']*(Cons['chi2chi2']/2 + 3*Vars['chi_a_ell']*Vars['chi_s_ell']) + Cons['nu']*(Cons['chi1chi2'] + 6*Vars['chi_a_ell']**2) + (Cons['chi1chi1'] + Cons['chi2chi2'])*(Cons['delta'] - 2*Cons['nu'] + 1)/4
    Vars['E_SO_3'] = (14*Vars['S_ell']/3 + 2*Vars['Sigma_ell']*Cons['delta'])/Cons['M']**2
    Vars['E_SO_5'] = (Vars['S_ell']*(11 - 61*Cons['nu']/9) + Vars['Sigma_ell']*Cons['delta']*(3 - 10*Cons['nu']/3))/Cons['M']**2
    return Vars

@njit
def OmegaVec_chiVec_1_2p5(Cons,Vars):
    Omega1_coeff = Vars['v']**5/Cons['M']
    return Omega1_coeff*(Vars['ellHat']*(-0.75*Cons['delta'] + 0.5*Cons['nu'] + Vars['v']**2*(Cons['delta']*(0.625*Cons['nu'] - 0.5625) + Cons['nu']*(1.25 - 0.0416666666666667*Cons['nu']) + Vars['v']**2*(Cons['delta']*(Cons['nu']*(4.875 - 0.15625*Cons['nu']) - 0.84375) + Cons['nu']*(Cons['nu']*(-0.0208333333333333*Cons['nu'] - 3.28125) + 0.1875) + 0.84375) + 0.5625) + 0.75) + Vars['nHat']*Vars['v']*(3.0*Vars['chi1_n']*Cons['nu'] + 3.0*Cons['M2']**2*Vars['chi2_n']/Cons['M']**2) - Cons['M2']**2*Vars['chiVec2']*Vars['v']/Cons['M']**2)

@njit
def OmegaVec_chiVec_2_2p5(Cons,Vars):
    Omega2_coeff = Vars['v']**5/Cons['M']
    return Omega2_coeff*(Vars['ellHat']*(0.75*Cons['delta'] + 0.5*Cons['nu'] + Vars['v']**2*(Cons['delta']*(0.5625 - 0.625*Cons['nu']) + Cons['nu']*(1.25 - 0.0416666666666667*Cons['nu']) + Vars['v']**2*(Cons['delta']*(Cons['nu']*(0.15625*Cons['nu'] - 4.875) + 0.84375) + Cons['nu']*(Cons['nu']*(-0.0208333333333333*Cons['nu'] - 3.28125) + 0.1875) + 0.84375) + 0.5625) + 0.75) + Vars['nHat']*Vars['v']*(3.0*Vars['chi2_n']*Cons['nu'] + 3.0*Cons['M1']**2*Vars['chi1_n']/Cons['M']**2) - Cons['M1']**2*Vars['chiVec1']*Vars['v']/Cons['M']**2)

@njit
def OmegaVec_2p5(Cons,Vars):
    gamma_PN_5 = (Vars['S_ell']*(0.888888888888889*Cons['nu'] + 3.33333333333333) + 2.0*Vars['Sigma_ell']*Cons['delta'])/Cons['M']**2
    a_ell_0 = 7.0*Vars['S_n'] + 3.0*Vars['Sigma_n']*Cons['delta']
    gamma_PN_3 = (1.66666666666667*Vars['S_ell'] + Vars['Sigma_ell']*Cons['delta'])/Cons['M']**2
    gamma_PN_4 = 1.0 - 5.41666666666667*Cons['nu']
    gamma_PN_2 = 1.0 - 0.333333333333333*Cons['nu']
    gamma_PN_0 = 1.00000000000000
    a_ell_2 = Vars['S_n']*(-9.66666666666667*Cons['nu'] - 10.0) + Vars['Sigma_n']*Cons['delta']*(-4.5*Cons['nu'] - 6.0)
    a_ell_4 = Vars['S_n']*(5.77777777777778*Cons['nu']**2 + 14.75*Cons['nu'] + 1.5) + Vars['Sigma_n']*Cons['delta']*(2.83333333333333*Cons['nu']**2 + 9.125*Cons['nu'] + 1.5)
    return Vars['ellHat']*Vars['v']**3/Cons['M'] + Vars['nHat']*Vars['v']**6*(a_ell_0 + Vars['v']**2*(a_ell_2 + a_ell_4*Vars['v']**2))*(gamma_PN_0 + Vars['v']**2*(gamma_PN_2 + Vars['v']*(gamma_PN_3 + Vars['v']*(gamma_PN_4 + gamma_PN_5*Vars['v']))))/Cons['M']**3


@njit
def TaylorT1_2p5(Cons,Vars):
    Flux = Vars['Fcal_coeff']*(Cons['Fcal_0'] + Vars['v']**2*(Cons['Fcal_2'] + Vars['v']*(Cons['Fcal_3'] + Vars['Fcal_SO_3'] + Vars['v']*(Cons['Fcal_4'] + Vars['Fcal_SQ_4'] + Vars['v']*(Cons['Fcal_5'] + Vars['Fcal_SO_5'])))))
    dEdV = -0.5*Cons['M']*Cons['nu']*Vars['v']*(2.0*Cons['E_0'] + Vars['v']**2*(4.0*Cons['E_2'] + Vars['v']*(5.0*Vars['E_SO_3'] + Vars['v']*(6.0*Cons['E_4'] + 7.0*Vars['E_SO_5']*Vars['v'] + 6.0*Vars['E_SQ_4']))))
    Absorption = 0
    dvdt = (-Absorption - Flux)/dEdV
    dydt=np.zeros(9)
    rfrak_frame=[Vars['rfrak_frame_x'][0],Vars['rfrak_frame_y'][0],Vars['rfrak_frame_z'][0]]
    [dydt[5],dydt[6],dydt[7]] = FrameFromAngularVelocityIntegrand(rfrak_frame, OmegaVec_2p5(Cons,Vars)[1:])
    dydt[0] = dvdt[0]
    if(Cons['EvolveSpin1']):
        dydt[1], dydt[2]=FrameFromAngularVelocity_2D_Integrand(Vars['rfrak_chi1_x'][0], Vars['rfrak_chi1_y'][0],(mul(mul(inverse(Cons['S_chi1']),OmegaVec_chiVec_1_2p5(Cons,Vars)),Cons['S_chi1']))[1:])
    else:
        dydt[1] = 0.0
        dydt[2] = 0.0
    if(Cons['EvolveSpin2']):
        dydt[3], dydt[4]=FrameFromAngularVelocity_2D_Integrand(Vars['rfrak_chi2_x'][0], Vars['rfrak_chi2_y'][0],(mul(mul(inverse(Cons['S_chi2']),OmegaVec_chiVec_2_2p5(Cons,Vars)),Cons['S_chi2']))[1:])
    else:
        dydt[3] = 0.0
        dydt[4] = 0.0
    dydt[8] = (Vars['v']*Vars['v']*Vars['v']/Cons['M'])[0]
    return dydt      

@njit
def Recalculate_3p0(Cons,y):
    Vars=Dict()
    Vars['v'] = np.array([y[0]])
    Vars['rfrak_chi1_x'] = np.array([y[1]])
    Vars['rfrak_chi1_y'] = np.array([y[2]])
    Vars['rfrak_chi2_x'] = np.array([y[3]])
    Vars['rfrak_chi2_y'] = np.array([y[4]])
    Vars['rfrak_frame_x'] = np.array([y[5]])
    Vars['rfrak_frame_y'] = np.array([y[6]])
    Vars['rfrak_frame_z'] = np.array([y[7]])
    Vars['Phi'] = np.array([y[8]])
    Vars['R'] = exp(Vars['rfrak_frame_x']*Cons['xHat'] + Vars['rfrak_frame_y']*Cons['yHat'] + Vars['rfrak_frame_z']*Cons['zHat'])
    Vars['nHat'] = mul(mul(Vars['R'],Cons['xHat']),conjugate(Vars['R']))
    Vars['lambdaHat'] = mul(mul(Vars['R'],Cons['yHat']),conjugate(Vars['R']))
    Vars['ellHat'] = mul(mul(Vars['R'],Cons['zHat']),conjugate(Vars['R']))
    Vars['R_S1'] = exp(Vars['rfrak_chi1_x']*Cons['xHat'] + Vars['rfrak_chi1_y']*Cons['yHat'])
    Vars['R_S2'] = exp(Vars['rfrak_chi2_x']*Cons['xHat'] + Vars['rfrak_chi2_y']*Cons['yHat'])
    Vars['chiVec1'] = mul(mul(mul(Cons['S_chi1'],Vars['R_S1']),Cons['zHat']),mul(conjugate(Vars['R_S1']),conjugate(Cons['S_chi1'])))
    Vars['chiVec2'] = mul(mul(mul(Cons['S_chi2'],Vars['R_S2']),Cons['zHat']),mul(conjugate(Vars['R_S2']),conjugate(Cons['S_chi2'])))
    Vars['chi1_n'] = np.array([dot(Vars['chiVec1'][1:],Vars['nHat'][1:])])
    Vars['chi1_lambda'] = np.array([dot(Vars['chiVec1'][1:],Vars['lambdaHat'][1:])])
    Vars['chi1_ell'] = np.array([dot(Vars['chiVec1'][1:],Vars['ellHat'][1:])])
    Vars['chi2_n'] = np.array([dot(Vars['chiVec2'][1:],Vars['nHat'][1:])])
    Vars['chi2_lambda'] = np.array([dot(Vars['chiVec2'][1:],Vars['lambdaHat'][1:])])
    Vars['chi2_ell'] = np.array([dot(Vars['chiVec2'][1:],Vars['ellHat'][1:])])
    Vars['S_ell'] = Cons['M1']**2*Vars['chi1_ell'] + Cons['M2']**2*Vars['chi2_ell']
    Vars['S_n'] = Cons['M1']**2*Vars['chi1_n'] + Cons['M2']**2*Vars['chi2_n']
    Vars['S_lambda'] = Cons['M1']**2*Vars['chi1_lambda'] + Cons['M2']**2*Vars['chi2_lambda']
    Vars['Sigma_ell'] = Cons['M']*(-Cons['M1']*Vars['chi1_ell'] + Cons['M2']*Vars['chi2_ell'])
    Vars['Sigma_n'] = Cons['M']*(-Cons['M1']*Vars['chi1_n'] + Cons['M2']*Vars['chi2_n'])
    Vars['Sigma_lambda'] = Cons['M']*(-Cons['M1']*Vars['chi1_lambda'] + Cons['M2']*Vars['chi2_lambda'])
    Vars['chi_s_ell'] = Vars['chi1_ell']/2 + Vars['chi2_ell']/2
    Vars['chi_a_ell'] = Vars['chi1_ell']/2 - Vars['chi2_ell']/2
    Vars['logv'] = log(Vars['v'])
    Vars['Fcal_coeff'] = 32*Cons['nu']**2*Vars['v']**10/5
    Vars['Fcal_SQ_4'] = Cons['chi1chi1']*(-89*Cons['delta']/192 + 89*Cons['nu']/96 - 89/192) - 103*Cons['chi1chi2']*Cons['nu']/48 + Cons['chi2chi2']*(89*Cons['delta']/192 + 89*Cons['nu']/96 - 89/192) + Vars['chi_a_ell']*(Vars['chi_a_ell']*(287/96 - 12*Cons['nu']) + 287*Vars['chi_s_ell']*Cons['delta']/48) + Vars['chi_s_ell']**2*(Cons['nu']/24 + 287/96)
    Vars['Fcal_SO_3'] = (-4*Vars['S_ell'] - 5*Vars['Sigma_ell']*Cons['delta']/4)/Cons['M']**2
    Vars['Fcal_SO_5'] = (Vars['S_ell']*(272*Cons['nu']/9 - 9/2) + Vars['Sigma_ell']*Cons['delta']*(43*Cons['nu']/4 - 13/16))/Cons['M']**2
    Vars['Fcal_SO_6'] = (-16*Vars['S_ell']*pi - 31*Vars['Sigma_ell']*Cons['delta']*pi/6)/Cons['M']**2
    Vars['E_SQ_4'] = -3*Vars['chi_a_ell']**2/2 - 3*Vars['chi_s_ell']**2/2 - Cons['delta']*(Cons['chi2chi2']/2 + 3*Vars['chi_a_ell']*Vars['chi_s_ell']) + Cons['nu']*(Cons['chi1chi2'] + 6*Vars['chi_a_ell']**2) + (Cons['chi1chi1'] + Cons['chi2chi2'])*(Cons['delta'] - 2*Cons['nu'] + 1)/4
    Vars['E_SO_3'] = (14*Vars['S_ell']/3 + 2*Vars['Sigma_ell']*Cons['delta'])/Cons['M']**2
    Vars['E_SO_5'] = (Vars['S_ell']*(11 - 61*Cons['nu']/9) + Vars['Sigma_ell']*Cons['delta']*(3 - 10*Cons['nu']/3))/Cons['M']**2
    return Vars

@njit
def OmegaVec_chiVec_1_3p0(Cons,Vars):
    Omega1_coeff = Vars['v']**5/Cons['M']
    return Omega1_coeff*(Vars['ellHat']*(-0.75*Cons['delta'] + 0.5*Cons['nu'] + Vars['v']**2*(Cons['delta']*(0.625*Cons['nu'] - 0.5625) + Cons['nu']*(1.25 - 0.0416666666666667*Cons['nu']) + Vars['v']**2*(Cons['delta']*(Cons['nu']*(4.875 - 0.15625*Cons['nu']) - 0.84375) + Cons['nu']*(Cons['nu']*(-0.0208333333333333*Cons['nu'] - 3.28125) + 0.1875) + 0.84375) + 0.5625) + 0.75) + Vars['nHat']*Vars['v']*(3.0*Vars['chi1_n']*Cons['nu'] + 3.0*Cons['M2']**2*Vars['chi2_n']/Cons['M']**2) - Cons['M2']**2*Vars['chiVec2']*Vars['v']/Cons['M']**2)

@njit
def OmegaVec_chiVec_2_3p0(Cons,Vars):
    Omega2_coeff = Vars['v']**5/Cons['M']
    return Omega2_coeff*(Vars['ellHat']*(0.75*Cons['delta'] + 0.5*Cons['nu'] + Vars['v']**2*(Cons['delta']*(0.5625 - 0.625*Cons['nu']) + Cons['nu']*(1.25 - 0.0416666666666667*Cons['nu']) + Vars['v']**2*(Cons['delta']*(Cons['nu']*(0.15625*Cons['nu'] - 4.875) + 0.84375) + Cons['nu']*(Cons['nu']*(-0.0208333333333333*Cons['nu'] - 3.28125) + 0.1875) + 0.84375) + 0.5625) + 0.75) + Vars['nHat']*Vars['v']*(3.0*Vars['chi2_n']*Cons['nu'] + 3.0*Cons['M1']**2*Vars['chi1_n']/Cons['M']**2) - Cons['M1']**2*Vars['chiVec1']*Vars['v']/Cons['M']**2)

@njit
def OmegaVec_3p0(Cons,Vars):
    gamma_PN_5 = (Vars['S_ell']*(0.888888888888889*Cons['nu'] + 3.33333333333333) + 2.0*Vars['Sigma_ell']*Cons['delta'])/Cons['M']**2
    a_ell_0 = 7.0*Vars['S_n'] + 3.0*Vars['Sigma_n']*Cons['delta']
    gamma_PN_3 = (1.66666666666667*Vars['S_ell'] + Vars['Sigma_ell']*Cons['delta'])/Cons['M']**2
    gamma_PN_4 = 1.0 - 5.41666666666667*Cons['nu']
    gamma_PN_2 = 1.0 - 0.333333333333333*Cons['nu']
    gamma_PN_0 = 1.00000000000000
    a_ell_2 = Vars['S_n']*(-9.66666666666667*Cons['nu'] - 10.0) + Vars['Sigma_n']*Cons['delta']*(-4.5*Cons['nu'] - 6.0)
    gamma_PN_6 = 0.0123456790123457*Cons['nu']**3 + 6.36111111111111*Cons['nu']**2 - 2.98177812235564*Cons['nu'] + 1.0
    a_ell_4 = Vars['S_n']*(5.77777777777778*Cons['nu']**2 + 14.75*Cons['nu'] + 1.5) + Vars['Sigma_n']*Cons['delta']*(2.83333333333333*Cons['nu']**2 + 9.125*Cons['nu'] + 1.5)
    return Vars['ellHat']*Vars['v']**3/Cons['M'] + Vars['nHat']*Vars['v']**6*(a_ell_0 + Vars['v']**2*(a_ell_2 + a_ell_4*Vars['v']**2))*(gamma_PN_0 + Vars['v']**2*(gamma_PN_2 + Vars['v']*(gamma_PN_3 + Vars['v']*(gamma_PN_4 + Vars['v']*(gamma_PN_5 + gamma_PN_6*Vars['v'])))))/Cons['M']**3


@njit
def TaylorT1_3p0(Cons,Vars):
    Flux = Vars['Fcal_coeff']*(Cons['Fcal_0'] + Vars['v']**2*(Cons['Fcal_2'] + Vars['v']*(Cons['Fcal_3'] + Vars['Fcal_SO_3'] + Vars['v']*(Cons['Fcal_4'] + Vars['Fcal_SQ_4'] + Vars['v']*(Cons['Fcal_5'] + Vars['Fcal_SO_5'] + Vars['v']*(Cons['Fcal_6'] + Vars['Fcal_SO_6'] + Cons['Fcal_lnv_6']*Vars['logv']))))))
    dEdV = -0.5*Cons['M']*Cons['nu']*Vars['v']*(2.0*Cons['E_0'] + Vars['v']**2*(4.0*Cons['E_2'] + Vars['v']*(5.0*Vars['E_SO_3'] + Vars['v']*(6.0*Cons['E_4'] + 6.0*Vars['E_SQ_4'] + Vars['v']*(8.0*Cons['E_6']*Vars['v'] + 7.0*Vars['E_SO_5'])))))
    Absorption = 0
    dvdt = (-Absorption - Flux)/dEdV
    dydt=np.zeros(9)
    rfrak_frame=[Vars['rfrak_frame_x'][0],Vars['rfrak_frame_y'][0],Vars['rfrak_frame_z'][0]]
    [dydt[5],dydt[6],dydt[7]] = FrameFromAngularVelocityIntegrand(rfrak_frame, OmegaVec_3p0(Cons,Vars)[1:])
    dydt[0] = dvdt[0]
    if(Cons['EvolveSpin1']):
        dydt[1], dydt[2]=FrameFromAngularVelocity_2D_Integrand(Vars['rfrak_chi1_x'][0], Vars['rfrak_chi1_y'][0],(mul(mul(inverse(Cons['S_chi1']),OmegaVec_chiVec_1_3p0(Cons,Vars)),Cons['S_chi1']))[1:])
    else:
        dydt[1] = 0.0
        dydt[2] = 0.0
    if(Cons['EvolveSpin2']):
        dydt[3], dydt[4]=FrameFromAngularVelocity_2D_Integrand(Vars['rfrak_chi2_x'][0], Vars['rfrak_chi2_y'][0],(mul(mul(inverse(Cons['S_chi2']),OmegaVec_chiVec_2_3p0(Cons,Vars)),Cons['S_chi2']))[1:])
    else:
        dydt[3] = 0.0
        dydt[4] = 0.0
    dydt[8] = (Vars['v']*Vars['v']*Vars['v']/Cons['M'])[0]
    return dydt      

@njit
def Recalculate_3p5(Cons,y):
    Vars=Dict()
    Vars['v'] = np.array([y[0]])
    Vars['rfrak_chi1_x'] = np.array([y[1]])
    Vars['rfrak_chi1_y'] = np.array([y[2]])
    Vars['rfrak_chi2_x'] = np.array([y[3]])
    Vars['rfrak_chi2_y'] = np.array([y[4]])
    Vars['rfrak_frame_x'] = np.array([y[5]])
    Vars['rfrak_frame_y'] = np.array([y[6]])
    Vars['rfrak_frame_z'] = np.array([y[7]])
    Vars['Phi'] = np.array([y[8]])
    Vars['R'] = exp(Vars['rfrak_frame_x']*Cons['xHat'] + Vars['rfrak_frame_y']*Cons['yHat'] + Vars['rfrak_frame_z']*Cons['zHat'])
    Vars['nHat'] = mul(mul(Vars['R'],Cons['xHat']),conjugate(Vars['R']))
    Vars['lambdaHat'] = mul(mul(Vars['R'],Cons['yHat']),conjugate(Vars['R']))
    Vars['ellHat'] = mul(mul(Vars['R'],Cons['zHat']),conjugate(Vars['R']))
    Vars['R_S1'] = exp(Vars['rfrak_chi1_x']*Cons['xHat'] + Vars['rfrak_chi1_y']*Cons['yHat'])
    Vars['R_S2'] = exp(Vars['rfrak_chi2_x']*Cons['xHat'] + Vars['rfrak_chi2_y']*Cons['yHat'])
    Vars['chiVec1'] = mul(mul(mul(Cons['S_chi1'],Vars['R_S1']),Cons['zHat']),mul(conjugate(Vars['R_S1']),conjugate(Cons['S_chi1'])))
    Vars['chiVec2'] = mul(mul(mul(Cons['S_chi2'],Vars['R_S2']),Cons['zHat']),mul(conjugate(Vars['R_S2']),conjugate(Cons['S_chi2'])))
    Vars['chi1_n'] = np.array([dot(Vars['chiVec1'][1:],Vars['nHat'][1:])])
    Vars['chi1_lambda'] = np.array([dot(Vars['chiVec1'][1:],Vars['lambdaHat'][1:])])
    Vars['chi1_ell'] = np.array([dot(Vars['chiVec1'][1:],Vars['ellHat'][1:])])
    Vars['chi2_n'] = np.array([dot(Vars['chiVec2'][1:],Vars['nHat'][1:])])
    Vars['chi2_lambda'] = np.array([dot(Vars['chiVec2'][1:],Vars['lambdaHat'][1:])])
    Vars['chi2_ell'] = np.array([dot(Vars['chiVec2'][1:],Vars['ellHat'][1:])])
    Vars['S_ell'] = Cons['M1']**2*Vars['chi1_ell'] + Cons['M2']**2*Vars['chi2_ell']
    Vars['S_n'] = Cons['M1']**2*Vars['chi1_n'] + Cons['M2']**2*Vars['chi2_n']
    Vars['S_lambda'] = Cons['M1']**2*Vars['chi1_lambda'] + Cons['M2']**2*Vars['chi2_lambda']
    Vars['Sigma_ell'] = Cons['M']*(-Cons['M1']*Vars['chi1_ell'] + Cons['M2']*Vars['chi2_ell'])
    Vars['Sigma_n'] = Cons['M']*(-Cons['M1']*Vars['chi1_n'] + Cons['M2']*Vars['chi2_n'])
    Vars['Sigma_lambda'] = Cons['M']*(-Cons['M1']*Vars['chi1_lambda'] + Cons['M2']*Vars['chi2_lambda'])
    Vars['chi_s_ell'] = Vars['chi1_ell']/2 + Vars['chi2_ell']/2
    Vars['chi_a_ell'] = Vars['chi1_ell']/2 - Vars['chi2_ell']/2
    Vars['logv'] = log(Vars['v'])
    Vars['Fcal_coeff'] = 32*Cons['nu']**2*Vars['v']**10/5
    Vars['Fcal_SQ_4'] = Cons['chi1chi1']*(-89*Cons['delta']/192 + 89*Cons['nu']/96 - 89/192) - 103*Cons['chi1chi2']*Cons['nu']/48 + Cons['chi2chi2']*(89*Cons['delta']/192 + 89*Cons['nu']/96 - 89/192) + Vars['chi_a_ell']*(Vars['chi_a_ell']*(287/96 - 12*Cons['nu']) + 287*Vars['chi_s_ell']*Cons['delta']/48) + Vars['chi_s_ell']**2*(Cons['nu']/24 + 287/96)
    Vars['Fcal_SO_3'] = (-4*Vars['S_ell'] - 5*Vars['Sigma_ell']*Cons['delta']/4)/Cons['M']**2
    Vars['Fcal_SO_5'] = (Vars['S_ell']*(272*Cons['nu']/9 - 9/2) + Vars['Sigma_ell']*Cons['delta']*(43*Cons['nu']/4 - 13/16))/Cons['M']**2
    Vars['Fcal_SO_6'] = (-16*Vars['S_ell']*pi - 31*Vars['Sigma_ell']*Cons['delta']*pi/6)/Cons['M']**2
    Vars['Fcal_SO_7'] = (Vars['S_ell']*(-2810*Cons['nu']**2/27 + 6172*Cons['nu']/189 + 476645/6804) + Vars['Sigma_ell']*Cons['delta']*(-1501*Cons['nu']**2/36 + 1849*Cons['nu']/126 + 9535/336))/Cons['M']**2
    Vars['E_SQ_4'] = -3*Vars['chi_a_ell']**2/2 - 3*Vars['chi_s_ell']**2/2 - Cons['delta']*(Cons['chi2chi2']/2 + 3*Vars['chi_a_ell']*Vars['chi_s_ell']) + Cons['nu']*(Cons['chi1chi2'] + 6*Vars['chi_a_ell']**2) + (Cons['chi1chi1'] + Cons['chi2chi2'])*(Cons['delta'] - 2*Cons['nu'] + 1)/4
    Vars['E_SO_3'] = (14*Vars['S_ell']/3 + 2*Vars['Sigma_ell']*Cons['delta'])/Cons['M']**2
    Vars['E_SO_5'] = (Vars['S_ell']*(11 - 61*Cons['nu']/9) + Vars['Sigma_ell']*Cons['delta']*(3 - 10*Cons['nu']/3))/Cons['M']**2
    Vars['E_SO_7'] = (Vars['S_ell']*(29*Cons['nu']**2/12 - 367*Cons['nu']/4 + 135/4) + Vars['Sigma_ell']*Cons['delta']*(5*Cons['nu']**2/4 - 39*Cons['nu'] + 27/4))/Cons['M']**2
    return Vars

@njit
def OmegaVec_chiVec_1_3p5(Cons,Vars):
    Omega1_coeff = Vars['v']**5/Cons['M']
    return Omega1_coeff*(Vars['ellHat']*(-0.75*Cons['delta'] + 0.5*Cons['nu'] + Vars['v']**2*(Cons['delta']*(0.625*Cons['nu'] - 0.5625) + Cons['nu']*(1.25 - 0.0416666666666667*Cons['nu']) + Vars['v']**2*(Cons['delta']*(Cons['nu']*(4.875 - 0.15625*Cons['nu']) - 0.84375) + Cons['nu']*(Cons['nu']*(-0.0208333333333333*Cons['nu'] - 3.28125) + 0.1875) + 0.84375) + 0.5625) + 0.75) + Vars['nHat']*Vars['v']*(3.0*Vars['chi1_n']*Cons['nu'] + 3.0*Cons['M2']**2*Vars['chi2_n']/Cons['M']**2) - Cons['M2']**2*Vars['chiVec2']*Vars['v']/Cons['M']**2)

@njit
def OmegaVec_chiVec_2_3p5(Cons,Vars):
    Omega2_coeff = Vars['v']**5/Cons['M']
    return Omega2_coeff*(Vars['ellHat']*(0.75*Cons['delta'] + 0.5*Cons['nu'] + Vars['v']**2*(Cons['delta']*(0.5625 - 0.625*Cons['nu']) + Cons['nu']*(1.25 - 0.0416666666666667*Cons['nu']) + Vars['v']**2*(Cons['delta']*(Cons['nu']*(0.15625*Cons['nu'] - 4.875) + 0.84375) + Cons['nu']*(Cons['nu']*(-0.0208333333333333*Cons['nu'] - 3.28125) + 0.1875) + 0.84375) + 0.5625) + 0.75) + Vars['nHat']*Vars['v']*(3.0*Vars['chi2_n']*Cons['nu'] + 3.0*Cons['M1']**2*Vars['chi1_n']/Cons['M']**2) - Cons['M1']**2*Vars['chiVec1']*Vars['v']/Cons['M']**2)

@njit
def OmegaVec_3p5(Cons,Vars):
    gamma_PN_5 = (Vars['S_ell']*(0.888888888888889*Cons['nu'] + 3.33333333333333) + 2.0*Vars['Sigma_ell']*Cons['delta'])/Cons['M']**2
    a_ell_0 = 7.0*Vars['S_n'] + 3.0*Vars['Sigma_n']*Cons['delta']
    gamma_PN_3 = (1.66666666666667*Vars['S_ell'] + Vars['Sigma_ell']*Cons['delta'])/Cons['M']**2
    gamma_PN_4 = 1.0 - 5.41666666666667*Cons['nu']
    gamma_PN_2 = 1.0 - 0.333333333333333*Cons['nu']
    gamma_PN_0 = 1.00000000000000
    a_ell_2 = Vars['S_n']*(-9.66666666666667*Cons['nu'] - 10.0) + Vars['Sigma_n']*Cons['delta']*(-4.5*Cons['nu'] - 6.0)
    gamma_PN_7 = (Vars['S_ell']*(-6.0*Cons['nu']**2 - 10.5833333333333*Cons['nu'] + 5.0) - 2.66666666666667*Vars['Sigma_ell']*Cons['delta']*Cons['nu']**2 + Vars['Sigma_ell']*Cons['delta']*(3.0 - 10.1666666666667*Cons['nu']))/Cons['M']**2
    gamma_PN_6 = 0.0123456790123457*Cons['nu']**3 + 6.36111111111111*Cons['nu']**2 - 2.98177812235564*Cons['nu'] + 1.0
    a_ell_4 = Vars['S_n']*(5.77777777777778*Cons['nu']**2 + 14.75*Cons['nu'] + 1.5) + Vars['Sigma_n']*Cons['delta']*(2.83333333333333*Cons['nu']**2 + 9.125*Cons['nu'] + 1.5)
    return Vars['ellHat']*Vars['v']**3/Cons['M'] + Vars['nHat']*Vars['v']**6*(a_ell_0 + Vars['v']**2*(a_ell_2 + a_ell_4*Vars['v']**2))*(gamma_PN_0 + Vars['v']**2*(gamma_PN_2 + Vars['v']*(gamma_PN_3 + Vars['v']*(gamma_PN_4 + Vars['v']*(gamma_PN_5 + Vars['v']*(gamma_PN_6 + gamma_PN_7*Vars['v']))))))/Cons['M']**3


@njit
def TaylorT1_3p5(Cons,Vars):
    Flux = Vars['Fcal_coeff']*(Cons['Fcal_0'] + Vars['v']**2*(Cons['Fcal_2'] + Vars['v']*(Cons['Fcal_3'] + Vars['Fcal_SO_3'] + Vars['v']*(Cons['Fcal_4'] + Vars['Fcal_SQ_4'] + Vars['v']*(Cons['Fcal_5'] + Vars['Fcal_SO_5'] + Vars['v']*(Cons['Fcal_6'] + Vars['Fcal_SO_6'] + Cons['Fcal_lnv_6']*Vars['logv'] + Vars['v']*(Cons['Fcal_7'] + Vars['Fcal_SO_7'])))))))
    dEdV = -0.5*Cons['M']*Cons['nu']*Vars['v']*(2.0*Cons['E_0'] + Vars['v']**2*(4.0*Cons['E_2'] + Vars['v']*(5.0*Vars['E_SO_3'] + Vars['v']*(6.0*Cons['E_4'] + 6.0*Vars['E_SQ_4'] + Vars['v']*(7.0*Vars['E_SO_5'] + Vars['v']*(8.0*Cons['E_6'] + 9.0*Vars['E_SO_7']*Vars['v']))))))
    Absorption = 0
    dvdt = (-Absorption - Flux)/dEdV
    dydt=np.zeros(9)
    rfrak_frame=[Vars['rfrak_frame_x'][0],Vars['rfrak_frame_y'][0],Vars['rfrak_frame_z'][0]]
    [dydt[5],dydt[6],dydt[7]] = FrameFromAngularVelocityIntegrand(rfrak_frame, OmegaVec_3p5(Cons,Vars)[1:])
    dydt[0] = dvdt[0]
    if(Cons['EvolveSpin1']):
        dydt[1], dydt[2]=FrameFromAngularVelocity_2D_Integrand(Vars['rfrak_chi1_x'][0], Vars['rfrak_chi1_y'][0],(mul(mul(inverse(Cons['S_chi1']),OmegaVec_chiVec_1_3p5(Cons,Vars)),Cons['S_chi1']))[1:])
    else:
        dydt[1] = 0.0
        dydt[2] = 0.0
    if(Cons['EvolveSpin2']):
        dydt[3], dydt[4]=FrameFromAngularVelocity_2D_Integrand(Vars['rfrak_chi2_x'][0], Vars['rfrak_chi2_y'][0],(mul(mul(inverse(Cons['S_chi2']),OmegaVec_chiVec_2_3p5(Cons,Vars)),Cons['S_chi2']))[1:])
    else:
        dydt[3] = 0.0
        dydt[4] = 0.0
    dydt[8] = (Vars['v']*Vars['v']*Vars['v']/Cons['M'])[0]
    return dydt      

@njit
def Recalculate_4p0(Cons,y):
    Vars=Dict()
    Vars['v'] = np.array([y[0]])
    Vars['rfrak_chi1_x'] = np.array([y[1]])
    Vars['rfrak_chi1_y'] = np.array([y[2]])
    Vars['rfrak_chi2_x'] = np.array([y[3]])
    Vars['rfrak_chi2_y'] = np.array([y[4]])
    Vars['rfrak_frame_x'] = np.array([y[5]])
    Vars['rfrak_frame_y'] = np.array([y[6]])
    Vars['rfrak_frame_z'] = np.array([y[7]])
    Vars['Phi'] = np.array([y[8]])
    Vars['R'] = exp(Vars['rfrak_frame_x']*Cons['xHat'] + Vars['rfrak_frame_y']*Cons['yHat'] + Vars['rfrak_frame_z']*Cons['zHat'])
    Vars['nHat'] = mul(mul(Vars['R'],Cons['xHat']),conjugate(Vars['R']))
    Vars['lambdaHat'] = mul(mul(Vars['R'],Cons['yHat']),conjugate(Vars['R']))
    Vars['ellHat'] = mul(mul(Vars['R'],Cons['zHat']),conjugate(Vars['R']))
    Vars['R_S1'] = exp(Vars['rfrak_chi1_x']*Cons['xHat'] + Vars['rfrak_chi1_y']*Cons['yHat'])
    Vars['R_S2'] = exp(Vars['rfrak_chi2_x']*Cons['xHat'] + Vars['rfrak_chi2_y']*Cons['yHat'])
    Vars['chiVec1'] = mul(mul(mul(Cons['S_chi1'],Vars['R_S1']),Cons['zHat']),mul(conjugate(Vars['R_S1']),conjugate(Cons['S_chi1'])))
    Vars['chiVec2'] = mul(mul(mul(Cons['S_chi2'],Vars['R_S2']),Cons['zHat']),mul(conjugate(Vars['R_S2']),conjugate(Cons['S_chi2'])))
    Vars['chi1_n'] = np.array([dot(Vars['chiVec1'][1:],Vars['nHat'][1:])])
    Vars['chi1_lambda'] = np.array([dot(Vars['chiVec1'][1:],Vars['lambdaHat'][1:])])
    Vars['chi1_ell'] = np.array([dot(Vars['chiVec1'][1:],Vars['ellHat'][1:])])
    Vars['chi2_n'] = np.array([dot(Vars['chiVec2'][1:],Vars['nHat'][1:])])
    Vars['chi2_lambda'] = np.array([dot(Vars['chiVec2'][1:],Vars['lambdaHat'][1:])])
    Vars['chi2_ell'] = np.array([dot(Vars['chiVec2'][1:],Vars['ellHat'][1:])])
    Vars['S_ell'] = Cons['M1']**2*Vars['chi1_ell'] + Cons['M2']**2*Vars['chi2_ell']
    Vars['S_n'] = Cons['M1']**2*Vars['chi1_n'] + Cons['M2']**2*Vars['chi2_n']
    Vars['S_lambda'] = Cons['M1']**2*Vars['chi1_lambda'] + Cons['M2']**2*Vars['chi2_lambda']
    Vars['Sigma_ell'] = Cons['M']*(-Cons['M1']*Vars['chi1_ell'] + Cons['M2']*Vars['chi2_ell'])
    Vars['Sigma_n'] = Cons['M']*(-Cons['M1']*Vars['chi1_n'] + Cons['M2']*Vars['chi2_n'])
    Vars['Sigma_lambda'] = Cons['M']*(-Cons['M1']*Vars['chi1_lambda'] + Cons['M2']*Vars['chi2_lambda'])
    Vars['chi_s_ell'] = Vars['chi1_ell']/2 + Vars['chi2_ell']/2
    Vars['chi_a_ell'] = Vars['chi1_ell']/2 - Vars['chi2_ell']/2
    Vars['logv'] = log(Vars['v'])
    Vars['Fcal_coeff'] = 32*Cons['nu']**2*Vars['v']**10/5
    Vars['Fcal_SQ_4'] = Cons['chi1chi1']*(-89*Cons['delta']/192 + 89*Cons['nu']/96 - 89/192) - 103*Cons['chi1chi2']*Cons['nu']/48 + Cons['chi2chi2']*(89*Cons['delta']/192 + 89*Cons['nu']/96 - 89/192) + Vars['chi_a_ell']*(Vars['chi_a_ell']*(287/96 - 12*Cons['nu']) + 287*Vars['chi_s_ell']*Cons['delta']/48) + Vars['chi_s_ell']**2*(Cons['nu']/24 + 287/96)
    Vars['Fcal_SO_3'] = (-4*Vars['S_ell'] - 5*Vars['Sigma_ell']*Cons['delta']/4)/Cons['M']**2
    Vars['Fcal_SO_5'] = (Vars['S_ell']*(272*Cons['nu']/9 - 9/2) + Vars['Sigma_ell']*Cons['delta']*(43*Cons['nu']/4 - 13/16))/Cons['M']**2
    Vars['Fcal_SO_6'] = (-16*Vars['S_ell']*pi - 31*Vars['Sigma_ell']*Cons['delta']*pi/6)/Cons['M']**2
    Vars['Fcal_SO_7'] = (Vars['S_ell']*(-2810*Cons['nu']**2/27 + 6172*Cons['nu']/189 + 476645/6804) + Vars['Sigma_ell']*Cons['delta']*(-1501*Cons['nu']**2/36 + 1849*Cons['nu']/126 + 9535/336))/Cons['M']**2
    Vars['Fcal_SO_8'] = (Vars['S_ell']*pi*(13879*Cons['nu']/72 - 3485/96) + Vars['Sigma_ell']*Cons['delta']*pi*(130583*Cons['nu']/2016 - 7163/672))/Cons['M']**2
    Vars['E_SQ_4'] = -3*Vars['chi_a_ell']**2/2 - 3*Vars['chi_s_ell']**2/2 - Cons['delta']*(Cons['chi2chi2']/2 + 3*Vars['chi_a_ell']*Vars['chi_s_ell']) + Cons['nu']*(Cons['chi1chi2'] + 6*Vars['chi_a_ell']**2) + (Cons['chi1chi1'] + Cons['chi2chi2'])*(Cons['delta'] - 2*Cons['nu'] + 1)/4
    Vars['E_SO_3'] = (14*Vars['S_ell']/3 + 2*Vars['Sigma_ell']*Cons['delta'])/Cons['M']**2
    Vars['E_SO_5'] = (Vars['S_ell']*(11 - 61*Cons['nu']/9) + Vars['Sigma_ell']*Cons['delta']*(3 - 10*Cons['nu']/3))/Cons['M']**2
    Vars['E_SO_7'] = (Vars['S_ell']*(29*Cons['nu']**2/12 - 367*Cons['nu']/4 + 135/4) + Vars['Sigma_ell']*Cons['delta']*(5*Cons['nu']**2/4 - 39*Cons['nu'] + 27/4))/Cons['M']**2
    return Vars

@njit
def OmegaVec_chiVec_1_4p0(Cons,Vars):
    Omega1_coeff = Vars['v']**5/Cons['M']
    return Omega1_coeff*(Vars['ellHat']*(-0.75*Cons['delta'] + 0.5*Cons['nu'] + Vars['v']**2*(Cons['delta']*(0.625*Cons['nu'] - 0.5625) + Cons['nu']*(1.25 - 0.0416666666666667*Cons['nu']) + Vars['v']**2*(Cons['delta']*(Cons['nu']*(4.875 - 0.15625*Cons['nu']) - 0.84375) + Cons['nu']*(Cons['nu']*(-0.0208333333333333*Cons['nu'] - 3.28125) + 0.1875) + 0.84375) + 0.5625) + 0.75) + Vars['nHat']*Vars['v']*(3.0*Vars['chi1_n']*Cons['nu'] + 3.0*Cons['M2']**2*Vars['chi2_n']/Cons['M']**2) - Cons['M2']**2*Vars['chiVec2']*Vars['v']/Cons['M']**2)

@njit
def OmegaVec_chiVec_2_4p0(Cons,Vars):
    Omega2_coeff = Vars['v']**5/Cons['M']
    return Omega2_coeff*(Vars['ellHat']*(0.75*Cons['delta'] + 0.5*Cons['nu'] + Vars['v']**2*(Cons['delta']*(0.5625 - 0.625*Cons['nu']) + Cons['nu']*(1.25 - 0.0416666666666667*Cons['nu']) + Vars['v']**2*(Cons['delta']*(Cons['nu']*(0.15625*Cons['nu'] - 4.875) + 0.84375) + Cons['nu']*(Cons['nu']*(-0.0208333333333333*Cons['nu'] - 3.28125) + 0.1875) + 0.84375) + 0.5625) + 0.75) + Vars['nHat']*Vars['v']*(3.0*Vars['chi2_n']*Cons['nu'] + 3.0*Cons['M1']**2*Vars['chi1_n']/Cons['M']**2) - Cons['M1']**2*Vars['chiVec1']*Vars['v']/Cons['M']**2)

@njit
def OmegaVec_4p0(Cons,Vars):
    gamma_PN_5 = (Vars['S_ell']*(0.888888888888889*Cons['nu'] + 3.33333333333333) + 2.0*Vars['Sigma_ell']*Cons['delta'])/Cons['M']**2
    a_ell_0 = 7.0*Vars['S_n'] + 3.0*Vars['Sigma_n']*Cons['delta']
    gamma_PN_3 = (1.66666666666667*Vars['S_ell'] + Vars['Sigma_ell']*Cons['delta'])/Cons['M']**2
    gamma_PN_4 = 1.0 - 5.41666666666667*Cons['nu']
    gamma_PN_2 = 1.0 - 0.333333333333333*Cons['nu']
    gamma_PN_0 = 1.00000000000000
    a_ell_2 = Vars['S_n']*(-9.66666666666667*Cons['nu'] - 10.0) + Vars['Sigma_n']*Cons['delta']*(-4.5*Cons['nu'] - 6.0)
    gamma_PN_7 = (Vars['S_ell']*(-6.0*Cons['nu']**2 - 10.5833333333333*Cons['nu'] + 5.0) - 2.66666666666667*Vars['Sigma_ell']*Cons['delta']*Cons['nu']**2 + Vars['Sigma_ell']*Cons['delta']*(3.0 - 10.1666666666667*Cons['nu']))/Cons['M']**2
    gamma_PN_6 = 0.0123456790123457*Cons['nu']**3 + 6.36111111111111*Cons['nu']**2 - 2.98177812235564*Cons['nu'] + 1.0
    a_ell_4 = Vars['S_n']*(5.77777777777778*Cons['nu']**2 + 14.75*Cons['nu'] + 1.5) + Vars['Sigma_n']*Cons['delta']*(2.83333333333333*Cons['nu']**2 + 9.125*Cons['nu'] + 1.5)
    return Vars['ellHat']*Vars['v']**3/Cons['M'] + Vars['nHat']*Vars['v']**6*(a_ell_0 + Vars['v']**2*(a_ell_2 + a_ell_4*Vars['v']**2))*(gamma_PN_0 + Vars['v']**2*(gamma_PN_2 + Vars['v']*(gamma_PN_3 + Vars['v']*(gamma_PN_4 + Vars['v']*(gamma_PN_5 + Vars['v']*(gamma_PN_6 + gamma_PN_7*Vars['v']))))))/Cons['M']**3


@njit
def TaylorT1_4p0(Cons,Vars):
    Flux = Vars['Fcal_coeff']*(Cons['Fcal_0'] + Vars['v']**2*(Cons['Fcal_2'] + Vars['v']*(Cons['Fcal_3'] + Vars['Fcal_SO_3'] + Vars['v']*(Cons['Fcal_4'] + Vars['Fcal_SQ_4'] + Vars['v']*(Cons['Fcal_5'] + Vars['Fcal_SO_5'] + Vars['v']*(Cons['Fcal_6'] + Vars['Fcal_SO_6'] + Cons['Fcal_lnv_6']*Vars['logv'] + Vars['v']*(Cons['Fcal_7'] + Vars['Fcal_SO_7'] + Vars['v']*(Cons['Fcal_8'] + Vars['Fcal_SO_8'] + Cons['Fcal_lnv_8']*Vars['logv']))))))))
    dEdV = -0.5*Cons['M']*Cons['nu']*Vars['v']*(2.0*Cons['E_0'] + Vars['v']**2*(4.0*Cons['E_2'] + Vars['v']*(5.0*Vars['E_SO_3'] + Vars['v']*(6.0*Cons['E_4'] + 6.0*Vars['E_SQ_4'] + Vars['v']*(7.0*Vars['E_SO_5'] + Vars['v']*(8.0*Cons['E_6'] + Vars['v']*(9.0*Vars['E_SO_7'] + Vars['v']*(10.0*Cons['E_8'] + Cons['E_lnv_8']*(10.0*Vars['logv'] + 1.0)))))))))
    Absorption = 0
    dvdt = (-Absorption - Flux)/dEdV
    dydt=np.zeros(9)
    rfrak_frame=[Vars['rfrak_frame_x'][0],Vars['rfrak_frame_y'][0],Vars['rfrak_frame_z'][0]]
    [dydt[5],dydt[6],dydt[7]] = FrameFromAngularVelocityIntegrand(rfrak_frame, OmegaVec_4p0(Cons,Vars)[1:])
    dydt[0] = dvdt[0]
    if(Cons['EvolveSpin1']):
        dydt[1], dydt[2]=FrameFromAngularVelocity_2D_Integrand(Vars['rfrak_chi1_x'][0], Vars['rfrak_chi1_y'][0],(mul(mul(inverse(Cons['S_chi1']),OmegaVec_chiVec_1_4p0(Cons,Vars)),Cons['S_chi1']))[1:])
    else:
        dydt[1] = 0.0
        dydt[2] = 0.0
    if(Cons['EvolveSpin2']):
        dydt[3], dydt[4]=FrameFromAngularVelocity_2D_Integrand(Vars['rfrak_chi2_x'][0], Vars['rfrak_chi2_y'][0],(mul(mul(inverse(Cons['S_chi2']),OmegaVec_chiVec_2_4p0(Cons,Vars)),Cons['S_chi2']))[1:])
    else:
        dydt[3] = 0.0
        dydt[4] = 0.0
    dydt[8] = (Vars['v']*Vars['v']*Vars['v']/Cons['M'])[0]
    return dydt      

class PNEv:
    def Integrand(t,y):
        Cons=Dict()
        for i, j in PNEv.Cons.items():
            Cons[i]=j
        Y=List()
        [Y.append(i) for i in y]
        Vars=PNEv.SwitcherOrder.get(2*PNEv.PNEvolutionOrder)(PNEv.Cons,Y)
        dydt=PNEv.SwitcherTn.get(PNEv.TaylorTn+20*PNEv.PNEvolutionOrder)(PNEv.Cons,Vars)
        if Vars['v']>=1.0 and PNEv.NotForward:
            print("Beyond domain of PN validity, this is a good way to terminate.")
            PNEv.terminal1=False
        if dydt[0]<1.0e-12 and PNEv.NotForward:
            print("v is decreasing, which is not an uncommon way to stop.")
            PNEv.terminal2=False
        return dydt
        
    def Evolution(xHat_i, yHat_i, zHat_i, M1_i, M2_i, v_i, S_chi1_i, S_chi2_i, rfrak_chi1_x_i, rfrak_chi1_y_i, rfrak_chi2_x_i, rfrak_chi2_y_i, rfrak_frame_x_i, rfrak_frame_y_i, rfrak_frame_z_i, t_PNStart=False, t_PNEnd=False, PNEvolutionOrder=1, TaylorTn=1, StepsPerOrbit=32, ForwardInTime=False, tol=1e-12, MinStep=1e-7): 
        # Initialization of constants
        PNEv.terminal1=True
        PNEv.terminal2=True
        PNEv.NotForward=True
        PNEv.PNEvolutionOrder=PNEvolutionOrder
        PNEv.TaylorTn=TaylorTn
        PNEv.SwitcherOrder= {0:Recalculate_0,
            1:Recalculate_0p50,
            2:Recalculate_1p0,
            3:Recalculate_1p5,
            4:Recalculate_2p0,
            5:Recalculate_2p5,
            6:Recalculate_3p0,
            7:Recalculate_3p5,
            8:Recalculate_4p0}
        PNEv.SwitcherTn= {1:TaylorT1_0,
            11:TaylorT1_0p50,
            21:TaylorT1_1p0,
            31:TaylorT1_1p5,
            41:TaylorT1_2p0,
            51:TaylorT1_2p5,
            61:TaylorT1_3p0,
            71:TaylorT1_3p5,
            81:TaylorT1_4p0}
        PNEv.Cons=Initialization(xHat_i, yHat_i, zHat_i, M1_i, M2_i, v_i, S_chi1_i, S_chi2_i, rfrak_chi1_x_i, rfrak_chi1_y_i, rfrak_chi2_x_i, rfrak_chi2_y_i, rfrak_frame_x_i, rfrak_frame_y_i, rfrak_frame_z_i)
        PNEv.Cons['EvolveSpin1']=np.array([(np.linalg.norm(mul(PNEv.Cons['S_chi1'],conjugate(PNEv.Cons['S_chi1'])))>1e-12)*1.0])
        PNEv.Cons['EvolveSpin2']=np.array([(np.linalg.norm(mul(PNEv.Cons['S_chi2'],conjugate(PNEv.Cons['S_chi2'])))>1e-12)*1.0])
    
        def terminate(t,y):
            return 1.0*PNEv.terminal1*PNEv.terminal2
        terminate.terminal=True
        TMerger=5.0/(256.0*PNEv.Cons['nu']*v_i**8)
        TEnd=TMerger
        if t_PNEnd:
            TEnd=t_PNEnd
        time=[0.0]
        while time[-1]<TEnd and 2*PNEv.Cons['M']*(256*PNEv.Cons['nu']*(TMerger-time[-1])/5)**(3/8)/StepsPerOrbit>MinStep:
            time.append(time[-1]+(2*PNEv.Cons['M']*(256*PNEv.Cons['nu']*(TMerger-time[-1])/5)**(3/8)/StepsPerOrbit)[0])
        time=np.delete(time, -1)
       
        # Integrate
        yy=solve_ivp(PNEv.Integrand, [time[0],time[-1]], [v_i,rfrak_chi1_x_i,\
            rfrak_chi1_y_i,rfrak_chi2_x_i,rfrak_chi2_y_i,rfrak_frame_x_i,\
            rfrak_frame_y_i,rfrak_frame_z_i,0.0], method='DOP853',\
            t_eval=time, dense_output=True, events=terminate, rtol=tol, atol=tol)           
        if ForwardInTime:
            PNEv.NotForward=False
            time=[0.0]
            TStart=-3*TMerger
            if t_PNStart:
                TStart=t_PNStart
            while time[-1]>TStart:
                time.append(time[-1]-(2*PNEv.Cons['M']*(256*PNEv.Cons['nu']*(TMerger-time[-1])/5)**(3/8)/StepsPerOrbit)[0])
            yyForward=solve_ivp(PNEv.Integrand, [time[0],time[-1]], [v_i,rfrak_chi1_x_i,\
                rfrak_chi1_y_i,rfrak_chi2_x_i,rfrak_chi2_y_i,rfrak_frame_x_i,\
                rfrak_frame_y_i,rfrak_frame_z_i,0.0], method='DOP853',\
                t_eval=time, dense_output=True, rtol=tol, atol=tol)
            yy.t=np.append(yyForward.t[1:][::-1],yy.t)
            data=np.empty((9,len(yy.t)))
            for i in range(9):
                data[i]=np.append(yyForward.y[i][1:][::-1],yy.y[i])
            yy.y=data
             
        return yy
