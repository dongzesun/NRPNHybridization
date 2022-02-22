# File produced automatically by PNCodeGen.ipynb
from scipy.integrate import solve_ivp
import numpy as np
from numpy import conjugate, dot, cross, exp, log, sqrt, pi
from numpy import euler_gamma as EulerGamma
from numba import jit, njit
from numba.typed import List, Dict
from scipy.interpolate import InterpolatedUnivariateSpline as Spline
from scipy.special import zeta
import quaternionic

qmul = njit(quaternionic.algebra.multiply)
@njit
def mul(A,B):
    C=np.empty(4)
    qmul(A,B,C)
    return C

@njit
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

@njit
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

@njit
def Initialization(xHat_i, yHat_i, zHat_i, M1_i, M2_i, v_i, S_chi1_i, S_chi2_i, rfrak_chi1_x_i, rfrak_chi1_y_i, rfrak_chi2_x_i, rfrak_chi2_y_i, rfrak_frame_x_i, rfrak_frame_y_i, rfrak_frame_z_i): 
    Cons={'M1':M1_i, 'xHat':xHat_i}
    Cons['xHat']=xHat_i
    Cons['yHat']=yHat_i
    Cons['zHat']=zHat_i
    Cons['M1']=M1_i
    Cons['M2']=M2_i
    v=v_i
    Cons['S_chi1']=S_chi1_i
    Cons['S_chi2']=S_chi2_i
    rfrak_chi1_x=rfrak_chi1_x_i
    rfrak_chi1_y=rfrak_chi1_y_i
    rfrak_chi2_x=rfrak_chi2_x_i
    rfrak_chi2_y=rfrak_chi2_y_i
    rfrak_frame_x=rfrak_frame_x_i
    rfrak_frame_y=rfrak_frame_y_i
    rfrak_frame_z=rfrak_frame_z_i
    Cons['M']=Cons['M1'] + Cons['M2']
    Cons['delta']=(Cons['M1'] - Cons['M2'])/Cons['M']
    Cons['nu']=Cons['M1']*Cons['M2']/Cons['M']**2
    R=exp(rfrak_frame_x*Cons['xHat'] + rfrak_frame_y*Cons['yHat'] + rfrak_frame_z*Cons['zHat'])
    nHat=R*Cons['xHat']*conjugate(R)
    lambdaHat=R*Cons['yHat']*conjugate(R)
    ellHat=R*Cons['zHat']*conjugate(R)
    R_S1=exp(rfrak_chi1_x*Cons['xHat'] + rfrak_chi1_y*Cons['yHat'])
    R_S2=exp(rfrak_chi2_x*Cons['xHat'] + rfrak_chi2_y*Cons['yHat'])
    chiVec1=Cons['S_chi1']*R_S1*Cons['zHat']*conjugate(R_S1)*conjugate(Cons['S_chi1'])
    chiVec2=Cons['S_chi2']*R_S2*Cons['zHat']*conjugate(R_S2)*conjugate(Cons['S_chi2'])
    Cons['chi1chi1']=dot(chiVec1.vector,chiVec1.vector)
    Cons['chi1chi2']=dot(chiVec1.vector,chiVec2.vector)
    Cons['chi2chi2']=dot(chiVec2.vector,chiVec2.vector)
    chi1_n=dot(chiVec1.vector,nHat.vector)
    chi1_lambda=dot(chiVec1.vector,lambdaHat.vector)
    chi1_ell=dot(chiVec1.vector,ellHat.vector)
    chi2_n=dot(chiVec2.vector,nHat.vector)
    chi2_lambda=dot(chiVec2.vector,lambdaHat.vector)
    chi2_ell=dot(chiVec2.vector,ellHat.vector)
    S_ell=Cons['M1']**2*chi1_ell + Cons['M2']**2*chi2_ell
    S_n=Cons['M1']**2*chi1_n + Cons['M2']**2*chi2_n
    S_lambda=Cons['M1']**2*chi1_lambda + Cons['M2']**2*chi2_lambda
    Sigma_ell=Cons['M']*(-Cons['M1']*chi1_ell + Cons['M2']*chi2_ell)
    Sigma_n=Cons['M']*(-Cons['M1']*chi1_n + Cons['M2']*chi2_n)
    Sigma_lambda=Cons['M']*(-Cons['M1']*chi1_lambda + Cons['M2']*chi2_lambda)
    chi_s_ell=chi1_ell/2 + chi2_ell/2
    chi_a_ell=chi1_ell/2 - chi2_ell/2
    logv=log(v)
    Fcal_coeff=32*Cons['nu']**2*v**10/5
    Cons['Fcal_0']=1
    Cons['Fcal_2']=-35*Cons['nu']/12 - 1247/336
    Cons['Fcal_3']=4*pi
    Cons['Fcal_4']=65*Cons['nu']**2/18 + 9271*Cons['nu']/504 - 44711/9072
    Cons['Fcal_5']=pi*(-583*Cons['nu']/24 - 8191/672)
    Cons['Fcal_6']=-775*Cons['nu']**3/324 - 94403*Cons['nu']**2/3024 + Cons['nu']*(-134543/7776 + 41*pi**2/48) - 1712*log(4)/105 - 1712*EulerGamma/105 + 16*pi**2/3 + 6643739519/69854400
    Cons['Fcal_lnv_6']=-1712/105
    Cons['Fcal_7']=pi*(193385*Cons['nu']**2/3024 + 214745*Cons['nu']/1728 - 16285/504)
    Cons['Fcal_8']=-1369*pi**2/126 - 323105549467/3178375200 - 47385*log(3)/1568 + 232597*EulerGamma/4410 + 39931*log(2)/294
    Cons['Fcal_lnv_8']=232597/4410
    Fcal_SQ_4=Cons['chi1chi1']*(-89*Cons['delta']/192 + 89*Cons['nu']/96 - 89/192) - 103*Cons['chi1chi2']*Cons['nu']/48 + Cons['chi2chi2']*(89*Cons['delta']/192 + 89*Cons['nu']/96 - 89/192) + chi_a_ell*(chi_a_ell*(287/96 - 12*Cons['nu']) + 287*chi_s_ell*Cons['delta']/48) + chi_s_ell**2*(Cons['nu']/24 + 287/96)
    Fcal_SO_3=(-4*S_ell - 5*Sigma_ell*Cons['delta']/4)/Cons['M']**2
    Fcal_SO_5=(S_ell*(272*Cons['nu']/9 - 9/2) + Sigma_ell*Cons['delta']*(43*Cons['nu']/4 - 13/16))/Cons['M']**2
    Fcal_SO_6=(-16*S_ell*pi - 31*Sigma_ell*Cons['delta']*pi/6)/Cons['M']**2
    Fcal_SO_7=(S_ell*(-2810*Cons['nu']**2/27 + 6172*Cons['nu']/189 + 476645/6804) + Sigma_ell*Cons['delta']*(-1501*Cons['nu']**2/36 + 1849*Cons['nu']/126 + 9535/336))/Cons['M']**2
    Fcal_SO_8=(S_ell*pi*(13879*Cons['nu']/72 - 3485/96) + Sigma_ell*Cons['delta']*pi*(130583*Cons['nu']/2016 - 7163/672))/Cons['M']**2
    Cons['E_0']=1
    Cons['E_2']=-Cons['nu']/12 - 3/4
    Cons['E_4']=-Cons['nu']**2/24 + 19*Cons['nu']/8 - 27/8
    Cons['E_6']=-35*Cons['nu']**3/5184 - 155*Cons['nu']**2/96 + Cons['nu']*(34445/576 - 205*pi**2/96) - 675/64
    Cons['E_8']=77*Cons['nu']**4/31104 + 301*Cons['nu']**3/1728 + Cons['nu']**2*(-498449/3456 + 3157*pi**2/576) + Cons['nu']*(-123671/5760 + 896*EulerGamma/15 + 9037*pi**2/1536 + 1792*log(2)/15) - 3969/128
    Cons['E_lnv_8']=896*Cons['nu']/15
    E_SQ_4=-3*chi_a_ell**2/2 - 3*chi_s_ell**2/2 - Cons['delta']*(Cons['chi2chi2']/2 + 3*chi_a_ell*chi_s_ell) + Cons['nu']*(Cons['chi1chi2'] + 6*chi_a_ell**2) + (Cons['chi1chi1'] + Cons['chi2chi2'])*(Cons['delta'] - 2*Cons['nu'] + 1)/4
    E_SO_3=(14*S_ell/3 + 2*Sigma_ell*Cons['delta'])/Cons['M']**2
    E_SO_5=(S_ell*(11 - 61*Cons['nu']/9) + Sigma_ell*Cons['delta']*(3 - 10*Cons['nu']/3))/Cons['M']**2
    E_SO_7=(S_ell*(29*Cons['nu']**2/12 - 367*Cons['nu']/4 + 135/4) + Sigma_ell*Cons['delta']*(5*Cons['nu']**2/4 - 39*Cons['nu'] + 27/4))/Cons['M']**2
    return Cons

@njit
def Recalculate_0(Cons,y):
    Vars={'v':y[0],'ellHat':[0,0,0,0]}
    Vars['v'] = y[0]
    Vars['rfrak_chi1_x'] = y[1]
    Vars['rfrak_chi1_y'] = y[2]
    Vars['rfrak_chi2_x'] = y[3]
    Vars['rfrak_chi2_y'] = y[4]
    Vars['rfrak_frame_x'] = y[5]
    Vars['rfrak_frame_y'] = y[6]
    Vars['rfrak_frame_z'] = y[7]
    Vars['Phi'] = y[8]
    Vars['R'] = exp(Vars['rfrak_frame_x']*Cons['xHat'] + Vars['rfrak_frame_y']*Cons['yHat'] + Vars['rfrak_frame_z']*Cons['zHat'])
    Vars['nHat'] = Vars['R']*Cons['xHat']*conjugate(Vars['R'])
    Vars['ellHat'] = Vars['R']*Cons['zHat']*conjugate(Vars['R'])
    Vars['R_S1'] = exp(Vars['rfrak_chi1_x']*Cons['xHat'] + Vars['rfrak_chi1_y']*Cons['yHat'])
    Vars['R_S2'] = exp(Vars['rfrak_chi2_x']*Cons['xHat'] + Vars['rfrak_chi2_y']*Cons['yHat'])
    Vars['chiVec1'] = Cons['S_chi1']*Vars['R_S1']*Cons['zHat']*conjugate(Vars['R_S1'])*conjugate(Cons['S_chi1'])
    Vars['chiVec2'] = Cons['S_chi2']*Vars['R_S2']*Cons['zHat']*conjugate(Vars['R_S2'])*conjugate(Cons['S_chi2'])
    Vars['chi1_n'] = dot(Vars['chiVec1'].vector,Vars['nHat'].vector)
    Vars['chi1_ell'] = dot(Vars['chiVec1'].vector,Vars['ellHat'].vector)
    Vars['chi2_n'] = dot(Vars['chiVec2'].vector,Vars['nHat'].vector)
    Vars['chi2_ell'] = dot(Vars['chiVec2'].vector,Vars['ellHat'].vector)
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
    gamma_PN_0 = 1.00000000000000
    a_ell_0 = 7.0*Vars['S_n'] + 3.0*Vars['Sigma_n']*Cons['delta']
    return Vars['ellHat']*Vars['v']**3/Cons['M'] + a_ell_0*gamma_PN_0*Vars['nHat']*Vars['v']**6/Cons['M']**3


@njit
def TaylorT1_0(Cons,Vars):
    Flux = Cons['Fcal_0']*Vars['Fcal_coeff']
    dEdV = -Cons['E_0']*Cons['M']*Cons['nu']*Vars['v']
    Absorption = 0
    dvdt = (-Absorption - Flux)/dEdV
    dydt=np.zeros(9)
    rfrak_frame=[Vars['rfrak_frame_x'],Vars['rfrak_frame_y'],Vars['rfrak_frame_z']]
    [dydt[5],dydt[6],dydt[7]] = FrameFromAngularVelocityIntegrand(rfrak_frame, OmegaVec_0(Cons,Vars).vector)
    dydt[0] = dvdt
    if(Cons['EvolveSpin1']):
        dydt[1], dydt[2]= FrameFromAngularVelocity_2D_Integrand(Vars['rfrak_chi1_x'], Vars['rfrak_chi1_y'],(Cons['S_chi1'].inverse*OmegaVec_chiVec_1_0(Cons,Vars)*Cons['S_chi1']).vector)
    else:
        dydt[1] = 0.0
        dydt[2] = 0.0
    if(Cons['EvolveSpin2']):
        dydt[3], dydt[4] = FrameFromAngularVelocity_2D_Integrand(Vars['rfrak_chi2_x'], Vars['rfrak_chi2_y'],(Cons['S_chi2'].inverse*OmegaVec_chiVec_2_0(Cons,Vars)*Cons['S_chi2']).vector)
    else:
        dydt[3] = 0.0
        dydt[4] = 0.0
    dydt[8] = Vars['v']*Vars['v']*Vars['v']/Cons['M']
    return dydt      

@njit
def Recalculate_0p50(Cons,y):
    Vars={'v':y[0],'ellHat':[0,0,0,0]}
    Vars['v'] = y[0]
    Vars['rfrak_chi1_x'] = y[1]
    Vars['rfrak_chi1_y'] = y[2]
    Vars['rfrak_chi2_x'] = y[3]
    Vars['rfrak_chi2_y'] = y[4]
    Vars['rfrak_frame_x'] = y[5]
    Vars['rfrak_frame_y'] = y[6]
    Vars['rfrak_frame_z'] = y[7]
    Vars['Phi'] = y[8]
    Vars['R'] = exp(Vars['rfrak_frame_x']*Cons['xHat'] + Vars['rfrak_frame_y']*Cons['yHat'] + Vars['rfrak_frame_z']*Cons['zHat'])
    Vars['nHat'] = Vars['R']*Cons['xHat']*conjugate(Vars['R'])
    Vars['ellHat'] = Vars['R']*Cons['zHat']*conjugate(Vars['R'])
    Vars['R_S1'] = exp(Vars['rfrak_chi1_x']*Cons['xHat'] + Vars['rfrak_chi1_y']*Cons['yHat'])
    Vars['R_S2'] = exp(Vars['rfrak_chi2_x']*Cons['xHat'] + Vars['rfrak_chi2_y']*Cons['yHat'])
    Vars['chiVec1'] = Cons['S_chi1']*Vars['R_S1']*Cons['zHat']*conjugate(Vars['R_S1'])*conjugate(Cons['S_chi1'])
    Vars['chiVec2'] = Cons['S_chi2']*Vars['R_S2']*Cons['zHat']*conjugate(Vars['R_S2'])*conjugate(Cons['S_chi2'])
    Vars['chi1_n'] = dot(Vars['chiVec1'].vector,Vars['nHat'].vector)
    Vars['chi1_ell'] = dot(Vars['chiVec1'].vector,Vars['ellHat'].vector)
    Vars['chi2_n'] = dot(Vars['chiVec2'].vector,Vars['nHat'].vector)
    Vars['chi2_ell'] = dot(Vars['chiVec2'].vector,Vars['ellHat'].vector)
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
    gamma_PN_0 = 1.00000000000000
    a_ell_0 = 7.0*Vars['S_n'] + 3.0*Vars['Sigma_n']*Cons['delta']
    return Vars['ellHat']*Vars['v']**3/Cons['M'] + a_ell_0*gamma_PN_0*Vars['nHat']*Vars['v']**6/Cons['M']**3


@njit
def TaylorT1_0p50(Cons,Vars):
    Flux = Cons['Fcal_0']*Vars['Fcal_coeff']
    dEdV = -Cons['E_0']*Cons['M']*Cons['nu']*Vars['v']
    Absorption = 0
    dvdt = (-Absorption - Flux)/dEdV
    dydt=np.zeros(9)
    rfrak_frame=[Vars['rfrak_frame_x'],Vars['rfrak_frame_y'],Vars['rfrak_frame_z']]
    [dydt[5],dydt[6],dydt[7]] = FrameFromAngularVelocityIntegrand(rfrak_frame, OmegaVec_0p50(Cons,Vars).vector)
    dydt[0] = dvdt
    if(Cons['EvolveSpin1']):
        dydt[1], dydt[2]= FrameFromAngularVelocity_2D_Integrand(Vars['rfrak_chi1_x'], Vars['rfrak_chi1_y'],(Cons['S_chi1'].inverse*OmegaVec_chiVec_1_0p50(Cons,Vars)*Cons['S_chi1']).vector)
    else:
        dydt[1] = 0.0
        dydt[2] = 0.0
    if(Cons['EvolveSpin2']):
        dydt[3], dydt[4] = FrameFromAngularVelocity_2D_Integrand(Vars['rfrak_chi2_x'], Vars['rfrak_chi2_y'],(Cons['S_chi2'].inverse*OmegaVec_chiVec_2_0p50(Cons,Vars)*Cons['S_chi2']).vector)
    else:
        dydt[3] = 0.0
        dydt[4] = 0.0
    dydt[8] = Vars['v']*Vars['v']*Vars['v']/Cons['M']
    return dydt      

@njit
def Recalculate_1p0(Cons,y):
    Vars={'v':y[0],'ellHat':[0,0,0,0]}
    Vars['v'] = y[0]
    Vars['rfrak_chi1_x'] = y[1]
    Vars['rfrak_chi1_y'] = y[2]
    Vars['rfrak_chi2_x'] = y[3]
    Vars['rfrak_chi2_y'] = y[4]
    Vars['rfrak_frame_x'] = y[5]
    Vars['rfrak_frame_y'] = y[6]
    Vars['rfrak_frame_z'] = y[7]
    Vars['Phi'] = y[8]
    Vars['R'] = exp(Vars['rfrak_frame_x']*Cons['xHat'] + Vars['rfrak_frame_y']*Cons['yHat'] + Vars['rfrak_frame_z']*Cons['zHat'])
    Vars['nHat'] = Vars['R']*Cons['xHat']*conjugate(Vars['R'])
    Vars['ellHat'] = Vars['R']*Cons['zHat']*conjugate(Vars['R'])
    Vars['R_S1'] = exp(Vars['rfrak_chi1_x']*Cons['xHat'] + Vars['rfrak_chi1_y']*Cons['yHat'])
    Vars['R_S2'] = exp(Vars['rfrak_chi2_x']*Cons['xHat'] + Vars['rfrak_chi2_y']*Cons['yHat'])
    Vars['chiVec1'] = Cons['S_chi1']*Vars['R_S1']*Cons['zHat']*conjugate(Vars['R_S1'])*conjugate(Cons['S_chi1'])
    Vars['chiVec2'] = Cons['S_chi2']*Vars['R_S2']*Cons['zHat']*conjugate(Vars['R_S2'])*conjugate(Cons['S_chi2'])
    Vars['chi1_n'] = dot(Vars['chiVec1'].vector,Vars['nHat'].vector)
    Vars['chi1_ell'] = dot(Vars['chiVec1'].vector,Vars['ellHat'].vector)
    Vars['chi2_n'] = dot(Vars['chiVec2'].vector,Vars['nHat'].vector)
    Vars['chi2_ell'] = dot(Vars['chiVec2'].vector,Vars['ellHat'].vector)
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
    gamma_PN_0 = 1.00000000000000
    a_ell_0 = 7.0*Vars['S_n'] + 3.0*Vars['Sigma_n']*Cons['delta']
    a_ell_2 = Vars['S_n']*(-9.66666666666667*Cons['nu'] - 10.0) + Vars['Sigma_n']*Cons['delta']*(-4.5*Cons['nu'] - 6.0)
    gamma_PN_2 = 1.0 - 0.333333333333333*Cons['nu']
    return Vars['ellHat']*Vars['v']**3/Cons['M'] + Vars['nHat']*Vars['v']**6*(a_ell_0 + a_ell_2*Vars['v']**2)*(gamma_PN_0 + gamma_PN_2*Vars['v']**2)/Cons['M']**3


@njit
def TaylorT1_1p0(Cons,Vars):
    Flux = Vars['Fcal_coeff']*(Cons['Fcal_0'] + Cons['Fcal_2']*Vars['v']**2)
    dEdV = -Cons['M']*Cons['nu']*Vars['v']*(Cons['E_0'] + 2.0*Cons['E_2']*Vars['v']**2)
    Absorption = 0
    dvdt = (-Absorption - Flux)/dEdV
    dydt=np.zeros(9)
    rfrak_frame=[Vars['rfrak_frame_x'],Vars['rfrak_frame_y'],Vars['rfrak_frame_z']]
    [dydt[5],dydt[6],dydt[7]] = FrameFromAngularVelocityIntegrand(rfrak_frame, OmegaVec_1p0(Cons,Vars).vector)
    dydt[0] = dvdt
    if(Cons['EvolveSpin1']):
        dydt[1], dydt[2]= FrameFromAngularVelocity_2D_Integrand(Vars['rfrak_chi1_x'], Vars['rfrak_chi1_y'],(Cons['S_chi1'].inverse*OmegaVec_chiVec_1_1p0(Cons,Vars)*Cons['S_chi1']).vector)
    else:
        dydt[1] = 0.0
        dydt[2] = 0.0
    if(Cons['EvolveSpin2']):
        dydt[3], dydt[4] = FrameFromAngularVelocity_2D_Integrand(Vars['rfrak_chi2_x'], Vars['rfrak_chi2_y'],(Cons['S_chi2'].inverse*OmegaVec_chiVec_2_1p0(Cons,Vars)*Cons['S_chi2']).vector)
    else:
        dydt[3] = 0.0
        dydt[4] = 0.0
    dydt[8] = Vars['v']*Vars['v']*Vars['v']/Cons['M']
    return dydt      

@njit
def Recalculate_1p5(Cons,y):
    Vars={'v':y[0],'ellHat':[0,0,0,0]}
    Vars['v'] = y[0]
    Vars['rfrak_chi1_x'] = y[1]
    Vars['rfrak_chi1_y'] = y[2]
    Vars['rfrak_chi2_x'] = y[3]
    Vars['rfrak_chi2_y'] = y[4]
    Vars['rfrak_frame_x'] = y[5]
    Vars['rfrak_frame_y'] = y[6]
    Vars['rfrak_frame_z'] = y[7]
    Vars['Phi'] = y[8]
    Vars['R'] = exp(Vars['rfrak_frame_x']*Cons['xHat'] + Vars['rfrak_frame_y']*Cons['yHat'] + Vars['rfrak_frame_z']*Cons['zHat'])
    Vars['nHat'] = Vars['R']*Cons['xHat']*conjugate(Vars['R'])
    Vars['lambdaHat'] = Vars['R']*Cons['yHat']*conjugate(Vars['R'])
    Vars['ellHat'] = Vars['R']*Cons['zHat']*conjugate(Vars['R'])
    Vars['R_S1'] = exp(Vars['rfrak_chi1_x']*Cons['xHat'] + Vars['rfrak_chi1_y']*Cons['yHat'])
    Vars['R_S2'] = exp(Vars['rfrak_chi2_x']*Cons['xHat'] + Vars['rfrak_chi2_y']*Cons['yHat'])
    Vars['chiVec1'] = Cons['S_chi1']*Vars['R_S1']*Cons['zHat']*conjugate(Vars['R_S1'])*conjugate(Cons['S_chi1'])
    Vars['chiVec2'] = Cons['S_chi2']*Vars['R_S2']*Cons['zHat']*conjugate(Vars['R_S2'])*conjugate(Cons['S_chi2'])
    Vars['chi1_n'] = dot(Vars['chiVec1'].vector,Vars['nHat'].vector)
    Vars['chi1_lambda'] = dot(Vars['chiVec1'].vector,Vars['lambdaHat'].vector)
    Vars['chi1_ell'] = dot(Vars['chiVec1'].vector,Vars['ellHat'].vector)
    Vars['chi2_n'] = dot(Vars['chiVec2'].vector,Vars['nHat'].vector)
    Vars['chi2_lambda'] = dot(Vars['chiVec2'].vector,Vars['lambdaHat'].vector)
    Vars['chi2_ell'] = dot(Vars['chiVec2'].vector,Vars['ellHat'].vector)
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
    gamma_PN_0 = 1.00000000000000
    a_ell_2 = Vars['S_n']*(-9.66666666666667*Cons['nu'] - 10.0) + Vars['Sigma_n']*Cons['delta']*(-4.5*Cons['nu'] - 6.0)
    gamma_PN_2 = 1.0 - 0.333333333333333*Cons['nu']
    a_ell_0 = 7.0*Vars['S_n'] + 3.0*Vars['Sigma_n']*Cons['delta']
    gamma_PN_3 = (1.66666666666667*Vars['S_ell'] + Vars['Sigma_ell']*Cons['delta'])/Cons['M']**2
    return Vars['ellHat']*Vars['v']**3/Cons['M'] + Vars['nHat']*Vars['v']**6*(a_ell_0 + a_ell_2*Vars['v']**2)*(gamma_PN_0 + Vars['v']**2*(gamma_PN_2 + gamma_PN_3*Vars['v']))/Cons['M']**3


@njit
def TaylorT1_1p5(Cons,Vars):
    Flux = Vars['Fcal_coeff']*(Cons['Fcal_0'] + Vars['v']**2*(Cons['Fcal_2'] + Vars['v']*(Cons['Fcal_3'] + Vars['Fcal_SO_3'])))
    dEdV = -0.5*Cons['M']*Cons['nu']*Vars['v']*(2.0*Cons['E_0'] + Vars['v']**2*(4.0*Cons['E_2'] + 5.0*Vars['E_SO_3']*Vars['v']))
    Absorption = 0
    dvdt = (-Absorption - Flux)/dEdV
    dydt=np.zeros(9)
    rfrak_frame=[Vars['rfrak_frame_x'],Vars['rfrak_frame_y'],Vars['rfrak_frame_z']]
    [dydt[5],dydt[6],dydt[7]] = FrameFromAngularVelocityIntegrand(rfrak_frame, OmegaVec_1p5(Cons,Vars).vector)
    dydt[0] = dvdt
    if(Cons['EvolveSpin1']):
        dydt[1], dydt[2]= FrameFromAngularVelocity_2D_Integrand(Vars['rfrak_chi1_x'], Vars['rfrak_chi1_y'],(Cons['S_chi1'].inverse*OmegaVec_chiVec_1_1p5(Cons,Vars)*Cons['S_chi1']).vector)
    else:
        dydt[1] = 0.0
        dydt[2] = 0.0
    if(Cons['EvolveSpin2']):
        dydt[3], dydt[4] = FrameFromAngularVelocity_2D_Integrand(Vars['rfrak_chi2_x'], Vars['rfrak_chi2_y'],(Cons['S_chi2'].inverse*OmegaVec_chiVec_2_1p5(Cons,Vars)*Cons['S_chi2']).vector)
    else:
        dydt[3] = 0.0
        dydt[4] = 0.0
    dydt[8] = Vars['v']*Vars['v']*Vars['v']/Cons['M']
    return dydt      

@njit
def Recalculate_2p0(Cons,y):
    Vars={'v':y[0],'ellHat':[0,0,0,0]}
    Vars['v'] = y[0]
    Vars['rfrak_chi1_x'] = y[1]
    Vars['rfrak_chi1_y'] = y[2]
    Vars['rfrak_chi2_x'] = y[3]
    Vars['rfrak_chi2_y'] = y[4]
    Vars['rfrak_frame_x'] = y[5]
    Vars['rfrak_frame_y'] = y[6]
    Vars['rfrak_frame_z'] = y[7]
    Vars['Phi'] = y[8]
    Vars['R'] = exp(Vars['rfrak_frame_x']*Cons['xHat'] + Vars['rfrak_frame_y']*Cons['yHat'] + Vars['rfrak_frame_z']*Cons['zHat'])
    Vars['nHat'] = Vars['R']*Cons['xHat']*conjugate(Vars['R'])
    Vars['lambdaHat'] = Vars['R']*Cons['yHat']*conjugate(Vars['R'])
    Vars['ellHat'] = Vars['R']*Cons['zHat']*conjugate(Vars['R'])
    Vars['R_S1'] = exp(Vars['rfrak_chi1_x']*Cons['xHat'] + Vars['rfrak_chi1_y']*Cons['yHat'])
    Vars['R_S2'] = exp(Vars['rfrak_chi2_x']*Cons['xHat'] + Vars['rfrak_chi2_y']*Cons['yHat'])
    Vars['chiVec1'] = Cons['S_chi1']*Vars['R_S1']*Cons['zHat']*conjugate(Vars['R_S1'])*conjugate(Cons['S_chi1'])
    Vars['chiVec2'] = Cons['S_chi2']*Vars['R_S2']*Cons['zHat']*conjugate(Vars['R_S2'])*conjugate(Cons['S_chi2'])
    Vars['chi1_n'] = dot(Vars['chiVec1'].vector,Vars['nHat'].vector)
    Vars['chi1_lambda'] = dot(Vars['chiVec1'].vector,Vars['lambdaHat'].vector)
    Vars['chi1_ell'] = dot(Vars['chiVec1'].vector,Vars['ellHat'].vector)
    Vars['chi2_n'] = dot(Vars['chiVec2'].vector,Vars['nHat'].vector)
    Vars['chi2_lambda'] = dot(Vars['chiVec2'].vector,Vars['lambdaHat'].vector)
    Vars['chi2_ell'] = dot(Vars['chiVec2'].vector,Vars['ellHat'].vector)
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
    gamma_PN_0 = 1.00000000000000
    gamma_PN_4 = 1.0 - 5.41666666666667*Cons['nu']
    a_ell_2 = Vars['S_n']*(-9.66666666666667*Cons['nu'] - 10.0) + Vars['Sigma_n']*Cons['delta']*(-4.5*Cons['nu'] - 6.0)
    a_ell_0 = 7.0*Vars['S_n'] + 3.0*Vars['Sigma_n']*Cons['delta']
    gamma_PN_2 = 1.0 - 0.333333333333333*Cons['nu']
    a_ell_4 = Vars['S_n']*(5.77777777777778*Cons['nu']**2 + 14.75*Cons['nu'] + 1.5) + Vars['Sigma_n']*Cons['delta']*(2.83333333333333*Cons['nu']**2 + 9.125*Cons['nu'] + 1.5)
    gamma_PN_3 = (1.66666666666667*Vars['S_ell'] + Vars['Sigma_ell']*Cons['delta'])/Cons['M']**2
    return Vars['ellHat']*Vars['v']**3/Cons['M'] + Vars['nHat']*Vars['v']**6*(a_ell_0 + Vars['v']**2*(a_ell_2 + a_ell_4*Vars['v']**2))*(gamma_PN_0 + Vars['v']**2*(gamma_PN_2 + Vars['v']*(gamma_PN_3 + gamma_PN_4*Vars['v'])))/Cons['M']**3


@njit
def TaylorT1_2p0(Cons,Vars):
    Flux = Vars['Fcal_coeff']*(Cons['Fcal_0'] + Vars['v']**2*(Cons['Fcal_2'] + Vars['v']*(Cons['Fcal_3'] + Vars['Fcal_SO_3'] + Vars['v']*(Cons['Fcal_4'] + Vars['Fcal_SQ_4']))))
    dEdV = -0.5*Cons['M']*Cons['nu']*Vars['v']*(2.0*Cons['E_0'] + Vars['v']**2*(4.0*Cons['E_2'] + Vars['v']*(5.0*Vars['E_SO_3'] + 6.0*Vars['v']*(Cons['E_4'] + Vars['E_SQ_4']))))
    Absorption = 0
    dvdt = (-Absorption - Flux)/dEdV
    dydt=np.zeros(9)
    rfrak_frame=[Vars['rfrak_frame_x'],Vars['rfrak_frame_y'],Vars['rfrak_frame_z']]
    [dydt[5],dydt[6],dydt[7]] = FrameFromAngularVelocityIntegrand(rfrak_frame, OmegaVec_2p0(Cons,Vars).vector)
    dydt[0] = dvdt
    if(Cons['EvolveSpin1']):
        dydt[1], dydt[2]= FrameFromAngularVelocity_2D_Integrand(Vars['rfrak_chi1_x'], Vars['rfrak_chi1_y'],(Cons['S_chi1'].inverse*OmegaVec_chiVec_1_2p0(Cons,Vars)*Cons['S_chi1']).vector)
    else:
        dydt[1] = 0.0
        dydt[2] = 0.0
    if(Cons['EvolveSpin2']):
        dydt[3], dydt[4] = FrameFromAngularVelocity_2D_Integrand(Vars['rfrak_chi2_x'], Vars['rfrak_chi2_y'],(Cons['S_chi2'].inverse*OmegaVec_chiVec_2_2p0(Cons,Vars)*Cons['S_chi2']).vector)
    else:
        dydt[3] = 0.0
        dydt[4] = 0.0
    dydt[8] = Vars['v']*Vars['v']*Vars['v']/Cons['M']
    return dydt      

@njit
def Recalculate_2p5(Cons,y):
    Vars={'v':y[0],'ellHat':[0,0,0,0]}
    Vars['v'] = y[0]
    Vars['rfrak_chi1_x'] = y[1]
    Vars['rfrak_chi1_y'] = y[2]
    Vars['rfrak_chi2_x'] = y[3]
    Vars['rfrak_chi2_y'] = y[4]
    Vars['rfrak_frame_x'] = y[5]
    Vars['rfrak_frame_y'] = y[6]
    Vars['rfrak_frame_z'] = y[7]
    Vars['Phi'] = y[8]
    Vars['R'] = exp(Vars['rfrak_frame_x']*Cons['xHat'] + Vars['rfrak_frame_y']*Cons['yHat'] + Vars['rfrak_frame_z']*Cons['zHat'])
    Vars['nHat'] = Vars['R']*Cons['xHat']*conjugate(Vars['R'])
    Vars['lambdaHat'] = Vars['R']*Cons['yHat']*conjugate(Vars['R'])
    Vars['ellHat'] = Vars['R']*Cons['zHat']*conjugate(Vars['R'])
    Vars['R_S1'] = exp(Vars['rfrak_chi1_x']*Cons['xHat'] + Vars['rfrak_chi1_y']*Cons['yHat'])
    Vars['R_S2'] = exp(Vars['rfrak_chi2_x']*Cons['xHat'] + Vars['rfrak_chi2_y']*Cons['yHat'])
    Vars['chiVec1'] = Cons['S_chi1']*Vars['R_S1']*Cons['zHat']*conjugate(Vars['R_S1'])*conjugate(Cons['S_chi1'])
    Vars['chiVec2'] = Cons['S_chi2']*Vars['R_S2']*Cons['zHat']*conjugate(Vars['R_S2'])*conjugate(Cons['S_chi2'])
    Vars['chi1_n'] = dot(Vars['chiVec1'].vector,Vars['nHat'].vector)
    Vars['chi1_lambda'] = dot(Vars['chiVec1'].vector,Vars['lambdaHat'].vector)
    Vars['chi1_ell'] = dot(Vars['chiVec1'].vector,Vars['ellHat'].vector)
    Vars['chi2_n'] = dot(Vars['chiVec2'].vector,Vars['nHat'].vector)
    Vars['chi2_lambda'] = dot(Vars['chiVec2'].vector,Vars['lambdaHat'].vector)
    Vars['chi2_ell'] = dot(Vars['chiVec2'].vector,Vars['ellHat'].vector)
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
    gamma_PN_0 = 1.00000000000000
    gamma_PN_4 = 1.0 - 5.41666666666667*Cons['nu']
    a_ell_2 = Vars['S_n']*(-9.66666666666667*Cons['nu'] - 10.0) + Vars['Sigma_n']*Cons['delta']*(-4.5*Cons['nu'] - 6.0)
    a_ell_0 = 7.0*Vars['S_n'] + 3.0*Vars['Sigma_n']*Cons['delta']
    gamma_PN_2 = 1.0 - 0.333333333333333*Cons['nu']
    gamma_PN_5 = (Vars['S_ell']*(0.888888888888889*Cons['nu'] + 3.33333333333333) + 2.0*Vars['Sigma_ell']*Cons['delta'])/Cons['M']**2
    a_ell_4 = Vars['S_n']*(5.77777777777778*Cons['nu']**2 + 14.75*Cons['nu'] + 1.5) + Vars['Sigma_n']*Cons['delta']*(2.83333333333333*Cons['nu']**2 + 9.125*Cons['nu'] + 1.5)
    gamma_PN_3 = (1.66666666666667*Vars['S_ell'] + Vars['Sigma_ell']*Cons['delta'])/Cons['M']**2
    return Vars['ellHat']*Vars['v']**3/Cons['M'] + Vars['nHat']*Vars['v']**6*(a_ell_0 + Vars['v']**2*(a_ell_2 + a_ell_4*Vars['v']**2))*(gamma_PN_0 + Vars['v']**2*(gamma_PN_2 + Vars['v']*(gamma_PN_3 + Vars['v']*(gamma_PN_4 + gamma_PN_5*Vars['v']))))/Cons['M']**3


@njit
def TaylorT1_2p5(Cons,Vars):
    Flux = Vars['Fcal_coeff']*(Cons['Fcal_0'] + Vars['v']**2*(Cons['Fcal_2'] + Vars['v']*(Cons['Fcal_3'] + Vars['Fcal_SO_3'] + Vars['v']*(Cons['Fcal_4'] + Vars['Fcal_SQ_4'] + Vars['v']*(Cons['Fcal_5'] + Vars['Fcal_SO_5'])))))
    dEdV = -0.5*Cons['M']*Cons['nu']*Vars['v']*(2.0*Cons['E_0'] + Vars['v']**2*(4.0*Cons['E_2'] + Vars['v']*(5.0*Vars['E_SO_3'] + Vars['v']*(6.0*Cons['E_4'] + 7.0*Vars['E_SO_5']*Vars['v'] + 6.0*Vars['E_SQ_4']))))
    Absorption = 0
    dvdt = (-Absorption - Flux)/dEdV
    dydt=np.zeros(9)
    rfrak_frame=[Vars['rfrak_frame_x'],Vars['rfrak_frame_y'],Vars['rfrak_frame_z']]
    [dydt[5],dydt[6],dydt[7]] = FrameFromAngularVelocityIntegrand(rfrak_frame, OmegaVec_2p5(Cons,Vars).vector)
    dydt[0] = dvdt
    if(Cons['EvolveSpin1']):
        dydt[1], dydt[2]= FrameFromAngularVelocity_2D_Integrand(Vars['rfrak_chi1_x'], Vars['rfrak_chi1_y'],(Cons['S_chi1'].inverse*OmegaVec_chiVec_1_2p5(Cons,Vars)*Cons['S_chi1']).vector)
    else:
        dydt[1] = 0.0
        dydt[2] = 0.0
    if(Cons['EvolveSpin2']):
        dydt[3], dydt[4] = FrameFromAngularVelocity_2D_Integrand(Vars['rfrak_chi2_x'], Vars['rfrak_chi2_y'],(Cons['S_chi2'].inverse*OmegaVec_chiVec_2_2p5(Cons,Vars)*Cons['S_chi2']).vector)
    else:
        dydt[3] = 0.0
        dydt[4] = 0.0
    dydt[8] = Vars['v']*Vars['v']*Vars['v']/Cons['M']
    return dydt      

@njit
def Recalculate_3p0(Cons,y):
    Vars={'v':y[0],'ellHat':[0,0,0,0]}
    Vars['v'] = y[0]
    Vars['rfrak_chi1_x'] = y[1]
    Vars['rfrak_chi1_y'] = y[2]
    Vars['rfrak_chi2_x'] = y[3]
    Vars['rfrak_chi2_y'] = y[4]
    Vars['rfrak_frame_x'] = y[5]
    Vars['rfrak_frame_y'] = y[6]
    Vars['rfrak_frame_z'] = y[7]
    Vars['Phi'] = y[8]
    Vars['R'] = exp(Vars['rfrak_frame_x']*Cons['xHat'] + Vars['rfrak_frame_y']*Cons['yHat'] + Vars['rfrak_frame_z']*Cons['zHat'])
    Vars['nHat'] = Vars['R']*Cons['xHat']*conjugate(Vars['R'])
    Vars['lambdaHat'] = Vars['R']*Cons['yHat']*conjugate(Vars['R'])
    Vars['ellHat'] = Vars['R']*Cons['zHat']*conjugate(Vars['R'])
    Vars['R_S1'] = exp(Vars['rfrak_chi1_x']*Cons['xHat'] + Vars['rfrak_chi1_y']*Cons['yHat'])
    Vars['R_S2'] = exp(Vars['rfrak_chi2_x']*Cons['xHat'] + Vars['rfrak_chi2_y']*Cons['yHat'])
    Vars['chiVec1'] = Cons['S_chi1']*Vars['R_S1']*Cons['zHat']*conjugate(Vars['R_S1'])*conjugate(Cons['S_chi1'])
    Vars['chiVec2'] = Cons['S_chi2']*Vars['R_S2']*Cons['zHat']*conjugate(Vars['R_S2'])*conjugate(Cons['S_chi2'])
    Vars['chi1_n'] = dot(Vars['chiVec1'].vector,Vars['nHat'].vector)
    Vars['chi1_lambda'] = dot(Vars['chiVec1'].vector,Vars['lambdaHat'].vector)
    Vars['chi1_ell'] = dot(Vars['chiVec1'].vector,Vars['ellHat'].vector)
    Vars['chi2_n'] = dot(Vars['chiVec2'].vector,Vars['nHat'].vector)
    Vars['chi2_lambda'] = dot(Vars['chiVec2'].vector,Vars['lambdaHat'].vector)
    Vars['chi2_ell'] = dot(Vars['chiVec2'].vector,Vars['ellHat'].vector)
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
    gamma_PN_6 = 0.0123456790123457*Cons['nu']**3 + 6.36111111111111*Cons['nu']**2 - 2.98177812235564*Cons['nu'] + 1.0
    gamma_PN_0 = 1.00000000000000
    gamma_PN_4 = 1.0 - 5.41666666666667*Cons['nu']
    a_ell_2 = Vars['S_n']*(-9.66666666666667*Cons['nu'] - 10.0) + Vars['Sigma_n']*Cons['delta']*(-4.5*Cons['nu'] - 6.0)
    a_ell_0 = 7.0*Vars['S_n'] + 3.0*Vars['Sigma_n']*Cons['delta']
    gamma_PN_2 = 1.0 - 0.333333333333333*Cons['nu']
    gamma_PN_5 = (Vars['S_ell']*(0.888888888888889*Cons['nu'] + 3.33333333333333) + 2.0*Vars['Sigma_ell']*Cons['delta'])/Cons['M']**2
    a_ell_4 = Vars['S_n']*(5.77777777777778*Cons['nu']**2 + 14.75*Cons['nu'] + 1.5) + Vars['Sigma_n']*Cons['delta']*(2.83333333333333*Cons['nu']**2 + 9.125*Cons['nu'] + 1.5)
    gamma_PN_3 = (1.66666666666667*Vars['S_ell'] + Vars['Sigma_ell']*Cons['delta'])/Cons['M']**2
    return Vars['ellHat']*Vars['v']**3/Cons['M'] + Vars['nHat']*Vars['v']**6*(a_ell_0 + Vars['v']**2*(a_ell_2 + a_ell_4*Vars['v']**2))*(gamma_PN_0 + Vars['v']**2*(gamma_PN_2 + Vars['v']*(gamma_PN_3 + Vars['v']*(gamma_PN_4 + Vars['v']*(gamma_PN_5 + gamma_PN_6*Vars['v'])))))/Cons['M']**3


@njit
def TaylorT1_3p0(Cons,Vars):
    Flux = Vars['Fcal_coeff']*(Cons['Fcal_0'] + Vars['v']**2*(Cons['Fcal_2'] + Vars['v']*(Cons['Fcal_3'] + Vars['Fcal_SO_3'] + Vars['v']*(Cons['Fcal_4'] + Vars['Fcal_SQ_4'] + Vars['v']*(Cons['Fcal_5'] + Vars['Fcal_SO_5'] + Vars['v']*(Cons['Fcal_6'] + Vars['Fcal_SO_6'] + Cons['Fcal_lnv_6']*Vars['logv']))))))
    dEdV = -0.5*Cons['M']*Cons['nu']*Vars['v']*(2.0*Cons['E_0'] + Vars['v']**2*(4.0*Cons['E_2'] + Vars['v']*(5.0*Vars['E_SO_3'] + Vars['v']*(6.0*Cons['E_4'] + 6.0*Vars['E_SQ_4'] + Vars['v']*(8.0*Cons['E_6']*Vars['v'] + 7.0*Vars['E_SO_5'])))))
    Absorption = 0
    dvdt = (-Absorption - Flux)/dEdV
    dydt=np.zeros(9)
    rfrak_frame=[Vars['rfrak_frame_x'],Vars['rfrak_frame_y'],Vars['rfrak_frame_z']]
    [dydt[5],dydt[6],dydt[7]] = FrameFromAngularVelocityIntegrand(rfrak_frame, OmegaVec_3p0(Cons,Vars).vector)
    dydt[0] = dvdt
    if(Cons['EvolveSpin1']):
        dydt[1], dydt[2]= FrameFromAngularVelocity_2D_Integrand(Vars['rfrak_chi1_x'], Vars['rfrak_chi1_y'],(Cons['S_chi1'].inverse*OmegaVec_chiVec_1_3p0(Cons,Vars)*Cons['S_chi1']).vector)
    else:
        dydt[1] = 0.0
        dydt[2] = 0.0
    if(Cons['EvolveSpin2']):
        dydt[3], dydt[4] = FrameFromAngularVelocity_2D_Integrand(Vars['rfrak_chi2_x'], Vars['rfrak_chi2_y'],(Cons['S_chi2'].inverse*OmegaVec_chiVec_2_3p0(Cons,Vars)*Cons['S_chi2']).vector)
    else:
        dydt[3] = 0.0
        dydt[4] = 0.0
    dydt[8] = Vars['v']*Vars['v']*Vars['v']/Cons['M']
    return dydt      

@njit
def Recalculate_3p5(Cons,y):
    Vars={'v':y[0],'ellHat':[0,0,0,0]}
    Vars['v'] = y[0]
    Vars['rfrak_chi1_x'] = y[1]
    Vars['rfrak_chi1_y'] = y[2]
    Vars['rfrak_chi2_x'] = y[3]
    Vars['rfrak_chi2_y'] = y[4]
    Vars['rfrak_frame_x'] = y[5]
    Vars['rfrak_frame_y'] = y[6]
    Vars['rfrak_frame_z'] = y[7]
    Vars['Phi'] = y[8]
    Vars['R'] = exp(Vars['rfrak_frame_x']*Cons['xHat'] + Vars['rfrak_frame_y']*Cons['yHat'] + Vars['rfrak_frame_z']*Cons['zHat'])
    Vars['nHat'] = Vars['R']*Cons['xHat']*conjugate(Vars['R'])
    Vars['lambdaHat'] = Vars['R']*Cons['yHat']*conjugate(Vars['R'])
    Vars['ellHat'] = Vars['R']*Cons['zHat']*conjugate(Vars['R'])
    Vars['R_S1'] = exp(Vars['rfrak_chi1_x']*Cons['xHat'] + Vars['rfrak_chi1_y']*Cons['yHat'])
    Vars['R_S2'] = exp(Vars['rfrak_chi2_x']*Cons['xHat'] + Vars['rfrak_chi2_y']*Cons['yHat'])
    Vars['chiVec1'] = Cons['S_chi1']*Vars['R_S1']*Cons['zHat']*conjugate(Vars['R_S1'])*conjugate(Cons['S_chi1'])
    Vars['chiVec2'] = Cons['S_chi2']*Vars['R_S2']*Cons['zHat']*conjugate(Vars['R_S2'])*conjugate(Cons['S_chi2'])
    Vars['chi1_n'] = dot(Vars['chiVec1'].vector,Vars['nHat'].vector)
    Vars['chi1_lambda'] = dot(Vars['chiVec1'].vector,Vars['lambdaHat'].vector)
    Vars['chi1_ell'] = dot(Vars['chiVec1'].vector,Vars['ellHat'].vector)
    Vars['chi2_n'] = dot(Vars['chiVec2'].vector,Vars['nHat'].vector)
    Vars['chi2_lambda'] = dot(Vars['chiVec2'].vector,Vars['lambdaHat'].vector)
    Vars['chi2_ell'] = dot(Vars['chiVec2'].vector,Vars['ellHat'].vector)
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
    gamma_PN_6 = 0.0123456790123457*Cons['nu']**3 + 6.36111111111111*Cons['nu']**2 - 2.98177812235564*Cons['nu'] + 1.0
    gamma_PN_0 = 1.00000000000000
    gamma_PN_4 = 1.0 - 5.41666666666667*Cons['nu']
    gamma_PN_7 = (Vars['S_ell']*(-6.0*Cons['nu']**2 - 10.5833333333333*Cons['nu'] + 5.0) - 2.66666666666667*Vars['Sigma_ell']*Cons['delta']*Cons['nu']**2 + Vars['Sigma_ell']*Cons['delta']*(3.0 - 10.1666666666667*Cons['nu']))/Cons['M']**2
    a_ell_2 = Vars['S_n']*(-9.66666666666667*Cons['nu'] - 10.0) + Vars['Sigma_n']*Cons['delta']*(-4.5*Cons['nu'] - 6.0)
    a_ell_0 = 7.0*Vars['S_n'] + 3.0*Vars['Sigma_n']*Cons['delta']
    gamma_PN_2 = 1.0 - 0.333333333333333*Cons['nu']
    gamma_PN_5 = (Vars['S_ell']*(0.888888888888889*Cons['nu'] + 3.33333333333333) + 2.0*Vars['Sigma_ell']*Cons['delta'])/Cons['M']**2
    a_ell_4 = Vars['S_n']*(5.77777777777778*Cons['nu']**2 + 14.75*Cons['nu'] + 1.5) + Vars['Sigma_n']*Cons['delta']*(2.83333333333333*Cons['nu']**2 + 9.125*Cons['nu'] + 1.5)
    gamma_PN_3 = (1.66666666666667*Vars['S_ell'] + Vars['Sigma_ell']*Cons['delta'])/Cons['M']**2
    return Vars['ellHat']*Vars['v']**3/Cons['M'] + Vars['nHat']*Vars['v']**6*(a_ell_0 + Vars['v']**2*(a_ell_2 + a_ell_4*Vars['v']**2))*(gamma_PN_0 + Vars['v']**2*(gamma_PN_2 + Vars['v']*(gamma_PN_3 + Vars['v']*(gamma_PN_4 + Vars['v']*(gamma_PN_5 + Vars['v']*(gamma_PN_6 + gamma_PN_7*Vars['v']))))))/Cons['M']**3


@njit
def TaylorT1_3p5(Cons,Vars):
    Flux = Vars['Fcal_coeff']*(Cons['Fcal_0'] + Vars['v']**2*(Cons['Fcal_2'] + Vars['v']*(Cons['Fcal_3'] + Vars['Fcal_SO_3'] + Vars['v']*(Cons['Fcal_4'] + Vars['Fcal_SQ_4'] + Vars['v']*(Cons['Fcal_5'] + Vars['Fcal_SO_5'] + Vars['v']*(Cons['Fcal_6'] + Vars['Fcal_SO_6'] + Cons['Fcal_lnv_6']*Vars['logv'] + Vars['v']*(Cons['Fcal_7'] + Vars['Fcal_SO_7'])))))))
    dEdV = -0.5*Cons['M']*Cons['nu']*Vars['v']*(2.0*Cons['E_0'] + Vars['v']**2*(4.0*Cons['E_2'] + Vars['v']*(5.0*Vars['E_SO_3'] + Vars['v']*(6.0*Cons['E_4'] + 6.0*Vars['E_SQ_4'] + Vars['v']*(7.0*Vars['E_SO_5'] + Vars['v']*(8.0*Cons['E_6'] + 9.0*Vars['E_SO_7']*Vars['v']))))))
    Absorption = 0
    dvdt = (-Absorption - Flux)/dEdV
    dydt=np.zeros(9)
    rfrak_frame=[Vars['rfrak_frame_x'],Vars['rfrak_frame_y'],Vars['rfrak_frame_z']]
    [dydt[5],dydt[6],dydt[7]] = FrameFromAngularVelocityIntegrand(rfrak_frame, OmegaVec_3p5(Cons,Vars).vector)
    dydt[0] = dvdt
    if(Cons['EvolveSpin1']):
        dydt[1], dydt[2]= FrameFromAngularVelocity_2D_Integrand(Vars['rfrak_chi1_x'], Vars['rfrak_chi1_y'],(Cons['S_chi1'].inverse*OmegaVec_chiVec_1_3p5(Cons,Vars)*Cons['S_chi1']).vector)
    else:
        dydt[1] = 0.0
        dydt[2] = 0.0
    if(Cons['EvolveSpin2']):
        dydt[3], dydt[4] = FrameFromAngularVelocity_2D_Integrand(Vars['rfrak_chi2_x'], Vars['rfrak_chi2_y'],(Cons['S_chi2'].inverse*OmegaVec_chiVec_2_3p5(Cons,Vars)*Cons['S_chi2']).vector)
    else:
        dydt[3] = 0.0
        dydt[4] = 0.0
    dydt[8] = Vars['v']*Vars['v']*Vars['v']/Cons['M']
    return dydt      

@njit
def Recalculate_4p0(Cons,y):
    Vars={'v':y[0],'ellHat':[0,0,0,0]}
    Vars['v'] = y[0]
    Vars['rfrak_chi1_x'] = y[1]
    Vars['rfrak_chi1_y'] = y[2]
    Vars['rfrak_chi2_x'] = y[3]
    Vars['rfrak_chi2_y'] = y[4]
    Vars['rfrak_frame_x'] = y[5]
    Vars['rfrak_frame_y'] = y[6]
    Vars['rfrak_frame_z'] = y[7]
    Vars['Phi'] = y[8]
    Vars['R'] = exp(Vars['rfrak_frame_x']*Cons['xHat'] + Vars['rfrak_frame_y']*Cons['yHat'] + Vars['rfrak_frame_z']*Cons['zHat'])
    Vars['nHat'] = Vars['R']*Cons['xHat']*conjugate(Vars['R'])
    Vars['lambdaHat'] = Vars['R']*Cons['yHat']*conjugate(Vars['R'])
    Vars['ellHat'] = Vars['R']*Cons['zHat']*conjugate(Vars['R'])
    Vars['R_S1'] = exp(Vars['rfrak_chi1_x']*Cons['xHat'] + Vars['rfrak_chi1_y']*Cons['yHat'])
    Vars['R_S2'] = exp(Vars['rfrak_chi2_x']*Cons['xHat'] + Vars['rfrak_chi2_y']*Cons['yHat'])
    Vars['chiVec1'] = Cons['S_chi1']*Vars['R_S1']*Cons['zHat']*conjugate(Vars['R_S1'])*conjugate(Cons['S_chi1'])
    Vars['chiVec2'] = Cons['S_chi2']*Vars['R_S2']*Cons['zHat']*conjugate(Vars['R_S2'])*conjugate(Cons['S_chi2'])
    Vars['chi1_n'] = dot(Vars['chiVec1'].vector,Vars['nHat'].vector)
    Vars['chi1_lambda'] = dot(Vars['chiVec1'].vector,Vars['lambdaHat'].vector)
    Vars['chi1_ell'] = dot(Vars['chiVec1'].vector,Vars['ellHat'].vector)
    Vars['chi2_n'] = dot(Vars['chiVec2'].vector,Vars['nHat'].vector)
    Vars['chi2_lambda'] = dot(Vars['chiVec2'].vector,Vars['lambdaHat'].vector)
    Vars['chi2_ell'] = dot(Vars['chiVec2'].vector,Vars['ellHat'].vector)
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
    gamma_PN_6 = 0.0123456790123457*Cons['nu']**3 + 6.36111111111111*Cons['nu']**2 - 2.98177812235564*Cons['nu'] + 1.0
    gamma_PN_0 = 1.00000000000000
    gamma_PN_4 = 1.0 - 5.41666666666667*Cons['nu']
    gamma_PN_7 = (Vars['S_ell']*(-6.0*Cons['nu']**2 - 10.5833333333333*Cons['nu'] + 5.0) - 2.66666666666667*Vars['Sigma_ell']*Cons['delta']*Cons['nu']**2 + Vars['Sigma_ell']*Cons['delta']*(3.0 - 10.1666666666667*Cons['nu']))/Cons['M']**2
    a_ell_2 = Vars['S_n']*(-9.66666666666667*Cons['nu'] - 10.0) + Vars['Sigma_n']*Cons['delta']*(-4.5*Cons['nu'] - 6.0)
    a_ell_0 = 7.0*Vars['S_n'] + 3.0*Vars['Sigma_n']*Cons['delta']
    gamma_PN_2 = 1.0 - 0.333333333333333*Cons['nu']
    gamma_PN_5 = (Vars['S_ell']*(0.888888888888889*Cons['nu'] + 3.33333333333333) + 2.0*Vars['Sigma_ell']*Cons['delta'])/Cons['M']**2
    a_ell_4 = Vars['S_n']*(5.77777777777778*Cons['nu']**2 + 14.75*Cons['nu'] + 1.5) + Vars['Sigma_n']*Cons['delta']*(2.83333333333333*Cons['nu']**2 + 9.125*Cons['nu'] + 1.5)
    gamma_PN_3 = (1.66666666666667*Vars['S_ell'] + Vars['Sigma_ell']*Cons['delta'])/Cons['M']**2
    return Vars['ellHat']*Vars['v']**3/Cons['M'] + Vars['nHat']*Vars['v']**6*(a_ell_0 + Vars['v']**2*(a_ell_2 + a_ell_4*Vars['v']**2))*(gamma_PN_0 + Vars['v']**2*(gamma_PN_2 + Vars['v']*(gamma_PN_3 + Vars['v']*(gamma_PN_4 + Vars['v']*(gamma_PN_5 + Vars['v']*(gamma_PN_6 + gamma_PN_7*Vars['v']))))))/Cons['M']**3


@njit
def TaylorT1_4p0(Cons,Vars):
    Flux = Vars['Fcal_coeff']*(Cons['Fcal_0'] + Vars['v']**2*(Cons['Fcal_2'] + Vars['v']*(Cons['Fcal_3'] + Vars['Fcal_SO_3'] + Vars['v']*(Cons['Fcal_4'] + Vars['Fcal_SQ_4'] + Vars['v']*(Cons['Fcal_5'] + Vars['Fcal_SO_5'] + Vars['v']*(Cons['Fcal_6'] + Vars['Fcal_SO_6'] + Cons['Fcal_lnv_6']*Vars['logv'] + Vars['v']*(Cons['Fcal_7'] + Vars['Fcal_SO_7'] + Vars['v']*(Cons['Fcal_8'] + Vars['Fcal_SO_8'] + Cons['Fcal_lnv_8']*Vars['logv']))))))))
    dEdV = -0.5*Cons['M']*Cons['nu']*Vars['v']*(2.0*Cons['E_0'] + Vars['v']**2*(4.0*Cons['E_2'] + Vars['v']*(5.0*Vars['E_SO_3'] + Vars['v']*(6.0*Cons['E_4'] + 6.0*Vars['E_SQ_4'] + Vars['v']*(7.0*Vars['E_SO_5'] + Vars['v']*(8.0*Cons['E_6'] + Vars['v']*(9.0*Vars['E_SO_7'] + Vars['v']*(10.0*Cons['E_8'] + Cons['E_lnv_8']*(10.0*Vars['logv'] + 1.0)))))))))
    Absorption = 0
    dvdt = (-Absorption - Flux)/dEdV
    dydt=np.zeros(9)
    rfrak_frame=[Vars['rfrak_frame_x'],Vars['rfrak_frame_y'],Vars['rfrak_frame_z']]
    [dydt[5],dydt[6],dydt[7]] = FrameFromAngularVelocityIntegrand(rfrak_frame, OmegaVec_4p0(Cons,Vars).vector)
    dydt[0] = dvdt
    if(Cons['EvolveSpin1']):
        dydt[1], dydt[2]= FrameFromAngularVelocity_2D_Integrand(Vars['rfrak_chi1_x'], Vars['rfrak_chi1_y'],(Cons['S_chi1'].inverse*OmegaVec_chiVec_1_4p0(Cons,Vars)*Cons['S_chi1']).vector)
    else:
        dydt[1] = 0.0
        dydt[2] = 0.0
    if(Cons['EvolveSpin2']):
        dydt[3], dydt[4] = FrameFromAngularVelocity_2D_Integrand(Vars['rfrak_chi2_x'], Vars['rfrak_chi2_y'],(Cons['S_chi2'].inverse*OmegaVec_chiVec_2_4p0(Cons,Vars)*Cons['S_chi2']).vector)
    else:
        dydt[3] = 0.0
        dydt[4] = 0.0
    dydt[8] = Vars['v']*Vars['v']*Vars['v']/Cons['M']
    return dydt      

class PNEv:
    def Integrand(t,y):
        #Cons=Dict()
        #for i, j in PNEv.Cons.items():
        #    Cons[i]=j
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
        
    def Evolution(xHat_i, yHat_i, zHat_i, M1_i, M2_i, v_i, S_chi1_i, S_chi2_i, rfrak_chi1_x_i, rfrak_chi1_y_i, rfrak_chi2_x_i, rfrak_chi2_y_i, rfrak_frame_x_i, rfrak_frame_y_i, rfrak_frame_z_i, t_PNStart=False, t_PNEnd=False, PNEvolutionOrder=3.5, TaylorTn=1, StepsPerOrbit=32, ForwardInTime=True, tol=1e-12, MinStep=1e-7): 
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
        PNEv.Cons['EvolveSpin1']=(PNEv.Cons['S_chi1']*conjugate(PNEv.Cons['S_chi1'])).norm>1e-12
        PNEv.Cons['EvolveSpin2']=(PNEv.Cons['S_chi2']*conjugate(PNEv.Cons['S_chi2'])).norm>1e-12
    
        def terminate(t,y):
            return 1.0*PNEv.terminal1*PNEv.terminal2
        terminate.terminal=True
        TMerger=5.0/(256.0*PNEv.Cons['nu']*v_i**8)
        TEnd=TMerger
        if t_PNEnd:
            TEnd=t_PNEnd
        time=[0.0]
        while time[-1]<TEnd and 2*PNEv.Cons['M']*(256*PNEv.Cons['nu']*(TMerger-time[-1])/5)**(3/8)/StepsPerOrbit>MinStep:
            time.append(time[-1]+2*PNEv.Cons['M']*(256*PNEv.Cons['nu']*(TMerger-time[-1])/5)**(3/8)/StepsPerOrbit)
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
                time.append(time[-1]-2*PNEv.Cons['M']*(256*PNEv.Cons['nu']*(TMerger-time[-1])/5)**(3/8)/StepsPerOrbit)
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
