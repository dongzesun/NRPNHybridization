# File produced automatically by PNCodeGen.ipynb
from scipy.integrate import solve_ivp
import numpy as np
from numpy import dot, cross, log, sqrt, pi
from numpy import euler_gamma as EulerGamma
from numba import jit, njit, float64, boolean
from numba.experimental import jitclass
from scipy.interpolate import InterpolatedUnivariateSpline as Spline
from scipy.special import zeta
import quaternionic

qmul = njit(quaternionic.algebra.multiply)
qexp=njit(quaternionic.algebra.exp)
qconj=njit(quaternionic.algebra.conj)
qinverse=njit(quaternionic.algebra.reciprocal)

@njit(cache=True)
def mul(A,B):
    C=np.empty(4)
    qmul(A,B,C)
    return C
    
@njit(cache=True)
def exp(A):
    B=np.empty(4)
    qexp(A,B)
    return B
    
@njit(cache=True)
def conjugate(A):
    B=np.empty(4)
    qconj(A,B)
    return B
    
@njit(cache=True)
def inverse(A):
    B=np.empty(4)
    qinverse(A,B)
    return B
    
@njit(cache=True)
def normalized(A):
    return A/np.linalg.norm(A)

@njit(cache=True)
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

@njit(cache=True)
def FrameFromAngularVelocityIntegrand(R, Omega):
    return 0.5*mul(np.append(0.0,Omega),R)  

ConsSpec=[('wHat', float64[:]),('xHat', float64[:]),('yHat', float64[:]),('zHat', float64[:]),('M1', float64[:]),('M2', float64[:]),('S_chi1', float64[:]),('S_chi2', float64[:]),('lambda_1', float64[:]),('lambda_2', float64[:]),('kappa_1', float64[:]),('kappa_2', float64[:]),('M', float64[:]),('delta', float64[:]),('nu', float64[:]),('kappa_p', float64[:]),('kappa_m', float64[:]),('lambda_p', float64[:]),('lambda_m', float64[:]),('chi1chi1', float64[:]),('chi1chi2', float64[:]),('chi2chi2', float64[:]),('Fcal_0', float64[:]),('Fcal_2', float64[:]),('Fcal_3', float64[:]),('Fcal_4', float64[:]),('Fcal_5', float64[:]),('Fcal_6', float64[:]),('Fcal_lnv_6', float64[:]),('Fcal_7', float64[:]),('Fcal_8', float64[:]),('Fcal_lnv_8', float64[:]),('Fcal_9', float64[:]),('Fcal_lnv_9', float64[:]),('Fcal_10', float64[:]),('Fcal_lnv_10', float64[:]),('Fcal_11', float64[:]),('Fcal_lnv_11', float64[:]),('Fcal_12', float64[:]),('Fcal_lnv_12', float64[:]),('Fcal_lnv2_12', float64[:]),('E_0', float64[:]),('E_2', float64[:]),('E_4', float64[:]),('E_6', float64[:]),('E_8', float64[:]),('E_lnv_8', float64[:]),('E_10', float64[:]),('E_lnv_10', float64[:]),('E_11', float64[:]),('E_12', float64[:]),('E_lnv_12', float64[:]),('MDot_Alvi_8', float64[:]),('EvolveSpin1',boolean),('EvolveSpin2',boolean)]
@jitclass(ConsSpec)
class Cons:
    def __init__(self,wHat,xHat,yHat,zHat,M1,M2,S_chi1,S_chi2,lambda_1,lambda_2,kappa_1,kappa_2,M,delta,nu,kappa_p,kappa_m,lambda_p,lambda_m,chi1chi1,chi1chi2,chi2chi2,Fcal_0,Fcal_2,Fcal_3,Fcal_4,Fcal_5,Fcal_6,Fcal_lnv_6,Fcal_7,Fcal_8,Fcal_lnv_8,Fcal_9,Fcal_lnv_9,Fcal_10,Fcal_lnv_10,Fcal_11,Fcal_lnv_11,Fcal_12,Fcal_lnv_12,Fcal_lnv2_12,E_0,E_2,E_4,E_6,E_8,E_lnv_8,E_10,E_lnv_10,E_11,E_12,E_lnv_12,MDot_Alvi_8,EvolveSpin1,EvolveSpin2):
        self.wHat=wHat
        self.xHat=xHat
        self.yHat=yHat
        self.zHat=zHat
        self.M1=M1
        self.M2=M2
        self.S_chi1=S_chi1
        self.S_chi2=S_chi2
        self.lambda_1=lambda_1
        self.lambda_2=lambda_2
        self.kappa_1=kappa_1
        self.kappa_2=kappa_2
        self.M=M
        self.delta=delta
        self.nu=nu
        self.kappa_p=kappa_p
        self.kappa_m=kappa_m
        self.lambda_p=lambda_p
        self.lambda_m=lambda_m
        self.chi1chi1=chi1chi1
        self.chi1chi2=chi1chi2
        self.chi2chi2=chi2chi2
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
        self.MDot_Alvi_8=MDot_Alvi_8
        self.EvolveSpin1=EvolveSpin1
        self.EvolveSpin2=EvolveSpin2

VarsSpec=[('v', float64[:]),('rfrak_chi1', float64[:]),('rfrak_chi2', float64[:]),('R', float64[:]),('nHat', float64[:]),('lambdaHat', float64[:]),('ellHat', float64[:]),('R_S1', float64[:]),('R_S2', float64[:]),('chiVec1', float64[:]),('chiVec2', float64[:]),('chi1_n', float64[:]),('chi1_lambda', float64[:]),('chi1_ell', float64[:]),('chi2_n', float64[:]),('chi2_lambda', float64[:]),('chi2_ell', float64[:]),('S', float64[:]),('S_ell', float64[:]),('S_n', float64[:]),('S_lambda', float64[:]),('Sigma', float64[:]),('Sigma_ell', float64[:]),('Sigma_n', float64[:]),('Sigma_lambda', float64[:]),('SS', float64[:]),('SigmaSigma', float64[:]),('SSigma', float64[:]),('S1', float64[:]),('S1_n', float64[:]),('S1_lambda', float64[:]),('S2', float64[:]),('S2_n', float64[:]),('S2_lambda', float64[:]),('chi_s_ell', float64[:]),('chi_a_ell', float64[:]),('logv', float64[:]),('Fcal_coeff', float64[:]),('Fcal_SQ_4', float64[:]),('Fcal_SQ_6', float64[:]),('Fcal_SQ_7', float64[:]),('Fcal_SO_3', float64[:]),('Fcal_SO_5', float64[:]),('Fcal_SO_6', float64[:]),('Fcal_SO_7', float64[:]),('Fcal_SO_8', float64[:]),('E_SQ_4', float64[:]),('E_SQ_6', float64[:]),('E_SQ_7', float64[:]),('E_SO_3', float64[:]),('E_SO_5', float64[:]),('E_SO_7', float64[:]),('MDot_Alvi_5', float64[:])]
@jitclass(VarsSpec)
class Vars:
    def __init__(self,v,rfrak_chi1,rfrak_chi2,R,nHat,lambdaHat,ellHat,R_S1,R_S2,chiVec1,chiVec2,chi1_n,chi1_lambda,chi1_ell,chi2_n,chi2_lambda,chi2_ell,S,S_ell,S_n,S_lambda,Sigma,Sigma_ell,Sigma_n,Sigma_lambda,SS,SigmaSigma,SSigma,S1,S1_n,S1_lambda,S2,S2_n,S2_lambda,chi_s_ell,chi_a_ell,logv,Fcal_coeff,Fcal_SQ_4,Fcal_SQ_6,Fcal_SQ_7,Fcal_SO_3,Fcal_SO_5,Fcal_SO_6,Fcal_SO_7,Fcal_SO_8,E_SQ_4,E_SQ_6,E_SQ_7,E_SO_3,E_SO_5,E_SO_7,MDot_Alvi_5):
        self.v=v
        self.rfrak_chi1=rfrak_chi1
        self.rfrak_chi2=rfrak_chi2
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
        self.S=S
        self.S_ell=S_ell
        self.S_n=S_n
        self.S_lambda=S_lambda
        self.Sigma=Sigma
        self.Sigma_ell=Sigma_ell
        self.Sigma_n=Sigma_n
        self.Sigma_lambda=Sigma_lambda
        self.SS=SS
        self.SigmaSigma=SigmaSigma
        self.SSigma=SSigma
        self.S1=S1
        self.S1_n=S1_n
        self.S1_lambda=S1_lambda
        self.S2=S2
        self.S2_n=S2_n
        self.S2_lambda=S2_lambda
        self.chi_s_ell=chi_s_ell
        self.chi_a_ell=chi_a_ell
        self.logv=logv
        self.Fcal_coeff=Fcal_coeff
        self.Fcal_SQ_4=Fcal_SQ_4
        self.Fcal_SQ_6=Fcal_SQ_6
        self.Fcal_SQ_7=Fcal_SQ_7
        self.Fcal_SO_3=Fcal_SO_3
        self.Fcal_SO_5=Fcal_SO_5
        self.Fcal_SO_6=Fcal_SO_6
        self.Fcal_SO_7=Fcal_SO_7
        self.Fcal_SO_8=Fcal_SO_8
        self.E_SQ_4=E_SQ_4
        self.E_SQ_6=E_SQ_6
        self.E_SQ_7=E_SQ_7
        self.E_SO_3=E_SO_3
        self.E_SO_5=E_SO_5
        self.E_SO_7=E_SO_7
        self.MDot_Alvi_5=MDot_Alvi_5

@njit(cache=True)
def Initialization(Cons, wHat_i, xHat_i, yHat_i, zHat_i, M1_i, M2_i, v_i, S_chi1_i, S_chi2_i, lambda_1_i, lambda_2_i, kappa_1_i, kappa_2_i): 
    Cons.wHat=wHat_i
    Cons.xHat=xHat_i
    Cons.yHat=yHat_i
    Cons.zHat=zHat_i
    Cons.M1=np.array([M1_i])
    Cons.M2=np.array([M2_i])
    Cons.S_chi1=S_chi1_i
    Cons.S_chi2=S_chi2_i
    rfrak_chi1=np.array([0.0,0.0])
    rfrak_chi2=np.array([0.0,0.0])
    Cons.lambda_1=np.array([lambda_1_i])
    Cons.lambda_2=np.array([lambda_2_i])
    Cons.kappa_1=np.array([kappa_1_i])
    Cons.kappa_2=np.array([kappa_2_i])
    Cons.M=Cons.M1 + Cons.M2
    Cons.delta=(Cons.M1 - Cons.M2)/Cons.M
    Cons.nu=Cons.M1*Cons.M2/Cons.M**2
    Cons.kappa_p=Cons.kappa_1 + Cons.kappa_2
    Cons.kappa_m=Cons.kappa_1 - Cons.kappa_2
    Cons.lambda_p=Cons.lambda_1 + Cons.lambda_2
    Cons.lambda_m=Cons.lambda_1 - Cons.lambda_2
    R_S1=exp(rfrak_chi1[0]*Cons.xHat + rfrak_chi1[1]*Cons.yHat)
    R_S2=exp(rfrak_chi2[0]*Cons.xHat + rfrak_chi2[1]*Cons.yHat)
    chiVec1=mul(mul(mul(Cons.S_chi1,R_S1),Cons.zHat),mul(conjugate(R_S1),conjugate(Cons.S_chi1)))
    chiVec2=mul(mul(mul(Cons.S_chi2,R_S2),Cons.zHat),mul(conjugate(R_S2),conjugate(Cons.S_chi2)))
    Cons.chi1chi1=np.array([dot(chiVec1[1:],chiVec1[1:])])
    Cons.chi1chi2=np.array([dot(chiVec1[1:],chiVec2[1:])])
    Cons.chi2chi2=np.array([dot(chiVec2[1:],chiVec2[1:])])
    Cons.Fcal_0=np.array([1.0])
    Cons.Fcal_2=-35*Cons.nu/12 - 1247/336
    Cons.Fcal_3=np.array([4*pi])
    Cons.Fcal_4=65*Cons.nu**2/18 + 9271*Cons.nu/504 - 44711/9072
    Cons.Fcal_5=pi*(-583*Cons.nu/24 - 8191/672)
    Cons.Fcal_6=-775*Cons.nu**3/324 - 94403*Cons.nu**2/3024 + Cons.nu*(-134543/7776 + 41*pi**2/48) - 1712*log(4)/105 - 1712*EulerGamma/105 + 16*pi**2/3 + 6643739519/69854400
    Cons.Fcal_lnv_6=np.array([-1712/105])
    Cons.Fcal_7=pi*(193385*Cons.nu**2/3024 + 214745*Cons.nu/1728 - 16285/504)
    Cons.Fcal_8=5*Cons.nu**4/6 + 6875*Cons.nu**3/504 + Cons.nu**2*(1607125/6804 - 3157*pi**2/384) + Cons.nu*(-1452202403629/1466942400 - 267127*pi**2/4608 + 41478*EulerGamma/245 + 47385*log(3)/392 + 479062*log(2)/2205) - 1369*pi**2/126 - 323105549467/3178375200 - 47385*log(3)/1568 + 232597*EulerGamma/4410 + 39931*log(2)/294
    Cons.Fcal_lnv_8=41478*Cons.nu/245 + 232597/4410
    Cons.Fcal_9=pi*(-3719141*Cons.nu**3/38016 - 133112905*Cons.nu**2/290304 + Cons.nu*(41*pi**2/12 + 2062241/22176) - 13696*log(2)/105 - 6848*EulerGamma/105 + 265978667519/745113600)
    Cons.Fcal_lnv_9=np.array([-6848*pi/105])
    Cons.Fcal_10=np.array([-2500861660823683/2831932303200 - 424223*pi**2/6804 - 83217611*log(2)/1122660 + 916628467*EulerGamma/7858620 + 47385*log(3)/196])
    Cons.Fcal_lnv_10=np.array([916628467/7858620])
    Cons.Fcal_11=np.array([-142155*pi*log(3)/784 + 8399309750401*pi/101708006400 + 177293*EulerGamma*pi/1176 + 8521283*pi*log(2)/17640])
    Cons.Fcal_lnv_11=np.array([177293*pi/1176])
    Cons.Fcal_12=np.array([-271272899815409*log(2)/157329572400 - 54784*pi**2*log(2)/315 - 246137536815857*EulerGamma/157329572400 - 437114506833*log(3)/789268480 - 256*pi**4/45 - 27392*EulerGamma*pi**2/315 - 37744140625*log(5)/260941824 + 1465472*EulerGamma**2/11025 + 5861888*EulerGamma*log(2)/11025 + 5861888*log(2)**2/11025 + 3118.73176332243 + 3803225263*pi**2/10478160])
    Cons.Fcal_lnv_12=np.array([-246137536815857/157329572400 - 27392*pi**2/315 + 2930944*EulerGamma/11025 + 5861888*log(2)/11025])
    Cons.Fcal_lnv2_12=np.array([1465472/11025])
    Cons.E_0=np.array([1.0])
    Cons.E_2=-Cons.nu/12 - 3/4
    Cons.E_4=-Cons.nu**2/24 + 19*Cons.nu/8 - 27/8
    Cons.E_6=-35*Cons.nu**3/5184 - 155*Cons.nu**2/96 + Cons.nu*(34445/576 - 205*pi**2/96) - 675/64
    Cons.E_8=77*Cons.nu**4/31104 + 301*Cons.nu**3/1728 + Cons.nu**2*(-498449/3456 + 3157*pi**2/576) + Cons.nu*(-123671/5760 + 896*EulerGamma/15 + 9037*pi**2/1536 + 1792*log(2)/15) - 3969/128
    Cons.E_lnv_8=896*Cons.nu/15
    Cons.E_10=Cons.nu**5/512 + 55*Cons.nu**4/512 + Cons.nu**3*(69423/512 - 1353*pi**2/256) + Cons.nu**2*(-21337*pi**2/1024 - 896*log(2)/5 - 448*EulerGamma/5 + 893429/2880) + Cons.nu*(-228916843/115200 - 23672*log(2)/35 - 9976*EulerGamma/35 + 729*log(3)/7 + 126779*pi**2/512) - 45927/512
    Cons.E_lnv_10=-1312*Cons.nu**2/5 - 9976*Cons.nu/35
    Cons.E_11=27392*Cons.nu*pi/315
    Cons.E_12=2717*Cons.nu**6/6718464 + 5159*Cons.nu**5/248832 + Cons.nu**4*(-20543435/373248 + 272855*pi**2/124416) + Cons.nu**3*(-71700787/51840 + 1232*EulerGamma/27 + 2464*log(2)/27 + 6634243*pi**2/110592) + Cons.nu**2*(-86017789*pi**2/110592 - 2673*log(3)/14 + 112772*EulerGamma/105 + 18491*pi**4/2304 + 246004*log(2)/105 + 113176680983/14515200) + Cons.nu*(-389727504721/43545600 - 30809603*pi**4/786432 - 7128*log(3)/7 - 3934568*EulerGamma/8505 + 74888*log(2)/243 + 9118627045*pi**2/5308416) - 264627/1024
    Cons.E_lnv_12=48928*Cons.nu**3/135 + 79508*Cons.nu**2/105 - 3934568*Cons.nu/8505
    Cons.MDot_Alvi_8=Cons.M1**4*(3*Cons.chi1chi1 + 1)*(sqrt(1 - Cons.chi1chi1) + 1)/(2*Cons.M**4) + Cons.M2**4*(3*Cons.chi2chi2 + 1)*(sqrt(1 - Cons.chi2chi2) + 1)/(2*Cons.M**4)
    Cons.EvolveSpin1=np.linalg.norm(mul(Cons.S_chi1,conjugate(Cons.S_chi1)))>1e-8
    Cons.EvolveSpin2=np.linalg.norm(mul(Cons.S_chi2,conjugate(Cons.S_chi2)))>1e-8

@njit(cache=True)
def Recalculate_0(Cons,Vars,y):
    Vars.v = np.array([y[0]])
    Vars.rfrak_chi1 = np.array([y[1],y[2]])
    Vars.rfrak_chi2 = np.array([y[3],y[4]])
    Vars.R = y[5]*Cons.wHat + y[6]*Cons.xHat + y[7]*Cons.yHat + y[8]*Cons.zHat
    Vars.nHat = mul(mul(Vars.R,Cons.xHat),inverse(Vars.R))
    Vars.lambdaHat = mul(mul(Vars.R,Cons.yHat),inverse(Vars.R))
    Vars.ellHat = mul(mul(Vars.R,Cons.zHat),inverse(Vars.R))
    Vars.R_S1 = exp(Vars.rfrak_chi1[0]*Cons.xHat + Vars.rfrak_chi1[1]*Cons.yHat)
    Vars.R_S2 = exp(Vars.rfrak_chi2[0]*Cons.xHat + Vars.rfrak_chi2[1]*Cons.yHat)
    Vars.chiVec1 = mul(mul(mul(Cons.S_chi1,Vars.R_S1),Cons.zHat),mul(conjugate(Vars.R_S1),conjugate(Cons.S_chi1)))
    Vars.chiVec2 = mul(mul(mul(Cons.S_chi2,Vars.R_S2),Cons.zHat),mul(conjugate(Vars.R_S2),conjugate(Cons.S_chi2)))
    Vars.chi1_n = np.array([dot(Vars.chiVec1[1:],Vars.nHat[1:])])
    Vars.chi1_lambda = np.array([dot(Vars.chiVec1[1:],Vars.lambdaHat[1:])])
    Vars.chi1_ell = np.array([dot(Vars.chiVec1[1:],Vars.ellHat[1:])])
    Vars.chi2_n = np.array([dot(Vars.chiVec2[1:],Vars.nHat[1:])])
    Vars.chi2_lambda = np.array([dot(Vars.chiVec2[1:],Vars.lambdaHat[1:])])
    Vars.chi2_ell = np.array([dot(Vars.chiVec2[1:],Vars.ellHat[1:])])
    Vars.S = Cons.M1**2*Vars.chiVec1 + Cons.M2**2*Vars.chiVec2
    Vars.S_ell = Cons.M1**2*Vars.chi1_ell + Cons.M2**2*Vars.chi2_ell
    Vars.S_n = Cons.M1**2*Vars.chi1_n + Cons.M2**2*Vars.chi2_n
    Vars.S_lambda = Cons.M1**2*Vars.chi1_lambda + Cons.M2**2*Vars.chi2_lambda
    Vars.Sigma = Cons.M*(-Cons.M1*Vars.chiVec1 + Cons.M2*Vars.chiVec2)
    Vars.Sigma_ell = Cons.M*(-Cons.M1*Vars.chi1_ell + Cons.M2*Vars.chi2_ell)
    Vars.Sigma_n = Cons.M*(-Cons.M1*Vars.chi1_n + Cons.M2*Vars.chi2_n)
    Vars.Sigma_lambda = Cons.M*(-Cons.M1*Vars.chi1_lambda + Cons.M2*Vars.chi2_lambda)
    Vars.SS = np.array([dot(Vars.S[1:],Vars.S[1:])])
    Vars.SigmaSigma = np.array([dot(Vars.Sigma[1:],Vars.Sigma[1:])])
    Vars.SSigma = np.array([dot(Vars.S[1:],Vars.Sigma[1:])])
    Vars.S1 = Cons.M1**2*Vars.chiVec1
    Vars.S1_n = Cons.M1**2*Vars.chi1_n
    Vars.S1_lambda = Cons.M1**2*Vars.chi1_lambda
    Vars.S2 = Cons.M2**2*Vars.chiVec2
    Vars.S2_n = Cons.M2**2*Vars.chi2_n
    Vars.S2_lambda = Cons.M2**2*Vars.chi2_lambda
    Vars.chi_s_ell = Vars.chi1_ell/2 + Vars.chi2_ell/2
    Vars.chi_a_ell = Vars.chi1_ell/2 - Vars.chi2_ell/2
    Vars.Fcal_coeff = 32*Cons.nu**2*Vars.v**10/5

@njit
def OmegaVec_chiVec_1_0(Cons,Vars):
    Omega1_coeff = Vars.v**5/Cons.M
    return Omega1_coeff*Vars.ellHat*(-0.75*Cons.delta + 0.5*Cons.nu + 0.75)

@njit
def OmegaVec_chiVec_2_0(Cons,Vars):
    Omega2_coeff = Vars.v**5/Cons.M
    return Omega2_coeff*Vars.ellHat*(0.75*Cons.delta + 0.5*Cons.nu + 0.75)

@njit
def OmegaVec_0(Cons,Vars):
    gamma_PN_0 = 1.00000000000000
    a_ell_0 = 7.0*Vars.S_n + 3.0*Vars.Sigma_n*Cons.delta
    return Vars.ellHat*Vars.v**3/Cons.M + a_ell_0*gamma_PN_0*Vars.nHat*Vars.v**6/Cons.M**3


@njit(cache=True)
def TaylorT1_0(Cons,Vars):
    Flux = Cons.Fcal_0*Vars.Fcal_coeff
    dEdV = -Cons.E_0*Cons.M*Cons.nu*Vars.v
    Absorption = 0
    dvdt_T1 = (-Absorption - Flux)/dEdV
    dydt=np.zeros(9)
    [dydt[5],dydt[6],dydt[7],dydt[8]] = FrameFromAngularVelocityIntegrand(Vars.R, OmegaVec_0(Cons,Vars)[1:])
    dydt[0] = dvdt_T1[0] 
    if(Cons.EvolveSpin1):
        dydt[1], dydt[2]=FrameFromAngularVelocity_2D_Integrand(Vars.rfrak_chi1[0], Vars.rfrak_chi1[1],(mul(mul(inverse(Cons.S_chi1),OmegaVec_chiVec_1_0(Cons,Vars)),Cons.S_chi1))[1:])
    else:
        dydt[1], dydt[2] = 0.0, 0.0
    if(Cons.EvolveSpin2):
        dydt[3], dydt[4]=FrameFromAngularVelocity_2D_Integrand(Vars.rfrak_chi2[0], Vars.rfrak_chi2[1],(mul(mul(inverse(Cons.S_chi2),OmegaVec_chiVec_2_0(Cons,Vars)),Cons.S_chi2))[1:])
    else:
        dydt[3], dydt[4] = 0.0, 0.0
    return dydt

@njit(cache=True)
def TaylorT4_0(Cons,Vars):
    Flux = Cons.Fcal_0*Vars.Fcal_coeff
    dEdV = -Cons.E_0*Cons.M*Cons.nu*Vars.v
    Absorption = 0
    dvdt_T1 = (-Absorption - Flux)/dEdV
    dydt=np.zeros(9)
    [dydt[5],dydt[6],dydt[7],dydt[8]] = FrameFromAngularVelocityIntegrand(Vars.R, OmegaVec_0(Cons,Vars)[1:])
    dydt[0] = dvdt_T4[0]
    if(Cons.EvolveSpin1):
        dydt[1], dydt[2]=FrameFromAngularVelocity_2D_Integrand(Vars.rfrak_chi1[0], Vars.rfrak_chi1[1],(mul(mul(inverse(Cons.S_chi1),OmegaVec_chiVec_1_0(Cons,Vars)),Cons.S_chi1))[1:])
    else:
        dydt[1], dydt[2] = 0.0, 0.0
    if(Cons.EvolveSpin2):
        dydt[3], dydt[4]=FrameFromAngularVelocity_2D_Integrand(Vars.rfrak_chi2[0], Vars.rfrak_chi2[1],(mul(mul(inverse(Cons.S_chi2),OmegaVec_chiVec_2_0(Cons,Vars)),Cons.S_chi2))[1:])
    else:
        dydt[3], dydt[4] = 0.0, 0.0
    return dydt      

@njit(cache=True)
def TaylorT5_0(Cons,Vars):
    Flux = Cons.Fcal_0*Vars.Fcal_coeff
    dEdV = -Cons.E_0*Cons.M*Cons.nu*Vars.v
    Absorption = 0
    dvdt_T1 = (-Absorption - Flux)/dEdV
    dydt=np.zeros(8)
    [dydt[5],dydt[6],dydt[7],dydt[8]] = FrameFromAngularVelocityIntegrand(Vars.R, OmegaVec_0(Cons,Vars)[1:])
    dydt[0] = dvdt_T5[0]
    if(Cons.EvolveSpin1):
        dydt[1], dydt[2]=FrameFromAngularVelocity_2D_Integrand(Vars.rfrak_chi1[0], Vars.rfrak_chi1[1],(mul(mul(inverse(Cons.S_chi1),OmegaVec_chiVec_1_0(Cons,Vars)),Cons.S_chi1))[1:])
    else:
        dydt[1], dydt[2] = 0.0, 0.0
    if(Cons.EvolveSpin2):
        dydt[3], dydt[4]=FrameFromAngularVelocity_2D_Integrand(Vars.rfrak_chi2[0], Vars.rfrak_chi2[1],(mul(mul(inverse(Cons.S_chi2),OmegaVec_chiVec_2_0(Cons,Vars)),Cons.S_chi2))[1:])
    else:
        dydt[3], dydt[4] = 0.0, 0.0
    return dydt      

@njit(cache=True)
def Recalculate_0p50(Cons,Vars,y):
    Vars.v = np.array([y[0]])
    Vars.rfrak_chi1 = np.array([y[1],y[2]])
    Vars.rfrak_chi2 = np.array([y[3],y[4]])
    Vars.R = y[5]*Cons.wHat + y[6]*Cons.xHat + y[7]*Cons.yHat + y[8]*Cons.zHat
    Vars.nHat = mul(mul(Vars.R,Cons.xHat),inverse(Vars.R))
    Vars.lambdaHat = mul(mul(Vars.R,Cons.yHat),inverse(Vars.R))
    Vars.ellHat = mul(mul(Vars.R,Cons.zHat),inverse(Vars.R))
    Vars.R_S1 = exp(Vars.rfrak_chi1[0]*Cons.xHat + Vars.rfrak_chi1[1]*Cons.yHat)
    Vars.R_S2 = exp(Vars.rfrak_chi2[0]*Cons.xHat + Vars.rfrak_chi2[1]*Cons.yHat)
    Vars.chiVec1 = mul(mul(mul(Cons.S_chi1,Vars.R_S1),Cons.zHat),mul(conjugate(Vars.R_S1),conjugate(Cons.S_chi1)))
    Vars.chiVec2 = mul(mul(mul(Cons.S_chi2,Vars.R_S2),Cons.zHat),mul(conjugate(Vars.R_S2),conjugate(Cons.S_chi2)))
    Vars.chi1_n = np.array([dot(Vars.chiVec1[1:],Vars.nHat[1:])])
    Vars.chi1_lambda = np.array([dot(Vars.chiVec1[1:],Vars.lambdaHat[1:])])
    Vars.chi1_ell = np.array([dot(Vars.chiVec1[1:],Vars.ellHat[1:])])
    Vars.chi2_n = np.array([dot(Vars.chiVec2[1:],Vars.nHat[1:])])
    Vars.chi2_lambda = np.array([dot(Vars.chiVec2[1:],Vars.lambdaHat[1:])])
    Vars.chi2_ell = np.array([dot(Vars.chiVec2[1:],Vars.ellHat[1:])])
    Vars.S = Cons.M1**2*Vars.chiVec1 + Cons.M2**2*Vars.chiVec2
    Vars.S_ell = Cons.M1**2*Vars.chi1_ell + Cons.M2**2*Vars.chi2_ell
    Vars.S_n = Cons.M1**2*Vars.chi1_n + Cons.M2**2*Vars.chi2_n
    Vars.S_lambda = Cons.M1**2*Vars.chi1_lambda + Cons.M2**2*Vars.chi2_lambda
    Vars.Sigma = Cons.M*(-Cons.M1*Vars.chiVec1 + Cons.M2*Vars.chiVec2)
    Vars.Sigma_ell = Cons.M*(-Cons.M1*Vars.chi1_ell + Cons.M2*Vars.chi2_ell)
    Vars.Sigma_n = Cons.M*(-Cons.M1*Vars.chi1_n + Cons.M2*Vars.chi2_n)
    Vars.Sigma_lambda = Cons.M*(-Cons.M1*Vars.chi1_lambda + Cons.M2*Vars.chi2_lambda)
    Vars.SS = np.array([dot(Vars.S[1:],Vars.S[1:])])
    Vars.SigmaSigma = np.array([dot(Vars.Sigma[1:],Vars.Sigma[1:])])
    Vars.SSigma = np.array([dot(Vars.S[1:],Vars.Sigma[1:])])
    Vars.S1 = Cons.M1**2*Vars.chiVec1
    Vars.S1_n = Cons.M1**2*Vars.chi1_n
    Vars.S1_lambda = Cons.M1**2*Vars.chi1_lambda
    Vars.S2 = Cons.M2**2*Vars.chiVec2
    Vars.S2_n = Cons.M2**2*Vars.chi2_n
    Vars.S2_lambda = Cons.M2**2*Vars.chi2_lambda
    Vars.chi_s_ell = Vars.chi1_ell/2 + Vars.chi2_ell/2
    Vars.chi_a_ell = Vars.chi1_ell/2 - Vars.chi2_ell/2
    Vars.Fcal_coeff = 32*Cons.nu**2*Vars.v**10/5

@njit
def OmegaVec_chiVec_1_0p50(Cons,Vars):
    Omega1_coeff = Vars.v**5/Cons.M
    return Omega1_coeff*(Vars.ellHat*(-0.75*Cons.delta + 0.5*Cons.nu + 0.75) + Vars.nHat*Vars.v*(3.0*Vars.chi1_n*Cons.nu + 3.0*Cons.M2**2*Vars.chi2_n/Cons.M**2) - Cons.M2**2*Vars.chiVec2*Vars.v/Cons.M**2)

@njit
def OmegaVec_chiVec_2_0p50(Cons,Vars):
    Omega2_coeff = Vars.v**5/Cons.M
    return Omega2_coeff*(Vars.ellHat*(0.75*Cons.delta + 0.5*Cons.nu + 0.75) + Vars.nHat*Vars.v*(3.0*Vars.chi2_n*Cons.nu + 3.0*Cons.M1**2*Vars.chi1_n/Cons.M**2) - Cons.M1**2*Vars.chiVec1*Vars.v/Cons.M**2)

@njit
def OmegaVec_0p50(Cons,Vars):
    gamma_PN_0 = 1.00000000000000
    a_ell_1 = 3.0*Vars.Sigma_ell*Vars.Sigma_n*Cons.nu*(Cons.kappa_p + 2.0)/Cons.M**2 - 3.0*Cons.delta*(Vars.S_ell*Vars.Sigma_n + Vars.S_n*Vars.Sigma_ell)/Cons.M**2 - 0.5*(3.0*Vars.S_n*(2.0*Vars.S_ell*(Cons.kappa_p + 2.0) + Vars.Sigma_ell*(Cons.delta*Cons.kappa_p - Cons.kappa_m)) + 3.0*Vars.Sigma_n*(Vars.S_ell*(Cons.delta*Cons.kappa_p - Cons.kappa_m) + Vars.Sigma_ell*(-Cons.delta*Cons.kappa_m + Cons.kappa_p)))/Cons.M**2
    a_ell_0 = 7.0*Vars.S_n + 3.0*Vars.Sigma_n*Cons.delta
    return Vars.ellHat*Vars.v**3/Cons.M + gamma_PN_0*Vars.nHat*Vars.v**6*(a_ell_0 + a_ell_1*Vars.v)/Cons.M**3


@njit(cache=True)
def TaylorT1_0p50(Cons,Vars):
    Flux = Cons.Fcal_0*Vars.Fcal_coeff
    dEdV = -Cons.E_0*Cons.M*Cons.nu*Vars.v
    Absorption = 0
    dvdt_T1 = (-Absorption - Flux)/dEdV
    dydt=np.zeros(9)
    [dydt[5],dydt[6],dydt[7],dydt[8]] = FrameFromAngularVelocityIntegrand(Vars.R, OmegaVec_0p50(Cons,Vars)[1:])
    dydt[0] = dvdt_T1[0] 
    if(Cons.EvolveSpin1):
        dydt[1], dydt[2]=FrameFromAngularVelocity_2D_Integrand(Vars.rfrak_chi1[0], Vars.rfrak_chi1[1],(mul(mul(inverse(Cons.S_chi1),OmegaVec_chiVec_1_0p50(Cons,Vars)),Cons.S_chi1))[1:])
    else:
        dydt[1], dydt[2] = 0.0, 0.0
    if(Cons.EvolveSpin2):
        dydt[3], dydt[4]=FrameFromAngularVelocity_2D_Integrand(Vars.rfrak_chi2[0], Vars.rfrak_chi2[1],(mul(mul(inverse(Cons.S_chi2),OmegaVec_chiVec_2_0p50(Cons,Vars)),Cons.S_chi2))[1:])
    else:
        dydt[3], dydt[4] = 0.0, 0.0
    return dydt

@njit(cache=True)
def TaylorT4_0p50(Cons,Vars):
    Flux = Cons.Fcal_0*Vars.Fcal_coeff
    dEdV = -Cons.E_0*Cons.M*Cons.nu*Vars.v
    Absorption = 0
    dvdt_T1 = (-Absorption - Flux)/dEdV
    dydt=np.zeros(9)
    [dydt[5],dydt[6],dydt[7],dydt[8]] = FrameFromAngularVelocityIntegrand(Vars.R, OmegaVec_0p50(Cons,Vars)[1:])
    dydt[0] = dvdt_T4[0]
    if(Cons.EvolveSpin1):
        dydt[1], dydt[2]=FrameFromAngularVelocity_2D_Integrand(Vars.rfrak_chi1[0], Vars.rfrak_chi1[1],(mul(mul(inverse(Cons.S_chi1),OmegaVec_chiVec_1_0p50(Cons,Vars)),Cons.S_chi1))[1:])
    else:
        dydt[1], dydt[2] = 0.0, 0.0
    if(Cons.EvolveSpin2):
        dydt[3], dydt[4]=FrameFromAngularVelocity_2D_Integrand(Vars.rfrak_chi2[0], Vars.rfrak_chi2[1],(mul(mul(inverse(Cons.S_chi2),OmegaVec_chiVec_2_0p50(Cons,Vars)),Cons.S_chi2))[1:])
    else:
        dydt[3], dydt[4] = 0.0, 0.0
    return dydt      

@njit(cache=True)
def TaylorT5_0p50(Cons,Vars):
    Flux = Cons.Fcal_0*Vars.Fcal_coeff
    dEdV = -Cons.E_0*Cons.M*Cons.nu*Vars.v
    Absorption = 0
    dvdt_T1 = (-Absorption - Flux)/dEdV
    dydt=np.zeros(8)
    [dydt[5],dydt[6],dydt[7],dydt[8]] = FrameFromAngularVelocityIntegrand(Vars.R, OmegaVec_0p50(Cons,Vars)[1:])
    dydt[0] = dvdt_T5[0]
    if(Cons.EvolveSpin1):
        dydt[1], dydt[2]=FrameFromAngularVelocity_2D_Integrand(Vars.rfrak_chi1[0], Vars.rfrak_chi1[1],(mul(mul(inverse(Cons.S_chi1),OmegaVec_chiVec_1_0p50(Cons,Vars)),Cons.S_chi1))[1:])
    else:
        dydt[1], dydt[2] = 0.0, 0.0
    if(Cons.EvolveSpin2):
        dydt[3], dydt[4]=FrameFromAngularVelocity_2D_Integrand(Vars.rfrak_chi2[0], Vars.rfrak_chi2[1],(mul(mul(inverse(Cons.S_chi2),OmegaVec_chiVec_2_0p50(Cons,Vars)),Cons.S_chi2))[1:])
    else:
        dydt[3], dydt[4] = 0.0, 0.0
    return dydt      

@njit(cache=True)
def Recalculate_1p0(Cons,Vars,y):
    Vars.v = np.array([y[0]])
    Vars.rfrak_chi1 = np.array([y[1],y[2]])
    Vars.rfrak_chi2 = np.array([y[3],y[4]])
    Vars.R = y[5]*Cons.wHat + y[6]*Cons.xHat + y[7]*Cons.yHat + y[8]*Cons.zHat
    Vars.nHat = mul(mul(Vars.R,Cons.xHat),inverse(Vars.R))
    Vars.lambdaHat = mul(mul(Vars.R,Cons.yHat),inverse(Vars.R))
    Vars.ellHat = mul(mul(Vars.R,Cons.zHat),inverse(Vars.R))
    Vars.R_S1 = exp(Vars.rfrak_chi1[0]*Cons.xHat + Vars.rfrak_chi1[1]*Cons.yHat)
    Vars.R_S2 = exp(Vars.rfrak_chi2[0]*Cons.xHat + Vars.rfrak_chi2[1]*Cons.yHat)
    Vars.chiVec1 = mul(mul(mul(Cons.S_chi1,Vars.R_S1),Cons.zHat),mul(conjugate(Vars.R_S1),conjugate(Cons.S_chi1)))
    Vars.chiVec2 = mul(mul(mul(Cons.S_chi2,Vars.R_S2),Cons.zHat),mul(conjugate(Vars.R_S2),conjugate(Cons.S_chi2)))
    Vars.chi1_n = np.array([dot(Vars.chiVec1[1:],Vars.nHat[1:])])
    Vars.chi1_lambda = np.array([dot(Vars.chiVec1[1:],Vars.lambdaHat[1:])])
    Vars.chi1_ell = np.array([dot(Vars.chiVec1[1:],Vars.ellHat[1:])])
    Vars.chi2_n = np.array([dot(Vars.chiVec2[1:],Vars.nHat[1:])])
    Vars.chi2_lambda = np.array([dot(Vars.chiVec2[1:],Vars.lambdaHat[1:])])
    Vars.chi2_ell = np.array([dot(Vars.chiVec2[1:],Vars.ellHat[1:])])
    Vars.S = Cons.M1**2*Vars.chiVec1 + Cons.M2**2*Vars.chiVec2
    Vars.S_ell = Cons.M1**2*Vars.chi1_ell + Cons.M2**2*Vars.chi2_ell
    Vars.S_n = Cons.M1**2*Vars.chi1_n + Cons.M2**2*Vars.chi2_n
    Vars.S_lambda = Cons.M1**2*Vars.chi1_lambda + Cons.M2**2*Vars.chi2_lambda
    Vars.Sigma = Cons.M*(-Cons.M1*Vars.chiVec1 + Cons.M2*Vars.chiVec2)
    Vars.Sigma_ell = Cons.M*(-Cons.M1*Vars.chi1_ell + Cons.M2*Vars.chi2_ell)
    Vars.Sigma_n = Cons.M*(-Cons.M1*Vars.chi1_n + Cons.M2*Vars.chi2_n)
    Vars.Sigma_lambda = Cons.M*(-Cons.M1*Vars.chi1_lambda + Cons.M2*Vars.chi2_lambda)
    Vars.SS = np.array([dot(Vars.S[1:],Vars.S[1:])])
    Vars.SigmaSigma = np.array([dot(Vars.Sigma[1:],Vars.Sigma[1:])])
    Vars.SSigma = np.array([dot(Vars.S[1:],Vars.Sigma[1:])])
    Vars.S1 = Cons.M1**2*Vars.chiVec1
    Vars.S1_n = Cons.M1**2*Vars.chi1_n
    Vars.S1_lambda = Cons.M1**2*Vars.chi1_lambda
    Vars.S2 = Cons.M2**2*Vars.chiVec2
    Vars.S2_n = Cons.M2**2*Vars.chi2_n
    Vars.S2_lambda = Cons.M2**2*Vars.chi2_lambda
    Vars.chi_s_ell = Vars.chi1_ell/2 + Vars.chi2_ell/2
    Vars.chi_a_ell = Vars.chi1_ell/2 - Vars.chi2_ell/2
    Vars.Fcal_coeff = 32*Cons.nu**2*Vars.v**10/5

@njit
def OmegaVec_chiVec_1_1p0(Cons,Vars):
    Omega1_coeff = Vars.v**5/Cons.M
    return Omega1_coeff*(Vars.ellHat*(-0.75*Cons.delta + 0.5*Cons.nu + Vars.v**2*(Cons.delta*(0.625*Cons.nu - 0.5625) + Cons.nu*(1.25 - 0.0416666666666667*Cons.nu) + 0.5625) + 0.75) + Vars.nHat*Vars.v*(3.0*Vars.chi1_n*Cons.nu + 3.0*Cons.M2**2*Vars.chi2_n/Cons.M**2) - Cons.M2**2*Vars.chiVec2*Vars.v/Cons.M**2)

@njit
def OmegaVec_chiVec_2_1p0(Cons,Vars):
    Omega2_coeff = Vars.v**5/Cons.M
    return Omega2_coeff*(Vars.ellHat*(0.75*Cons.delta + 0.5*Cons.nu + Vars.v**2*(Cons.delta*(0.5625 - 0.625*Cons.nu) + Cons.nu*(1.25 - 0.0416666666666667*Cons.nu) + 0.5625) + 0.75) + Vars.nHat*Vars.v*(3.0*Vars.chi2_n*Cons.nu + 3.0*Cons.M1**2*Vars.chi1_n/Cons.M**2) - Cons.M1**2*Vars.chiVec1*Vars.v/Cons.M**2)

@njit
def OmegaVec_1p0(Cons,Vars):
    gamma_PN_0 = 1.00000000000000
    a_ell_1 = 3.0*Vars.Sigma_ell*Vars.Sigma_n*Cons.nu*(Cons.kappa_p + 2.0)/Cons.M**2 - 3.0*Cons.delta*(Vars.S_ell*Vars.Sigma_n + Vars.S_n*Vars.Sigma_ell)/Cons.M**2 - 0.5*(3.0*Vars.S_n*(2.0*Vars.S_ell*(Cons.kappa_p + 2.0) + Vars.Sigma_ell*(Cons.delta*Cons.kappa_p - Cons.kappa_m)) + 3.0*Vars.Sigma_n*(Vars.S_ell*(Cons.delta*Cons.kappa_p - Cons.kappa_m) + Vars.Sigma_ell*(-Cons.delta*Cons.kappa_m + Cons.kappa_p)))/Cons.M**2
    a_ell_2 = Vars.S_n*(-9.66666666666667*Cons.nu - 10.0) + Vars.Sigma_n*Cons.delta*(-4.5*Cons.nu - 6.0)
    gamma_PN_2 = 1.0 - 0.333333333333333*Cons.nu
    a_ell_0 = 7.0*Vars.S_n + 3.0*Vars.Sigma_n*Cons.delta
    return Vars.ellHat*Vars.v**3/Cons.M + Vars.nHat*Vars.v**6*(a_ell_0 + Vars.v*(a_ell_1 + a_ell_2*Vars.v))*(gamma_PN_0 + gamma_PN_2*Vars.v**2)/Cons.M**3


@njit(cache=True)
def TaylorT1_1p0(Cons,Vars):
    Flux = Vars.Fcal_coeff*(Cons.Fcal_0 + Cons.Fcal_2*Vars.v**2)
    dEdV = -Cons.M*Cons.nu*Vars.v*(Cons.E_0 + 2.0*Cons.E_2*Vars.v**2)
    Absorption = 0
    dvdt_T1 = (-Absorption - Flux)/dEdV
    dydt=np.zeros(9)
    [dydt[5],dydt[6],dydt[7],dydt[8]] = FrameFromAngularVelocityIntegrand(Vars.R, OmegaVec_1p0(Cons,Vars)[1:])
    dydt[0] = dvdt_T1[0] 
    if(Cons.EvolveSpin1):
        dydt[1], dydt[2]=FrameFromAngularVelocity_2D_Integrand(Vars.rfrak_chi1[0], Vars.rfrak_chi1[1],(mul(mul(inverse(Cons.S_chi1),OmegaVec_chiVec_1_1p0(Cons,Vars)),Cons.S_chi1))[1:])
    else:
        dydt[1], dydt[2] = 0.0, 0.0
    if(Cons.EvolveSpin2):
        dydt[3], dydt[4]=FrameFromAngularVelocity_2D_Integrand(Vars.rfrak_chi2[0], Vars.rfrak_chi2[1],(mul(mul(inverse(Cons.S_chi2),OmegaVec_chiVec_2_1p0(Cons,Vars)),Cons.S_chi2))[1:])
    else:
        dydt[3], dydt[4] = 0.0, 0.0
    return dydt

@njit(cache=True)
def TaylorT4_1p0(Cons,Vars):
    Flux = Vars.Fcal_coeff*(Cons.Fcal_0 + Cons.Fcal_2*Vars.v**2)
    dEdV = -Cons.M*Cons.nu*Vars.v*(Cons.E_0 + 2.0*Cons.E_2*Vars.v**2)
    Absorption = 0
    dvdt_T1 = (-Absorption - Flux)/dEdV
    dydt=np.zeros(9)
    [dydt[5],dydt[6],dydt[7],dydt[8]] = FrameFromAngularVelocityIntegrand(Vars.R, OmegaVec_1p0(Cons,Vars)[1:])
    dydt[0] = dvdt_T4[0]
    if(Cons.EvolveSpin1):
        dydt[1], dydt[2]=FrameFromAngularVelocity_2D_Integrand(Vars.rfrak_chi1[0], Vars.rfrak_chi1[1],(mul(mul(inverse(Cons.S_chi1),OmegaVec_chiVec_1_1p0(Cons,Vars)),Cons.S_chi1))[1:])
    else:
        dydt[1], dydt[2] = 0.0, 0.0
    if(Cons.EvolveSpin2):
        dydt[3], dydt[4]=FrameFromAngularVelocity_2D_Integrand(Vars.rfrak_chi2[0], Vars.rfrak_chi2[1],(mul(mul(inverse(Cons.S_chi2),OmegaVec_chiVec_2_1p0(Cons,Vars)),Cons.S_chi2))[1:])
    else:
        dydt[3], dydt[4] = 0.0, 0.0
    return dydt      

@njit(cache=True)
def TaylorT5_1p0(Cons,Vars):
    Flux = Vars.Fcal_coeff*(Cons.Fcal_0 + Cons.Fcal_2*Vars.v**2)
    dEdV = -Cons.M*Cons.nu*Vars.v*(Cons.E_0 + 2.0*Cons.E_2*Vars.v**2)
    Absorption = 0
    dvdt_T1 = (-Absorption - Flux)/dEdV
    dydt=np.zeros(8)
    [dydt[5],dydt[6],dydt[7],dydt[8]] = FrameFromAngularVelocityIntegrand(Vars.R, OmegaVec_1p0(Cons,Vars)[1:])
    dydt[0] = dvdt_T5[0]
    if(Cons.EvolveSpin1):
        dydt[1], dydt[2]=FrameFromAngularVelocity_2D_Integrand(Vars.rfrak_chi1[0], Vars.rfrak_chi1[1],(mul(mul(inverse(Cons.S_chi1),OmegaVec_chiVec_1_1p0(Cons,Vars)),Cons.S_chi1))[1:])
    else:
        dydt[1], dydt[2] = 0.0, 0.0
    if(Cons.EvolveSpin2):
        dydt[3], dydt[4]=FrameFromAngularVelocity_2D_Integrand(Vars.rfrak_chi2[0], Vars.rfrak_chi2[1],(mul(mul(inverse(Cons.S_chi2),OmegaVec_chiVec_2_1p0(Cons,Vars)),Cons.S_chi2))[1:])
    else:
        dydt[3], dydt[4] = 0.0, 0.0
    return dydt      

@njit(cache=True)
def Recalculate_1p5(Cons,Vars,y):
    Vars.v = np.array([y[0]])
    Vars.rfrak_chi1 = np.array([y[1],y[2]])
    Vars.rfrak_chi2 = np.array([y[3],y[4]])
    Vars.R = y[5]*Cons.wHat + y[6]*Cons.xHat + y[7]*Cons.yHat + y[8]*Cons.zHat
    Vars.nHat = mul(mul(Vars.R,Cons.xHat),inverse(Vars.R))
    Vars.lambdaHat = mul(mul(Vars.R,Cons.yHat),inverse(Vars.R))
    Vars.ellHat = mul(mul(Vars.R,Cons.zHat),inverse(Vars.R))
    Vars.R_S1 = exp(Vars.rfrak_chi1[0]*Cons.xHat + Vars.rfrak_chi1[1]*Cons.yHat)
    Vars.R_S2 = exp(Vars.rfrak_chi2[0]*Cons.xHat + Vars.rfrak_chi2[1]*Cons.yHat)
    Vars.chiVec1 = mul(mul(mul(Cons.S_chi1,Vars.R_S1),Cons.zHat),mul(conjugate(Vars.R_S1),conjugate(Cons.S_chi1)))
    Vars.chiVec2 = mul(mul(mul(Cons.S_chi2,Vars.R_S2),Cons.zHat),mul(conjugate(Vars.R_S2),conjugate(Cons.S_chi2)))
    Vars.chi1_n = np.array([dot(Vars.chiVec1[1:],Vars.nHat[1:])])
    Vars.chi1_lambda = np.array([dot(Vars.chiVec1[1:],Vars.lambdaHat[1:])])
    Vars.chi1_ell = np.array([dot(Vars.chiVec1[1:],Vars.ellHat[1:])])
    Vars.chi2_n = np.array([dot(Vars.chiVec2[1:],Vars.nHat[1:])])
    Vars.chi2_lambda = np.array([dot(Vars.chiVec2[1:],Vars.lambdaHat[1:])])
    Vars.chi2_ell = np.array([dot(Vars.chiVec2[1:],Vars.ellHat[1:])])
    Vars.S = Cons.M1**2*Vars.chiVec1 + Cons.M2**2*Vars.chiVec2
    Vars.S_ell = Cons.M1**2*Vars.chi1_ell + Cons.M2**2*Vars.chi2_ell
    Vars.S_n = Cons.M1**2*Vars.chi1_n + Cons.M2**2*Vars.chi2_n
    Vars.S_lambda = Cons.M1**2*Vars.chi1_lambda + Cons.M2**2*Vars.chi2_lambda
    Vars.Sigma = Cons.M*(-Cons.M1*Vars.chiVec1 + Cons.M2*Vars.chiVec2)
    Vars.Sigma_ell = Cons.M*(-Cons.M1*Vars.chi1_ell + Cons.M2*Vars.chi2_ell)
    Vars.Sigma_n = Cons.M*(-Cons.M1*Vars.chi1_n + Cons.M2*Vars.chi2_n)
    Vars.Sigma_lambda = Cons.M*(-Cons.M1*Vars.chi1_lambda + Cons.M2*Vars.chi2_lambda)
    Vars.SS = np.array([dot(Vars.S[1:],Vars.S[1:])])
    Vars.SigmaSigma = np.array([dot(Vars.Sigma[1:],Vars.Sigma[1:])])
    Vars.SSigma = np.array([dot(Vars.S[1:],Vars.Sigma[1:])])
    Vars.S1 = Cons.M1**2*Vars.chiVec1
    Vars.S1_n = Cons.M1**2*Vars.chi1_n
    Vars.S1_lambda = Cons.M1**2*Vars.chi1_lambda
    Vars.S2 = Cons.M2**2*Vars.chiVec2
    Vars.S2_n = Cons.M2**2*Vars.chi2_n
    Vars.S2_lambda = Cons.M2**2*Vars.chi2_lambda
    Vars.chi_s_ell = Vars.chi1_ell/2 + Vars.chi2_ell/2
    Vars.chi_a_ell = Vars.chi1_ell/2 - Vars.chi2_ell/2
    Vars.Fcal_coeff = 32*Cons.nu**2*Vars.v**10/5
    Vars.Fcal_SO_3 = (-4*Vars.S_ell - 5*Vars.Sigma_ell*Cons.delta/4)/Cons.M**2
    Vars.E_SO_3 = (14*Vars.S_ell/3 + 2*Vars.Sigma_ell*Cons.delta)/Cons.M**2

@njit
def OmegaVec_chiVec_1_1p5(Cons,Vars):
    Omega1_coeff = Vars.v**5/Cons.M
    return Omega1_coeff*(Vars.ellHat*(-0.75*Cons.delta + 0.5*Cons.nu + Vars.v**2*(Cons.delta*(0.625*Cons.nu - 0.5625) + Cons.nu*(1.25 - 0.0416666666666667*Cons.nu) + 0.5625 + Vars.v*(Vars.S_ell*(0.833333333333333*Cons.nu + 1.25) + (-1.25*Vars.S_ell*Cons.delta + Vars.Sigma_ell*Cons.delta*(0.5*Cons.nu + 0.75) - 0.75*Vars.Sigma_ell*Cons.delta**2/Cons.M)/Cons.M)/Cons.M**2) + 0.75) + Vars.nHat*Vars.v*(3.0*Vars.chi1_n*Cons.nu + 3.0*Cons.M2**2*Vars.chi2_n/Cons.M**2 + Vars.v**2*(Vars.S1_n*(Cons.delta*(-2.25*Cons.kappa_1 - 0.25) + Cons.kappa_1*(1.5*Cons.nu + 0.75) + 0.5*Cons.nu + 3.25) + Vars.S1_n*(Cons.delta*(1.5 - 0.75*Cons.kappa_1) + 0.75*Cons.kappa_1 - 1.5)/Cons.nu + Vars.S2_n*(1.5*Cons.delta - 1.0*Cons.nu + 3.0))/Cons.M**2) - Cons.M2**2*Vars.chiVec2*Vars.v/Cons.M**2 + Vars.S2*Vars.v**3*(2.0 - 0.5*Cons.delta)/Cons.M**2 + Vars.S2_lambda*Vars.lambdaHat*Vars.v**3*(0.5*Cons.delta + 2.0)/Cons.M**2)

@njit
def OmegaVec_chiVec_2_1p5(Cons,Vars):
    Omega2_coeff = Vars.v**5/Cons.M
    return Omega2_coeff*(Vars.ellHat*(0.75*Cons.delta + 0.5*Cons.nu + Vars.v**2*(Cons.delta*(0.5625 - 0.625*Cons.nu) + Cons.nu*(1.25 - 0.0416666666666667*Cons.nu) + 0.5625 + Vars.v*(Vars.S_ell*(0.833333333333333*Cons.nu + 1.25) + (1.25*Vars.S_ell*Cons.delta + Vars.Sigma_ell*Cons.delta*(-0.5*Cons.nu - 0.75) - 0.75*Vars.Sigma_ell*Cons.delta**2/Cons.M)/Cons.M)/Cons.M**2) + 0.75) + Vars.nHat*Vars.v*(3.0*Vars.chi2_n*Cons.nu + 3.0*Cons.M1**2*Vars.chi1_n/Cons.M**2 + Vars.v**2*(Vars.S1_n*(-1.5*Cons.delta - 1.0*Cons.nu + 3.0) + Vars.S2_n*(Cons.delta*(2.25*Cons.kappa_2 + 0.25) + Cons.kappa_2*(1.5*Cons.nu + 0.75) + 0.5*Cons.nu + 3.25) + Vars.S2_n*(Cons.delta*(0.75*Cons.kappa_2 - 1.5) + 0.75*Cons.kappa_2 - 1.5)/Cons.nu)/Cons.M**2) - Cons.M1**2*Vars.chiVec1*Vars.v/Cons.M**2 + Vars.S1*Vars.v**3*(0.5*Cons.delta + 2.0)/Cons.M**2 + Vars.S1_lambda*Vars.lambdaHat*Vars.v**3*(2.0 - 0.5*Cons.delta)/Cons.M**2)

@njit
def OmegaVec_1p5(Cons,Vars):
    gamma_PN_0 = 1.00000000000000
    a_ell_1 = 3.0*Vars.Sigma_ell*Vars.Sigma_n*Cons.nu*(Cons.kappa_p + 2.0)/Cons.M**2 - 3.0*Cons.delta*(Vars.S_ell*Vars.Sigma_n + Vars.S_n*Vars.Sigma_ell)/Cons.M**2 - 0.5*(3.0*Vars.S_n*(2.0*Vars.S_ell*(Cons.kappa_p + 2.0) + Vars.Sigma_ell*(Cons.delta*Cons.kappa_p - Cons.kappa_m)) + 3.0*Vars.Sigma_n*(Vars.S_ell*(Cons.delta*Cons.kappa_p - Cons.kappa_m) + Vars.Sigma_ell*(-Cons.delta*Cons.kappa_m + Cons.kappa_p)))/Cons.M**2
    a_ell_2 = Vars.S_n*(-9.66666666666667*Cons.nu - 10.0) + Vars.Sigma_n*Cons.delta*(-4.5*Cons.nu - 6.0)
    gamma_PN_2 = 1.0 - 0.333333333333333*Cons.nu
    a_ell_3 = Cons.nu*(4.0*Cons.delta*(Vars.S_ell*Vars.Sigma_n + Vars.S_n*Vars.Sigma_ell)/Cons.M**2 + 0.5*(2.0*Vars.S_n*(4.0*Vars.S_ell*(Cons.kappa_p + 2.0) + Vars.Sigma_ell*(2.0*Cons.delta*Cons.kappa_p + 7.0*Cons.kappa_m)) + Vars.Sigma_n*(2.0*Vars.S_ell*(2.0*Cons.delta*Cons.kappa_p + 7.0*Cons.kappa_m) + Vars.Sigma_ell*(5.0*Cons.delta*Cons.kappa_m - 17.0*Cons.kappa_p - 26.0)))/Cons.M**2) - 4.0*Vars.Sigma_ell*Vars.Sigma_n*Cons.nu**2*(Cons.kappa_p + 2.0)/Cons.M**2 + 6.0*Cons.delta*(Vars.S_ell*Vars.Sigma_n + Vars.S_n*Vars.Sigma_ell)/Cons.M**2 + 0.5*(Vars.S_n*(Vars.S_ell*(-9.0*Cons.delta*Cons.kappa_m + 3.0*Cons.kappa_p + 22.0) + Vars.Sigma_ell*(6.0*Cons.delta*Cons.kappa_p - 6.0*Cons.kappa_m)) + 6.0*Vars.Sigma_n*(Vars.S_ell*(Cons.delta*Cons.kappa_p - Cons.kappa_m) + Vars.Sigma_ell*(-Cons.delta*Cons.kappa_m + Cons.kappa_p + 1.0)))/Cons.M**2
    gamma_PN_3 = (1.66666666666667*Vars.S_ell + Vars.Sigma_ell*Cons.delta)/Cons.M**2
    a_ell_0 = 7.0*Vars.S_n + 3.0*Vars.Sigma_n*Cons.delta
    return Vars.ellHat*Vars.v**3/Cons.M + Vars.nHat*Vars.v**6*(a_ell_0 + Vars.v*(a_ell_1 + Vars.v*(a_ell_2 + a_ell_3*Vars.v)))*(gamma_PN_0 + Vars.v**2*(gamma_PN_2 + gamma_PN_3*Vars.v))/Cons.M**3


@njit(cache=True)
def TaylorT1_1p5(Cons,Vars):
    Flux = Vars.Fcal_coeff*(Cons.Fcal_0 + Vars.v**2*(Cons.Fcal_2 + Vars.v*(Cons.Fcal_3 + Vars.Fcal_SO_3)))
    dEdV = -0.5*Cons.M*Cons.nu*Vars.v*(2.0*Cons.E_0 + Vars.v**2*(4.0*Cons.E_2 + 5.0*Vars.E_SO_3*Vars.v))
    Absorption = 0
    dvdt_T1 = (-Absorption - Flux)/dEdV
    dydt=np.zeros(9)
    [dydt[5],dydt[6],dydt[7],dydt[8]] = FrameFromAngularVelocityIntegrand(Vars.R, OmegaVec_1p5(Cons,Vars)[1:])
    dydt[0] = dvdt_T1[0] 
    if(Cons.EvolveSpin1):
        dydt[1], dydt[2]=FrameFromAngularVelocity_2D_Integrand(Vars.rfrak_chi1[0], Vars.rfrak_chi1[1],(mul(mul(inverse(Cons.S_chi1),OmegaVec_chiVec_1_1p5(Cons,Vars)),Cons.S_chi1))[1:])
    else:
        dydt[1], dydt[2] = 0.0, 0.0
    if(Cons.EvolveSpin2):
        dydt[3], dydt[4]=FrameFromAngularVelocity_2D_Integrand(Vars.rfrak_chi2[0], Vars.rfrak_chi2[1],(mul(mul(inverse(Cons.S_chi2),OmegaVec_chiVec_2_1p5(Cons,Vars)),Cons.S_chi2))[1:])
    else:
        dydt[3], dydt[4] = 0.0, 0.0
    return dydt

@njit(cache=True)
def TaylorT4_1p5(Cons,Vars):
    Flux = Vars.Fcal_coeff*(Cons.Fcal_0 + Vars.v**2*(Cons.Fcal_2 + Vars.v*(Cons.Fcal_3 + Vars.Fcal_SO_3)))
    dEdV = -0.5*Cons.M*Cons.nu*Vars.v*(2.0*Cons.E_0 + Vars.v**2*(4.0*Cons.E_2 + 5.0*Vars.E_SO_3*Vars.v))
    Absorption = 0
    dvdt_T1 = (-Absorption - Flux)/dEdV
    dydt=np.zeros(9)
    [dydt[5],dydt[6],dydt[7],dydt[8]] = FrameFromAngularVelocityIntegrand(Vars.R, OmegaVec_1p5(Cons,Vars)[1:])
    dydt[0] = dvdt_T4[0]
    if(Cons.EvolveSpin1):
        dydt[1], dydt[2]=FrameFromAngularVelocity_2D_Integrand(Vars.rfrak_chi1[0], Vars.rfrak_chi1[1],(mul(mul(inverse(Cons.S_chi1),OmegaVec_chiVec_1_1p5(Cons,Vars)),Cons.S_chi1))[1:])
    else:
        dydt[1], dydt[2] = 0.0, 0.0
    if(Cons.EvolveSpin2):
        dydt[3], dydt[4]=FrameFromAngularVelocity_2D_Integrand(Vars.rfrak_chi2[0], Vars.rfrak_chi2[1],(mul(mul(inverse(Cons.S_chi2),OmegaVec_chiVec_2_1p5(Cons,Vars)),Cons.S_chi2))[1:])
    else:
        dydt[3], dydt[4] = 0.0, 0.0
    return dydt      

@njit(cache=True)
def TaylorT5_1p5(Cons,Vars):
    Flux = Vars.Fcal_coeff*(Cons.Fcal_0 + Vars.v**2*(Cons.Fcal_2 + Vars.v*(Cons.Fcal_3 + Vars.Fcal_SO_3)))
    dEdV = -0.5*Cons.M*Cons.nu*Vars.v*(2.0*Cons.E_0 + Vars.v**2*(4.0*Cons.E_2 + 5.0*Vars.E_SO_3*Vars.v))
    Absorption = 0
    dvdt_T1 = (-Absorption - Flux)/dEdV
    dydt=np.zeros(8)
    [dydt[5],dydt[6],dydt[7],dydt[8]] = FrameFromAngularVelocityIntegrand(Vars.R, OmegaVec_1p5(Cons,Vars)[1:])
    dydt[0] = dvdt_T5[0]
    if(Cons.EvolveSpin1):
        dydt[1], dydt[2]=FrameFromAngularVelocity_2D_Integrand(Vars.rfrak_chi1[0], Vars.rfrak_chi1[1],(mul(mul(inverse(Cons.S_chi1),OmegaVec_chiVec_1_1p5(Cons,Vars)),Cons.S_chi1))[1:])
    else:
        dydt[1], dydt[2] = 0.0, 0.0
    if(Cons.EvolveSpin2):
        dydt[3], dydt[4]=FrameFromAngularVelocity_2D_Integrand(Vars.rfrak_chi2[0], Vars.rfrak_chi2[1],(mul(mul(inverse(Cons.S_chi2),OmegaVec_chiVec_2_1p5(Cons,Vars)),Cons.S_chi2))[1:])
    else:
        dydt[3], dydt[4] = 0.0, 0.0
    return dydt      

@njit(cache=True)
def Recalculate_2p0(Cons,Vars,y):
    Vars.v = np.array([y[0]])
    Vars.rfrak_chi1 = np.array([y[1],y[2]])
    Vars.rfrak_chi2 = np.array([y[3],y[4]])
    Vars.R = y[5]*Cons.wHat + y[6]*Cons.xHat + y[7]*Cons.yHat + y[8]*Cons.zHat
    Vars.nHat = mul(mul(Vars.R,Cons.xHat),inverse(Vars.R))
    Vars.lambdaHat = mul(mul(Vars.R,Cons.yHat),inverse(Vars.R))
    Vars.ellHat = mul(mul(Vars.R,Cons.zHat),inverse(Vars.R))
    Vars.R_S1 = exp(Vars.rfrak_chi1[0]*Cons.xHat + Vars.rfrak_chi1[1]*Cons.yHat)
    Vars.R_S2 = exp(Vars.rfrak_chi2[0]*Cons.xHat + Vars.rfrak_chi2[1]*Cons.yHat)
    Vars.chiVec1 = mul(mul(mul(Cons.S_chi1,Vars.R_S1),Cons.zHat),mul(conjugate(Vars.R_S1),conjugate(Cons.S_chi1)))
    Vars.chiVec2 = mul(mul(mul(Cons.S_chi2,Vars.R_S2),Cons.zHat),mul(conjugate(Vars.R_S2),conjugate(Cons.S_chi2)))
    Vars.chi1_n = np.array([dot(Vars.chiVec1[1:],Vars.nHat[1:])])
    Vars.chi1_lambda = np.array([dot(Vars.chiVec1[1:],Vars.lambdaHat[1:])])
    Vars.chi1_ell = np.array([dot(Vars.chiVec1[1:],Vars.ellHat[1:])])
    Vars.chi2_n = np.array([dot(Vars.chiVec2[1:],Vars.nHat[1:])])
    Vars.chi2_lambda = np.array([dot(Vars.chiVec2[1:],Vars.lambdaHat[1:])])
    Vars.chi2_ell = np.array([dot(Vars.chiVec2[1:],Vars.ellHat[1:])])
    Vars.S = Cons.M1**2*Vars.chiVec1 + Cons.M2**2*Vars.chiVec2
    Vars.S_ell = Cons.M1**2*Vars.chi1_ell + Cons.M2**2*Vars.chi2_ell
    Vars.S_n = Cons.M1**2*Vars.chi1_n + Cons.M2**2*Vars.chi2_n
    Vars.S_lambda = Cons.M1**2*Vars.chi1_lambda + Cons.M2**2*Vars.chi2_lambda
    Vars.Sigma = Cons.M*(-Cons.M1*Vars.chiVec1 + Cons.M2*Vars.chiVec2)
    Vars.Sigma_ell = Cons.M*(-Cons.M1*Vars.chi1_ell + Cons.M2*Vars.chi2_ell)
    Vars.Sigma_n = Cons.M*(-Cons.M1*Vars.chi1_n + Cons.M2*Vars.chi2_n)
    Vars.Sigma_lambda = Cons.M*(-Cons.M1*Vars.chi1_lambda + Cons.M2*Vars.chi2_lambda)
    Vars.SS = np.array([dot(Vars.S[1:],Vars.S[1:])])
    Vars.SigmaSigma = np.array([dot(Vars.Sigma[1:],Vars.Sigma[1:])])
    Vars.SSigma = np.array([dot(Vars.S[1:],Vars.Sigma[1:])])
    Vars.S1 = Cons.M1**2*Vars.chiVec1
    Vars.S1_n = Cons.M1**2*Vars.chi1_n
    Vars.S1_lambda = Cons.M1**2*Vars.chi1_lambda
    Vars.S2 = Cons.M2**2*Vars.chiVec2
    Vars.S2_n = Cons.M2**2*Vars.chi2_n
    Vars.S2_lambda = Cons.M2**2*Vars.chi2_lambda
    Vars.chi_s_ell = Vars.chi1_ell/2 + Vars.chi2_ell/2
    Vars.chi_a_ell = Vars.chi1_ell/2 - Vars.chi2_ell/2
    Vars.Fcal_coeff = 32*Cons.nu**2*Vars.v**10/5
    Vars.Fcal_SQ_4 = Cons.chi1chi1*(-89*Cons.delta/192 + 89*Cons.nu/96 - 89/192) - 103*Cons.chi1chi2*Cons.nu/48 + Cons.chi2chi2*(89*Cons.delta/192 + 89*Cons.nu/96 - 89/192) + Vars.chi_a_ell*(Vars.chi_a_ell*(287/96 - 12*Cons.nu) + 287*Vars.chi_s_ell*Cons.delta/48) + Vars.chi_s_ell**2*(Cons.nu/24 + 287/96)
    Vars.Fcal_SO_3 = (-4*Vars.S_ell - 5*Vars.Sigma_ell*Cons.delta/4)/Cons.M**2
    Vars.E_SQ_4 = -3*Vars.chi_a_ell**2/2 - 3*Vars.chi_s_ell**2/2 - Cons.delta*(Cons.chi2chi2/2 + 3*Vars.chi_a_ell*Vars.chi_s_ell) + Cons.nu*(Cons.chi1chi2 + 6*Vars.chi_a_ell**2) + (Cons.chi1chi1 + Cons.chi2chi2)*(Cons.delta - 2*Cons.nu + 1)/4
    Vars.E_SO_3 = (14*Vars.S_ell/3 + 2*Vars.Sigma_ell*Cons.delta)/Cons.M**2

@njit
def OmegaVec_chiVec_1_2p0(Cons,Vars):
    Omega1_coeff = Vars.v**5/Cons.M
    return Omega1_coeff*(Vars.S2*Vars.v**3*((2.0 - 0.5*Cons.delta)/Cons.M**2 + Vars.v*(-5.0*Vars.S_ell - 3.0*Vars.Sigma_ell*Cons.delta/Cons.M)/Cons.M**4) + Vars.ellHat*(-0.75*Cons.delta + 0.5*Cons.nu + Vars.v**2*(Cons.delta*(0.625*Cons.nu - 0.5625) + Cons.nu*(1.25 - 0.0416666666666667*Cons.nu) + Vars.v*(Vars.v*(Cons.nu*(Cons.nu*(0.229166666666667*Cons.nu - 4.94791666666667) + 2.89583333333333) + 0.96875 + (Cons.delta*(Cons.nu*(5.91666666666667 - 0.489583333333333*Cons.nu) - 0.96875) + (Vars.S_ell*(Vars.S_ell*(Cons.kappa_p*(-0.25*Cons.nu - 0.375) - 0.5*Cons.nu - 0.75) + Vars.Sigma_ell*(Cons.delta*(Cons.kappa_p*(-0.25*Cons.nu - 0.375) - 0.5*Cons.nu - 0.75) + Cons.kappa_m*(0.25*Cons.nu + 0.375))) + Vars.Sigma_ell**2*(Cons.delta*Cons.kappa_m*(0.125*Cons.nu + 0.1875) + Cons.kappa_p*(Cons.nu*(0.25*Cons.nu + 0.25) - 0.1875) + Cons.nu*(0.5*Cons.nu + 0.75)) + (Vars.S_ell*(Vars.S_ell*Cons.delta*(0.375*Cons.kappa_p + 0.75) + Vars.Sigma_ell*Cons.delta*(Cons.delta*(0.375*Cons.kappa_p + 0.75) - 0.375*Cons.kappa_m)) + Vars.Sigma_ell**2*Cons.delta*(-0.1875*Cons.delta*Cons.kappa_m + Cons.kappa_p*(0.1875 - 0.375*Cons.nu) - 0.75*Cons.nu))/Cons.M)/Cons.M**7)/Cons.M) + (Vars.S_ell*(0.833333333333333*Cons.nu + 1.25) + (-1.25*Vars.S_ell*Cons.delta + Vars.Sigma_ell*Cons.delta*(0.5*Cons.nu + 0.75) - 0.75*Vars.Sigma_ell*Cons.delta**2/Cons.M)/Cons.M)/Cons.M**2) + 0.5625) + 0.75) + Vars.lambdaHat*Vars.v**3*(Vars.S2_lambda*(0.5*Cons.delta + 2.0)/Cons.M**2 + Vars.v*(Vars.S1_lambda*(-0.75*Cons.delta - 3.0*Cons.kappa_1 + 2.25) + Vars.S1_lambda*(Cons.delta*(0.75 - 1.5*Cons.kappa_1) + 1.5*Cons.kappa_1 - 0.75)/Cons.nu)/Cons.M**2) + Vars.nHat*Vars.v*(3.0*Vars.chi1_n*Cons.nu + Vars.v**2*((Vars.S1_n*(Cons.delta*(-2.25*Cons.kappa_1 - 0.25) + Cons.kappa_1*(1.5*Cons.nu + 0.75) + 0.5*Cons.nu + 3.25) + Vars.S1_n*(Cons.delta*(1.5 - 0.75*Cons.kappa_1) + 0.75*Cons.kappa_1 - 1.5)/Cons.nu + Vars.S2_n*(1.5*Cons.delta - 1.0*Cons.nu + 3.0))/Cons.M**2 + Vars.v*(-15.0*Vars.S1_n*Vars.S_ell*Cons.kappa_1 + Vars.S1_n*Vars.S_ell*(-7.5*Cons.delta*Cons.kappa_1 + 7.5*Cons.kappa_1)/Cons.nu + 15.0*Vars.S2_n*Vars.S_ell + (-9.0*Vars.S1_n*Vars.Sigma_ell*Cons.delta*Cons.kappa_1 + Vars.S1_n*Vars.Sigma_ell*Cons.delta*(-4.5*Cons.delta*Cons.kappa_1 + 4.5*Cons.kappa_1)/Cons.nu + 9.0*Vars.S2_n*Vars.Sigma_ell*Cons.delta)/Cons.M)/Cons.M**4) + 3.0*Cons.M2**2*Vars.chi2_n/Cons.M**2) - Cons.M2**2*Vars.chiVec2*Vars.v/Cons.M**2)

@njit
def OmegaVec_chiVec_2_2p0(Cons,Vars):
    Omega2_coeff = Vars.v**5/Cons.M
    return Omega2_coeff*(Vars.S1*Vars.v**3*((0.5*Cons.delta + 2.0)/Cons.M**2 + Vars.v*(-5.0*Vars.S_ell + 3.0*Vars.Sigma_ell*Cons.delta/Cons.M)/Cons.M**4) + Vars.ellHat*(0.75*Cons.delta + 0.5*Cons.nu + Vars.v**2*(Cons.delta*(0.5625 - 0.625*Cons.nu) + Cons.nu*(1.25 - 0.0416666666666667*Cons.nu) + Vars.v*(Vars.v*(Cons.nu*(Cons.nu*(0.229166666666667*Cons.nu - 4.94791666666667) + 2.89583333333333) + 0.96875 + (Cons.delta*(Cons.nu*(0.489583333333333*Cons.nu - 5.91666666666667) + 0.96875) + (Vars.S_ell*(Vars.S_ell*(Cons.kappa_p*(-0.25*Cons.nu - 0.375) - 0.5*Cons.nu - 0.75) + Vars.Sigma_ell*(Cons.delta*(Cons.kappa_p*(0.25*Cons.nu + 0.375) + 0.5*Cons.nu + 0.75) + Cons.kappa_m*(-0.25*Cons.nu - 0.375))) + Vars.Sigma_ell**2*(Cons.delta*Cons.kappa_m*(0.125*Cons.nu + 0.1875) + Cons.kappa_p*(Cons.nu*(0.25*Cons.nu + 0.25) - 0.1875) + Cons.nu*(0.5*Cons.nu + 0.75)) + (Vars.S_ell*(Vars.S_ell*Cons.delta*(-0.375*Cons.kappa_p - 0.75) + Vars.Sigma_ell*Cons.delta*(Cons.delta*(0.375*Cons.kappa_p + 0.75) - 0.375*Cons.kappa_m)) + Vars.Sigma_ell**2*Cons.delta*(0.1875*Cons.delta*Cons.kappa_m + Cons.kappa_p*(0.375*Cons.nu - 0.1875) + 0.75*Cons.nu))/Cons.M)/Cons.M**7)/Cons.M) + (Vars.S_ell*(0.833333333333333*Cons.nu + 1.25) + (1.25*Vars.S_ell*Cons.delta + Vars.Sigma_ell*Cons.delta*(-0.5*Cons.nu - 0.75) - 0.75*Vars.Sigma_ell*Cons.delta**2/Cons.M)/Cons.M)/Cons.M**2) + 0.5625) + 0.75) + Vars.lambdaHat*Vars.v**3*(Vars.S1_lambda*(2.0 - 0.5*Cons.delta)/Cons.M**2 + Vars.v*(Vars.S2_lambda*(0.75*Cons.delta - 3.0*Cons.kappa_2 + 2.25) + Vars.S2_lambda*(Cons.delta*(1.5*Cons.kappa_2 - 0.75) + 1.5*Cons.kappa_2 - 0.75)/Cons.nu)/Cons.M**2) + Vars.nHat*Vars.v*(3.0*Vars.chi2_n*Cons.nu + Vars.v**2*((Vars.S1_n*(-1.5*Cons.delta - 1.0*Cons.nu + 3.0) + Vars.S2_n*(Cons.delta*(2.25*Cons.kappa_2 + 0.25) + Cons.kappa_2*(1.5*Cons.nu + 0.75) + 0.5*Cons.nu + 3.25) + Vars.S2_n*(Cons.delta*(0.75*Cons.kappa_2 - 1.5) + 0.75*Cons.kappa_2 - 1.5)/Cons.nu)/Cons.M**2 + Vars.v*(15.0*Vars.S1_n*Vars.S_ell - 15.0*Vars.S2_n*Vars.S_ell*Cons.kappa_2 + Vars.S2_n*Vars.S_ell*(7.5*Cons.delta*Cons.kappa_2 + 7.5*Cons.kappa_2)/Cons.nu + (-9.0*Vars.S1_n*Vars.Sigma_ell*Cons.delta + 9.0*Vars.S2_n*Vars.Sigma_ell*Cons.delta*Cons.kappa_2 + Vars.S2_n*Vars.Sigma_ell*Cons.delta*(-4.5*Cons.delta*Cons.kappa_2 - 4.5*Cons.kappa_2)/Cons.nu)/Cons.M)/Cons.M**4) + 3.0*Cons.M1**2*Vars.chi1_n/Cons.M**2) - Cons.M1**2*Vars.chiVec1*Vars.v/Cons.M**2)

@njit
def OmegaVec_2p0(Cons,Vars):
    gamma_PN_0 = 1.00000000000000
    a_ell_4 = Vars.S_n*(5.77777777777778*Cons.nu**2 + 14.75*Cons.nu + 1.5) + Vars.Sigma_n*Cons.delta*(2.83333333333333*Cons.nu**2 + 9.125*Cons.nu + 1.5)
    a_ell_1 = 3.0*Vars.Sigma_ell*Vars.Sigma_n*Cons.nu*(Cons.kappa_p + 2.0)/Cons.M**2 - 3.0*Cons.delta*(Vars.S_ell*Vars.Sigma_n + Vars.S_n*Vars.Sigma_ell)/Cons.M**2 - 0.5*(3.0*Vars.S_n*(2.0*Vars.S_ell*(Cons.kappa_p + 2.0) + Vars.Sigma_ell*(Cons.delta*Cons.kappa_p - Cons.kappa_m)) + 3.0*Vars.Sigma_n*(Vars.S_ell*(Cons.delta*Cons.kappa_p - Cons.kappa_m) + Vars.Sigma_ell*(-Cons.delta*Cons.kappa_m + Cons.kappa_p)))/Cons.M**2
    a_ell_2 = Vars.S_n*(-9.66666666666667*Cons.nu - 10.0) + Vars.Sigma_n*Cons.delta*(-4.5*Cons.nu - 6.0)
    gamma_PN_2 = 1.0 - 0.333333333333333*Cons.nu
    a_ell_3 = Cons.nu*(4.0*Cons.delta*(Vars.S_ell*Vars.Sigma_n + Vars.S_n*Vars.Sigma_ell)/Cons.M**2 + 0.5*(2.0*Vars.S_n*(4.0*Vars.S_ell*(Cons.kappa_p + 2.0) + Vars.Sigma_ell*(2.0*Cons.delta*Cons.kappa_p + 7.0*Cons.kappa_m)) + Vars.Sigma_n*(2.0*Vars.S_ell*(2.0*Cons.delta*Cons.kappa_p + 7.0*Cons.kappa_m) + Vars.Sigma_ell*(5.0*Cons.delta*Cons.kappa_m - 17.0*Cons.kappa_p - 26.0)))/Cons.M**2) - 4.0*Vars.Sigma_ell*Vars.Sigma_n*Cons.nu**2*(Cons.kappa_p + 2.0)/Cons.M**2 + 6.0*Cons.delta*(Vars.S_ell*Vars.Sigma_n + Vars.S_n*Vars.Sigma_ell)/Cons.M**2 + 0.5*(Vars.S_n*(Vars.S_ell*(-9.0*Cons.delta*Cons.kappa_m + 3.0*Cons.kappa_p + 22.0) + Vars.Sigma_ell*(6.0*Cons.delta*Cons.kappa_p - 6.0*Cons.kappa_m)) + 6.0*Vars.Sigma_n*(Vars.S_ell*(Cons.delta*Cons.kappa_p - Cons.kappa_m) + Vars.Sigma_ell*(-Cons.delta*Cons.kappa_m + Cons.kappa_p + 1.0)))/Cons.M**2
    gamma_PN_3 = (1.66666666666667*Vars.S_ell + Vars.Sigma_ell*Cons.delta)/Cons.M**2
    a_ell_0 = 7.0*Vars.S_n + 3.0*Vars.Sigma_n*Cons.delta
    gamma_PN_4 = -5.41666666666667*Cons.nu + 1.0 + (Vars.S_ell**2*(-0.5*Cons.kappa_p - 1.0) + Vars.S_ell*Vars.Sigma_ell*(-0.5*Cons.delta*Cons.kappa_p - Cons.delta + 0.5*Cons.kappa_m) + Vars.Sigma_ell**2*(0.25*Cons.delta*Cons.kappa_m - 0.25*Cons.kappa_p + Cons.nu*(0.5*Cons.kappa_p + 1.0)))/Cons.M**8
    return Vars.ellHat*Vars.v**3/Cons.M + Vars.nHat*Vars.v**6*(a_ell_0 + Vars.v*(a_ell_1 + Vars.v*(a_ell_2 + Vars.v*(a_ell_3 + a_ell_4*Vars.v))))*(gamma_PN_0 + Vars.v**2*(gamma_PN_2 + Vars.v*(gamma_PN_3 + gamma_PN_4*Vars.v)))/Cons.M**3


@njit(cache=True)
def TaylorT1_2p0(Cons,Vars):
    Flux = Vars.Fcal_coeff*(Cons.Fcal_0 + Vars.v**2*(Cons.Fcal_2 + Vars.v*(Cons.Fcal_3 + Vars.Fcal_SO_3 + Vars.v*(Cons.Fcal_4 + Vars.Fcal_SQ_4))))
    dEdV = -0.5*Cons.M*Cons.nu*Vars.v*(2.0*Cons.E_0 + Vars.v**2*(4.0*Cons.E_2 + Vars.v*(5.0*Vars.E_SO_3 + 6.0*Vars.v*(Cons.E_4 + Vars.E_SQ_4))))
    Absorption = 0
    dvdt_T1 = (-Absorption - Flux)/dEdV
    dydt=np.zeros(9)
    [dydt[5],dydt[6],dydt[7],dydt[8]] = FrameFromAngularVelocityIntegrand(Vars.R, OmegaVec_2p0(Cons,Vars)[1:])
    dydt[0] = dvdt_T1[0] 
    if(Cons.EvolveSpin1):
        dydt[1], dydt[2]=FrameFromAngularVelocity_2D_Integrand(Vars.rfrak_chi1[0], Vars.rfrak_chi1[1],(mul(mul(inverse(Cons.S_chi1),OmegaVec_chiVec_1_2p0(Cons,Vars)),Cons.S_chi1))[1:])
    else:
        dydt[1], dydt[2] = 0.0, 0.0
    if(Cons.EvolveSpin2):
        dydt[3], dydt[4]=FrameFromAngularVelocity_2D_Integrand(Vars.rfrak_chi2[0], Vars.rfrak_chi2[1],(mul(mul(inverse(Cons.S_chi2),OmegaVec_chiVec_2_2p0(Cons,Vars)),Cons.S_chi2))[1:])
    else:
        dydt[3], dydt[4] = 0.0, 0.0
    return dydt

@njit(cache=True)
def TaylorT4_2p0(Cons,Vars):
    Flux = Vars.Fcal_coeff*(Cons.Fcal_0 + Vars.v**2*(Cons.Fcal_2 + Vars.v*(Cons.Fcal_3 + Vars.Fcal_SO_3 + Vars.v*(Cons.Fcal_4 + Vars.Fcal_SQ_4))))
    dEdV = -0.5*Cons.M*Cons.nu*Vars.v*(2.0*Cons.E_0 + Vars.v**2*(4.0*Cons.E_2 + Vars.v*(5.0*Vars.E_SO_3 + 6.0*Vars.v*(Cons.E_4 + Vars.E_SQ_4))))
    Absorption = 0
    dvdt_T1 = (-Absorption - Flux)/dEdV
    dydt=np.zeros(9)
    [dydt[5],dydt[6],dydt[7],dydt[8]] = FrameFromAngularVelocityIntegrand(Vars.R, OmegaVec_2p0(Cons,Vars)[1:])
    dydt[0] = dvdt_T4[0]
    if(Cons.EvolveSpin1):
        dydt[1], dydt[2]=FrameFromAngularVelocity_2D_Integrand(Vars.rfrak_chi1[0], Vars.rfrak_chi1[1],(mul(mul(inverse(Cons.S_chi1),OmegaVec_chiVec_1_2p0(Cons,Vars)),Cons.S_chi1))[1:])
    else:
        dydt[1], dydt[2] = 0.0, 0.0
    if(Cons.EvolveSpin2):
        dydt[3], dydt[4]=FrameFromAngularVelocity_2D_Integrand(Vars.rfrak_chi2[0], Vars.rfrak_chi2[1],(mul(mul(inverse(Cons.S_chi2),OmegaVec_chiVec_2_2p0(Cons,Vars)),Cons.S_chi2))[1:])
    else:
        dydt[3], dydt[4] = 0.0, 0.0
    return dydt      

@njit(cache=True)
def TaylorT5_2p0(Cons,Vars):
    Flux = Vars.Fcal_coeff*(Cons.Fcal_0 + Vars.v**2*(Cons.Fcal_2 + Vars.v*(Cons.Fcal_3 + Vars.Fcal_SO_3 + Vars.v*(Cons.Fcal_4 + Vars.Fcal_SQ_4))))
    dEdV = -0.5*Cons.M*Cons.nu*Vars.v*(2.0*Cons.E_0 + Vars.v**2*(4.0*Cons.E_2 + Vars.v*(5.0*Vars.E_SO_3 + 6.0*Vars.v*(Cons.E_4 + Vars.E_SQ_4))))
    Absorption = 0
    dvdt_T1 = (-Absorption - Flux)/dEdV
    dydt=np.zeros(8)
    [dydt[5],dydt[6],dydt[7],dydt[8]] = FrameFromAngularVelocityIntegrand(Vars.R, OmegaVec_2p0(Cons,Vars)[1:])
    dydt[0] = dvdt_T5[0]
    if(Cons.EvolveSpin1):
        dydt[1], dydt[2]=FrameFromAngularVelocity_2D_Integrand(Vars.rfrak_chi1[0], Vars.rfrak_chi1[1],(mul(mul(inverse(Cons.S_chi1),OmegaVec_chiVec_1_2p0(Cons,Vars)),Cons.S_chi1))[1:])
    else:
        dydt[1], dydt[2] = 0.0, 0.0
    if(Cons.EvolveSpin2):
        dydt[3], dydt[4]=FrameFromAngularVelocity_2D_Integrand(Vars.rfrak_chi2[0], Vars.rfrak_chi2[1],(mul(mul(inverse(Cons.S_chi2),OmegaVec_chiVec_2_2p0(Cons,Vars)),Cons.S_chi2))[1:])
    else:
        dydt[3], dydt[4] = 0.0, 0.0
    return dydt      

@njit(cache=True)
def Recalculate_2p5(Cons,Vars,y):
    Vars.v = np.array([y[0]])
    Vars.rfrak_chi1 = np.array([y[1],y[2]])
    Vars.rfrak_chi2 = np.array([y[3],y[4]])
    Vars.R = y[5]*Cons.wHat + y[6]*Cons.xHat + y[7]*Cons.yHat + y[8]*Cons.zHat
    Vars.nHat = mul(mul(Vars.R,Cons.xHat),inverse(Vars.R))
    Vars.lambdaHat = mul(mul(Vars.R,Cons.yHat),inverse(Vars.R))
    Vars.ellHat = mul(mul(Vars.R,Cons.zHat),inverse(Vars.R))
    Vars.R_S1 = exp(Vars.rfrak_chi1[0]*Cons.xHat + Vars.rfrak_chi1[1]*Cons.yHat)
    Vars.R_S2 = exp(Vars.rfrak_chi2[0]*Cons.xHat + Vars.rfrak_chi2[1]*Cons.yHat)
    Vars.chiVec1 = mul(mul(mul(Cons.S_chi1,Vars.R_S1),Cons.zHat),mul(conjugate(Vars.R_S1),conjugate(Cons.S_chi1)))
    Vars.chiVec2 = mul(mul(mul(Cons.S_chi2,Vars.R_S2),Cons.zHat),mul(conjugate(Vars.R_S2),conjugate(Cons.S_chi2)))
    Vars.chi1_n = np.array([dot(Vars.chiVec1[1:],Vars.nHat[1:])])
    Vars.chi1_lambda = np.array([dot(Vars.chiVec1[1:],Vars.lambdaHat[1:])])
    Vars.chi1_ell = np.array([dot(Vars.chiVec1[1:],Vars.ellHat[1:])])
    Vars.chi2_n = np.array([dot(Vars.chiVec2[1:],Vars.nHat[1:])])
    Vars.chi2_lambda = np.array([dot(Vars.chiVec2[1:],Vars.lambdaHat[1:])])
    Vars.chi2_ell = np.array([dot(Vars.chiVec2[1:],Vars.ellHat[1:])])
    Vars.S = Cons.M1**2*Vars.chiVec1 + Cons.M2**2*Vars.chiVec2
    Vars.S_ell = Cons.M1**2*Vars.chi1_ell + Cons.M2**2*Vars.chi2_ell
    Vars.S_n = Cons.M1**2*Vars.chi1_n + Cons.M2**2*Vars.chi2_n
    Vars.S_lambda = Cons.M1**2*Vars.chi1_lambda + Cons.M2**2*Vars.chi2_lambda
    Vars.Sigma = Cons.M*(-Cons.M1*Vars.chiVec1 + Cons.M2*Vars.chiVec2)
    Vars.Sigma_ell = Cons.M*(-Cons.M1*Vars.chi1_ell + Cons.M2*Vars.chi2_ell)
    Vars.Sigma_n = Cons.M*(-Cons.M1*Vars.chi1_n + Cons.M2*Vars.chi2_n)
    Vars.Sigma_lambda = Cons.M*(-Cons.M1*Vars.chi1_lambda + Cons.M2*Vars.chi2_lambda)
    Vars.SS = np.array([dot(Vars.S[1:],Vars.S[1:])])
    Vars.SigmaSigma = np.array([dot(Vars.Sigma[1:],Vars.Sigma[1:])])
    Vars.SSigma = np.array([dot(Vars.S[1:],Vars.Sigma[1:])])
    Vars.S1 = Cons.M1**2*Vars.chiVec1
    Vars.S1_n = Cons.M1**2*Vars.chi1_n
    Vars.S1_lambda = Cons.M1**2*Vars.chi1_lambda
    Vars.S2 = Cons.M2**2*Vars.chiVec2
    Vars.S2_n = Cons.M2**2*Vars.chi2_n
    Vars.S2_lambda = Cons.M2**2*Vars.chi2_lambda
    Vars.chi_s_ell = Vars.chi1_ell/2 + Vars.chi2_ell/2
    Vars.chi_a_ell = Vars.chi1_ell/2 - Vars.chi2_ell/2
    Vars.Fcal_coeff = 32*Cons.nu**2*Vars.v**10/5
    Vars.Fcal_SQ_4 = Cons.chi1chi1*(-89*Cons.delta/192 + 89*Cons.nu/96 - 89/192) - 103*Cons.chi1chi2*Cons.nu/48 + Cons.chi2chi2*(89*Cons.delta/192 + 89*Cons.nu/96 - 89/192) + Vars.chi_a_ell*(Vars.chi_a_ell*(287/96 - 12*Cons.nu) + 287*Vars.chi_s_ell*Cons.delta/48) + Vars.chi_s_ell**2*(Cons.nu/24 + 287/96)
    Vars.Fcal_SO_3 = (-4*Vars.S_ell - 5*Vars.Sigma_ell*Cons.delta/4)/Cons.M**2
    Vars.Fcal_SO_5 = (Vars.S_ell*(272*Cons.nu/9 - 9/2) + Vars.Sigma_ell*Cons.delta*(43*Cons.nu/4 - 13/16))/Cons.M**2
    Vars.E_SQ_4 = -3*Vars.chi_a_ell**2/2 - 3*Vars.chi_s_ell**2/2 - Cons.delta*(Cons.chi2chi2/2 + 3*Vars.chi_a_ell*Vars.chi_s_ell) + Cons.nu*(Cons.chi1chi2 + 6*Vars.chi_a_ell**2) + (Cons.chi1chi1 + Cons.chi2chi2)*(Cons.delta - 2*Cons.nu + 1)/4
    Vars.E_SO_3 = (14*Vars.S_ell/3 + 2*Vars.Sigma_ell*Cons.delta)/Cons.M**2
    Vars.E_SO_5 = (Vars.S_ell*(11 - 61*Cons.nu/9) + Vars.Sigma_ell*Cons.delta*(3 - 10*Cons.nu/3))/Cons.M**2
    Vars.MDot_Alvi_5 = -Cons.M1**3*Vars.chi1_ell*(3*Cons.chi1chi1 + 1)/(4*Cons.M**3) - Cons.M2**3*Vars.chi2_ell*(3*Cons.chi2chi2 + 1)/(4*Cons.M**3)

@njit
def OmegaVec_chiVec_1_2p5(Cons,Vars):
    Omega1_coeff = Vars.v**5/Cons.M
    return Omega1_coeff*(Vars.S2*Vars.v**3*((2.0 - 0.5*Cons.delta)/Cons.M**2 + Vars.v*(-5.0*Vars.S_ell - 3.0*Vars.Sigma_ell*Cons.delta/Cons.M)/Cons.M**4) + Vars.ellHat*(-0.75*Cons.delta + 0.5*Cons.nu + Vars.v**2*(Cons.delta*(0.625*Cons.nu - 0.5625) + Cons.nu*(1.25 - 0.0416666666666667*Cons.nu) + Vars.v*(Vars.v*(Cons.nu*(Cons.nu*(0.229166666666667*Cons.nu - 4.94791666666667) + 2.89583333333333) + 0.96875 + (Cons.delta*(Cons.nu*(5.91666666666667 - 0.489583333333333*Cons.nu) - 0.96875) + (Vars.S_ell*(Vars.S_ell*(Cons.kappa_p*(-0.25*Cons.nu - 0.375) - 0.5*Cons.nu - 0.75) + Vars.Sigma_ell*(Cons.delta*(Cons.kappa_p*(-0.25*Cons.nu - 0.375) - 0.5*Cons.nu - 0.75) + Cons.kappa_m*(0.25*Cons.nu + 0.375))) + Vars.Sigma_ell**2*(Cons.delta*Cons.kappa_m*(0.125*Cons.nu + 0.1875) + Cons.kappa_p*(Cons.nu*(0.25*Cons.nu + 0.25) - 0.1875) + Cons.nu*(0.5*Cons.nu + 0.75)) + (Vars.S_ell*(Vars.S_ell*Cons.delta*(0.375*Cons.kappa_p + 0.75) + Vars.Sigma_ell*Cons.delta*(Cons.delta*(0.375*Cons.kappa_p + 0.75) - 0.375*Cons.kappa_m)) + Vars.Sigma_ell**2*Cons.delta*(-0.1875*Cons.delta*Cons.kappa_m + Cons.kappa_p*(0.1875 - 0.375*Cons.nu) - 0.75*Cons.nu))/Cons.M)/Cons.M**7)/Cons.M) + (Vars.S_ell*(0.833333333333333*Cons.nu + 1.25) + (-1.25*Vars.S_ell*Cons.delta + Vars.Sigma_ell*Cons.delta*(0.5*Cons.nu + 0.75) - 0.75*Vars.Sigma_ell*Cons.delta**2/Cons.M)/Cons.M)/Cons.M**2) + 0.5625) + 0.75) + Vars.lambdaHat*Vars.v**3*(Vars.S2_lambda*(0.5*Cons.delta + 2.0)/Cons.M**2 + Vars.v*(Vars.S1_lambda*(-0.75*Cons.delta - 3.0*Cons.kappa_1 + 2.25) + Vars.S1_lambda*(Cons.delta*(0.75 - 1.5*Cons.kappa_1) + 1.5*Cons.kappa_1 - 0.75)/Cons.nu)/Cons.M**2) + Vars.nHat*Vars.v*(3.0*Vars.chi1_n*Cons.nu + Vars.v**2*((Vars.S1_n*(Cons.delta*(-2.25*Cons.kappa_1 - 0.25) + Cons.kappa_1*(1.5*Cons.nu + 0.75) + 0.5*Cons.nu + 3.25) + Vars.S1_n*(Cons.delta*(1.5 - 0.75*Cons.kappa_1) + 0.75*Cons.kappa_1 - 1.5)/Cons.nu + Vars.S2_n*(1.5*Cons.delta - 1.0*Cons.nu + 3.0))/Cons.M**2 + Vars.v*(-15.0*Vars.S1_n*Vars.S_ell*Cons.kappa_1 + Vars.S1_n*Vars.S_ell*(-7.5*Cons.delta*Cons.kappa_1 + 7.5*Cons.kappa_1)/Cons.nu + 15.0*Vars.S2_n*Vars.S_ell + (-9.0*Vars.S1_n*Vars.Sigma_ell*Cons.delta*Cons.kappa_1 + Vars.S1_n*Vars.Sigma_ell*Cons.delta*(-4.5*Cons.delta*Cons.kappa_1 + 4.5*Cons.kappa_1)/Cons.nu + 9.0*Vars.S2_n*Vars.Sigma_ell*Cons.delta)/Cons.M)/Cons.M**4) + 3.0*Cons.M2**2*Vars.chi2_n/Cons.M**2) - Cons.M2**2*Vars.chiVec2*Vars.v/Cons.M**2)

@njit
def OmegaVec_chiVec_2_2p5(Cons,Vars):
    Omega2_coeff = Vars.v**5/Cons.M
    return Omega2_coeff*(Vars.S1*Vars.v**3*((0.5*Cons.delta + 2.0)/Cons.M**2 + Vars.v*(-5.0*Vars.S_ell + 3.0*Vars.Sigma_ell*Cons.delta/Cons.M)/Cons.M**4) + Vars.ellHat*(0.75*Cons.delta + 0.5*Cons.nu + Vars.v**2*(Cons.delta*(0.5625 - 0.625*Cons.nu) + Cons.nu*(1.25 - 0.0416666666666667*Cons.nu) + Vars.v*(Vars.v*(Cons.nu*(Cons.nu*(0.229166666666667*Cons.nu - 4.94791666666667) + 2.89583333333333) + 0.96875 + (Cons.delta*(Cons.nu*(0.489583333333333*Cons.nu - 5.91666666666667) + 0.96875) + (Vars.S_ell*(Vars.S_ell*(Cons.kappa_p*(-0.25*Cons.nu - 0.375) - 0.5*Cons.nu - 0.75) + Vars.Sigma_ell*(Cons.delta*(Cons.kappa_p*(0.25*Cons.nu + 0.375) + 0.5*Cons.nu + 0.75) + Cons.kappa_m*(-0.25*Cons.nu - 0.375))) + Vars.Sigma_ell**2*(Cons.delta*Cons.kappa_m*(0.125*Cons.nu + 0.1875) + Cons.kappa_p*(Cons.nu*(0.25*Cons.nu + 0.25) - 0.1875) + Cons.nu*(0.5*Cons.nu + 0.75)) + (Vars.S_ell*(Vars.S_ell*Cons.delta*(-0.375*Cons.kappa_p - 0.75) + Vars.Sigma_ell*Cons.delta*(Cons.delta*(0.375*Cons.kappa_p + 0.75) - 0.375*Cons.kappa_m)) + Vars.Sigma_ell**2*Cons.delta*(0.1875*Cons.delta*Cons.kappa_m + Cons.kappa_p*(0.375*Cons.nu - 0.1875) + 0.75*Cons.nu))/Cons.M)/Cons.M**7)/Cons.M) + (Vars.S_ell*(0.833333333333333*Cons.nu + 1.25) + (1.25*Vars.S_ell*Cons.delta + Vars.Sigma_ell*Cons.delta*(-0.5*Cons.nu - 0.75) - 0.75*Vars.Sigma_ell*Cons.delta**2/Cons.M)/Cons.M)/Cons.M**2) + 0.5625) + 0.75) + Vars.lambdaHat*Vars.v**3*(Vars.S1_lambda*(2.0 - 0.5*Cons.delta)/Cons.M**2 + Vars.v*(Vars.S2_lambda*(0.75*Cons.delta - 3.0*Cons.kappa_2 + 2.25) + Vars.S2_lambda*(Cons.delta*(1.5*Cons.kappa_2 - 0.75) + 1.5*Cons.kappa_2 - 0.75)/Cons.nu)/Cons.M**2) + Vars.nHat*Vars.v*(3.0*Vars.chi2_n*Cons.nu + Vars.v**2*((Vars.S1_n*(-1.5*Cons.delta - 1.0*Cons.nu + 3.0) + Vars.S2_n*(Cons.delta*(2.25*Cons.kappa_2 + 0.25) + Cons.kappa_2*(1.5*Cons.nu + 0.75) + 0.5*Cons.nu + 3.25) + Vars.S2_n*(Cons.delta*(0.75*Cons.kappa_2 - 1.5) + 0.75*Cons.kappa_2 - 1.5)/Cons.nu)/Cons.M**2 + Vars.v*(15.0*Vars.S1_n*Vars.S_ell - 15.0*Vars.S2_n*Vars.S_ell*Cons.kappa_2 + Vars.S2_n*Vars.S_ell*(7.5*Cons.delta*Cons.kappa_2 + 7.5*Cons.kappa_2)/Cons.nu + (-9.0*Vars.S1_n*Vars.Sigma_ell*Cons.delta + 9.0*Vars.S2_n*Vars.Sigma_ell*Cons.delta*Cons.kappa_2 + Vars.S2_n*Vars.Sigma_ell*Cons.delta*(-4.5*Cons.delta*Cons.kappa_2 - 4.5*Cons.kappa_2)/Cons.nu)/Cons.M)/Cons.M**4) + 3.0*Cons.M1**2*Vars.chi1_n/Cons.M**2) - Cons.M1**2*Vars.chiVec1*Vars.v/Cons.M**2)

@njit
def OmegaVec_2p5(Cons,Vars):
    gamma_PN_5 = (0.888888888888889*Vars.S_ell*Cons.nu + 3.33333333333333*Vars.S_ell)/Cons.M**2 + 2.0*Vars.Sigma_ell*Cons.delta/Cons.M**3
    gamma_PN_0 = 1.00000000000000
    a_ell_4 = Vars.S_n*(5.77777777777778*Cons.nu**2 + 14.75*Cons.nu + 1.5) + Vars.Sigma_n*Cons.delta*(2.83333333333333*Cons.nu**2 + 9.125*Cons.nu + 1.5)
    a_ell_1 = 3.0*Vars.Sigma_ell*Vars.Sigma_n*Cons.nu*(Cons.kappa_p + 2.0)/Cons.M**2 - 3.0*Cons.delta*(Vars.S_ell*Vars.Sigma_n + Vars.S_n*Vars.Sigma_ell)/Cons.M**2 - 0.5*(3.0*Vars.S_n*(2.0*Vars.S_ell*(Cons.kappa_p + 2.0) + Vars.Sigma_ell*(Cons.delta*Cons.kappa_p - Cons.kappa_m)) + 3.0*Vars.Sigma_n*(Vars.S_ell*(Cons.delta*Cons.kappa_p - Cons.kappa_m) + Vars.Sigma_ell*(-Cons.delta*Cons.kappa_m + Cons.kappa_p)))/Cons.M**2
    a_ell_2 = Vars.S_n*(-9.66666666666667*Cons.nu - 10.0) + Vars.Sigma_n*Cons.delta*(-4.5*Cons.nu - 6.0)
    gamma_PN_2 = 1.0 - 0.333333333333333*Cons.nu
    a_ell_3 = Cons.nu*(4.0*Cons.delta*(Vars.S_ell*Vars.Sigma_n + Vars.S_n*Vars.Sigma_ell)/Cons.M**2 + 0.5*(2.0*Vars.S_n*(4.0*Vars.S_ell*(Cons.kappa_p + 2.0) + Vars.Sigma_ell*(2.0*Cons.delta*Cons.kappa_p + 7.0*Cons.kappa_m)) + Vars.Sigma_n*(2.0*Vars.S_ell*(2.0*Cons.delta*Cons.kappa_p + 7.0*Cons.kappa_m) + Vars.Sigma_ell*(5.0*Cons.delta*Cons.kappa_m - 17.0*Cons.kappa_p - 26.0)))/Cons.M**2) - 4.0*Vars.Sigma_ell*Vars.Sigma_n*Cons.nu**2*(Cons.kappa_p + 2.0)/Cons.M**2 + 6.0*Cons.delta*(Vars.S_ell*Vars.Sigma_n + Vars.S_n*Vars.Sigma_ell)/Cons.M**2 + 0.5*(Vars.S_n*(Vars.S_ell*(-9.0*Cons.delta*Cons.kappa_m + 3.0*Cons.kappa_p + 22.0) + Vars.Sigma_ell*(6.0*Cons.delta*Cons.kappa_p - 6.0*Cons.kappa_m)) + 6.0*Vars.Sigma_n*(Vars.S_ell*(Cons.delta*Cons.kappa_p - Cons.kappa_m) + Vars.Sigma_ell*(-Cons.delta*Cons.kappa_m + Cons.kappa_p + 1.0)))/Cons.M**2
    gamma_PN_3 = (1.66666666666667*Vars.S_ell + Vars.Sigma_ell*Cons.delta)/Cons.M**2
    a_ell_0 = 7.0*Vars.S_n + 3.0*Vars.Sigma_n*Cons.delta
    gamma_PN_4 = -5.41666666666667*Cons.nu + 1.0 + (Vars.S_ell**2*(-0.5*Cons.kappa_p - 1.0) + Vars.S_ell*Vars.Sigma_ell*(-0.5*Cons.delta*Cons.kappa_p - Cons.delta + 0.5*Cons.kappa_m) + Vars.Sigma_ell**2*(0.25*Cons.delta*Cons.kappa_m - 0.25*Cons.kappa_p + Cons.nu*(0.5*Cons.kappa_p + 1.0)))/Cons.M**8
    return Vars.ellHat*Vars.v**3/Cons.M + Vars.nHat*Vars.v**6*(a_ell_0 + Vars.v*(a_ell_1 + Vars.v*(a_ell_2 + Vars.v*(a_ell_3 + a_ell_4*Vars.v))))*(gamma_PN_0 + Vars.v**2*(gamma_PN_2 + Vars.v*(gamma_PN_3 + Vars.v*(gamma_PN_4 + gamma_PN_5*Vars.v))))/Cons.M**3


@njit(cache=True)
def TaylorT1_2p5(Cons,Vars):
    Flux = Vars.Fcal_coeff*(Cons.Fcal_0 + Vars.v**2*(Cons.Fcal_2 + Vars.v*(Cons.Fcal_3 + Vars.Fcal_SO_3 + Vars.v*(Cons.Fcal_4 + Vars.Fcal_SQ_4 + Vars.v*(Cons.Fcal_5 + Vars.Fcal_SO_5)))))
    dEdV = -0.5*Cons.M*Cons.nu*Vars.v*(2.0*Cons.E_0 + Vars.v**2*(4.0*Cons.E_2 + Vars.v*(5.0*Vars.E_SO_3 + Vars.v*(6.0*Cons.E_4 + 7.0*Vars.E_SO_5*Vars.v + 6.0*Vars.E_SQ_4))))
    Absorption = Vars.Fcal_coeff*Vars.MDot_Alvi_5*Vars.v**5
    dvdt_T1 = (-Absorption - Flux)/dEdV
    dydt=np.zeros(9)
    [dydt[5],dydt[6],dydt[7],dydt[8]] = FrameFromAngularVelocityIntegrand(Vars.R, OmegaVec_2p5(Cons,Vars)[1:])
    dydt[0] = dvdt_T1[0] 
    if(Cons.EvolveSpin1):
        dydt[1], dydt[2]=FrameFromAngularVelocity_2D_Integrand(Vars.rfrak_chi1[0], Vars.rfrak_chi1[1],(mul(mul(inverse(Cons.S_chi1),OmegaVec_chiVec_1_2p5(Cons,Vars)),Cons.S_chi1))[1:])
    else:
        dydt[1], dydt[2] = 0.0, 0.0
    if(Cons.EvolveSpin2):
        dydt[3], dydt[4]=FrameFromAngularVelocity_2D_Integrand(Vars.rfrak_chi2[0], Vars.rfrak_chi2[1],(mul(mul(inverse(Cons.S_chi2),OmegaVec_chiVec_2_2p5(Cons,Vars)),Cons.S_chi2))[1:])
    else:
        dydt[3], dydt[4] = 0.0, 0.0
    return dydt

@njit(cache=True)
def TaylorT4_2p5(Cons,Vars):
    Flux = Vars.Fcal_coeff*(Cons.Fcal_0 + Vars.v**2*(Cons.Fcal_2 + Vars.v*(Cons.Fcal_3 + Vars.Fcal_SO_3 + Vars.v*(Cons.Fcal_4 + Vars.Fcal_SQ_4 + Vars.v*(Cons.Fcal_5 + Vars.Fcal_SO_5)))))
    dEdV = -0.5*Cons.M*Cons.nu*Vars.v*(2.0*Cons.E_0 + Vars.v**2*(4.0*Cons.E_2 + Vars.v*(5.0*Vars.E_SO_3 + Vars.v*(6.0*Cons.E_4 + 7.0*Vars.E_SO_5*Vars.v + 6.0*Vars.E_SQ_4))))
    Absorption = Vars.Fcal_coeff*Vars.MDot_Alvi_5*Vars.v**5
    dvdt_T1 = (-Absorption - Flux)/dEdV
    dydt=np.zeros(9)
    [dydt[5],dydt[6],dydt[7],dydt[8]] = FrameFromAngularVelocityIntegrand(Vars.R, OmegaVec_2p5(Cons,Vars)[1:])
    dydt[0] = dvdt_T4[0]
    if(Cons.EvolveSpin1):
        dydt[1], dydt[2]=FrameFromAngularVelocity_2D_Integrand(Vars.rfrak_chi1[0], Vars.rfrak_chi1[1],(mul(mul(inverse(Cons.S_chi1),OmegaVec_chiVec_1_2p5(Cons,Vars)),Cons.S_chi1))[1:])
    else:
        dydt[1], dydt[2] = 0.0, 0.0
    if(Cons.EvolveSpin2):
        dydt[3], dydt[4]=FrameFromAngularVelocity_2D_Integrand(Vars.rfrak_chi2[0], Vars.rfrak_chi2[1],(mul(mul(inverse(Cons.S_chi2),OmegaVec_chiVec_2_2p5(Cons,Vars)),Cons.S_chi2))[1:])
    else:
        dydt[3], dydt[4] = 0.0, 0.0
    return dydt      

@njit(cache=True)
def TaylorT5_2p5(Cons,Vars):
    Flux = Vars.Fcal_coeff*(Cons.Fcal_0 + Vars.v**2*(Cons.Fcal_2 + Vars.v*(Cons.Fcal_3 + Vars.Fcal_SO_3 + Vars.v*(Cons.Fcal_4 + Vars.Fcal_SQ_4 + Vars.v*(Cons.Fcal_5 + Vars.Fcal_SO_5)))))
    dEdV = -0.5*Cons.M*Cons.nu*Vars.v*(2.0*Cons.E_0 + Vars.v**2*(4.0*Cons.E_2 + Vars.v*(5.0*Vars.E_SO_3 + Vars.v*(6.0*Cons.E_4 + 7.0*Vars.E_SO_5*Vars.v + 6.0*Vars.E_SQ_4))))
    Absorption = Vars.Fcal_coeff*Vars.MDot_Alvi_5*Vars.v**5
    dvdt_T1 = (-Absorption - Flux)/dEdV
    dydt=np.zeros(8)
    [dydt[5],dydt[6],dydt[7],dydt[8]] = FrameFromAngularVelocityIntegrand(Vars.R, OmegaVec_2p5(Cons,Vars)[1:])
    dydt[0] = dvdt_T5[0]
    if(Cons.EvolveSpin1):
        dydt[1], dydt[2]=FrameFromAngularVelocity_2D_Integrand(Vars.rfrak_chi1[0], Vars.rfrak_chi1[1],(mul(mul(inverse(Cons.S_chi1),OmegaVec_chiVec_1_2p5(Cons,Vars)),Cons.S_chi1))[1:])
    else:
        dydt[1], dydt[2] = 0.0, 0.0
    if(Cons.EvolveSpin2):
        dydt[3], dydt[4]=FrameFromAngularVelocity_2D_Integrand(Vars.rfrak_chi2[0], Vars.rfrak_chi2[1],(mul(mul(inverse(Cons.S_chi2),OmegaVec_chiVec_2_2p5(Cons,Vars)),Cons.S_chi2))[1:])
    else:
        dydt[3], dydt[4] = 0.0, 0.0
    return dydt      

@njit(cache=True)
def Recalculate_3p0(Cons,Vars,y):
    Vars.v = np.array([y[0]])
    Vars.rfrak_chi1 = np.array([y[1],y[2]])
    Vars.rfrak_chi2 = np.array([y[3],y[4]])
    Vars.R = y[5]*Cons.wHat + y[6]*Cons.xHat + y[7]*Cons.yHat + y[8]*Cons.zHat
    Vars.nHat = mul(mul(Vars.R,Cons.xHat),inverse(Vars.R))
    Vars.lambdaHat = mul(mul(Vars.R,Cons.yHat),inverse(Vars.R))
    Vars.ellHat = mul(mul(Vars.R,Cons.zHat),inverse(Vars.R))
    Vars.R_S1 = exp(Vars.rfrak_chi1[0]*Cons.xHat + Vars.rfrak_chi1[1]*Cons.yHat)
    Vars.R_S2 = exp(Vars.rfrak_chi2[0]*Cons.xHat + Vars.rfrak_chi2[1]*Cons.yHat)
    Vars.chiVec1 = mul(mul(mul(Cons.S_chi1,Vars.R_S1),Cons.zHat),mul(conjugate(Vars.R_S1),conjugate(Cons.S_chi1)))
    Vars.chiVec2 = mul(mul(mul(Cons.S_chi2,Vars.R_S2),Cons.zHat),mul(conjugate(Vars.R_S2),conjugate(Cons.S_chi2)))
    Vars.chi1_n = np.array([dot(Vars.chiVec1[1:],Vars.nHat[1:])])
    Vars.chi1_lambda = np.array([dot(Vars.chiVec1[1:],Vars.lambdaHat[1:])])
    Vars.chi1_ell = np.array([dot(Vars.chiVec1[1:],Vars.ellHat[1:])])
    Vars.chi2_n = np.array([dot(Vars.chiVec2[1:],Vars.nHat[1:])])
    Vars.chi2_lambda = np.array([dot(Vars.chiVec2[1:],Vars.lambdaHat[1:])])
    Vars.chi2_ell = np.array([dot(Vars.chiVec2[1:],Vars.ellHat[1:])])
    Vars.S = Cons.M1**2*Vars.chiVec1 + Cons.M2**2*Vars.chiVec2
    Vars.S_ell = Cons.M1**2*Vars.chi1_ell + Cons.M2**2*Vars.chi2_ell
    Vars.S_n = Cons.M1**2*Vars.chi1_n + Cons.M2**2*Vars.chi2_n
    Vars.S_lambda = Cons.M1**2*Vars.chi1_lambda + Cons.M2**2*Vars.chi2_lambda
    Vars.Sigma = Cons.M*(-Cons.M1*Vars.chiVec1 + Cons.M2*Vars.chiVec2)
    Vars.Sigma_ell = Cons.M*(-Cons.M1*Vars.chi1_ell + Cons.M2*Vars.chi2_ell)
    Vars.Sigma_n = Cons.M*(-Cons.M1*Vars.chi1_n + Cons.M2*Vars.chi2_n)
    Vars.Sigma_lambda = Cons.M*(-Cons.M1*Vars.chi1_lambda + Cons.M2*Vars.chi2_lambda)
    Vars.SS = np.array([dot(Vars.S[1:],Vars.S[1:])])
    Vars.SigmaSigma = np.array([dot(Vars.Sigma[1:],Vars.Sigma[1:])])
    Vars.SSigma = np.array([dot(Vars.S[1:],Vars.Sigma[1:])])
    Vars.S1 = Cons.M1**2*Vars.chiVec1
    Vars.S1_n = Cons.M1**2*Vars.chi1_n
    Vars.S1_lambda = Cons.M1**2*Vars.chi1_lambda
    Vars.S2 = Cons.M2**2*Vars.chiVec2
    Vars.S2_n = Cons.M2**2*Vars.chi2_n
    Vars.S2_lambda = Cons.M2**2*Vars.chi2_lambda
    Vars.chi_s_ell = Vars.chi1_ell/2 + Vars.chi2_ell/2
    Vars.chi_a_ell = Vars.chi1_ell/2 - Vars.chi2_ell/2
    Vars.logv = log(Vars.v)
    Vars.Fcal_coeff = 32*Cons.nu**2*Vars.v**10/5
    Vars.Fcal_SQ_4 = Cons.chi1chi1*(-89*Cons.delta/192 + 89*Cons.nu/96 - 89/192) - 103*Cons.chi1chi2*Cons.nu/48 + Cons.chi2chi2*(89*Cons.delta/192 + 89*Cons.nu/96 - 89/192) + Vars.chi_a_ell*(Vars.chi_a_ell*(287/96 - 12*Cons.nu) + 287*Vars.chi_s_ell*Cons.delta/48) + Vars.chi_s_ell**2*(Cons.nu/24 + 287/96)
    Vars.Fcal_SQ_6 = Cons.nu*(-0.0714285714285714*Cons.delta*(353*Vars.SSigma + 86*Vars.S_lambda*Vars.Sigma_lambda - 1145*Vars.S_n*Vars.Sigma_n)/Cons.M**4 - 0.00297619047619048*(4236*Vars.SS*Cons.kappa_p + 8472*Vars.SS + 4236*Vars.SSigma*Cons.delta*Cons.kappa_p + 1672*Vars.SSigma*Cons.kappa_m + 1032*Vars.S_lambda**2*(Cons.kappa_p + 2) + 8*Vars.S_lambda*Vars.Sigma_lambda*(129*Cons.delta*Cons.kappa_p + 620*Cons.kappa_m) - 13740*Vars.S_n**2*(Cons.kappa_p + 2) - 4*Vars.S_n*Vars.Sigma_n*(3435*Cons.delta*Cons.kappa_p + 2494*Cons.kappa_m) - 641*Vars.SigmaSigma*Cons.delta*Cons.kappa_m - 2713*Vars.SigmaSigma*Cons.kappa_p - 4146*Vars.SigmaSigma + 982*Vars.Sigma_lambda**2*Cons.delta*Cons.kappa_m - 6944*Vars.Sigma_lambda**2*Cons.kappa_p - 3736*Vars.Sigma_lambda**2 + 941*Vars.Sigma_n**2*Cons.delta*Cons.kappa_m + 15083*Vars.Sigma_n**2*Cons.kappa_p + 44524*Vars.Sigma_n**2)/Cons.M**4) - 0.00595238095238095*Cons.delta*(911*Vars.SSigma + 3643*Vars.S_lambda*Vars.Sigma_lambda - 22399*Vars.S_n*Vars.Sigma_n)/Cons.M**4 + 0.0357142857142857*Cons.nu**2*(Cons.kappa_p + 2)*(353*Vars.SigmaSigma + 86*Vars.Sigma_lambda**2 - 1145*Vars.Sigma_n**2)/Cons.M**4 + 0.00297619047619048*(1477*Vars.SS*Cons.delta*Cons.kappa_m - 1877*Vars.SS*Cons.kappa_p + 502*Vars.SS - 3354*Vars.SSigma*Cons.delta*Cons.kappa_p + 3354*Vars.SSigma*Cons.kappa_m + Vars.S_lambda**2*(1498*Cons.delta*Cons.kappa_m - 4464*Cons.kappa_p - 10811) - 5962*Vars.S_lambda*Vars.Sigma_lambda*(Cons.delta*Cons.kappa_p - Cons.kappa_m) + Vars.S_n**2*(-5929*Cons.delta*Cons.kappa_m + 10095*Cons.kappa_p + 45299) + 16024*Vars.S_n*Vars.Sigma_n*(Cons.delta*Cons.kappa_p - Cons.kappa_m) + 1677*Vars.SigmaSigma*Cons.delta*Cons.kappa_m - 1677*Vars.SigmaSigma*Cons.kappa_p - 168*Vars.SigmaSigma + 2981*Vars.Sigma_lambda**2*Cons.delta*Cons.kappa_m - 2981*Vars.Sigma_lambda**2*Cons.kappa_p + 1313*Vars.Sigma_lambda**2 - 8012*Vars.Sigma_n**2*Cons.delta*Cons.kappa_m + 8012*Vars.Sigma_n**2*Cons.kappa_p + 6387*Vars.Sigma_n**2)/Cons.M**4
    Vars.Fcal_SO_3 = (-4*Vars.S_ell - 5*Vars.Sigma_ell*Cons.delta/4)/Cons.M**2
    Vars.Fcal_SO_5 = (Vars.S_ell*(272*Cons.nu/9 - 9/2) + Vars.Sigma_ell*Cons.delta*(43*Cons.nu/4 - 13/16))/Cons.M**2
    Vars.Fcal_SO_6 = (-16*Vars.S_ell*pi - 31*Vars.Sigma_ell*Cons.delta*pi/6)/Cons.M**2
    Vars.E_SQ_4 = -3*Vars.chi_a_ell**2/2 - 3*Vars.chi_s_ell**2/2 - Cons.delta*(Cons.chi2chi2/2 + 3*Vars.chi_a_ell*Vars.chi_s_ell) + Cons.nu*(Cons.chi1chi2 + 6*Vars.chi_a_ell**2) + (Cons.chi1chi1 + Cons.chi2chi2)*(Cons.delta - 2*Cons.nu + 1)/4
    Vars.E_SQ_6 = Cons.nu*(-3*Cons.delta*(Vars.SSigma - 3*Vars.S_n*Vars.Sigma_n)/Cons.M**4 + 0.25*(-6*Vars.SS*Cons.kappa_p - 12*Vars.SS - 6*Vars.SSigma*Cons.delta*Cons.kappa_p - 26*Vars.SSigma*Cons.kappa_m + 24*Vars.S_lambda*Vars.Sigma_lambda*Cons.kappa_m + 18*Vars.S_n**2*(Cons.kappa_p + 2) + 18*Vars.S_n*Vars.Sigma_n*(Cons.delta*Cons.kappa_p + 3*Cons.kappa_m) + 19*Vars.SigmaSigma*Cons.kappa_p - 80*Vars.SigmaSigma - 18*Vars.Sigma_lambda**2*Cons.kappa_p + 76*Vars.Sigma_lambda**2 - 39*Vars.Sigma_n**2*Cons.kappa_p + 168*Vars.Sigma_n**2 + Cons.delta*Cons.kappa_m*(-5*Vars.SigmaSigma + 6*Vars.Sigma_lambda**2 + 9*Vars.Sigma_n**2))/Cons.M**4) + Cons.delta*(15*Vars.SSigma - 13*Vars.S_lambda*Vars.Sigma_lambda - 29*Vars.S_n*Vars.Sigma_n)/Cons.M**4 + 1.5*Cons.nu**2*(Vars.SigmaSigma - 3*Vars.Sigma_n**2)*(Cons.kappa_p + 2)/Cons.M**4 + 0.25*(8*Vars.SS*Cons.delta*Cons.kappa_m - 6*Vars.SS*Cons.kappa_p + 36*Vars.SS - 14*Vars.SSigma*Cons.delta*Cons.kappa_p + 14*Vars.SSigma*Cons.kappa_m + Vars.S_lambda**2*(-6*Cons.delta*Cons.kappa_m + 6*Cons.kappa_p - 28) + 12*Vars.S_lambda*Vars.Sigma_lambda*(Cons.delta*Cons.kappa_p - Cons.kappa_m) - 6*Vars.S_n**2*(3*Cons.delta*Cons.kappa_m - 2*Cons.kappa_p + 10) + 30*Vars.S_n*Vars.Sigma_n*(Cons.delta*Cons.kappa_p - Cons.kappa_m) + 7*Vars.SigmaSigma*Cons.delta*Cons.kappa_m - 7*Vars.SigmaSigma*Cons.kappa_p + 24*Vars.SigmaSigma - 6*Vars.Sigma_lambda**2*Cons.delta*Cons.kappa_m + 6*Vars.Sigma_lambda**2*Cons.kappa_p - 24*Vars.Sigma_lambda**2 - 15*Vars.Sigma_n**2*Cons.delta*Cons.kappa_m + 15*Vars.Sigma_n**2*Cons.kappa_p - 48*Vars.Sigma_n**2)/Cons.M**4
    Vars.E_SO_3 = (14*Vars.S_ell/3 + 2*Vars.Sigma_ell*Cons.delta)/Cons.M**2
    Vars.E_SO_5 = (Vars.S_ell*(11 - 61*Cons.nu/9) + Vars.Sigma_ell*Cons.delta*(3 - 10*Cons.nu/3))/Cons.M**2
    Vars.MDot_Alvi_5 = -Cons.M1**3*Vars.chi1_ell*(3*Cons.chi1chi1 + 1)/(4*Cons.M**3) - Cons.M2**3*Vars.chi2_ell*(3*Cons.chi2chi2 + 1)/(4*Cons.M**3)

@njit
def OmegaVec_chiVec_1_3p0(Cons,Vars):
    Omega1_coeff = Vars.v**5/Cons.M
    return Omega1_coeff*(Vars.S2*Vars.v**3*((2.0 - 0.5*Cons.delta)/Cons.M**2 + Vars.v*(-5.0*Vars.S_ell - 3.0*Vars.Sigma_ell*Cons.delta/Cons.M)/Cons.M**4) + Vars.ellHat*(-0.75*Cons.delta + 0.5*Cons.nu + Vars.v**2*(Cons.delta*(0.625*Cons.nu - 0.5625) + Cons.nu*(1.25 - 0.0416666666666667*Cons.nu) + Vars.v*(Vars.v*(Cons.nu*(Cons.nu*(0.229166666666667*Cons.nu - 4.94791666666667) + 2.89583333333333) + 0.96875 + (Cons.delta*(Cons.nu*(5.91666666666667 - 0.489583333333333*Cons.nu) - 0.96875) + (Vars.S_ell*(Vars.S_ell*(Cons.kappa_p*(-0.25*Cons.nu - 0.375) - 0.5*Cons.nu - 0.75) + Vars.Sigma_ell*(Cons.delta*(Cons.kappa_p*(-0.25*Cons.nu - 0.375) - 0.5*Cons.nu - 0.75) + Cons.kappa_m*(0.25*Cons.nu + 0.375))) + Vars.Sigma_ell**2*(Cons.delta*Cons.kappa_m*(0.125*Cons.nu + 0.1875) + Cons.kappa_p*(Cons.nu*(0.25*Cons.nu + 0.25) - 0.1875) + Cons.nu*(0.5*Cons.nu + 0.75)) + (Vars.S_ell*(Vars.S_ell*Cons.delta*(0.375*Cons.kappa_p + 0.75) + Vars.Sigma_ell*Cons.delta*(Cons.delta*(0.375*Cons.kappa_p + 0.75) - 0.375*Cons.kappa_m)) + Vars.Sigma_ell**2*Cons.delta*(-0.1875*Cons.delta*Cons.kappa_m + Cons.kappa_p*(0.1875 - 0.375*Cons.nu) - 0.75*Cons.nu))/Cons.M)/Cons.M**7)/Cons.M) + (Vars.S_ell*(0.833333333333333*Cons.nu + 1.25) + (-1.25*Vars.S_ell*Cons.delta + Vars.Sigma_ell*Cons.delta*(0.5*Cons.nu + 0.75) - 0.75*Vars.Sigma_ell*Cons.delta**2/Cons.M)/Cons.M)/Cons.M**2) + 0.5625) + 0.75) + Vars.lambdaHat*Vars.v**3*(Vars.S2_lambda*(0.5*Cons.delta + 2.0)/Cons.M**2 + Vars.v*(Vars.S1_lambda*(-0.75*Cons.delta - 3.0*Cons.kappa_1 + 2.25) + Vars.S1_lambda*(Cons.delta*(0.75 - 1.5*Cons.kappa_1) + 1.5*Cons.kappa_1 - 0.75)/Cons.nu)/Cons.M**2) + Vars.nHat*Vars.v*(3.0*Vars.chi1_n*Cons.nu + Vars.v**2*((Vars.S1_n*(Cons.delta*(-2.25*Cons.kappa_1 - 0.25) + Cons.kappa_1*(1.5*Cons.nu + 0.75) + 0.5*Cons.nu + 3.25) + Vars.S1_n*(Cons.delta*(1.5 - 0.75*Cons.kappa_1) + 0.75*Cons.kappa_1 - 1.5)/Cons.nu + Vars.S2_n*(1.5*Cons.delta - 1.0*Cons.nu + 3.0))/Cons.M**2 + Vars.v*(-15.0*Vars.S1_n*Vars.S_ell*Cons.kappa_1 + Vars.S1_n*Vars.S_ell*(-7.5*Cons.delta*Cons.kappa_1 + 7.5*Cons.kappa_1)/Cons.nu + 15.0*Vars.S2_n*Vars.S_ell + (-9.0*Vars.S1_n*Vars.Sigma_ell*Cons.delta*Cons.kappa_1 + Vars.S1_n*Vars.Sigma_ell*Cons.delta*(-4.5*Cons.delta*Cons.kappa_1 + 4.5*Cons.kappa_1)/Cons.nu + 9.0*Vars.S2_n*Vars.Sigma_ell*Cons.delta)/Cons.M)/Cons.M**4) + 3.0*Cons.M2**2*Vars.chi2_n/Cons.M**2) - Cons.M2**2*Vars.chiVec2*Vars.v/Cons.M**2)

@njit
def OmegaVec_chiVec_2_3p0(Cons,Vars):
    Omega2_coeff = Vars.v**5/Cons.M
    return Omega2_coeff*(Vars.S1*Vars.v**3*((0.5*Cons.delta + 2.0)/Cons.M**2 + Vars.v*(-5.0*Vars.S_ell + 3.0*Vars.Sigma_ell*Cons.delta/Cons.M)/Cons.M**4) + Vars.ellHat*(0.75*Cons.delta + 0.5*Cons.nu + Vars.v**2*(Cons.delta*(0.5625 - 0.625*Cons.nu) + Cons.nu*(1.25 - 0.0416666666666667*Cons.nu) + Vars.v*(Vars.v*(Cons.nu*(Cons.nu*(0.229166666666667*Cons.nu - 4.94791666666667) + 2.89583333333333) + 0.96875 + (Cons.delta*(Cons.nu*(0.489583333333333*Cons.nu - 5.91666666666667) + 0.96875) + (Vars.S_ell*(Vars.S_ell*(Cons.kappa_p*(-0.25*Cons.nu - 0.375) - 0.5*Cons.nu - 0.75) + Vars.Sigma_ell*(Cons.delta*(Cons.kappa_p*(0.25*Cons.nu + 0.375) + 0.5*Cons.nu + 0.75) + Cons.kappa_m*(-0.25*Cons.nu - 0.375))) + Vars.Sigma_ell**2*(Cons.delta*Cons.kappa_m*(0.125*Cons.nu + 0.1875) + Cons.kappa_p*(Cons.nu*(0.25*Cons.nu + 0.25) - 0.1875) + Cons.nu*(0.5*Cons.nu + 0.75)) + (Vars.S_ell*(Vars.S_ell*Cons.delta*(-0.375*Cons.kappa_p - 0.75) + Vars.Sigma_ell*Cons.delta*(Cons.delta*(0.375*Cons.kappa_p + 0.75) - 0.375*Cons.kappa_m)) + Vars.Sigma_ell**2*Cons.delta*(0.1875*Cons.delta*Cons.kappa_m + Cons.kappa_p*(0.375*Cons.nu - 0.1875) + 0.75*Cons.nu))/Cons.M)/Cons.M**7)/Cons.M) + (Vars.S_ell*(0.833333333333333*Cons.nu + 1.25) + (1.25*Vars.S_ell*Cons.delta + Vars.Sigma_ell*Cons.delta*(-0.5*Cons.nu - 0.75) - 0.75*Vars.Sigma_ell*Cons.delta**2/Cons.M)/Cons.M)/Cons.M**2) + 0.5625) + 0.75) + Vars.lambdaHat*Vars.v**3*(Vars.S1_lambda*(2.0 - 0.5*Cons.delta)/Cons.M**2 + Vars.v*(Vars.S2_lambda*(0.75*Cons.delta - 3.0*Cons.kappa_2 + 2.25) + Vars.S2_lambda*(Cons.delta*(1.5*Cons.kappa_2 - 0.75) + 1.5*Cons.kappa_2 - 0.75)/Cons.nu)/Cons.M**2) + Vars.nHat*Vars.v*(3.0*Vars.chi2_n*Cons.nu + Vars.v**2*((Vars.S1_n*(-1.5*Cons.delta - 1.0*Cons.nu + 3.0) + Vars.S2_n*(Cons.delta*(2.25*Cons.kappa_2 + 0.25) + Cons.kappa_2*(1.5*Cons.nu + 0.75) + 0.5*Cons.nu + 3.25) + Vars.S2_n*(Cons.delta*(0.75*Cons.kappa_2 - 1.5) + 0.75*Cons.kappa_2 - 1.5)/Cons.nu)/Cons.M**2 + Vars.v*(15.0*Vars.S1_n*Vars.S_ell - 15.0*Vars.S2_n*Vars.S_ell*Cons.kappa_2 + Vars.S2_n*Vars.S_ell*(7.5*Cons.delta*Cons.kappa_2 + 7.5*Cons.kappa_2)/Cons.nu + (-9.0*Vars.S1_n*Vars.Sigma_ell*Cons.delta + 9.0*Vars.S2_n*Vars.Sigma_ell*Cons.delta*Cons.kappa_2 + Vars.S2_n*Vars.Sigma_ell*Cons.delta*(-4.5*Cons.delta*Cons.kappa_2 - 4.5*Cons.kappa_2)/Cons.nu)/Cons.M)/Cons.M**4) + 3.0*Cons.M1**2*Vars.chi1_n/Cons.M**2) - Cons.M1**2*Vars.chiVec1*Vars.v/Cons.M**2)

@njit
def OmegaVec_3p0(Cons,Vars):
    gamma_PN_5 = (0.888888888888889*Vars.S_ell*Cons.nu + 3.33333333333333*Vars.S_ell)/Cons.M**2 + 2.0*Vars.Sigma_ell*Cons.delta/Cons.M**3
    gamma_PN_6 = Vars.S_ell**2*(-0.916666666666667*Cons.delta*Cons.kappa_m - 0.916666666666667*Cons.kappa_p + Cons.nu*(-0.166666666666667*Cons.kappa_p - 0.333333333333333) + 1.55555555555556) + Vars.S_ell*Vars.Sigma_ell*(1.66666666666667*Cons.delta + Cons.nu*(-0.166666666666667*Cons.delta*Cons.kappa_p - 0.333333333333333*Cons.delta + 3.83333333333333*Cons.kappa_m)) + Vars.Sigma_ell**2*(Cons.nu**2*(0.166666666666667*Cons.kappa_p + 0.333333333333333) + Cons.nu*(Cons.delta*Cons.kappa_m - Cons.kappa_p - 2.0) + 1.0) + 0.0123456790123457*Cons.nu**3 + 6.36111111111111*Cons.nu**2 - 2.98177812235564*Cons.nu + 1.0
    gamma_PN_0 = 1.00000000000000
    a_ell_4 = Vars.S_n*(5.77777777777778*Cons.nu**2 + 14.75*Cons.nu + 1.5) + Vars.Sigma_n*Cons.delta*(2.83333333333333*Cons.nu**2 + 9.125*Cons.nu + 1.5)
    a_ell_1 = 3.0*Vars.Sigma_ell*Vars.Sigma_n*Cons.nu*(Cons.kappa_p + 2.0)/Cons.M**2 - 3.0*Cons.delta*(Vars.S_ell*Vars.Sigma_n + Vars.S_n*Vars.Sigma_ell)/Cons.M**2 - 0.5*(3.0*Vars.S_n*(2.0*Vars.S_ell*(Cons.kappa_p + 2.0) + Vars.Sigma_ell*(Cons.delta*Cons.kappa_p - Cons.kappa_m)) + 3.0*Vars.Sigma_n*(Vars.S_ell*(Cons.delta*Cons.kappa_p - Cons.kappa_m) + Vars.Sigma_ell*(-Cons.delta*Cons.kappa_m + Cons.kappa_p)))/Cons.M**2
    a_ell_2 = Vars.S_n*(-9.66666666666667*Cons.nu - 10.0) + Vars.Sigma_n*Cons.delta*(-4.5*Cons.nu - 6.0)
    gamma_PN_2 = 1.0 - 0.333333333333333*Cons.nu
    a_ell_3 = Cons.nu*(4.0*Cons.delta*(Vars.S_ell*Vars.Sigma_n + Vars.S_n*Vars.Sigma_ell)/Cons.M**2 + 0.5*(2.0*Vars.S_n*(4.0*Vars.S_ell*(Cons.kappa_p + 2.0) + Vars.Sigma_ell*(2.0*Cons.delta*Cons.kappa_p + 7.0*Cons.kappa_m)) + Vars.Sigma_n*(2.0*Vars.S_ell*(2.0*Cons.delta*Cons.kappa_p + 7.0*Cons.kappa_m) + Vars.Sigma_ell*(5.0*Cons.delta*Cons.kappa_m - 17.0*Cons.kappa_p - 26.0)))/Cons.M**2) - 4.0*Vars.Sigma_ell*Vars.Sigma_n*Cons.nu**2*(Cons.kappa_p + 2.0)/Cons.M**2 + 6.0*Cons.delta*(Vars.S_ell*Vars.Sigma_n + Vars.S_n*Vars.Sigma_ell)/Cons.M**2 + 0.5*(Vars.S_n*(Vars.S_ell*(-9.0*Cons.delta*Cons.kappa_m + 3.0*Cons.kappa_p + 22.0) + Vars.Sigma_ell*(6.0*Cons.delta*Cons.kappa_p - 6.0*Cons.kappa_m)) + 6.0*Vars.Sigma_n*(Vars.S_ell*(Cons.delta*Cons.kappa_p - Cons.kappa_m) + Vars.Sigma_ell*(-Cons.delta*Cons.kappa_m + Cons.kappa_p + 1.0)))/Cons.M**2
    gamma_PN_3 = (1.66666666666667*Vars.S_ell + Vars.Sigma_ell*Cons.delta)/Cons.M**2
    a_ell_0 = 7.0*Vars.S_n + 3.0*Vars.Sigma_n*Cons.delta
    gamma_PN_4 = -5.41666666666667*Cons.nu + 1.0 + (Vars.S_ell**2*(-0.5*Cons.kappa_p - 1.0) + Vars.S_ell*Vars.Sigma_ell*(-0.5*Cons.delta*Cons.kappa_p - Cons.delta + 0.5*Cons.kappa_m) + Vars.Sigma_ell**2*(0.25*Cons.delta*Cons.kappa_m - 0.25*Cons.kappa_p + Cons.nu*(0.5*Cons.kappa_p + 1.0)))/Cons.M**8
    return Vars.ellHat*Vars.v**3/Cons.M + Vars.nHat*Vars.v**6*(a_ell_0 + Vars.v*(a_ell_1 + Vars.v*(a_ell_2 + Vars.v*(a_ell_3 + a_ell_4*Vars.v))))*(gamma_PN_0 + Vars.v**2*(gamma_PN_2 + Vars.v*(gamma_PN_3 + Vars.v*(gamma_PN_4 + Vars.v*(gamma_PN_5 + gamma_PN_6*Vars.v)))))/Cons.M**3


@njit(cache=True)
def TaylorT1_3p0(Cons,Vars):
    Flux = Vars.Fcal_coeff*(Cons.Fcal_0 + Vars.v**2*(Cons.Fcal_2 + Vars.v*(Cons.Fcal_3 + Vars.Fcal_SO_3 + Vars.v*(Cons.Fcal_4 + Vars.Fcal_SQ_4 + Vars.v*(Cons.Fcal_5 + Vars.Fcal_SO_5 + Vars.v*(Cons.Fcal_6 + Vars.Fcal_SO_6 + Vars.Fcal_SQ_6 + Cons.Fcal_lnv_6*Vars.logv))))))
    dEdV = -0.5*Cons.M*Cons.nu*Vars.v*(2.0*Cons.E_0 + Vars.v**2*(4.0*Cons.E_2 + Vars.v*(5.0*Vars.E_SO_3 + Vars.v*(6.0*Cons.E_4 + 6.0*Vars.E_SQ_4 + Vars.v*(7.0*Vars.E_SO_5 + 8.0*Vars.v*(Cons.E_6 + Vars.E_SQ_6))))))
    Absorption = Vars.Fcal_coeff*Vars.MDot_Alvi_5*Vars.v**5
    dvdt_T1 = (-Absorption - Flux)/dEdV
    dydt=np.zeros(9)
    [dydt[5],dydt[6],dydt[7],dydt[8]] = FrameFromAngularVelocityIntegrand(Vars.R, OmegaVec_3p0(Cons,Vars)[1:])
    dydt[0] = dvdt_T1[0] 
    if(Cons.EvolveSpin1):
        dydt[1], dydt[2]=FrameFromAngularVelocity_2D_Integrand(Vars.rfrak_chi1[0], Vars.rfrak_chi1[1],(mul(mul(inverse(Cons.S_chi1),OmegaVec_chiVec_1_3p0(Cons,Vars)),Cons.S_chi1))[1:])
    else:
        dydt[1], dydt[2] = 0.0, 0.0
    if(Cons.EvolveSpin2):
        dydt[3], dydt[4]=FrameFromAngularVelocity_2D_Integrand(Vars.rfrak_chi2[0], Vars.rfrak_chi2[1],(mul(mul(inverse(Cons.S_chi2),OmegaVec_chiVec_2_3p0(Cons,Vars)),Cons.S_chi2))[1:])
    else:
        dydt[3], dydt[4] = 0.0, 0.0
    return dydt

@njit(cache=True)
def TaylorT4_3p0(Cons,Vars):
    Flux = Vars.Fcal_coeff*(Cons.Fcal_0 + Vars.v**2*(Cons.Fcal_2 + Vars.v*(Cons.Fcal_3 + Vars.Fcal_SO_3 + Vars.v*(Cons.Fcal_4 + Vars.Fcal_SQ_4 + Vars.v*(Cons.Fcal_5 + Vars.Fcal_SO_5 + Vars.v*(Cons.Fcal_6 + Vars.Fcal_SO_6 + Vars.Fcal_SQ_6 + Cons.Fcal_lnv_6*Vars.logv))))))
    dEdV = -0.5*Cons.M*Cons.nu*Vars.v*(2.0*Cons.E_0 + Vars.v**2*(4.0*Cons.E_2 + Vars.v*(5.0*Vars.E_SO_3 + Vars.v*(6.0*Cons.E_4 + 6.0*Vars.E_SQ_4 + Vars.v*(7.0*Vars.E_SO_5 + 8.0*Vars.v*(Cons.E_6 + Vars.E_SQ_6))))))
    Absorption = Vars.Fcal_coeff*Vars.MDot_Alvi_5*Vars.v**5
    dvdt_T1 = (-Absorption - Flux)/dEdV
    dydt=np.zeros(9)
    [dydt[5],dydt[6],dydt[7],dydt[8]] = FrameFromAngularVelocityIntegrand(Vars.R, OmegaVec_3p0(Cons,Vars)[1:])
    dydt[0] = dvdt_T4[0]
    if(Cons.EvolveSpin1):
        dydt[1], dydt[2]=FrameFromAngularVelocity_2D_Integrand(Vars.rfrak_chi1[0], Vars.rfrak_chi1[1],(mul(mul(inverse(Cons.S_chi1),OmegaVec_chiVec_1_3p0(Cons,Vars)),Cons.S_chi1))[1:])
    else:
        dydt[1], dydt[2] = 0.0, 0.0
    if(Cons.EvolveSpin2):
        dydt[3], dydt[4]=FrameFromAngularVelocity_2D_Integrand(Vars.rfrak_chi2[0], Vars.rfrak_chi2[1],(mul(mul(inverse(Cons.S_chi2),OmegaVec_chiVec_2_3p0(Cons,Vars)),Cons.S_chi2))[1:])
    else:
        dydt[3], dydt[4] = 0.0, 0.0
    return dydt      

@njit(cache=True)
def TaylorT5_3p0(Cons,Vars):
    Flux = Vars.Fcal_coeff*(Cons.Fcal_0 + Vars.v**2*(Cons.Fcal_2 + Vars.v*(Cons.Fcal_3 + Vars.Fcal_SO_3 + Vars.v*(Cons.Fcal_4 + Vars.Fcal_SQ_4 + Vars.v*(Cons.Fcal_5 + Vars.Fcal_SO_5 + Vars.v*(Cons.Fcal_6 + Vars.Fcal_SO_6 + Vars.Fcal_SQ_6 + Cons.Fcal_lnv_6*Vars.logv))))))
    dEdV = -0.5*Cons.M*Cons.nu*Vars.v*(2.0*Cons.E_0 + Vars.v**2*(4.0*Cons.E_2 + Vars.v*(5.0*Vars.E_SO_3 + Vars.v*(6.0*Cons.E_4 + 6.0*Vars.E_SQ_4 + Vars.v*(7.0*Vars.E_SO_5 + 8.0*Vars.v*(Cons.E_6 + Vars.E_SQ_6))))))
    Absorption = Vars.Fcal_coeff*Vars.MDot_Alvi_5*Vars.v**5
    dvdt_T1 = (-Absorption - Flux)/dEdV
    dydt=np.zeros(8)
    [dydt[5],dydt[6],dydt[7],dydt[8]] = FrameFromAngularVelocityIntegrand(Vars.R, OmegaVec_3p0(Cons,Vars)[1:])
    dydt[0] = dvdt_T5[0]
    if(Cons.EvolveSpin1):
        dydt[1], dydt[2]=FrameFromAngularVelocity_2D_Integrand(Vars.rfrak_chi1[0], Vars.rfrak_chi1[1],(mul(mul(inverse(Cons.S_chi1),OmegaVec_chiVec_1_3p0(Cons,Vars)),Cons.S_chi1))[1:])
    else:
        dydt[1], dydt[2] = 0.0, 0.0
    if(Cons.EvolveSpin2):
        dydt[3], dydt[4]=FrameFromAngularVelocity_2D_Integrand(Vars.rfrak_chi2[0], Vars.rfrak_chi2[1],(mul(mul(inverse(Cons.S_chi2),OmegaVec_chiVec_2_3p0(Cons,Vars)),Cons.S_chi2))[1:])
    else:
        dydt[3], dydt[4] = 0.0, 0.0
    return dydt      

@njit(cache=True)
def Recalculate_3p5(Cons,Vars,y):
    Vars.v = np.array([y[0]])
    Vars.rfrak_chi1 = np.array([y[1],y[2]])
    Vars.rfrak_chi2 = np.array([y[3],y[4]])
    Vars.R = y[5]*Cons.wHat + y[6]*Cons.xHat + y[7]*Cons.yHat + y[8]*Cons.zHat
    Vars.nHat = mul(mul(Vars.R,Cons.xHat),inverse(Vars.R))
    Vars.lambdaHat = mul(mul(Vars.R,Cons.yHat),inverse(Vars.R))
    Vars.ellHat = mul(mul(Vars.R,Cons.zHat),inverse(Vars.R))
    Vars.R_S1 = exp(Vars.rfrak_chi1[0]*Cons.xHat + Vars.rfrak_chi1[1]*Cons.yHat)
    Vars.R_S2 = exp(Vars.rfrak_chi2[0]*Cons.xHat + Vars.rfrak_chi2[1]*Cons.yHat)
    Vars.chiVec1 = mul(mul(mul(Cons.S_chi1,Vars.R_S1),Cons.zHat),mul(conjugate(Vars.R_S1),conjugate(Cons.S_chi1)))
    Vars.chiVec2 = mul(mul(mul(Cons.S_chi2,Vars.R_S2),Cons.zHat),mul(conjugate(Vars.R_S2),conjugate(Cons.S_chi2)))
    Vars.chi1_n = np.array([dot(Vars.chiVec1[1:],Vars.nHat[1:])])
    Vars.chi1_lambda = np.array([dot(Vars.chiVec1[1:],Vars.lambdaHat[1:])])
    Vars.chi1_ell = np.array([dot(Vars.chiVec1[1:],Vars.ellHat[1:])])
    Vars.chi2_n = np.array([dot(Vars.chiVec2[1:],Vars.nHat[1:])])
    Vars.chi2_lambda = np.array([dot(Vars.chiVec2[1:],Vars.lambdaHat[1:])])
    Vars.chi2_ell = np.array([dot(Vars.chiVec2[1:],Vars.ellHat[1:])])
    Vars.S = Cons.M1**2*Vars.chiVec1 + Cons.M2**2*Vars.chiVec2
    Vars.S_ell = Cons.M1**2*Vars.chi1_ell + Cons.M2**2*Vars.chi2_ell
    Vars.S_n = Cons.M1**2*Vars.chi1_n + Cons.M2**2*Vars.chi2_n
    Vars.S_lambda = Cons.M1**2*Vars.chi1_lambda + Cons.M2**2*Vars.chi2_lambda
    Vars.Sigma = Cons.M*(-Cons.M1*Vars.chiVec1 + Cons.M2*Vars.chiVec2)
    Vars.Sigma_ell = Cons.M*(-Cons.M1*Vars.chi1_ell + Cons.M2*Vars.chi2_ell)
    Vars.Sigma_n = Cons.M*(-Cons.M1*Vars.chi1_n + Cons.M2*Vars.chi2_n)
    Vars.Sigma_lambda = Cons.M*(-Cons.M1*Vars.chi1_lambda + Cons.M2*Vars.chi2_lambda)
    Vars.SS = np.array([dot(Vars.S[1:],Vars.S[1:])])
    Vars.SigmaSigma = np.array([dot(Vars.Sigma[1:],Vars.Sigma[1:])])
    Vars.SSigma = np.array([dot(Vars.S[1:],Vars.Sigma[1:])])
    Vars.S1 = Cons.M1**2*Vars.chiVec1
    Vars.S1_n = Cons.M1**2*Vars.chi1_n
    Vars.S1_lambda = Cons.M1**2*Vars.chi1_lambda
    Vars.S2 = Cons.M2**2*Vars.chiVec2
    Vars.S2_n = Cons.M2**2*Vars.chi2_n
    Vars.S2_lambda = Cons.M2**2*Vars.chi2_lambda
    Vars.chi_s_ell = Vars.chi1_ell/2 + Vars.chi2_ell/2
    Vars.chi_a_ell = Vars.chi1_ell/2 - Vars.chi2_ell/2
    Vars.logv = log(Vars.v)
    Vars.Fcal_coeff = 32*Cons.nu**2*Vars.v**10/5
    Vars.Fcal_SQ_4 = Cons.chi1chi1*(-89*Cons.delta/192 + 89*Cons.nu/96 - 89/192) - 103*Cons.chi1chi2*Cons.nu/48 + Cons.chi2chi2*(89*Cons.delta/192 + 89*Cons.nu/96 - 89/192) + Vars.chi_a_ell*(Vars.chi_a_ell*(287/96 - 12*Cons.nu) + 287*Vars.chi_s_ell*Cons.delta/48) + Vars.chi_s_ell**2*(Cons.nu/24 + 287/96)
    Vars.Fcal_SQ_6 = Cons.nu*(-0.0714285714285714*Cons.delta*(353*Vars.SSigma + 86*Vars.S_lambda*Vars.Sigma_lambda - 1145*Vars.S_n*Vars.Sigma_n)/Cons.M**4 - 0.00297619047619048*(4236*Vars.SS*Cons.kappa_p + 8472*Vars.SS + 4236*Vars.SSigma*Cons.delta*Cons.kappa_p + 1672*Vars.SSigma*Cons.kappa_m + 1032*Vars.S_lambda**2*(Cons.kappa_p + 2) + 8*Vars.S_lambda*Vars.Sigma_lambda*(129*Cons.delta*Cons.kappa_p + 620*Cons.kappa_m) - 13740*Vars.S_n**2*(Cons.kappa_p + 2) - 4*Vars.S_n*Vars.Sigma_n*(3435*Cons.delta*Cons.kappa_p + 2494*Cons.kappa_m) - 641*Vars.SigmaSigma*Cons.delta*Cons.kappa_m - 2713*Vars.SigmaSigma*Cons.kappa_p - 4146*Vars.SigmaSigma + 982*Vars.Sigma_lambda**2*Cons.delta*Cons.kappa_m - 6944*Vars.Sigma_lambda**2*Cons.kappa_p - 3736*Vars.Sigma_lambda**2 + 941*Vars.Sigma_n**2*Cons.delta*Cons.kappa_m + 15083*Vars.Sigma_n**2*Cons.kappa_p + 44524*Vars.Sigma_n**2)/Cons.M**4) - 0.00595238095238095*Cons.delta*(911*Vars.SSigma + 3643*Vars.S_lambda*Vars.Sigma_lambda - 22399*Vars.S_n*Vars.Sigma_n)/Cons.M**4 + 0.0357142857142857*Cons.nu**2*(Cons.kappa_p + 2)*(353*Vars.SigmaSigma + 86*Vars.Sigma_lambda**2 - 1145*Vars.Sigma_n**2)/Cons.M**4 + 0.00297619047619048*(1477*Vars.SS*Cons.delta*Cons.kappa_m - 1877*Vars.SS*Cons.kappa_p + 502*Vars.SS - 3354*Vars.SSigma*Cons.delta*Cons.kappa_p + 3354*Vars.SSigma*Cons.kappa_m + Vars.S_lambda**2*(1498*Cons.delta*Cons.kappa_m - 4464*Cons.kappa_p - 10811) - 5962*Vars.S_lambda*Vars.Sigma_lambda*(Cons.delta*Cons.kappa_p - Cons.kappa_m) + Vars.S_n**2*(-5929*Cons.delta*Cons.kappa_m + 10095*Cons.kappa_p + 45299) + 16024*Vars.S_n*Vars.Sigma_n*(Cons.delta*Cons.kappa_p - Cons.kappa_m) + 1677*Vars.SigmaSigma*Cons.delta*Cons.kappa_m - 1677*Vars.SigmaSigma*Cons.kappa_p - 168*Vars.SigmaSigma + 2981*Vars.Sigma_lambda**2*Cons.delta*Cons.kappa_m - 2981*Vars.Sigma_lambda**2*Cons.kappa_p + 1313*Vars.Sigma_lambda**2 - 8012*Vars.Sigma_n**2*Cons.delta*Cons.kappa_m + 8012*Vars.Sigma_n**2*Cons.kappa_p + 6387*Vars.Sigma_n**2)/Cons.M**4
    Vars.Fcal_SQ_7 = Vars.S_ell**3*(-16*Cons.kappa_p/3 - 4*Cons.lambda_p + 40/3) + Vars.S_ell**2*Vars.Sigma_ell*(-35*Cons.delta*Cons.kappa_p/6 - 6*Cons.delta*Cons.lambda_p + 73*Cons.delta/3 - 3*Cons.kappa_m/4 + 6*Cons.lambda_m) + Vars.S_ell*Vars.Sigma_ell**2*(-35*Cons.delta*Cons.kappa_m/12 + 6*Cons.delta*Cons.lambda_m + 35*Cons.kappa_p/12 - 6*Cons.lambda_p + Cons.nu*(22*Cons.kappa_p/3 + 12*Cons.lambda_p - 172/3) + 32/3) + Vars.Sigma_ell**3*(67*Cons.delta*Cons.kappa_p/24 - 2*Cons.delta*Cons.lambda_p - Cons.delta/8 - 67*Cons.kappa_m/24 + 2*Cons.lambda_m + Cons.nu*(Cons.delta*Cons.kappa_p/2 + 2*Cons.delta*Cons.lambda_p - 11*Cons.delta + 61*Cons.kappa_m/12 - 6*Cons.lambda_m))
    Vars.Fcal_SO_3 = (-4*Vars.S_ell - 5*Vars.Sigma_ell*Cons.delta/4)/Cons.M**2
    Vars.Fcal_SO_5 = (Vars.S_ell*(272*Cons.nu/9 - 9/2) + Vars.Sigma_ell*Cons.delta*(43*Cons.nu/4 - 13/16))/Cons.M**2
    Vars.Fcal_SO_6 = (-16*Vars.S_ell*pi - 31*Vars.Sigma_ell*Cons.delta*pi/6)/Cons.M**2
    Vars.Fcal_SO_7 = (Vars.S_ell*(-2810*Cons.nu**2/27 + 6172*Cons.nu/189 + 476645/6804) + Vars.Sigma_ell*Cons.delta*(-1501*Cons.nu**2/36 + 1849*Cons.nu/126 + 9535/336))/Cons.M**2
    Vars.E_SQ_4 = -3*Vars.chi_a_ell**2/2 - 3*Vars.chi_s_ell**2/2 - Cons.delta*(Cons.chi2chi2/2 + 3*Vars.chi_a_ell*Vars.chi_s_ell) + Cons.nu*(Cons.chi1chi2 + 6*Vars.chi_a_ell**2) + (Cons.chi1chi1 + Cons.chi2chi2)*(Cons.delta - 2*Cons.nu + 1)/4
    Vars.E_SQ_6 = Cons.nu*(-3*Cons.delta*(Vars.SSigma - 3*Vars.S_n*Vars.Sigma_n)/Cons.M**4 + 0.25*(-6*Vars.SS*Cons.kappa_p - 12*Vars.SS - 6*Vars.SSigma*Cons.delta*Cons.kappa_p - 26*Vars.SSigma*Cons.kappa_m + 24*Vars.S_lambda*Vars.Sigma_lambda*Cons.kappa_m + 18*Vars.S_n**2*(Cons.kappa_p + 2) + 18*Vars.S_n*Vars.Sigma_n*(Cons.delta*Cons.kappa_p + 3*Cons.kappa_m) + 19*Vars.SigmaSigma*Cons.kappa_p - 80*Vars.SigmaSigma - 18*Vars.Sigma_lambda**2*Cons.kappa_p + 76*Vars.Sigma_lambda**2 - 39*Vars.Sigma_n**2*Cons.kappa_p + 168*Vars.Sigma_n**2 + Cons.delta*Cons.kappa_m*(-5*Vars.SigmaSigma + 6*Vars.Sigma_lambda**2 + 9*Vars.Sigma_n**2))/Cons.M**4) + Cons.delta*(15*Vars.SSigma - 13*Vars.S_lambda*Vars.Sigma_lambda - 29*Vars.S_n*Vars.Sigma_n)/Cons.M**4 + 1.5*Cons.nu**2*(Vars.SigmaSigma - 3*Vars.Sigma_n**2)*(Cons.kappa_p + 2)/Cons.M**4 + 0.25*(8*Vars.SS*Cons.delta*Cons.kappa_m - 6*Vars.SS*Cons.kappa_p + 36*Vars.SS - 14*Vars.SSigma*Cons.delta*Cons.kappa_p + 14*Vars.SSigma*Cons.kappa_m + Vars.S_lambda**2*(-6*Cons.delta*Cons.kappa_m + 6*Cons.kappa_p - 28) + 12*Vars.S_lambda*Vars.Sigma_lambda*(Cons.delta*Cons.kappa_p - Cons.kappa_m) - 6*Vars.S_n**2*(3*Cons.delta*Cons.kappa_m - 2*Cons.kappa_p + 10) + 30*Vars.S_n*Vars.Sigma_n*(Cons.delta*Cons.kappa_p - Cons.kappa_m) + 7*Vars.SigmaSigma*Cons.delta*Cons.kappa_m - 7*Vars.SigmaSigma*Cons.kappa_p + 24*Vars.SigmaSigma - 6*Vars.Sigma_lambda**2*Cons.delta*Cons.kappa_m + 6*Vars.Sigma_lambda**2*Cons.kappa_p - 24*Vars.Sigma_lambda**2 - 15*Vars.Sigma_n**2*Cons.delta*Cons.kappa_m + 15*Vars.Sigma_n**2*Cons.kappa_p - 48*Vars.Sigma_n**2)/Cons.M**4
    Vars.E_SQ_7 = Vars.S_ell**3*(2*Cons.kappa_p + 4*Cons.lambda_p - 20) + Vars.S_ell**2*Vars.Sigma_ell*(2*Cons.delta*Cons.kappa_p + 6*Cons.delta*Cons.lambda_p - 32*Cons.delta + 4*Cons.kappa_m - 6*Cons.lambda_m) + Vars.S_ell*Vars.Sigma_ell**2*(5*Cons.delta*Cons.kappa_m - 6*Cons.delta*Cons.lambda_m - 5*Cons.kappa_p + 6*Cons.lambda_p + Cons.nu*(-2*Cons.kappa_p - 12*Cons.lambda_p + 68) - 12) + Vars.Sigma_ell**3*(-3*Cons.delta*Cons.kappa_p + 2*Cons.delta*Cons.lambda_p + 3*Cons.kappa_m - 2*Cons.lambda_m + Cons.nu*(-2*Cons.delta*Cons.lambda_p + 12*Cons.delta - 6*Cons.kappa_m + 6*Cons.lambda_m))
    Vars.E_SO_3 = (14*Vars.S_ell/3 + 2*Vars.Sigma_ell*Cons.delta)/Cons.M**2
    Vars.E_SO_5 = (Vars.S_ell*(11 - 61*Cons.nu/9) + Vars.Sigma_ell*Cons.delta*(3 - 10*Cons.nu/3))/Cons.M**2
    Vars.E_SO_7 = (Vars.S_ell*(29*Cons.nu**2/12 - 367*Cons.nu/4 + 135/4) + Vars.Sigma_ell*Cons.delta*(5*Cons.nu**2/4 - 39*Cons.nu + 27/4))/Cons.M**2
    Vars.MDot_Alvi_5 = -Cons.M1**3*Vars.chi1_ell*(3*Cons.chi1chi1 + 1)/(4*Cons.M**3) - Cons.M2**3*Vars.chi2_ell*(3*Cons.chi2chi2 + 1)/(4*Cons.M**3)

@njit
def OmegaVec_chiVec_1_3p5(Cons,Vars):
    Omega1_coeff = Vars.v**5/Cons.M
    return Omega1_coeff*(Vars.S2*Vars.v**3*((2.0 - 0.5*Cons.delta)/Cons.M**2 + Vars.v*(-5.0*Vars.S_ell - 3.0*Vars.Sigma_ell*Cons.delta/Cons.M)/Cons.M**4) + Vars.ellHat*(-0.75*Cons.delta + 0.5*Cons.nu + Vars.v**2*(Cons.delta*(0.625*Cons.nu - 0.5625) + Cons.nu*(1.25 - 0.0416666666666667*Cons.nu) + Vars.v*(Vars.v*(Cons.nu*(Cons.nu*(0.229166666666667*Cons.nu - 4.94791666666667) + 2.89583333333333) + 0.96875 + (Cons.delta*(Cons.nu*(5.91666666666667 - 0.489583333333333*Cons.nu) - 0.96875) + (Vars.S_ell*(Vars.S_ell*(Cons.kappa_p*(-0.25*Cons.nu - 0.375) - 0.5*Cons.nu - 0.75) + Vars.Sigma_ell*(Cons.delta*(Cons.kappa_p*(-0.25*Cons.nu - 0.375) - 0.5*Cons.nu - 0.75) + Cons.kappa_m*(0.25*Cons.nu + 0.375))) + Vars.Sigma_ell**2*(Cons.delta*Cons.kappa_m*(0.125*Cons.nu + 0.1875) + Cons.kappa_p*(Cons.nu*(0.25*Cons.nu + 0.25) - 0.1875) + Cons.nu*(0.5*Cons.nu + 0.75)) + (Vars.S_ell*(Vars.S_ell*Cons.delta*(0.375*Cons.kappa_p + 0.75) + Vars.Sigma_ell*Cons.delta*(Cons.delta*(0.375*Cons.kappa_p + 0.75) - 0.375*Cons.kappa_m)) + Vars.Sigma_ell**2*Cons.delta*(-0.1875*Cons.delta*Cons.kappa_m + Cons.kappa_p*(0.1875 - 0.375*Cons.nu) - 0.75*Cons.nu))/Cons.M)/Cons.M**7)/Cons.M) + (Vars.S_ell*(0.833333333333333*Cons.nu + 1.25) + (-1.25*Vars.S_ell*Cons.delta + Vars.Sigma_ell*Cons.delta*(0.5*Cons.nu + 0.75) - 0.75*Vars.Sigma_ell*Cons.delta**2/Cons.M)/Cons.M)/Cons.M**2) + 0.5625) + 0.75) + Vars.lambdaHat*Vars.v**3*(Vars.S2_lambda*(0.5*Cons.delta + 2.0)/Cons.M**2 + Vars.v*(Vars.S1_lambda*(-0.75*Cons.delta - 3.0*Cons.kappa_1 + 2.25) + Vars.S1_lambda*(Cons.delta*(0.75 - 1.5*Cons.kappa_1) + 1.5*Cons.kappa_1 - 0.75)/Cons.nu)/Cons.M**2) + Vars.nHat*Vars.v*(3.0*Vars.chi1_n*Cons.nu + Vars.v**2*((Vars.S1_n*(Cons.delta*(-2.25*Cons.kappa_1 - 0.25) + Cons.kappa_1*(1.5*Cons.nu + 0.75) + 0.5*Cons.nu + 3.25) + Vars.S1_n*(Cons.delta*(1.5 - 0.75*Cons.kappa_1) + 0.75*Cons.kappa_1 - 1.5)/Cons.nu + Vars.S2_n*(1.5*Cons.delta - 1.0*Cons.nu + 3.0))/Cons.M**2 + Vars.v*(-15.0*Vars.S1_n*Vars.S_ell*Cons.kappa_1 + Vars.S1_n*Vars.S_ell*(-7.5*Cons.delta*Cons.kappa_1 + 7.5*Cons.kappa_1)/Cons.nu + 15.0*Vars.S2_n*Vars.S_ell + (-9.0*Vars.S1_n*Vars.Sigma_ell*Cons.delta*Cons.kappa_1 + Vars.S1_n*Vars.Sigma_ell*Cons.delta*(-4.5*Cons.delta*Cons.kappa_1 + 4.5*Cons.kappa_1)/Cons.nu + 9.0*Vars.S2_n*Vars.Sigma_ell*Cons.delta)/Cons.M)/Cons.M**4) + 3.0*Cons.M2**2*Vars.chi2_n/Cons.M**2) - Cons.M2**2*Vars.chiVec2*Vars.v/Cons.M**2)

@njit
def OmegaVec_chiVec_2_3p5(Cons,Vars):
    Omega2_coeff = Vars.v**5/Cons.M
    return Omega2_coeff*(Vars.S1*Vars.v**3*((0.5*Cons.delta + 2.0)/Cons.M**2 + Vars.v*(-5.0*Vars.S_ell + 3.0*Vars.Sigma_ell*Cons.delta/Cons.M)/Cons.M**4) + Vars.ellHat*(0.75*Cons.delta + 0.5*Cons.nu + Vars.v**2*(Cons.delta*(0.5625 - 0.625*Cons.nu) + Cons.nu*(1.25 - 0.0416666666666667*Cons.nu) + Vars.v*(Vars.v*(Cons.nu*(Cons.nu*(0.229166666666667*Cons.nu - 4.94791666666667) + 2.89583333333333) + 0.96875 + (Cons.delta*(Cons.nu*(0.489583333333333*Cons.nu - 5.91666666666667) + 0.96875) + (Vars.S_ell*(Vars.S_ell*(Cons.kappa_p*(-0.25*Cons.nu - 0.375) - 0.5*Cons.nu - 0.75) + Vars.Sigma_ell*(Cons.delta*(Cons.kappa_p*(0.25*Cons.nu + 0.375) + 0.5*Cons.nu + 0.75) + Cons.kappa_m*(-0.25*Cons.nu - 0.375))) + Vars.Sigma_ell**2*(Cons.delta*Cons.kappa_m*(0.125*Cons.nu + 0.1875) + Cons.kappa_p*(Cons.nu*(0.25*Cons.nu + 0.25) - 0.1875) + Cons.nu*(0.5*Cons.nu + 0.75)) + (Vars.S_ell*(Vars.S_ell*Cons.delta*(-0.375*Cons.kappa_p - 0.75) + Vars.Sigma_ell*Cons.delta*(Cons.delta*(0.375*Cons.kappa_p + 0.75) - 0.375*Cons.kappa_m)) + Vars.Sigma_ell**2*Cons.delta*(0.1875*Cons.delta*Cons.kappa_m + Cons.kappa_p*(0.375*Cons.nu - 0.1875) + 0.75*Cons.nu))/Cons.M)/Cons.M**7)/Cons.M) + (Vars.S_ell*(0.833333333333333*Cons.nu + 1.25) + (1.25*Vars.S_ell*Cons.delta + Vars.Sigma_ell*Cons.delta*(-0.5*Cons.nu - 0.75) - 0.75*Vars.Sigma_ell*Cons.delta**2/Cons.M)/Cons.M)/Cons.M**2) + 0.5625) + 0.75) + Vars.lambdaHat*Vars.v**3*(Vars.S1_lambda*(2.0 - 0.5*Cons.delta)/Cons.M**2 + Vars.v*(Vars.S2_lambda*(0.75*Cons.delta - 3.0*Cons.kappa_2 + 2.25) + Vars.S2_lambda*(Cons.delta*(1.5*Cons.kappa_2 - 0.75) + 1.5*Cons.kappa_2 - 0.75)/Cons.nu)/Cons.M**2) + Vars.nHat*Vars.v*(3.0*Vars.chi2_n*Cons.nu + Vars.v**2*((Vars.S1_n*(-1.5*Cons.delta - 1.0*Cons.nu + 3.0) + Vars.S2_n*(Cons.delta*(2.25*Cons.kappa_2 + 0.25) + Cons.kappa_2*(1.5*Cons.nu + 0.75) + 0.5*Cons.nu + 3.25) + Vars.S2_n*(Cons.delta*(0.75*Cons.kappa_2 - 1.5) + 0.75*Cons.kappa_2 - 1.5)/Cons.nu)/Cons.M**2 + Vars.v*(15.0*Vars.S1_n*Vars.S_ell - 15.0*Vars.S2_n*Vars.S_ell*Cons.kappa_2 + Vars.S2_n*Vars.S_ell*(7.5*Cons.delta*Cons.kappa_2 + 7.5*Cons.kappa_2)/Cons.nu + (-9.0*Vars.S1_n*Vars.Sigma_ell*Cons.delta + 9.0*Vars.S2_n*Vars.Sigma_ell*Cons.delta*Cons.kappa_2 + Vars.S2_n*Vars.Sigma_ell*Cons.delta*(-4.5*Cons.delta*Cons.kappa_2 - 4.5*Cons.kappa_2)/Cons.nu)/Cons.M)/Cons.M**4) + 3.0*Cons.M1**2*Vars.chi1_n/Cons.M**2) - Cons.M1**2*Vars.chiVec1*Vars.v/Cons.M**2)

@njit
def OmegaVec_3p5(Cons,Vars):
    gamma_PN_5 = (0.888888888888889*Vars.S_ell*Cons.nu + 3.33333333333333*Vars.S_ell)/Cons.M**2 + 2.0*Vars.Sigma_ell*Cons.delta/Cons.M**3
    gamma_PN_6 = Vars.S_ell**2*(-0.916666666666667*Cons.delta*Cons.kappa_m - 0.916666666666667*Cons.kappa_p + Cons.nu*(-0.166666666666667*Cons.kappa_p - 0.333333333333333) + 1.55555555555556) + Vars.S_ell*Vars.Sigma_ell*(1.66666666666667*Cons.delta + Cons.nu*(-0.166666666666667*Cons.delta*Cons.kappa_p - 0.333333333333333*Cons.delta + 3.83333333333333*Cons.kappa_m)) + Vars.Sigma_ell**2*(Cons.nu**2*(0.166666666666667*Cons.kappa_p + 0.333333333333333) + Cons.nu*(Cons.delta*Cons.kappa_m - Cons.kappa_p - 2.0) + 1.0) + 0.0123456790123457*Cons.nu**3 + 6.36111111111111*Cons.nu**2 - 2.98177812235564*Cons.nu + 1.0
    gamma_PN_0 = 1.00000000000000
    a_ell_4 = Vars.S_n*(5.77777777777778*Cons.nu**2 + 14.75*Cons.nu + 1.5) + Vars.Sigma_n*Cons.delta*(2.83333333333333*Cons.nu**2 + 9.125*Cons.nu + 1.5)
    a_ell_1 = 3.0*Vars.Sigma_ell*Vars.Sigma_n*Cons.nu*(Cons.kappa_p + 2.0)/Cons.M**2 - 3.0*Cons.delta*(Vars.S_ell*Vars.Sigma_n + Vars.S_n*Vars.Sigma_ell)/Cons.M**2 - 0.5*(3.0*Vars.S_n*(2.0*Vars.S_ell*(Cons.kappa_p + 2.0) + Vars.Sigma_ell*(Cons.delta*Cons.kappa_p - Cons.kappa_m)) + 3.0*Vars.Sigma_n*(Vars.S_ell*(Cons.delta*Cons.kappa_p - Cons.kappa_m) + Vars.Sigma_ell*(-Cons.delta*Cons.kappa_m + Cons.kappa_p)))/Cons.M**2
    a_ell_2 = Vars.S_n*(-9.66666666666667*Cons.nu - 10.0) + Vars.Sigma_n*Cons.delta*(-4.5*Cons.nu - 6.0)
    gamma_PN_7 = (Vars.S_ell*(-6.0*Cons.nu**2 - 10.5833333333333*Cons.nu + 5.0) - 2.66666666666667*Vars.Sigma_ell*Cons.delta*Cons.nu**2 + Vars.Sigma_ell*Cons.delta*(3.0 - 10.1666666666667*Cons.nu))/Cons.M**2
    gamma_PN_2 = 1.0 - 0.333333333333333*Cons.nu
    a_ell_3 = Cons.nu*(4.0*Cons.delta*(Vars.S_ell*Vars.Sigma_n + Vars.S_n*Vars.Sigma_ell)/Cons.M**2 + 0.5*(2.0*Vars.S_n*(4.0*Vars.S_ell*(Cons.kappa_p + 2.0) + Vars.Sigma_ell*(2.0*Cons.delta*Cons.kappa_p + 7.0*Cons.kappa_m)) + Vars.Sigma_n*(2.0*Vars.S_ell*(2.0*Cons.delta*Cons.kappa_p + 7.0*Cons.kappa_m) + Vars.Sigma_ell*(5.0*Cons.delta*Cons.kappa_m - 17.0*Cons.kappa_p - 26.0)))/Cons.M**2) - 4.0*Vars.Sigma_ell*Vars.Sigma_n*Cons.nu**2*(Cons.kappa_p + 2.0)/Cons.M**2 + 6.0*Cons.delta*(Vars.S_ell*Vars.Sigma_n + Vars.S_n*Vars.Sigma_ell)/Cons.M**2 + 0.5*(Vars.S_n*(Vars.S_ell*(-9.0*Cons.delta*Cons.kappa_m + 3.0*Cons.kappa_p + 22.0) + Vars.Sigma_ell*(6.0*Cons.delta*Cons.kappa_p - 6.0*Cons.kappa_m)) + 6.0*Vars.Sigma_n*(Vars.S_ell*(Cons.delta*Cons.kappa_p - Cons.kappa_m) + Vars.Sigma_ell*(-Cons.delta*Cons.kappa_m + Cons.kappa_p + 1.0)))/Cons.M**2
    gamma_PN_3 = (1.66666666666667*Vars.S_ell + Vars.Sigma_ell*Cons.delta)/Cons.M**2
    a_ell_0 = 7.0*Vars.S_n + 3.0*Vars.Sigma_n*Cons.delta
    gamma_PN_4 = -5.41666666666667*Cons.nu + 1.0 + (Vars.S_ell**2*(-0.5*Cons.kappa_p - 1.0) + Vars.S_ell*Vars.Sigma_ell*(-0.5*Cons.delta*Cons.kappa_p - Cons.delta + 0.5*Cons.kappa_m) + Vars.Sigma_ell**2*(0.25*Cons.delta*Cons.kappa_m - 0.25*Cons.kappa_p + Cons.nu*(0.5*Cons.kappa_p + 1.0)))/Cons.M**8
    return Vars.ellHat*Vars.v**3/Cons.M + Vars.nHat*Vars.v**6*(a_ell_0 + Vars.v*(a_ell_1 + Vars.v*(a_ell_2 + Vars.v*(a_ell_3 + a_ell_4*Vars.v))))*(gamma_PN_0 + Vars.v**2*(gamma_PN_2 + Vars.v*(gamma_PN_3 + Vars.v*(gamma_PN_4 + Vars.v*(gamma_PN_5 + Vars.v*(gamma_PN_6 + gamma_PN_7*Vars.v))))))/Cons.M**3


@njit(cache=True)
def TaylorT1_3p5(Cons,Vars):
    Flux = Vars.Fcal_coeff*(Cons.Fcal_0 + Vars.v**2*(Cons.Fcal_2 + Vars.v*(Cons.Fcal_3 + Vars.Fcal_SO_3 + Vars.v*(Cons.Fcal_4 + Vars.Fcal_SQ_4 + Vars.v*(Cons.Fcal_5 + Vars.Fcal_SO_5 + Vars.v*(Cons.Fcal_6 + Vars.Fcal_SO_6 + Vars.Fcal_SQ_6 + Cons.Fcal_lnv_6*Vars.logv + Vars.v*(Cons.Fcal_7 + Vars.Fcal_SO_7 + Vars.Fcal_SQ_7)))))))
    dEdV = -0.5*Cons.M*Cons.nu*Vars.v*(2.0*Cons.E_0 + Vars.v**2*(4.0*Cons.E_2 + Vars.v*(5.0*Vars.E_SO_3 + Vars.v*(6.0*Cons.E_4 + 6.0*Vars.E_SQ_4 + Vars.v*(7.0*Vars.E_SO_5 + Vars.v*(8.0*Cons.E_6 + 8.0*Vars.E_SQ_6 + 9.0*Vars.v*(Vars.E_SO_7 + Vars.E_SQ_7)))))))
    Absorption = Vars.Fcal_coeff*Vars.MDot_Alvi_5*Vars.v**5
    dvdt_T1 = (-Absorption - Flux)/dEdV
    dydt=np.zeros(9)
    [dydt[5],dydt[6],dydt[7],dydt[8]] = FrameFromAngularVelocityIntegrand(Vars.R, OmegaVec_3p5(Cons,Vars)[1:])
    dydt[0] = dvdt_T1[0] 
    if(Cons.EvolveSpin1):
        dydt[1], dydt[2]=FrameFromAngularVelocity_2D_Integrand(Vars.rfrak_chi1[0], Vars.rfrak_chi1[1],(mul(mul(inverse(Cons.S_chi1),OmegaVec_chiVec_1_3p5(Cons,Vars)),Cons.S_chi1))[1:])
    else:
        dydt[1], dydt[2] = 0.0, 0.0
    if(Cons.EvolveSpin2):
        dydt[3], dydt[4]=FrameFromAngularVelocity_2D_Integrand(Vars.rfrak_chi2[0], Vars.rfrak_chi2[1],(mul(mul(inverse(Cons.S_chi2),OmegaVec_chiVec_2_3p5(Cons,Vars)),Cons.S_chi2))[1:])
    else:
        dydt[3], dydt[4] = 0.0, 0.0
    return dydt

@njit(cache=True)
def TaylorT4_3p5(Cons,Vars):
    Flux = Vars.Fcal_coeff*(Cons.Fcal_0 + Vars.v**2*(Cons.Fcal_2 + Vars.v*(Cons.Fcal_3 + Vars.Fcal_SO_3 + Vars.v*(Cons.Fcal_4 + Vars.Fcal_SQ_4 + Vars.v*(Cons.Fcal_5 + Vars.Fcal_SO_5 + Vars.v*(Cons.Fcal_6 + Vars.Fcal_SO_6 + Vars.Fcal_SQ_6 + Cons.Fcal_lnv_6*Vars.logv + Vars.v*(Cons.Fcal_7 + Vars.Fcal_SO_7 + Vars.Fcal_SQ_7)))))))
    dEdV = -0.5*Cons.M*Cons.nu*Vars.v*(2.0*Cons.E_0 + Vars.v**2*(4.0*Cons.E_2 + Vars.v*(5.0*Vars.E_SO_3 + Vars.v*(6.0*Cons.E_4 + 6.0*Vars.E_SQ_4 + Vars.v*(7.0*Vars.E_SO_5 + Vars.v*(8.0*Cons.E_6 + 8.0*Vars.E_SQ_6 + 9.0*Vars.v*(Vars.E_SO_7 + Vars.E_SQ_7)))))))
    Absorption = Vars.Fcal_coeff*Vars.MDot_Alvi_5*Vars.v**5
    dvdt_T1 = (-Absorption - Flux)/dEdV
    dydt=np.zeros(9)
    [dydt[5],dydt[6],dydt[7],dydt[8]] = FrameFromAngularVelocityIntegrand(Vars.R, OmegaVec_3p5(Cons,Vars)[1:])
    dydt[0] = dvdt_T4[0]
    if(Cons.EvolveSpin1):
        dydt[1], dydt[2]=FrameFromAngularVelocity_2D_Integrand(Vars.rfrak_chi1[0], Vars.rfrak_chi1[1],(mul(mul(inverse(Cons.S_chi1),OmegaVec_chiVec_1_3p5(Cons,Vars)),Cons.S_chi1))[1:])
    else:
        dydt[1], dydt[2] = 0.0, 0.0
    if(Cons.EvolveSpin2):
        dydt[3], dydt[4]=FrameFromAngularVelocity_2D_Integrand(Vars.rfrak_chi2[0], Vars.rfrak_chi2[1],(mul(mul(inverse(Cons.S_chi2),OmegaVec_chiVec_2_3p5(Cons,Vars)),Cons.S_chi2))[1:])
    else:
        dydt[3], dydt[4] = 0.0, 0.0
    return dydt      

@njit(cache=True)
def TaylorT5_3p5(Cons,Vars):
    Flux = Vars.Fcal_coeff*(Cons.Fcal_0 + Vars.v**2*(Cons.Fcal_2 + Vars.v*(Cons.Fcal_3 + Vars.Fcal_SO_3 + Vars.v*(Cons.Fcal_4 + Vars.Fcal_SQ_4 + Vars.v*(Cons.Fcal_5 + Vars.Fcal_SO_5 + Vars.v*(Cons.Fcal_6 + Vars.Fcal_SO_6 + Vars.Fcal_SQ_6 + Cons.Fcal_lnv_6*Vars.logv + Vars.v*(Cons.Fcal_7 + Vars.Fcal_SO_7 + Vars.Fcal_SQ_7)))))))
    dEdV = -0.5*Cons.M*Cons.nu*Vars.v*(2.0*Cons.E_0 + Vars.v**2*(4.0*Cons.E_2 + Vars.v*(5.0*Vars.E_SO_3 + Vars.v*(6.0*Cons.E_4 + 6.0*Vars.E_SQ_4 + Vars.v*(7.0*Vars.E_SO_5 + Vars.v*(8.0*Cons.E_6 + 8.0*Vars.E_SQ_6 + 9.0*Vars.v*(Vars.E_SO_7 + Vars.E_SQ_7)))))))
    Absorption = Vars.Fcal_coeff*Vars.MDot_Alvi_5*Vars.v**5
    dvdt_T1 = (-Absorption - Flux)/dEdV
    dydt=np.zeros(8)
    [dydt[5],dydt[6],dydt[7],dydt[8]] = FrameFromAngularVelocityIntegrand(Vars.R, OmegaVec_3p5(Cons,Vars)[1:])
    dydt[0] = dvdt_T5[0]
    if(Cons.EvolveSpin1):
        dydt[1], dydt[2]=FrameFromAngularVelocity_2D_Integrand(Vars.rfrak_chi1[0], Vars.rfrak_chi1[1],(mul(mul(inverse(Cons.S_chi1),OmegaVec_chiVec_1_3p5(Cons,Vars)),Cons.S_chi1))[1:])
    else:
        dydt[1], dydt[2] = 0.0, 0.0
    if(Cons.EvolveSpin2):
        dydt[3], dydt[4]=FrameFromAngularVelocity_2D_Integrand(Vars.rfrak_chi2[0], Vars.rfrak_chi2[1],(mul(mul(inverse(Cons.S_chi2),OmegaVec_chiVec_2_3p5(Cons,Vars)),Cons.S_chi2))[1:])
    else:
        dydt[3], dydt[4] = 0.0, 0.0
    return dydt      

@njit(cache=True)
def Recalculate_4p0(Cons,Vars,y):
    Vars.v = np.array([y[0]])
    Vars.rfrak_chi1 = np.array([y[1],y[2]])
    Vars.rfrak_chi2 = np.array([y[3],y[4]])
    Vars.R = y[5]*Cons.wHat + y[6]*Cons.xHat + y[7]*Cons.yHat + y[8]*Cons.zHat
    Vars.nHat = mul(mul(Vars.R,Cons.xHat),inverse(Vars.R))
    Vars.lambdaHat = mul(mul(Vars.R,Cons.yHat),inverse(Vars.R))
    Vars.ellHat = mul(mul(Vars.R,Cons.zHat),inverse(Vars.R))
    Vars.R_S1 = exp(Vars.rfrak_chi1[0]*Cons.xHat + Vars.rfrak_chi1[1]*Cons.yHat)
    Vars.R_S2 = exp(Vars.rfrak_chi2[0]*Cons.xHat + Vars.rfrak_chi2[1]*Cons.yHat)
    Vars.chiVec1 = mul(mul(mul(Cons.S_chi1,Vars.R_S1),Cons.zHat),mul(conjugate(Vars.R_S1),conjugate(Cons.S_chi1)))
    Vars.chiVec2 = mul(mul(mul(Cons.S_chi2,Vars.R_S2),Cons.zHat),mul(conjugate(Vars.R_S2),conjugate(Cons.S_chi2)))
    Vars.chi1_n = np.array([dot(Vars.chiVec1[1:],Vars.nHat[1:])])
    Vars.chi1_lambda = np.array([dot(Vars.chiVec1[1:],Vars.lambdaHat[1:])])
    Vars.chi1_ell = np.array([dot(Vars.chiVec1[1:],Vars.ellHat[1:])])
    Vars.chi2_n = np.array([dot(Vars.chiVec2[1:],Vars.nHat[1:])])
    Vars.chi2_lambda = np.array([dot(Vars.chiVec2[1:],Vars.lambdaHat[1:])])
    Vars.chi2_ell = np.array([dot(Vars.chiVec2[1:],Vars.ellHat[1:])])
    Vars.S = Cons.M1**2*Vars.chiVec1 + Cons.M2**2*Vars.chiVec2
    Vars.S_ell = Cons.M1**2*Vars.chi1_ell + Cons.M2**2*Vars.chi2_ell
    Vars.S_n = Cons.M1**2*Vars.chi1_n + Cons.M2**2*Vars.chi2_n
    Vars.S_lambda = Cons.M1**2*Vars.chi1_lambda + Cons.M2**2*Vars.chi2_lambda
    Vars.Sigma = Cons.M*(-Cons.M1*Vars.chiVec1 + Cons.M2*Vars.chiVec2)
    Vars.Sigma_ell = Cons.M*(-Cons.M1*Vars.chi1_ell + Cons.M2*Vars.chi2_ell)
    Vars.Sigma_n = Cons.M*(-Cons.M1*Vars.chi1_n + Cons.M2*Vars.chi2_n)
    Vars.Sigma_lambda = Cons.M*(-Cons.M1*Vars.chi1_lambda + Cons.M2*Vars.chi2_lambda)
    Vars.SS = np.array([dot(Vars.S[1:],Vars.S[1:])])
    Vars.SigmaSigma = np.array([dot(Vars.Sigma[1:],Vars.Sigma[1:])])
    Vars.SSigma = np.array([dot(Vars.S[1:],Vars.Sigma[1:])])
    Vars.S1 = Cons.M1**2*Vars.chiVec1
    Vars.S1_n = Cons.M1**2*Vars.chi1_n
    Vars.S1_lambda = Cons.M1**2*Vars.chi1_lambda
    Vars.S2 = Cons.M2**2*Vars.chiVec2
    Vars.S2_n = Cons.M2**2*Vars.chi2_n
    Vars.S2_lambda = Cons.M2**2*Vars.chi2_lambda
    Vars.chi_s_ell = Vars.chi1_ell/2 + Vars.chi2_ell/2
    Vars.chi_a_ell = Vars.chi1_ell/2 - Vars.chi2_ell/2
    Vars.logv = log(Vars.v)
    Vars.Fcal_coeff = 32*Cons.nu**2*Vars.v**10/5
    Vars.Fcal_SQ_4 = Cons.chi1chi1*(-89*Cons.delta/192 + 89*Cons.nu/96 - 89/192) - 103*Cons.chi1chi2*Cons.nu/48 + Cons.chi2chi2*(89*Cons.delta/192 + 89*Cons.nu/96 - 89/192) + Vars.chi_a_ell*(Vars.chi_a_ell*(287/96 - 12*Cons.nu) + 287*Vars.chi_s_ell*Cons.delta/48) + Vars.chi_s_ell**2*(Cons.nu/24 + 287/96)
    Vars.Fcal_SQ_6 = Cons.nu*(-0.0714285714285714*Cons.delta*(353*Vars.SSigma + 86*Vars.S_lambda*Vars.Sigma_lambda - 1145*Vars.S_n*Vars.Sigma_n)/Cons.M**4 - 0.00297619047619048*(4236*Vars.SS*Cons.kappa_p + 8472*Vars.SS + 4236*Vars.SSigma*Cons.delta*Cons.kappa_p + 1672*Vars.SSigma*Cons.kappa_m + 1032*Vars.S_lambda**2*(Cons.kappa_p + 2) + 8*Vars.S_lambda*Vars.Sigma_lambda*(129*Cons.delta*Cons.kappa_p + 620*Cons.kappa_m) - 13740*Vars.S_n**2*(Cons.kappa_p + 2) - 4*Vars.S_n*Vars.Sigma_n*(3435*Cons.delta*Cons.kappa_p + 2494*Cons.kappa_m) - 641*Vars.SigmaSigma*Cons.delta*Cons.kappa_m - 2713*Vars.SigmaSigma*Cons.kappa_p - 4146*Vars.SigmaSigma + 982*Vars.Sigma_lambda**2*Cons.delta*Cons.kappa_m - 6944*Vars.Sigma_lambda**2*Cons.kappa_p - 3736*Vars.Sigma_lambda**2 + 941*Vars.Sigma_n**2*Cons.delta*Cons.kappa_m + 15083*Vars.Sigma_n**2*Cons.kappa_p + 44524*Vars.Sigma_n**2)/Cons.M**4) - 0.00595238095238095*Cons.delta*(911*Vars.SSigma + 3643*Vars.S_lambda*Vars.Sigma_lambda - 22399*Vars.S_n*Vars.Sigma_n)/Cons.M**4 + 0.0357142857142857*Cons.nu**2*(Cons.kappa_p + 2)*(353*Vars.SigmaSigma + 86*Vars.Sigma_lambda**2 - 1145*Vars.Sigma_n**2)/Cons.M**4 + 0.00297619047619048*(1477*Vars.SS*Cons.delta*Cons.kappa_m - 1877*Vars.SS*Cons.kappa_p + 502*Vars.SS - 3354*Vars.SSigma*Cons.delta*Cons.kappa_p + 3354*Vars.SSigma*Cons.kappa_m + Vars.S_lambda**2*(1498*Cons.delta*Cons.kappa_m - 4464*Cons.kappa_p - 10811) - 5962*Vars.S_lambda*Vars.Sigma_lambda*(Cons.delta*Cons.kappa_p - Cons.kappa_m) + Vars.S_n**2*(-5929*Cons.delta*Cons.kappa_m + 10095*Cons.kappa_p + 45299) + 16024*Vars.S_n*Vars.Sigma_n*(Cons.delta*Cons.kappa_p - Cons.kappa_m) + 1677*Vars.SigmaSigma*Cons.delta*Cons.kappa_m - 1677*Vars.SigmaSigma*Cons.kappa_p - 168*Vars.SigmaSigma + 2981*Vars.Sigma_lambda**2*Cons.delta*Cons.kappa_m - 2981*Vars.Sigma_lambda**2*Cons.kappa_p + 1313*Vars.Sigma_lambda**2 - 8012*Vars.Sigma_n**2*Cons.delta*Cons.kappa_m + 8012*Vars.Sigma_n**2*Cons.kappa_p + 6387*Vars.Sigma_n**2)/Cons.M**4
    Vars.Fcal_SQ_7 = Vars.S_ell**3*(-16*Cons.kappa_p/3 - 4*Cons.lambda_p + 40/3) + Vars.S_ell**2*Vars.Sigma_ell*(-35*Cons.delta*Cons.kappa_p/6 - 6*Cons.delta*Cons.lambda_p + 73*Cons.delta/3 - 3*Cons.kappa_m/4 + 6*Cons.lambda_m) + Vars.S_ell*Vars.Sigma_ell**2*(-35*Cons.delta*Cons.kappa_m/12 + 6*Cons.delta*Cons.lambda_m + 35*Cons.kappa_p/12 - 6*Cons.lambda_p + Cons.nu*(22*Cons.kappa_p/3 + 12*Cons.lambda_p - 172/3) + 32/3) + Vars.Sigma_ell**3*(67*Cons.delta*Cons.kappa_p/24 - 2*Cons.delta*Cons.lambda_p - Cons.delta/8 - 67*Cons.kappa_m/24 + 2*Cons.lambda_m + Cons.nu*(Cons.delta*Cons.kappa_p/2 + 2*Cons.delta*Cons.lambda_p - 11*Cons.delta + 61*Cons.kappa_m/12 - 6*Cons.lambda_m))
    Vars.Fcal_SO_3 = (-4*Vars.S_ell - 5*Vars.Sigma_ell*Cons.delta/4)/Cons.M**2
    Vars.Fcal_SO_5 = (Vars.S_ell*(272*Cons.nu/9 - 9/2) + Vars.Sigma_ell*Cons.delta*(43*Cons.nu/4 - 13/16))/Cons.M**2
    Vars.Fcal_SO_6 = (-16*Vars.S_ell*pi - 31*Vars.Sigma_ell*Cons.delta*pi/6)/Cons.M**2
    Vars.Fcal_SO_7 = (Vars.S_ell*(-2810*Cons.nu**2/27 + 6172*Cons.nu/189 + 476645/6804) + Vars.Sigma_ell*Cons.delta*(-1501*Cons.nu**2/36 + 1849*Cons.nu/126 + 9535/336))/Cons.M**2
    Vars.Fcal_SO_8 = (Vars.S_ell*pi*(13879*Cons.nu/72 - 3485/96) + Vars.Sigma_ell*Cons.delta*pi*(130583*Cons.nu/2016 - 7163/672))/Cons.M**2
    Vars.E_SQ_4 = -3*Vars.chi_a_ell**2/2 - 3*Vars.chi_s_ell**2/2 - Cons.delta*(Cons.chi2chi2/2 + 3*Vars.chi_a_ell*Vars.chi_s_ell) + Cons.nu*(Cons.chi1chi2 + 6*Vars.chi_a_ell**2) + (Cons.chi1chi1 + Cons.chi2chi2)*(Cons.delta - 2*Cons.nu + 1)/4
    Vars.E_SQ_6 = Cons.nu*(-3*Cons.delta*(Vars.SSigma - 3*Vars.S_n*Vars.Sigma_n)/Cons.M**4 + 0.25*(-6*Vars.SS*Cons.kappa_p - 12*Vars.SS - 6*Vars.SSigma*Cons.delta*Cons.kappa_p - 26*Vars.SSigma*Cons.kappa_m + 24*Vars.S_lambda*Vars.Sigma_lambda*Cons.kappa_m + 18*Vars.S_n**2*(Cons.kappa_p + 2) + 18*Vars.S_n*Vars.Sigma_n*(Cons.delta*Cons.kappa_p + 3*Cons.kappa_m) + 19*Vars.SigmaSigma*Cons.kappa_p - 80*Vars.SigmaSigma - 18*Vars.Sigma_lambda**2*Cons.kappa_p + 76*Vars.Sigma_lambda**2 - 39*Vars.Sigma_n**2*Cons.kappa_p + 168*Vars.Sigma_n**2 + Cons.delta*Cons.kappa_m*(-5*Vars.SigmaSigma + 6*Vars.Sigma_lambda**2 + 9*Vars.Sigma_n**2))/Cons.M**4) + Cons.delta*(15*Vars.SSigma - 13*Vars.S_lambda*Vars.Sigma_lambda - 29*Vars.S_n*Vars.Sigma_n)/Cons.M**4 + 1.5*Cons.nu**2*(Vars.SigmaSigma - 3*Vars.Sigma_n**2)*(Cons.kappa_p + 2)/Cons.M**4 + 0.25*(8*Vars.SS*Cons.delta*Cons.kappa_m - 6*Vars.SS*Cons.kappa_p + 36*Vars.SS - 14*Vars.SSigma*Cons.delta*Cons.kappa_p + 14*Vars.SSigma*Cons.kappa_m + Vars.S_lambda**2*(-6*Cons.delta*Cons.kappa_m + 6*Cons.kappa_p - 28) + 12*Vars.S_lambda*Vars.Sigma_lambda*(Cons.delta*Cons.kappa_p - Cons.kappa_m) - 6*Vars.S_n**2*(3*Cons.delta*Cons.kappa_m - 2*Cons.kappa_p + 10) + 30*Vars.S_n*Vars.Sigma_n*(Cons.delta*Cons.kappa_p - Cons.kappa_m) + 7*Vars.SigmaSigma*Cons.delta*Cons.kappa_m - 7*Vars.SigmaSigma*Cons.kappa_p + 24*Vars.SigmaSigma - 6*Vars.Sigma_lambda**2*Cons.delta*Cons.kappa_m + 6*Vars.Sigma_lambda**2*Cons.kappa_p - 24*Vars.Sigma_lambda**2 - 15*Vars.Sigma_n**2*Cons.delta*Cons.kappa_m + 15*Vars.Sigma_n**2*Cons.kappa_p - 48*Vars.Sigma_n**2)/Cons.M**4
    Vars.E_SQ_7 = Vars.S_ell**3*(2*Cons.kappa_p + 4*Cons.lambda_p - 20) + Vars.S_ell**2*Vars.Sigma_ell*(2*Cons.delta*Cons.kappa_p + 6*Cons.delta*Cons.lambda_p - 32*Cons.delta + 4*Cons.kappa_m - 6*Cons.lambda_m) + Vars.S_ell*Vars.Sigma_ell**2*(5*Cons.delta*Cons.kappa_m - 6*Cons.delta*Cons.lambda_m - 5*Cons.kappa_p + 6*Cons.lambda_p + Cons.nu*(-2*Cons.kappa_p - 12*Cons.lambda_p + 68) - 12) + Vars.Sigma_ell**3*(-3*Cons.delta*Cons.kappa_p + 2*Cons.delta*Cons.lambda_p + 3*Cons.kappa_m - 2*Cons.lambda_m + Cons.nu*(-2*Cons.delta*Cons.lambda_p + 12*Cons.delta - 6*Cons.kappa_m + 6*Cons.lambda_m))
    Vars.E_SO_3 = (14*Vars.S_ell/3 + 2*Vars.Sigma_ell*Cons.delta)/Cons.M**2
    Vars.E_SO_5 = (Vars.S_ell*(11 - 61*Cons.nu/9) + Vars.Sigma_ell*Cons.delta*(3 - 10*Cons.nu/3))/Cons.M**2
    Vars.E_SO_7 = (Vars.S_ell*(29*Cons.nu**2/12 - 367*Cons.nu/4 + 135/4) + Vars.Sigma_ell*Cons.delta*(5*Cons.nu**2/4 - 39*Cons.nu + 27/4))/Cons.M**2
    Vars.MDot_Alvi_5 = -Cons.M1**3*Vars.chi1_ell*(3*Cons.chi1chi1 + 1)/(4*Cons.M**3) - Cons.M2**3*Vars.chi2_ell*(3*Cons.chi2chi2 + 1)/(4*Cons.M**3)

@njit
def OmegaVec_chiVec_1_4p0(Cons,Vars):
    Omega1_coeff = Vars.v**5/Cons.M
    return Omega1_coeff*(Vars.S2*Vars.v**3*((2.0 - 0.5*Cons.delta)/Cons.M**2 + Vars.v*(-5.0*Vars.S_ell - 3.0*Vars.Sigma_ell*Cons.delta/Cons.M)/Cons.M**4) + Vars.ellHat*(-0.75*Cons.delta + 0.5*Cons.nu + Vars.v**2*(Cons.delta*(0.625*Cons.nu - 0.5625) + Cons.nu*(1.25 - 0.0416666666666667*Cons.nu) + Vars.v*(Vars.v*(Cons.nu*(Cons.nu*(0.229166666666667*Cons.nu - 4.94791666666667) + 2.89583333333333) + 0.96875 + (Cons.delta*(Cons.nu*(5.91666666666667 - 0.489583333333333*Cons.nu) - 0.96875) + (Vars.S_ell*(Vars.S_ell*(Cons.kappa_p*(-0.25*Cons.nu - 0.375) - 0.5*Cons.nu - 0.75) + Vars.Sigma_ell*(Cons.delta*(Cons.kappa_p*(-0.25*Cons.nu - 0.375) - 0.5*Cons.nu - 0.75) + Cons.kappa_m*(0.25*Cons.nu + 0.375))) + Vars.Sigma_ell**2*(Cons.delta*Cons.kappa_m*(0.125*Cons.nu + 0.1875) + Cons.kappa_p*(Cons.nu*(0.25*Cons.nu + 0.25) - 0.1875) + Cons.nu*(0.5*Cons.nu + 0.75)) + (Vars.S_ell*(Vars.S_ell*Cons.delta*(0.375*Cons.kappa_p + 0.75) + Vars.Sigma_ell*Cons.delta*(Cons.delta*(0.375*Cons.kappa_p + 0.75) - 0.375*Cons.kappa_m)) + Vars.Sigma_ell**2*Cons.delta*(-0.1875*Cons.delta*Cons.kappa_m + Cons.kappa_p*(0.1875 - 0.375*Cons.nu) - 0.75*Cons.nu))/Cons.M)/Cons.M**7)/Cons.M) + (Vars.S_ell*(0.833333333333333*Cons.nu + 1.25) + (-1.25*Vars.S_ell*Cons.delta + Vars.Sigma_ell*Cons.delta*(0.5*Cons.nu + 0.75) - 0.75*Vars.Sigma_ell*Cons.delta**2/Cons.M)/Cons.M)/Cons.M**2) + 0.5625) + 0.75) + Vars.lambdaHat*Vars.v**3*(Vars.S2_lambda*(0.5*Cons.delta + 2.0)/Cons.M**2 + Vars.v*(Vars.S1_lambda*(-0.75*Cons.delta - 3.0*Cons.kappa_1 + 2.25) + Vars.S1_lambda*(Cons.delta*(0.75 - 1.5*Cons.kappa_1) + 1.5*Cons.kappa_1 - 0.75)/Cons.nu)/Cons.M**2) + Vars.nHat*Vars.v*(3.0*Vars.chi1_n*Cons.nu + Vars.v**2*((Vars.S1_n*(Cons.delta*(-2.25*Cons.kappa_1 - 0.25) + Cons.kappa_1*(1.5*Cons.nu + 0.75) + 0.5*Cons.nu + 3.25) + Vars.S1_n*(Cons.delta*(1.5 - 0.75*Cons.kappa_1) + 0.75*Cons.kappa_1 - 1.5)/Cons.nu + Vars.S2_n*(1.5*Cons.delta - 1.0*Cons.nu + 3.0))/Cons.M**2 + Vars.v*(-15.0*Vars.S1_n*Vars.S_ell*Cons.kappa_1 + Vars.S1_n*Vars.S_ell*(-7.5*Cons.delta*Cons.kappa_1 + 7.5*Cons.kappa_1)/Cons.nu + 15.0*Vars.S2_n*Vars.S_ell + (-9.0*Vars.S1_n*Vars.Sigma_ell*Cons.delta*Cons.kappa_1 + Vars.S1_n*Vars.Sigma_ell*Cons.delta*(-4.5*Cons.delta*Cons.kappa_1 + 4.5*Cons.kappa_1)/Cons.nu + 9.0*Vars.S2_n*Vars.Sigma_ell*Cons.delta)/Cons.M)/Cons.M**4) + 3.0*Cons.M2**2*Vars.chi2_n/Cons.M**2) - Cons.M2**2*Vars.chiVec2*Vars.v/Cons.M**2)

@njit
def OmegaVec_chiVec_2_4p0(Cons,Vars):
    Omega2_coeff = Vars.v**5/Cons.M
    return Omega2_coeff*(Vars.S1*Vars.v**3*((0.5*Cons.delta + 2.0)/Cons.M**2 + Vars.v*(-5.0*Vars.S_ell + 3.0*Vars.Sigma_ell*Cons.delta/Cons.M)/Cons.M**4) + Vars.ellHat*(0.75*Cons.delta + 0.5*Cons.nu + Vars.v**2*(Cons.delta*(0.5625 - 0.625*Cons.nu) + Cons.nu*(1.25 - 0.0416666666666667*Cons.nu) + Vars.v*(Vars.v*(Cons.nu*(Cons.nu*(0.229166666666667*Cons.nu - 4.94791666666667) + 2.89583333333333) + 0.96875 + (Cons.delta*(Cons.nu*(0.489583333333333*Cons.nu - 5.91666666666667) + 0.96875) + (Vars.S_ell*(Vars.S_ell*(Cons.kappa_p*(-0.25*Cons.nu - 0.375) - 0.5*Cons.nu - 0.75) + Vars.Sigma_ell*(Cons.delta*(Cons.kappa_p*(0.25*Cons.nu + 0.375) + 0.5*Cons.nu + 0.75) + Cons.kappa_m*(-0.25*Cons.nu - 0.375))) + Vars.Sigma_ell**2*(Cons.delta*Cons.kappa_m*(0.125*Cons.nu + 0.1875) + Cons.kappa_p*(Cons.nu*(0.25*Cons.nu + 0.25) - 0.1875) + Cons.nu*(0.5*Cons.nu + 0.75)) + (Vars.S_ell*(Vars.S_ell*Cons.delta*(-0.375*Cons.kappa_p - 0.75) + Vars.Sigma_ell*Cons.delta*(Cons.delta*(0.375*Cons.kappa_p + 0.75) - 0.375*Cons.kappa_m)) + Vars.Sigma_ell**2*Cons.delta*(0.1875*Cons.delta*Cons.kappa_m + Cons.kappa_p*(0.375*Cons.nu - 0.1875) + 0.75*Cons.nu))/Cons.M)/Cons.M**7)/Cons.M) + (Vars.S_ell*(0.833333333333333*Cons.nu + 1.25) + (1.25*Vars.S_ell*Cons.delta + Vars.Sigma_ell*Cons.delta*(-0.5*Cons.nu - 0.75) - 0.75*Vars.Sigma_ell*Cons.delta**2/Cons.M)/Cons.M)/Cons.M**2) + 0.5625) + 0.75) + Vars.lambdaHat*Vars.v**3*(Vars.S1_lambda*(2.0 - 0.5*Cons.delta)/Cons.M**2 + Vars.v*(Vars.S2_lambda*(0.75*Cons.delta - 3.0*Cons.kappa_2 + 2.25) + Vars.S2_lambda*(Cons.delta*(1.5*Cons.kappa_2 - 0.75) + 1.5*Cons.kappa_2 - 0.75)/Cons.nu)/Cons.M**2) + Vars.nHat*Vars.v*(3.0*Vars.chi2_n*Cons.nu + Vars.v**2*((Vars.S1_n*(-1.5*Cons.delta - 1.0*Cons.nu + 3.0) + Vars.S2_n*(Cons.delta*(2.25*Cons.kappa_2 + 0.25) + Cons.kappa_2*(1.5*Cons.nu + 0.75) + 0.5*Cons.nu + 3.25) + Vars.S2_n*(Cons.delta*(0.75*Cons.kappa_2 - 1.5) + 0.75*Cons.kappa_2 - 1.5)/Cons.nu)/Cons.M**2 + Vars.v*(15.0*Vars.S1_n*Vars.S_ell - 15.0*Vars.S2_n*Vars.S_ell*Cons.kappa_2 + Vars.S2_n*Vars.S_ell*(7.5*Cons.delta*Cons.kappa_2 + 7.5*Cons.kappa_2)/Cons.nu + (-9.0*Vars.S1_n*Vars.Sigma_ell*Cons.delta + 9.0*Vars.S2_n*Vars.Sigma_ell*Cons.delta*Cons.kappa_2 + Vars.S2_n*Vars.Sigma_ell*Cons.delta*(-4.5*Cons.delta*Cons.kappa_2 - 4.5*Cons.kappa_2)/Cons.nu)/Cons.M)/Cons.M**4) + 3.0*Cons.M1**2*Vars.chi1_n/Cons.M**2) - Cons.M1**2*Vars.chiVec1*Vars.v/Cons.M**2)

@njit
def OmegaVec_4p0(Cons,Vars):
    gamma_PN_5 = (0.888888888888889*Vars.S_ell*Cons.nu + 3.33333333333333*Vars.S_ell)/Cons.M**2 + 2.0*Vars.Sigma_ell*Cons.delta/Cons.M**3
    gamma_PN_6 = Vars.S_ell**2*(-0.916666666666667*Cons.delta*Cons.kappa_m - 0.916666666666667*Cons.kappa_p + Cons.nu*(-0.166666666666667*Cons.kappa_p - 0.333333333333333) + 1.55555555555556) + Vars.S_ell*Vars.Sigma_ell*(1.66666666666667*Cons.delta + Cons.nu*(-0.166666666666667*Cons.delta*Cons.kappa_p - 0.333333333333333*Cons.delta + 3.83333333333333*Cons.kappa_m)) + Vars.Sigma_ell**2*(Cons.nu**2*(0.166666666666667*Cons.kappa_p + 0.333333333333333) + Cons.nu*(Cons.delta*Cons.kappa_m - Cons.kappa_p - 2.0) + 1.0) + 0.0123456790123457*Cons.nu**3 + 6.36111111111111*Cons.nu**2 - 2.98177812235564*Cons.nu + 1.0
    gamma_PN_0 = 1.00000000000000
    a_ell_4 = Vars.S_n*(5.77777777777778*Cons.nu**2 + 14.75*Cons.nu + 1.5) + Vars.Sigma_n*Cons.delta*(2.83333333333333*Cons.nu**2 + 9.125*Cons.nu + 1.5)
    a_ell_1 = 3.0*Vars.Sigma_ell*Vars.Sigma_n*Cons.nu*(Cons.kappa_p + 2.0)/Cons.M**2 - 3.0*Cons.delta*(Vars.S_ell*Vars.Sigma_n + Vars.S_n*Vars.Sigma_ell)/Cons.M**2 - 0.5*(3.0*Vars.S_n*(2.0*Vars.S_ell*(Cons.kappa_p + 2.0) + Vars.Sigma_ell*(Cons.delta*Cons.kappa_p - Cons.kappa_m)) + 3.0*Vars.Sigma_n*(Vars.S_ell*(Cons.delta*Cons.kappa_p - Cons.kappa_m) + Vars.Sigma_ell*(-Cons.delta*Cons.kappa_m + Cons.kappa_p)))/Cons.M**2
    a_ell_2 = Vars.S_n*(-9.66666666666667*Cons.nu - 10.0) + Vars.Sigma_n*Cons.delta*(-4.5*Cons.nu - 6.0)
    gamma_PN_7 = (Vars.S_ell*(-6.0*Cons.nu**2 - 10.5833333333333*Cons.nu + 5.0) - 2.66666666666667*Vars.Sigma_ell*Cons.delta*Cons.nu**2 + Vars.Sigma_ell*Cons.delta*(3.0 - 10.1666666666667*Cons.nu))/Cons.M**2
    gamma_PN_2 = 1.0 - 0.333333333333333*Cons.nu
    a_ell_3 = Cons.nu*(4.0*Cons.delta*(Vars.S_ell*Vars.Sigma_n + Vars.S_n*Vars.Sigma_ell)/Cons.M**2 + 0.5*(2.0*Vars.S_n*(4.0*Vars.S_ell*(Cons.kappa_p + 2.0) + Vars.Sigma_ell*(2.0*Cons.delta*Cons.kappa_p + 7.0*Cons.kappa_m)) + Vars.Sigma_n*(2.0*Vars.S_ell*(2.0*Cons.delta*Cons.kappa_p + 7.0*Cons.kappa_m) + Vars.Sigma_ell*(5.0*Cons.delta*Cons.kappa_m - 17.0*Cons.kappa_p - 26.0)))/Cons.M**2) - 4.0*Vars.Sigma_ell*Vars.Sigma_n*Cons.nu**2*(Cons.kappa_p + 2.0)/Cons.M**2 + 6.0*Cons.delta*(Vars.S_ell*Vars.Sigma_n + Vars.S_n*Vars.Sigma_ell)/Cons.M**2 + 0.5*(Vars.S_n*(Vars.S_ell*(-9.0*Cons.delta*Cons.kappa_m + 3.0*Cons.kappa_p + 22.0) + Vars.Sigma_ell*(6.0*Cons.delta*Cons.kappa_p - 6.0*Cons.kappa_m)) + 6.0*Vars.Sigma_n*(Vars.S_ell*(Cons.delta*Cons.kappa_p - Cons.kappa_m) + Vars.Sigma_ell*(-Cons.delta*Cons.kappa_m + Cons.kappa_p + 1.0)))/Cons.M**2
    gamma_PN_3 = (1.66666666666667*Vars.S_ell + Vars.Sigma_ell*Cons.delta)/Cons.M**2
    a_ell_0 = 7.0*Vars.S_n + 3.0*Vars.Sigma_n*Cons.delta
    gamma_PN_4 = -5.41666666666667*Cons.nu + 1.0 + (Vars.S_ell**2*(-0.5*Cons.kappa_p - 1.0) + Vars.S_ell*Vars.Sigma_ell*(-0.5*Cons.delta*Cons.kappa_p - Cons.delta + 0.5*Cons.kappa_m) + Vars.Sigma_ell**2*(0.25*Cons.delta*Cons.kappa_m - 0.25*Cons.kappa_p + Cons.nu*(0.5*Cons.kappa_p + 1.0)))/Cons.M**8
    return Vars.ellHat*Vars.v**3/Cons.M + Vars.nHat*Vars.v**6*(a_ell_0 + Vars.v*(a_ell_1 + Vars.v*(a_ell_2 + Vars.v*(a_ell_3 + a_ell_4*Vars.v))))*(gamma_PN_0 + Vars.v**2*(gamma_PN_2 + Vars.v*(gamma_PN_3 + Vars.v*(gamma_PN_4 + Vars.v*(gamma_PN_5 + Vars.v*(gamma_PN_6 + gamma_PN_7*Vars.v))))))/Cons.M**3


@njit(cache=True)
def TaylorT1_4p0(Cons,Vars):
    Flux = Vars.Fcal_coeff*(Cons.Fcal_0 + Vars.v**2*(Cons.Fcal_2 + Vars.v*(Cons.Fcal_3 + Vars.Fcal_SO_3 + Vars.v*(Cons.Fcal_4 + Vars.Fcal_SQ_4 + Vars.v*(Cons.Fcal_5 + Vars.Fcal_SO_5 + Vars.v*(Cons.Fcal_6 + Vars.Fcal_SO_6 + Vars.Fcal_SQ_6 + Cons.Fcal_lnv_6*Vars.logv + Vars.v*(Cons.Fcal_7 + Vars.Fcal_SO_7 + Vars.Fcal_SQ_7 + Vars.v*(Cons.Fcal_8 + Vars.Fcal_SO_8 + Cons.Fcal_lnv_8*Vars.logv))))))))
    dEdV = -0.5*Cons.M*Cons.nu*Vars.v*(2.0*Cons.E_0 + Vars.v**2*(4.0*Cons.E_2 + Vars.v*(5.0*Vars.E_SO_3 + Vars.v*(6.0*Cons.E_4 + 6.0*Vars.E_SQ_4 + Vars.v*(7.0*Vars.E_SO_5 + Vars.v*(8.0*Cons.E_6 + 8.0*Vars.E_SQ_6 + Vars.v*(9.0*Vars.E_SO_7 + 9.0*Vars.E_SQ_7 + Vars.v*(10.0*Cons.E_8 + Cons.E_lnv_8*(10.0*Vars.logv + 1.0)))))))))
    Absorption = Vars.Fcal_coeff*Vars.v**5*(Vars.MDot_Alvi_5 + Cons.MDot_Alvi_8*Vars.v**3)
    dvdt_T1 = (-Absorption - Flux)/dEdV
    dydt=np.zeros(9)
    [dydt[5],dydt[6],dydt[7],dydt[8]] = FrameFromAngularVelocityIntegrand(Vars.R, OmegaVec_4p0(Cons,Vars)[1:])
    dydt[0] = dvdt_T1[0] 
    if(Cons.EvolveSpin1):
        dydt[1], dydt[2]=FrameFromAngularVelocity_2D_Integrand(Vars.rfrak_chi1[0], Vars.rfrak_chi1[1],(mul(mul(inverse(Cons.S_chi1),OmegaVec_chiVec_1_4p0(Cons,Vars)),Cons.S_chi1))[1:])
    else:
        dydt[1], dydt[2] = 0.0, 0.0
    if(Cons.EvolveSpin2):
        dydt[3], dydt[4]=FrameFromAngularVelocity_2D_Integrand(Vars.rfrak_chi2[0], Vars.rfrak_chi2[1],(mul(mul(inverse(Cons.S_chi2),OmegaVec_chiVec_2_4p0(Cons,Vars)),Cons.S_chi2))[1:])
    else:
        dydt[3], dydt[4] = 0.0, 0.0
    return dydt

@njit(cache=True)
def TaylorT4_4p0(Cons,Vars):
    Flux = Vars.Fcal_coeff*(Cons.Fcal_0 + Vars.v**2*(Cons.Fcal_2 + Vars.v*(Cons.Fcal_3 + Vars.Fcal_SO_3 + Vars.v*(Cons.Fcal_4 + Vars.Fcal_SQ_4 + Vars.v*(Cons.Fcal_5 + Vars.Fcal_SO_5 + Vars.v*(Cons.Fcal_6 + Vars.Fcal_SO_6 + Vars.Fcal_SQ_6 + Cons.Fcal_lnv_6*Vars.logv + Vars.v*(Cons.Fcal_7 + Vars.Fcal_SO_7 + Vars.Fcal_SQ_7 + Vars.v*(Cons.Fcal_8 + Vars.Fcal_SO_8 + Cons.Fcal_lnv_8*Vars.logv))))))))
    dEdV = -0.5*Cons.M*Cons.nu*Vars.v*(2.0*Cons.E_0 + Vars.v**2*(4.0*Cons.E_2 + Vars.v*(5.0*Vars.E_SO_3 + Vars.v*(6.0*Cons.E_4 + 6.0*Vars.E_SQ_4 + Vars.v*(7.0*Vars.E_SO_5 + Vars.v*(8.0*Cons.E_6 + 8.0*Vars.E_SQ_6 + Vars.v*(9.0*Vars.E_SO_7 + 9.0*Vars.E_SQ_7 + Vars.v*(10.0*Cons.E_8 + Cons.E_lnv_8*(10.0*Vars.logv + 1.0)))))))))
    Absorption = Vars.Fcal_coeff*Vars.v**5*(Vars.MDot_Alvi_5 + Cons.MDot_Alvi_8*Vars.v**3)
    dvdt_T1 = (-Absorption - Flux)/dEdV
    dydt=np.zeros(9)
    [dydt[5],dydt[6],dydt[7],dydt[8]] = FrameFromAngularVelocityIntegrand(Vars.R, OmegaVec_4p0(Cons,Vars)[1:])
    dydt[0] = dvdt_T4[0]
    if(Cons.EvolveSpin1):
        dydt[1], dydt[2]=FrameFromAngularVelocity_2D_Integrand(Vars.rfrak_chi1[0], Vars.rfrak_chi1[1],(mul(mul(inverse(Cons.S_chi1),OmegaVec_chiVec_1_4p0(Cons,Vars)),Cons.S_chi1))[1:])
    else:
        dydt[1], dydt[2] = 0.0, 0.0
    if(Cons.EvolveSpin2):
        dydt[3], dydt[4]=FrameFromAngularVelocity_2D_Integrand(Vars.rfrak_chi2[0], Vars.rfrak_chi2[1],(mul(mul(inverse(Cons.S_chi2),OmegaVec_chiVec_2_4p0(Cons,Vars)),Cons.S_chi2))[1:])
    else:
        dydt[3], dydt[4] = 0.0, 0.0
    return dydt      

@njit(cache=True)
def TaylorT5_4p0(Cons,Vars):
    Flux = Vars.Fcal_coeff*(Cons.Fcal_0 + Vars.v**2*(Cons.Fcal_2 + Vars.v*(Cons.Fcal_3 + Vars.Fcal_SO_3 + Vars.v*(Cons.Fcal_4 + Vars.Fcal_SQ_4 + Vars.v*(Cons.Fcal_5 + Vars.Fcal_SO_5 + Vars.v*(Cons.Fcal_6 + Vars.Fcal_SO_6 + Vars.Fcal_SQ_6 + Cons.Fcal_lnv_6*Vars.logv + Vars.v*(Cons.Fcal_7 + Vars.Fcal_SO_7 + Vars.Fcal_SQ_7 + Vars.v*(Cons.Fcal_8 + Vars.Fcal_SO_8 + Cons.Fcal_lnv_8*Vars.logv))))))))
    dEdV = -0.5*Cons.M*Cons.nu*Vars.v*(2.0*Cons.E_0 + Vars.v**2*(4.0*Cons.E_2 + Vars.v*(5.0*Vars.E_SO_3 + Vars.v*(6.0*Cons.E_4 + 6.0*Vars.E_SQ_4 + Vars.v*(7.0*Vars.E_SO_5 + Vars.v*(8.0*Cons.E_6 + 8.0*Vars.E_SQ_6 + Vars.v*(9.0*Vars.E_SO_7 + 9.0*Vars.E_SQ_7 + Vars.v*(10.0*Cons.E_8 + Cons.E_lnv_8*(10.0*Vars.logv + 1.0)))))))))
    Absorption = Vars.Fcal_coeff*Vars.v**5*(Vars.MDot_Alvi_5 + Cons.MDot_Alvi_8*Vars.v**3)
    dvdt_T1 = (-Absorption - Flux)/dEdV
    dydt=np.zeros(8)
    [dydt[5],dydt[6],dydt[7],dydt[8]] = FrameFromAngularVelocityIntegrand(Vars.R, OmegaVec_4p0(Cons,Vars)[1:])
    dydt[0] = dvdt_T5[0]
    if(Cons.EvolveSpin1):
        dydt[1], dydt[2]=FrameFromAngularVelocity_2D_Integrand(Vars.rfrak_chi1[0], Vars.rfrak_chi1[1],(mul(mul(inverse(Cons.S_chi1),OmegaVec_chiVec_1_4p0(Cons,Vars)),Cons.S_chi1))[1:])
    else:
        dydt[1], dydt[2] = 0.0, 0.0
    if(Cons.EvolveSpin2):
        dydt[3], dydt[4]=FrameFromAngularVelocity_2D_Integrand(Vars.rfrak_chi2[0], Vars.rfrak_chi2[1],(mul(mul(inverse(Cons.S_chi2),OmegaVec_chiVec_2_4p0(Cons,Vars)),Cons.S_chi2))[1:])
    else:
        dydt[3], dydt[4] = 0.0, 0.0
    return dydt      

@njit(cache=True)
def Recalculate_4p5(Cons,Vars,y):
    Vars.v = np.array([y[0]])
    Vars.rfrak_chi1 = np.array([y[1],y[2]])
    Vars.rfrak_chi2 = np.array([y[3],y[4]])
    Vars.R = y[5]*Cons.wHat + y[6]*Cons.xHat + y[7]*Cons.yHat + y[8]*Cons.zHat
    Vars.nHat = mul(mul(Vars.R,Cons.xHat),inverse(Vars.R))
    Vars.lambdaHat = mul(mul(Vars.R,Cons.yHat),inverse(Vars.R))
    Vars.ellHat = mul(mul(Vars.R,Cons.zHat),inverse(Vars.R))
    Vars.R_S1 = exp(Vars.rfrak_chi1[0]*Cons.xHat + Vars.rfrak_chi1[1]*Cons.yHat)
    Vars.R_S2 = exp(Vars.rfrak_chi2[0]*Cons.xHat + Vars.rfrak_chi2[1]*Cons.yHat)
    Vars.chiVec1 = mul(mul(mul(Cons.S_chi1,Vars.R_S1),Cons.zHat),mul(conjugate(Vars.R_S1),conjugate(Cons.S_chi1)))
    Vars.chiVec2 = mul(mul(mul(Cons.S_chi2,Vars.R_S2),Cons.zHat),mul(conjugate(Vars.R_S2),conjugate(Cons.S_chi2)))
    Vars.chi1_n = np.array([dot(Vars.chiVec1[1:],Vars.nHat[1:])])
    Vars.chi1_lambda = np.array([dot(Vars.chiVec1[1:],Vars.lambdaHat[1:])])
    Vars.chi1_ell = np.array([dot(Vars.chiVec1[1:],Vars.ellHat[1:])])
    Vars.chi2_n = np.array([dot(Vars.chiVec2[1:],Vars.nHat[1:])])
    Vars.chi2_lambda = np.array([dot(Vars.chiVec2[1:],Vars.lambdaHat[1:])])
    Vars.chi2_ell = np.array([dot(Vars.chiVec2[1:],Vars.ellHat[1:])])
    Vars.S = Cons.M1**2*Vars.chiVec1 + Cons.M2**2*Vars.chiVec2
    Vars.S_ell = Cons.M1**2*Vars.chi1_ell + Cons.M2**2*Vars.chi2_ell
    Vars.S_n = Cons.M1**2*Vars.chi1_n + Cons.M2**2*Vars.chi2_n
    Vars.S_lambda = Cons.M1**2*Vars.chi1_lambda + Cons.M2**2*Vars.chi2_lambda
    Vars.Sigma = Cons.M*(-Cons.M1*Vars.chiVec1 + Cons.M2*Vars.chiVec2)
    Vars.Sigma_ell = Cons.M*(-Cons.M1*Vars.chi1_ell + Cons.M2*Vars.chi2_ell)
    Vars.Sigma_n = Cons.M*(-Cons.M1*Vars.chi1_n + Cons.M2*Vars.chi2_n)
    Vars.Sigma_lambda = Cons.M*(-Cons.M1*Vars.chi1_lambda + Cons.M2*Vars.chi2_lambda)
    Vars.SS = np.array([dot(Vars.S[1:],Vars.S[1:])])
    Vars.SigmaSigma = np.array([dot(Vars.Sigma[1:],Vars.Sigma[1:])])
    Vars.SSigma = np.array([dot(Vars.S[1:],Vars.Sigma[1:])])
    Vars.S1 = Cons.M1**2*Vars.chiVec1
    Vars.S1_n = Cons.M1**2*Vars.chi1_n
    Vars.S1_lambda = Cons.M1**2*Vars.chi1_lambda
    Vars.S2 = Cons.M2**2*Vars.chiVec2
    Vars.S2_n = Cons.M2**2*Vars.chi2_n
    Vars.S2_lambda = Cons.M2**2*Vars.chi2_lambda
    Vars.chi_s_ell = Vars.chi1_ell/2 + Vars.chi2_ell/2
    Vars.chi_a_ell = Vars.chi1_ell/2 - Vars.chi2_ell/2
    Vars.logv = log(Vars.v)
    Vars.Fcal_coeff = 32*Cons.nu**2*Vars.v**10/5
    Vars.Fcal_SQ_4 = Cons.chi1chi1*(-89*Cons.delta/192 + 89*Cons.nu/96 - 89/192) - 103*Cons.chi1chi2*Cons.nu/48 + Cons.chi2chi2*(89*Cons.delta/192 + 89*Cons.nu/96 - 89/192) + Vars.chi_a_ell*(Vars.chi_a_ell*(287/96 - 12*Cons.nu) + 287*Vars.chi_s_ell*Cons.delta/48) + Vars.chi_s_ell**2*(Cons.nu/24 + 287/96)
    Vars.Fcal_SQ_6 = Cons.nu*(-0.0714285714285714*Cons.delta*(353*Vars.SSigma + 86*Vars.S_lambda*Vars.Sigma_lambda - 1145*Vars.S_n*Vars.Sigma_n)/Cons.M**4 - 0.00297619047619048*(4236*Vars.SS*Cons.kappa_p + 8472*Vars.SS + 4236*Vars.SSigma*Cons.delta*Cons.kappa_p + 1672*Vars.SSigma*Cons.kappa_m + 1032*Vars.S_lambda**2*(Cons.kappa_p + 2) + 8*Vars.S_lambda*Vars.Sigma_lambda*(129*Cons.delta*Cons.kappa_p + 620*Cons.kappa_m) - 13740*Vars.S_n**2*(Cons.kappa_p + 2) - 4*Vars.S_n*Vars.Sigma_n*(3435*Cons.delta*Cons.kappa_p + 2494*Cons.kappa_m) - 641*Vars.SigmaSigma*Cons.delta*Cons.kappa_m - 2713*Vars.SigmaSigma*Cons.kappa_p - 4146*Vars.SigmaSigma + 982*Vars.Sigma_lambda**2*Cons.delta*Cons.kappa_m - 6944*Vars.Sigma_lambda**2*Cons.kappa_p - 3736*Vars.Sigma_lambda**2 + 941*Vars.Sigma_n**2*Cons.delta*Cons.kappa_m + 15083*Vars.Sigma_n**2*Cons.kappa_p + 44524*Vars.Sigma_n**2)/Cons.M**4) - 0.00595238095238095*Cons.delta*(911*Vars.SSigma + 3643*Vars.S_lambda*Vars.Sigma_lambda - 22399*Vars.S_n*Vars.Sigma_n)/Cons.M**4 + 0.0357142857142857*Cons.nu**2*(Cons.kappa_p + 2)*(353*Vars.SigmaSigma + 86*Vars.Sigma_lambda**2 - 1145*Vars.Sigma_n**2)/Cons.M**4 + 0.00297619047619048*(1477*Vars.SS*Cons.delta*Cons.kappa_m - 1877*Vars.SS*Cons.kappa_p + 502*Vars.SS - 3354*Vars.SSigma*Cons.delta*Cons.kappa_p + 3354*Vars.SSigma*Cons.kappa_m + Vars.S_lambda**2*(1498*Cons.delta*Cons.kappa_m - 4464*Cons.kappa_p - 10811) - 5962*Vars.S_lambda*Vars.Sigma_lambda*(Cons.delta*Cons.kappa_p - Cons.kappa_m) + Vars.S_n**2*(-5929*Cons.delta*Cons.kappa_m + 10095*Cons.kappa_p + 45299) + 16024*Vars.S_n*Vars.Sigma_n*(Cons.delta*Cons.kappa_p - Cons.kappa_m) + 1677*Vars.SigmaSigma*Cons.delta*Cons.kappa_m - 1677*Vars.SigmaSigma*Cons.kappa_p - 168*Vars.SigmaSigma + 2981*Vars.Sigma_lambda**2*Cons.delta*Cons.kappa_m - 2981*Vars.Sigma_lambda**2*Cons.kappa_p + 1313*Vars.Sigma_lambda**2 - 8012*Vars.Sigma_n**2*Cons.delta*Cons.kappa_m + 8012*Vars.Sigma_n**2*Cons.kappa_p + 6387*Vars.Sigma_n**2)/Cons.M**4
    Vars.Fcal_SQ_7 = Vars.S_ell**3*(-16*Cons.kappa_p/3 - 4*Cons.lambda_p + 40/3) + Vars.S_ell**2*Vars.Sigma_ell*(-35*Cons.delta*Cons.kappa_p/6 - 6*Cons.delta*Cons.lambda_p + 73*Cons.delta/3 - 3*Cons.kappa_m/4 + 6*Cons.lambda_m) + Vars.S_ell*Vars.Sigma_ell**2*(-35*Cons.delta*Cons.kappa_m/12 + 6*Cons.delta*Cons.lambda_m + 35*Cons.kappa_p/12 - 6*Cons.lambda_p + Cons.nu*(22*Cons.kappa_p/3 + 12*Cons.lambda_p - 172/3) + 32/3) + Vars.Sigma_ell**3*(67*Cons.delta*Cons.kappa_p/24 - 2*Cons.delta*Cons.lambda_p - Cons.delta/8 - 67*Cons.kappa_m/24 + 2*Cons.lambda_m + Cons.nu*(Cons.delta*Cons.kappa_p/2 + 2*Cons.delta*Cons.lambda_p - 11*Cons.delta + 61*Cons.kappa_m/12 - 6*Cons.lambda_m))
    Vars.Fcal_SO_3 = (-4*Vars.S_ell - 5*Vars.Sigma_ell*Cons.delta/4)/Cons.M**2
    Vars.Fcal_SO_5 = (Vars.S_ell*(272*Cons.nu/9 - 9/2) + Vars.Sigma_ell*Cons.delta*(43*Cons.nu/4 - 13/16))/Cons.M**2
    Vars.Fcal_SO_6 = (-16*Vars.S_ell*pi - 31*Vars.Sigma_ell*Cons.delta*pi/6)/Cons.M**2
    Vars.Fcal_SO_7 = (Vars.S_ell*(-2810*Cons.nu**2/27 + 6172*Cons.nu/189 + 476645/6804) + Vars.Sigma_ell*Cons.delta*(-1501*Cons.nu**2/36 + 1849*Cons.nu/126 + 9535/336))/Cons.M**2
    Vars.Fcal_SO_8 = (Vars.S_ell*pi*(13879*Cons.nu/72 - 3485/96) + Vars.Sigma_ell*Cons.delta*pi*(130583*Cons.nu/2016 - 7163/672))/Cons.M**2
    Vars.E_SQ_4 = -3*Vars.chi_a_ell**2/2 - 3*Vars.chi_s_ell**2/2 - Cons.delta*(Cons.chi2chi2/2 + 3*Vars.chi_a_ell*Vars.chi_s_ell) + Cons.nu*(Cons.chi1chi2 + 6*Vars.chi_a_ell**2) + (Cons.chi1chi1 + Cons.chi2chi2)*(Cons.delta - 2*Cons.nu + 1)/4
    Vars.E_SQ_6 = Cons.nu*(-3*Cons.delta*(Vars.SSigma - 3*Vars.S_n*Vars.Sigma_n)/Cons.M**4 + 0.25*(-6*Vars.SS*Cons.kappa_p - 12*Vars.SS - 6*Vars.SSigma*Cons.delta*Cons.kappa_p - 26*Vars.SSigma*Cons.kappa_m + 24*Vars.S_lambda*Vars.Sigma_lambda*Cons.kappa_m + 18*Vars.S_n**2*(Cons.kappa_p + 2) + 18*Vars.S_n*Vars.Sigma_n*(Cons.delta*Cons.kappa_p + 3*Cons.kappa_m) + 19*Vars.SigmaSigma*Cons.kappa_p - 80*Vars.SigmaSigma - 18*Vars.Sigma_lambda**2*Cons.kappa_p + 76*Vars.Sigma_lambda**2 - 39*Vars.Sigma_n**2*Cons.kappa_p + 168*Vars.Sigma_n**2 + Cons.delta*Cons.kappa_m*(-5*Vars.SigmaSigma + 6*Vars.Sigma_lambda**2 + 9*Vars.Sigma_n**2))/Cons.M**4) + Cons.delta*(15*Vars.SSigma - 13*Vars.S_lambda*Vars.Sigma_lambda - 29*Vars.S_n*Vars.Sigma_n)/Cons.M**4 + 1.5*Cons.nu**2*(Vars.SigmaSigma - 3*Vars.Sigma_n**2)*(Cons.kappa_p + 2)/Cons.M**4 + 0.25*(8*Vars.SS*Cons.delta*Cons.kappa_m - 6*Vars.SS*Cons.kappa_p + 36*Vars.SS - 14*Vars.SSigma*Cons.delta*Cons.kappa_p + 14*Vars.SSigma*Cons.kappa_m + Vars.S_lambda**2*(-6*Cons.delta*Cons.kappa_m + 6*Cons.kappa_p - 28) + 12*Vars.S_lambda*Vars.Sigma_lambda*(Cons.delta*Cons.kappa_p - Cons.kappa_m) - 6*Vars.S_n**2*(3*Cons.delta*Cons.kappa_m - 2*Cons.kappa_p + 10) + 30*Vars.S_n*Vars.Sigma_n*(Cons.delta*Cons.kappa_p - Cons.kappa_m) + 7*Vars.SigmaSigma*Cons.delta*Cons.kappa_m - 7*Vars.SigmaSigma*Cons.kappa_p + 24*Vars.SigmaSigma - 6*Vars.Sigma_lambda**2*Cons.delta*Cons.kappa_m + 6*Vars.Sigma_lambda**2*Cons.kappa_p - 24*Vars.Sigma_lambda**2 - 15*Vars.Sigma_n**2*Cons.delta*Cons.kappa_m + 15*Vars.Sigma_n**2*Cons.kappa_p - 48*Vars.Sigma_n**2)/Cons.M**4
    Vars.E_SQ_7 = Vars.S_ell**3*(2*Cons.kappa_p + 4*Cons.lambda_p - 20) + Vars.S_ell**2*Vars.Sigma_ell*(2*Cons.delta*Cons.kappa_p + 6*Cons.delta*Cons.lambda_p - 32*Cons.delta + 4*Cons.kappa_m - 6*Cons.lambda_m) + Vars.S_ell*Vars.Sigma_ell**2*(5*Cons.delta*Cons.kappa_m - 6*Cons.delta*Cons.lambda_m - 5*Cons.kappa_p + 6*Cons.lambda_p + Cons.nu*(-2*Cons.kappa_p - 12*Cons.lambda_p + 68) - 12) + Vars.Sigma_ell**3*(-3*Cons.delta*Cons.kappa_p + 2*Cons.delta*Cons.lambda_p + 3*Cons.kappa_m - 2*Cons.lambda_m + Cons.nu*(-2*Cons.delta*Cons.lambda_p + 12*Cons.delta - 6*Cons.kappa_m + 6*Cons.lambda_m))
    Vars.E_SO_3 = (14*Vars.S_ell/3 + 2*Vars.Sigma_ell*Cons.delta)/Cons.M**2
    Vars.E_SO_5 = (Vars.S_ell*(11 - 61*Cons.nu/9) + Vars.Sigma_ell*Cons.delta*(3 - 10*Cons.nu/3))/Cons.M**2
    Vars.E_SO_7 = (Vars.S_ell*(29*Cons.nu**2/12 - 367*Cons.nu/4 + 135/4) + Vars.Sigma_ell*Cons.delta*(5*Cons.nu**2/4 - 39*Cons.nu + 27/4))/Cons.M**2
    Vars.MDot_Alvi_5 = -Cons.M1**3*Vars.chi1_ell*(3*Cons.chi1chi1 + 1)/(4*Cons.M**3) - Cons.M2**3*Vars.chi2_ell*(3*Cons.chi2chi2 + 1)/(4*Cons.M**3)

@njit
def OmegaVec_chiVec_1_4p5(Cons,Vars):
    Omega1_coeff = Vars.v**5/Cons.M
    return Omega1_coeff*(Vars.S2*Vars.v**3*((2.0 - 0.5*Cons.delta)/Cons.M**2 + Vars.v*(-5.0*Vars.S_ell - 3.0*Vars.Sigma_ell*Cons.delta/Cons.M)/Cons.M**4) + Vars.ellHat*(-0.75*Cons.delta + 0.5*Cons.nu + Vars.v**2*(Cons.delta*(0.625*Cons.nu - 0.5625) + Cons.nu*(1.25 - 0.0416666666666667*Cons.nu) + Vars.v*(Vars.v*(Cons.nu*(Cons.nu*(0.229166666666667*Cons.nu - 4.94791666666667) + 2.89583333333333) + 0.96875 + (Cons.delta*(Cons.nu*(5.91666666666667 - 0.489583333333333*Cons.nu) - 0.96875) + (Vars.S_ell*(Vars.S_ell*(Cons.kappa_p*(-0.25*Cons.nu - 0.375) - 0.5*Cons.nu - 0.75) + Vars.Sigma_ell*(Cons.delta*(Cons.kappa_p*(-0.25*Cons.nu - 0.375) - 0.5*Cons.nu - 0.75) + Cons.kappa_m*(0.25*Cons.nu + 0.375))) + Vars.Sigma_ell**2*(Cons.delta*Cons.kappa_m*(0.125*Cons.nu + 0.1875) + Cons.kappa_p*(Cons.nu*(0.25*Cons.nu + 0.25) - 0.1875) + Cons.nu*(0.5*Cons.nu + 0.75)) + (Vars.S_ell*(Vars.S_ell*Cons.delta*(0.375*Cons.kappa_p + 0.75) + Vars.Sigma_ell*Cons.delta*(Cons.delta*(0.375*Cons.kappa_p + 0.75) - 0.375*Cons.kappa_m)) + Vars.Sigma_ell**2*Cons.delta*(-0.1875*Cons.delta*Cons.kappa_m + Cons.kappa_p*(0.1875 - 0.375*Cons.nu) - 0.75*Cons.nu))/Cons.M)/Cons.M**7)/Cons.M) + (Vars.S_ell*(0.833333333333333*Cons.nu + 1.25) + (-1.25*Vars.S_ell*Cons.delta + Vars.Sigma_ell*Cons.delta*(0.5*Cons.nu + 0.75) - 0.75*Vars.Sigma_ell*Cons.delta**2/Cons.M)/Cons.M)/Cons.M**2) + 0.5625) + 0.75) + Vars.lambdaHat*Vars.v**3*(Vars.S2_lambda*(0.5*Cons.delta + 2.0)/Cons.M**2 + Vars.v*(Vars.S1_lambda*(-0.75*Cons.delta - 3.0*Cons.kappa_1 + 2.25) + Vars.S1_lambda*(Cons.delta*(0.75 - 1.5*Cons.kappa_1) + 1.5*Cons.kappa_1 - 0.75)/Cons.nu)/Cons.M**2) + Vars.nHat*Vars.v*(3.0*Vars.chi1_n*Cons.nu + Vars.v**2*((Vars.S1_n*(Cons.delta*(-2.25*Cons.kappa_1 - 0.25) + Cons.kappa_1*(1.5*Cons.nu + 0.75) + 0.5*Cons.nu + 3.25) + Vars.S1_n*(Cons.delta*(1.5 - 0.75*Cons.kappa_1) + 0.75*Cons.kappa_1 - 1.5)/Cons.nu + Vars.S2_n*(1.5*Cons.delta - 1.0*Cons.nu + 3.0))/Cons.M**2 + Vars.v*(-15.0*Vars.S1_n*Vars.S_ell*Cons.kappa_1 + Vars.S1_n*Vars.S_ell*(-7.5*Cons.delta*Cons.kappa_1 + 7.5*Cons.kappa_1)/Cons.nu + 15.0*Vars.S2_n*Vars.S_ell + (-9.0*Vars.S1_n*Vars.Sigma_ell*Cons.delta*Cons.kappa_1 + Vars.S1_n*Vars.Sigma_ell*Cons.delta*(-4.5*Cons.delta*Cons.kappa_1 + 4.5*Cons.kappa_1)/Cons.nu + 9.0*Vars.S2_n*Vars.Sigma_ell*Cons.delta)/Cons.M)/Cons.M**4) + 3.0*Cons.M2**2*Vars.chi2_n/Cons.M**2) - Cons.M2**2*Vars.chiVec2*Vars.v/Cons.M**2)

@njit
def OmegaVec_chiVec_2_4p5(Cons,Vars):
    Omega2_coeff = Vars.v**5/Cons.M
    return Omega2_coeff*(Vars.S1*Vars.v**3*((0.5*Cons.delta + 2.0)/Cons.M**2 + Vars.v*(-5.0*Vars.S_ell + 3.0*Vars.Sigma_ell*Cons.delta/Cons.M)/Cons.M**4) + Vars.ellHat*(0.75*Cons.delta + 0.5*Cons.nu + Vars.v**2*(Cons.delta*(0.5625 - 0.625*Cons.nu) + Cons.nu*(1.25 - 0.0416666666666667*Cons.nu) + Vars.v*(Vars.v*(Cons.nu*(Cons.nu*(0.229166666666667*Cons.nu - 4.94791666666667) + 2.89583333333333) + 0.96875 + (Cons.delta*(Cons.nu*(0.489583333333333*Cons.nu - 5.91666666666667) + 0.96875) + (Vars.S_ell*(Vars.S_ell*(Cons.kappa_p*(-0.25*Cons.nu - 0.375) - 0.5*Cons.nu - 0.75) + Vars.Sigma_ell*(Cons.delta*(Cons.kappa_p*(0.25*Cons.nu + 0.375) + 0.5*Cons.nu + 0.75) + Cons.kappa_m*(-0.25*Cons.nu - 0.375))) + Vars.Sigma_ell**2*(Cons.delta*Cons.kappa_m*(0.125*Cons.nu + 0.1875) + Cons.kappa_p*(Cons.nu*(0.25*Cons.nu + 0.25) - 0.1875) + Cons.nu*(0.5*Cons.nu + 0.75)) + (Vars.S_ell*(Vars.S_ell*Cons.delta*(-0.375*Cons.kappa_p - 0.75) + Vars.Sigma_ell*Cons.delta*(Cons.delta*(0.375*Cons.kappa_p + 0.75) - 0.375*Cons.kappa_m)) + Vars.Sigma_ell**2*Cons.delta*(0.1875*Cons.delta*Cons.kappa_m + Cons.kappa_p*(0.375*Cons.nu - 0.1875) + 0.75*Cons.nu))/Cons.M)/Cons.M**7)/Cons.M) + (Vars.S_ell*(0.833333333333333*Cons.nu + 1.25) + (1.25*Vars.S_ell*Cons.delta + Vars.Sigma_ell*Cons.delta*(-0.5*Cons.nu - 0.75) - 0.75*Vars.Sigma_ell*Cons.delta**2/Cons.M)/Cons.M)/Cons.M**2) + 0.5625) + 0.75) + Vars.lambdaHat*Vars.v**3*(Vars.S1_lambda*(2.0 - 0.5*Cons.delta)/Cons.M**2 + Vars.v*(Vars.S2_lambda*(0.75*Cons.delta - 3.0*Cons.kappa_2 + 2.25) + Vars.S2_lambda*(Cons.delta*(1.5*Cons.kappa_2 - 0.75) + 1.5*Cons.kappa_2 - 0.75)/Cons.nu)/Cons.M**2) + Vars.nHat*Vars.v*(3.0*Vars.chi2_n*Cons.nu + Vars.v**2*((Vars.S1_n*(-1.5*Cons.delta - 1.0*Cons.nu + 3.0) + Vars.S2_n*(Cons.delta*(2.25*Cons.kappa_2 + 0.25) + Cons.kappa_2*(1.5*Cons.nu + 0.75) + 0.5*Cons.nu + 3.25) + Vars.S2_n*(Cons.delta*(0.75*Cons.kappa_2 - 1.5) + 0.75*Cons.kappa_2 - 1.5)/Cons.nu)/Cons.M**2 + Vars.v*(15.0*Vars.S1_n*Vars.S_ell - 15.0*Vars.S2_n*Vars.S_ell*Cons.kappa_2 + Vars.S2_n*Vars.S_ell*(7.5*Cons.delta*Cons.kappa_2 + 7.5*Cons.kappa_2)/Cons.nu + (-9.0*Vars.S1_n*Vars.Sigma_ell*Cons.delta + 9.0*Vars.S2_n*Vars.Sigma_ell*Cons.delta*Cons.kappa_2 + Vars.S2_n*Vars.Sigma_ell*Cons.delta*(-4.5*Cons.delta*Cons.kappa_2 - 4.5*Cons.kappa_2)/Cons.nu)/Cons.M)/Cons.M**4) + 3.0*Cons.M1**2*Vars.chi1_n/Cons.M**2) - Cons.M1**2*Vars.chiVec1*Vars.v/Cons.M**2)

@njit
def OmegaVec_4p5(Cons,Vars):
    gamma_PN_5 = (0.888888888888889*Vars.S_ell*Cons.nu + 3.33333333333333*Vars.S_ell)/Cons.M**2 + 2.0*Vars.Sigma_ell*Cons.delta/Cons.M**3
    gamma_PN_6 = Vars.S_ell**2*(-0.916666666666667*Cons.delta*Cons.kappa_m - 0.916666666666667*Cons.kappa_p + Cons.nu*(-0.166666666666667*Cons.kappa_p - 0.333333333333333) + 1.55555555555556) + Vars.S_ell*Vars.Sigma_ell*(1.66666666666667*Cons.delta + Cons.nu*(-0.166666666666667*Cons.delta*Cons.kappa_p - 0.333333333333333*Cons.delta + 3.83333333333333*Cons.kappa_m)) + Vars.Sigma_ell**2*(Cons.nu**2*(0.166666666666667*Cons.kappa_p + 0.333333333333333) + Cons.nu*(Cons.delta*Cons.kappa_m - Cons.kappa_p - 2.0) + 1.0) + 0.0123456790123457*Cons.nu**3 + 6.36111111111111*Cons.nu**2 - 2.98177812235564*Cons.nu + 1.0
    gamma_PN_0 = 1.00000000000000
    a_ell_4 = Vars.S_n*(5.77777777777778*Cons.nu**2 + 14.75*Cons.nu + 1.5) + Vars.Sigma_n*Cons.delta*(2.83333333333333*Cons.nu**2 + 9.125*Cons.nu + 1.5)
    a_ell_1 = 3.0*Vars.Sigma_ell*Vars.Sigma_n*Cons.nu*(Cons.kappa_p + 2.0)/Cons.M**2 - 3.0*Cons.delta*(Vars.S_ell*Vars.Sigma_n + Vars.S_n*Vars.Sigma_ell)/Cons.M**2 - 0.5*(3.0*Vars.S_n*(2.0*Vars.S_ell*(Cons.kappa_p + 2.0) + Vars.Sigma_ell*(Cons.delta*Cons.kappa_p - Cons.kappa_m)) + 3.0*Vars.Sigma_n*(Vars.S_ell*(Cons.delta*Cons.kappa_p - Cons.kappa_m) + Vars.Sigma_ell*(-Cons.delta*Cons.kappa_m + Cons.kappa_p)))/Cons.M**2
    a_ell_2 = Vars.S_n*(-9.66666666666667*Cons.nu - 10.0) + Vars.Sigma_n*Cons.delta*(-4.5*Cons.nu - 6.0)
    gamma_PN_7 = (Vars.S_ell*(-6.0*Cons.nu**2 - 10.5833333333333*Cons.nu + 5.0) - 2.66666666666667*Vars.Sigma_ell*Cons.delta*Cons.nu**2 + Vars.Sigma_ell*Cons.delta*(3.0 - 10.1666666666667*Cons.nu))/Cons.M**2
    gamma_PN_2 = 1.0 - 0.333333333333333*Cons.nu
    a_ell_3 = Cons.nu*(4.0*Cons.delta*(Vars.S_ell*Vars.Sigma_n + Vars.S_n*Vars.Sigma_ell)/Cons.M**2 + 0.5*(2.0*Vars.S_n*(4.0*Vars.S_ell*(Cons.kappa_p + 2.0) + Vars.Sigma_ell*(2.0*Cons.delta*Cons.kappa_p + 7.0*Cons.kappa_m)) + Vars.Sigma_n*(2.0*Vars.S_ell*(2.0*Cons.delta*Cons.kappa_p + 7.0*Cons.kappa_m) + Vars.Sigma_ell*(5.0*Cons.delta*Cons.kappa_m - 17.0*Cons.kappa_p - 26.0)))/Cons.M**2) - 4.0*Vars.Sigma_ell*Vars.Sigma_n*Cons.nu**2*(Cons.kappa_p + 2.0)/Cons.M**2 + 6.0*Cons.delta*(Vars.S_ell*Vars.Sigma_n + Vars.S_n*Vars.Sigma_ell)/Cons.M**2 + 0.5*(Vars.S_n*(Vars.S_ell*(-9.0*Cons.delta*Cons.kappa_m + 3.0*Cons.kappa_p + 22.0) + Vars.Sigma_ell*(6.0*Cons.delta*Cons.kappa_p - 6.0*Cons.kappa_m)) + 6.0*Vars.Sigma_n*(Vars.S_ell*(Cons.delta*Cons.kappa_p - Cons.kappa_m) + Vars.Sigma_ell*(-Cons.delta*Cons.kappa_m + Cons.kappa_p + 1.0)))/Cons.M**2
    gamma_PN_3 = (1.66666666666667*Vars.S_ell + Vars.Sigma_ell*Cons.delta)/Cons.M**2
    a_ell_0 = 7.0*Vars.S_n + 3.0*Vars.Sigma_n*Cons.delta
    gamma_PN_4 = -5.41666666666667*Cons.nu + 1.0 + (Vars.S_ell**2*(-0.5*Cons.kappa_p - 1.0) + Vars.S_ell*Vars.Sigma_ell*(-0.5*Cons.delta*Cons.kappa_p - Cons.delta + 0.5*Cons.kappa_m) + Vars.Sigma_ell**2*(0.25*Cons.delta*Cons.kappa_m - 0.25*Cons.kappa_p + Cons.nu*(0.5*Cons.kappa_p + 1.0)))/Cons.M**8
    return Vars.ellHat*Vars.v**3/Cons.M + Vars.nHat*Vars.v**6*(a_ell_0 + Vars.v*(a_ell_1 + Vars.v*(a_ell_2 + Vars.v*(a_ell_3 + a_ell_4*Vars.v))))*(gamma_PN_0 + Vars.v**2*(gamma_PN_2 + Vars.v*(gamma_PN_3 + Vars.v*(gamma_PN_4 + Vars.v*(gamma_PN_5 + Vars.v*(gamma_PN_6 + gamma_PN_7*Vars.v))))))/Cons.M**3


@njit(cache=True)
def TaylorT1_4p5(Cons,Vars):
    Flux = Vars.Fcal_coeff*(Cons.Fcal_0 + Vars.v**2*(Cons.Fcal_2 + Vars.v*(Cons.Fcal_3 + Vars.Fcal_SO_3 + Vars.v*(Cons.Fcal_4 + Vars.Fcal_SQ_4 + Vars.v*(Cons.Fcal_5 + Vars.Fcal_SO_5 + Vars.v*(Cons.Fcal_6 + Vars.Fcal_SO_6 + Vars.Fcal_SQ_6 + Cons.Fcal_lnv_6*Vars.logv + Vars.v*(Cons.Fcal_7 + Vars.Fcal_SO_7 + Vars.Fcal_SQ_7 + Vars.v*(Cons.Fcal_8 + Vars.Fcal_SO_8 + Cons.Fcal_lnv_8*Vars.logv + Vars.v*(Cons.Fcal_9 + Cons.Fcal_lnv_9*Vars.logv)))))))))
    dEdV = -0.5*Cons.M*Cons.nu*Vars.v*(2.0*Cons.E_0 + Vars.v**2*(4.0*Cons.E_2 + Vars.v*(5.0*Vars.E_SO_3 + Vars.v*(6.0*Cons.E_4 + 6.0*Vars.E_SQ_4 + Vars.v*(7.0*Vars.E_SO_5 + Vars.v*(8.0*Cons.E_6 + 8.0*Vars.E_SQ_6 + Vars.v*(9.0*Vars.E_SO_7 + 9.0*Vars.E_SQ_7 + Vars.v*(10.0*Cons.E_8 + Cons.E_lnv_8*(10.0*Vars.logv + 1.0)))))))))
    Absorption = Vars.Fcal_coeff*Vars.v**5*(Vars.MDot_Alvi_5 + Cons.MDot_Alvi_8*Vars.v**3)
    dvdt_T1 = (-Absorption - Flux)/dEdV
    dydt=np.zeros(9)
    [dydt[5],dydt[6],dydt[7],dydt[8]] = FrameFromAngularVelocityIntegrand(Vars.R, OmegaVec_4p5(Cons,Vars)[1:])
    dydt[0] = dvdt_T1[0] 
    if(Cons.EvolveSpin1):
        dydt[1], dydt[2]=FrameFromAngularVelocity_2D_Integrand(Vars.rfrak_chi1[0], Vars.rfrak_chi1[1],(mul(mul(inverse(Cons.S_chi1),OmegaVec_chiVec_1_4p5(Cons,Vars)),Cons.S_chi1))[1:])
    else:
        dydt[1], dydt[2] = 0.0, 0.0
    if(Cons.EvolveSpin2):
        dydt[3], dydt[4]=FrameFromAngularVelocity_2D_Integrand(Vars.rfrak_chi2[0], Vars.rfrak_chi2[1],(mul(mul(inverse(Cons.S_chi2),OmegaVec_chiVec_2_4p5(Cons,Vars)),Cons.S_chi2))[1:])
    else:
        dydt[3], dydt[4] = 0.0, 0.0
    return dydt

@njit(cache=True)
def TaylorT4_4p5(Cons,Vars):
    Flux = Vars.Fcal_coeff*(Cons.Fcal_0 + Vars.v**2*(Cons.Fcal_2 + Vars.v*(Cons.Fcal_3 + Vars.Fcal_SO_3 + Vars.v*(Cons.Fcal_4 + Vars.Fcal_SQ_4 + Vars.v*(Cons.Fcal_5 + Vars.Fcal_SO_5 + Vars.v*(Cons.Fcal_6 + Vars.Fcal_SO_6 + Vars.Fcal_SQ_6 + Cons.Fcal_lnv_6*Vars.logv + Vars.v*(Cons.Fcal_7 + Vars.Fcal_SO_7 + Vars.Fcal_SQ_7 + Vars.v*(Cons.Fcal_8 + Vars.Fcal_SO_8 + Cons.Fcal_lnv_8*Vars.logv + Vars.v*(Cons.Fcal_9 + Cons.Fcal_lnv_9*Vars.logv)))))))))
    dEdV = -0.5*Cons.M*Cons.nu*Vars.v*(2.0*Cons.E_0 + Vars.v**2*(4.0*Cons.E_2 + Vars.v*(5.0*Vars.E_SO_3 + Vars.v*(6.0*Cons.E_4 + 6.0*Vars.E_SQ_4 + Vars.v*(7.0*Vars.E_SO_5 + Vars.v*(8.0*Cons.E_6 + 8.0*Vars.E_SQ_6 + Vars.v*(9.0*Vars.E_SO_7 + 9.0*Vars.E_SQ_7 + Vars.v*(10.0*Cons.E_8 + Cons.E_lnv_8*(10.0*Vars.logv + 1.0)))))))))
    Absorption = Vars.Fcal_coeff*Vars.v**5*(Vars.MDot_Alvi_5 + Cons.MDot_Alvi_8*Vars.v**3)
    dvdt_T1 = (-Absorption - Flux)/dEdV
    dydt=np.zeros(9)
    [dydt[5],dydt[6],dydt[7],dydt[8]] = FrameFromAngularVelocityIntegrand(Vars.R, OmegaVec_4p5(Cons,Vars)[1:])
    dydt[0] = dvdt_T4[0]
    if(Cons.EvolveSpin1):
        dydt[1], dydt[2]=FrameFromAngularVelocity_2D_Integrand(Vars.rfrak_chi1[0], Vars.rfrak_chi1[1],(mul(mul(inverse(Cons.S_chi1),OmegaVec_chiVec_1_4p5(Cons,Vars)),Cons.S_chi1))[1:])
    else:
        dydt[1], dydt[2] = 0.0, 0.0
    if(Cons.EvolveSpin2):
        dydt[3], dydt[4]=FrameFromAngularVelocity_2D_Integrand(Vars.rfrak_chi2[0], Vars.rfrak_chi2[1],(mul(mul(inverse(Cons.S_chi2),OmegaVec_chiVec_2_4p5(Cons,Vars)),Cons.S_chi2))[1:])
    else:
        dydt[3], dydt[4] = 0.0, 0.0
    return dydt      

@njit(cache=True)
def TaylorT5_4p5(Cons,Vars):
    Flux = Vars.Fcal_coeff*(Cons.Fcal_0 + Vars.v**2*(Cons.Fcal_2 + Vars.v*(Cons.Fcal_3 + Vars.Fcal_SO_3 + Vars.v*(Cons.Fcal_4 + Vars.Fcal_SQ_4 + Vars.v*(Cons.Fcal_5 + Vars.Fcal_SO_5 + Vars.v*(Cons.Fcal_6 + Vars.Fcal_SO_6 + Vars.Fcal_SQ_6 + Cons.Fcal_lnv_6*Vars.logv + Vars.v*(Cons.Fcal_7 + Vars.Fcal_SO_7 + Vars.Fcal_SQ_7 + Vars.v*(Cons.Fcal_8 + Vars.Fcal_SO_8 + Cons.Fcal_lnv_8*Vars.logv + Vars.v*(Cons.Fcal_9 + Cons.Fcal_lnv_9*Vars.logv)))))))))
    dEdV = -0.5*Cons.M*Cons.nu*Vars.v*(2.0*Cons.E_0 + Vars.v**2*(4.0*Cons.E_2 + Vars.v*(5.0*Vars.E_SO_3 + Vars.v*(6.0*Cons.E_4 + 6.0*Vars.E_SQ_4 + Vars.v*(7.0*Vars.E_SO_5 + Vars.v*(8.0*Cons.E_6 + 8.0*Vars.E_SQ_6 + Vars.v*(9.0*Vars.E_SO_7 + 9.0*Vars.E_SQ_7 + Vars.v*(10.0*Cons.E_8 + Cons.E_lnv_8*(10.0*Vars.logv + 1.0)))))))))
    Absorption = Vars.Fcal_coeff*Vars.v**5*(Vars.MDot_Alvi_5 + Cons.MDot_Alvi_8*Vars.v**3)
    dvdt_T1 = (-Absorption - Flux)/dEdV
    dydt=np.zeros(8)
    [dydt[5],dydt[6],dydt[7],dydt[8]] = FrameFromAngularVelocityIntegrand(Vars.R, OmegaVec_4p5(Cons,Vars)[1:])
    dydt[0] = dvdt_T5[0]
    if(Cons.EvolveSpin1):
        dydt[1], dydt[2]=FrameFromAngularVelocity_2D_Integrand(Vars.rfrak_chi1[0], Vars.rfrak_chi1[1],(mul(mul(inverse(Cons.S_chi1),OmegaVec_chiVec_1_4p5(Cons,Vars)),Cons.S_chi1))[1:])
    else:
        dydt[1], dydt[2] = 0.0, 0.0
    if(Cons.EvolveSpin2):
        dydt[3], dydt[4]=FrameFromAngularVelocity_2D_Integrand(Vars.rfrak_chi2[0], Vars.rfrak_chi2[1],(mul(mul(inverse(Cons.S_chi2),OmegaVec_chiVec_2_4p5(Cons,Vars)),Cons.S_chi2))[1:])
    else:
        dydt[3], dydt[4] = 0.0, 0.0
    return dydt      

@njit(cache=True)
def Recalculate_5p0(Cons,Vars,y):
    Vars.v = np.array([y[0]])
    Vars.rfrak_chi1 = np.array([y[1],y[2]])
    Vars.rfrak_chi2 = np.array([y[3],y[4]])
    Vars.R = y[5]*Cons.wHat + y[6]*Cons.xHat + y[7]*Cons.yHat + y[8]*Cons.zHat
    Vars.nHat = mul(mul(Vars.R,Cons.xHat),inverse(Vars.R))
    Vars.lambdaHat = mul(mul(Vars.R,Cons.yHat),inverse(Vars.R))
    Vars.ellHat = mul(mul(Vars.R,Cons.zHat),inverse(Vars.R))
    Vars.R_S1 = exp(Vars.rfrak_chi1[0]*Cons.xHat + Vars.rfrak_chi1[1]*Cons.yHat)
    Vars.R_S2 = exp(Vars.rfrak_chi2[0]*Cons.xHat + Vars.rfrak_chi2[1]*Cons.yHat)
    Vars.chiVec1 = mul(mul(mul(Cons.S_chi1,Vars.R_S1),Cons.zHat),mul(conjugate(Vars.R_S1),conjugate(Cons.S_chi1)))
    Vars.chiVec2 = mul(mul(mul(Cons.S_chi2,Vars.R_S2),Cons.zHat),mul(conjugate(Vars.R_S2),conjugate(Cons.S_chi2)))
    Vars.chi1_n = np.array([dot(Vars.chiVec1[1:],Vars.nHat[1:])])
    Vars.chi1_lambda = np.array([dot(Vars.chiVec1[1:],Vars.lambdaHat[1:])])
    Vars.chi1_ell = np.array([dot(Vars.chiVec1[1:],Vars.ellHat[1:])])
    Vars.chi2_n = np.array([dot(Vars.chiVec2[1:],Vars.nHat[1:])])
    Vars.chi2_lambda = np.array([dot(Vars.chiVec2[1:],Vars.lambdaHat[1:])])
    Vars.chi2_ell = np.array([dot(Vars.chiVec2[1:],Vars.ellHat[1:])])
    Vars.S = Cons.M1**2*Vars.chiVec1 + Cons.M2**2*Vars.chiVec2
    Vars.S_ell = Cons.M1**2*Vars.chi1_ell + Cons.M2**2*Vars.chi2_ell
    Vars.S_n = Cons.M1**2*Vars.chi1_n + Cons.M2**2*Vars.chi2_n
    Vars.S_lambda = Cons.M1**2*Vars.chi1_lambda + Cons.M2**2*Vars.chi2_lambda
    Vars.Sigma = Cons.M*(-Cons.M1*Vars.chiVec1 + Cons.M2*Vars.chiVec2)
    Vars.Sigma_ell = Cons.M*(-Cons.M1*Vars.chi1_ell + Cons.M2*Vars.chi2_ell)
    Vars.Sigma_n = Cons.M*(-Cons.M1*Vars.chi1_n + Cons.M2*Vars.chi2_n)
    Vars.Sigma_lambda = Cons.M*(-Cons.M1*Vars.chi1_lambda + Cons.M2*Vars.chi2_lambda)
    Vars.SS = np.array([dot(Vars.S[1:],Vars.S[1:])])
    Vars.SigmaSigma = np.array([dot(Vars.Sigma[1:],Vars.Sigma[1:])])
    Vars.SSigma = np.array([dot(Vars.S[1:],Vars.Sigma[1:])])
    Vars.S1 = Cons.M1**2*Vars.chiVec1
    Vars.S1_n = Cons.M1**2*Vars.chi1_n
    Vars.S1_lambda = Cons.M1**2*Vars.chi1_lambda
    Vars.S2 = Cons.M2**2*Vars.chiVec2
    Vars.S2_n = Cons.M2**2*Vars.chi2_n
    Vars.S2_lambda = Cons.M2**2*Vars.chi2_lambda
    Vars.chi_s_ell = Vars.chi1_ell/2 + Vars.chi2_ell/2
    Vars.chi_a_ell = Vars.chi1_ell/2 - Vars.chi2_ell/2
    Vars.logv = log(Vars.v)
    Vars.Fcal_coeff = 32*Cons.nu**2*Vars.v**10/5
    Vars.Fcal_SQ_4 = Cons.chi1chi1*(-89*Cons.delta/192 + 89*Cons.nu/96 - 89/192) - 103*Cons.chi1chi2*Cons.nu/48 + Cons.chi2chi2*(89*Cons.delta/192 + 89*Cons.nu/96 - 89/192) + Vars.chi_a_ell*(Vars.chi_a_ell*(287/96 - 12*Cons.nu) + 287*Vars.chi_s_ell*Cons.delta/48) + Vars.chi_s_ell**2*(Cons.nu/24 + 287/96)
    Vars.Fcal_SQ_6 = Cons.nu*(-0.0714285714285714*Cons.delta*(353*Vars.SSigma + 86*Vars.S_lambda*Vars.Sigma_lambda - 1145*Vars.S_n*Vars.Sigma_n)/Cons.M**4 - 0.00297619047619048*(4236*Vars.SS*Cons.kappa_p + 8472*Vars.SS + 4236*Vars.SSigma*Cons.delta*Cons.kappa_p + 1672*Vars.SSigma*Cons.kappa_m + 1032*Vars.S_lambda**2*(Cons.kappa_p + 2) + 8*Vars.S_lambda*Vars.Sigma_lambda*(129*Cons.delta*Cons.kappa_p + 620*Cons.kappa_m) - 13740*Vars.S_n**2*(Cons.kappa_p + 2) - 4*Vars.S_n*Vars.Sigma_n*(3435*Cons.delta*Cons.kappa_p + 2494*Cons.kappa_m) - 641*Vars.SigmaSigma*Cons.delta*Cons.kappa_m - 2713*Vars.SigmaSigma*Cons.kappa_p - 4146*Vars.SigmaSigma + 982*Vars.Sigma_lambda**2*Cons.delta*Cons.kappa_m - 6944*Vars.Sigma_lambda**2*Cons.kappa_p - 3736*Vars.Sigma_lambda**2 + 941*Vars.Sigma_n**2*Cons.delta*Cons.kappa_m + 15083*Vars.Sigma_n**2*Cons.kappa_p + 44524*Vars.Sigma_n**2)/Cons.M**4) - 0.00595238095238095*Cons.delta*(911*Vars.SSigma + 3643*Vars.S_lambda*Vars.Sigma_lambda - 22399*Vars.S_n*Vars.Sigma_n)/Cons.M**4 + 0.0357142857142857*Cons.nu**2*(Cons.kappa_p + 2)*(353*Vars.SigmaSigma + 86*Vars.Sigma_lambda**2 - 1145*Vars.Sigma_n**2)/Cons.M**4 + 0.00297619047619048*(1477*Vars.SS*Cons.delta*Cons.kappa_m - 1877*Vars.SS*Cons.kappa_p + 502*Vars.SS - 3354*Vars.SSigma*Cons.delta*Cons.kappa_p + 3354*Vars.SSigma*Cons.kappa_m + Vars.S_lambda**2*(1498*Cons.delta*Cons.kappa_m - 4464*Cons.kappa_p - 10811) - 5962*Vars.S_lambda*Vars.Sigma_lambda*(Cons.delta*Cons.kappa_p - Cons.kappa_m) + Vars.S_n**2*(-5929*Cons.delta*Cons.kappa_m + 10095*Cons.kappa_p + 45299) + 16024*Vars.S_n*Vars.Sigma_n*(Cons.delta*Cons.kappa_p - Cons.kappa_m) + 1677*Vars.SigmaSigma*Cons.delta*Cons.kappa_m - 1677*Vars.SigmaSigma*Cons.kappa_p - 168*Vars.SigmaSigma + 2981*Vars.Sigma_lambda**2*Cons.delta*Cons.kappa_m - 2981*Vars.Sigma_lambda**2*Cons.kappa_p + 1313*Vars.Sigma_lambda**2 - 8012*Vars.Sigma_n**2*Cons.delta*Cons.kappa_m + 8012*Vars.Sigma_n**2*Cons.kappa_p + 6387*Vars.Sigma_n**2)/Cons.M**4
    Vars.Fcal_SQ_7 = Vars.S_ell**3*(-16*Cons.kappa_p/3 - 4*Cons.lambda_p + 40/3) + Vars.S_ell**2*Vars.Sigma_ell*(-35*Cons.delta*Cons.kappa_p/6 - 6*Cons.delta*Cons.lambda_p + 73*Cons.delta/3 - 3*Cons.kappa_m/4 + 6*Cons.lambda_m) + Vars.S_ell*Vars.Sigma_ell**2*(-35*Cons.delta*Cons.kappa_m/12 + 6*Cons.delta*Cons.lambda_m + 35*Cons.kappa_p/12 - 6*Cons.lambda_p + Cons.nu*(22*Cons.kappa_p/3 + 12*Cons.lambda_p - 172/3) + 32/3) + Vars.Sigma_ell**3*(67*Cons.delta*Cons.kappa_p/24 - 2*Cons.delta*Cons.lambda_p - Cons.delta/8 - 67*Cons.kappa_m/24 + 2*Cons.lambda_m + Cons.nu*(Cons.delta*Cons.kappa_p/2 + 2*Cons.delta*Cons.lambda_p - 11*Cons.delta + 61*Cons.kappa_m/12 - 6*Cons.lambda_m))
    Vars.Fcal_SO_3 = (-4*Vars.S_ell - 5*Vars.Sigma_ell*Cons.delta/4)/Cons.M**2
    Vars.Fcal_SO_5 = (Vars.S_ell*(272*Cons.nu/9 - 9/2) + Vars.Sigma_ell*Cons.delta*(43*Cons.nu/4 - 13/16))/Cons.M**2
    Vars.Fcal_SO_6 = (-16*Vars.S_ell*pi - 31*Vars.Sigma_ell*Cons.delta*pi/6)/Cons.M**2
    Vars.Fcal_SO_7 = (Vars.S_ell*(-2810*Cons.nu**2/27 + 6172*Cons.nu/189 + 476645/6804) + Vars.Sigma_ell*Cons.delta*(-1501*Cons.nu**2/36 + 1849*Cons.nu/126 + 9535/336))/Cons.M**2
    Vars.Fcal_SO_8 = (Vars.S_ell*pi*(13879*Cons.nu/72 - 3485/96) + Vars.Sigma_ell*Cons.delta*pi*(130583*Cons.nu/2016 - 7163/672))/Cons.M**2
    Vars.E_SQ_4 = -3*Vars.chi_a_ell**2/2 - 3*Vars.chi_s_ell**2/2 - Cons.delta*(Cons.chi2chi2/2 + 3*Vars.chi_a_ell*Vars.chi_s_ell) + Cons.nu*(Cons.chi1chi2 + 6*Vars.chi_a_ell**2) + (Cons.chi1chi1 + Cons.chi2chi2)*(Cons.delta - 2*Cons.nu + 1)/4
    Vars.E_SQ_6 = Cons.nu*(-3*Cons.delta*(Vars.SSigma - 3*Vars.S_n*Vars.Sigma_n)/Cons.M**4 + 0.25*(-6*Vars.SS*Cons.kappa_p - 12*Vars.SS - 6*Vars.SSigma*Cons.delta*Cons.kappa_p - 26*Vars.SSigma*Cons.kappa_m + 24*Vars.S_lambda*Vars.Sigma_lambda*Cons.kappa_m + 18*Vars.S_n**2*(Cons.kappa_p + 2) + 18*Vars.S_n*Vars.Sigma_n*(Cons.delta*Cons.kappa_p + 3*Cons.kappa_m) + 19*Vars.SigmaSigma*Cons.kappa_p - 80*Vars.SigmaSigma - 18*Vars.Sigma_lambda**2*Cons.kappa_p + 76*Vars.Sigma_lambda**2 - 39*Vars.Sigma_n**2*Cons.kappa_p + 168*Vars.Sigma_n**2 + Cons.delta*Cons.kappa_m*(-5*Vars.SigmaSigma + 6*Vars.Sigma_lambda**2 + 9*Vars.Sigma_n**2))/Cons.M**4) + Cons.delta*(15*Vars.SSigma - 13*Vars.S_lambda*Vars.Sigma_lambda - 29*Vars.S_n*Vars.Sigma_n)/Cons.M**4 + 1.5*Cons.nu**2*(Vars.SigmaSigma - 3*Vars.Sigma_n**2)*(Cons.kappa_p + 2)/Cons.M**4 + 0.25*(8*Vars.SS*Cons.delta*Cons.kappa_m - 6*Vars.SS*Cons.kappa_p + 36*Vars.SS - 14*Vars.SSigma*Cons.delta*Cons.kappa_p + 14*Vars.SSigma*Cons.kappa_m + Vars.S_lambda**2*(-6*Cons.delta*Cons.kappa_m + 6*Cons.kappa_p - 28) + 12*Vars.S_lambda*Vars.Sigma_lambda*(Cons.delta*Cons.kappa_p - Cons.kappa_m) - 6*Vars.S_n**2*(3*Cons.delta*Cons.kappa_m - 2*Cons.kappa_p + 10) + 30*Vars.S_n*Vars.Sigma_n*(Cons.delta*Cons.kappa_p - Cons.kappa_m) + 7*Vars.SigmaSigma*Cons.delta*Cons.kappa_m - 7*Vars.SigmaSigma*Cons.kappa_p + 24*Vars.SigmaSigma - 6*Vars.Sigma_lambda**2*Cons.delta*Cons.kappa_m + 6*Vars.Sigma_lambda**2*Cons.kappa_p - 24*Vars.Sigma_lambda**2 - 15*Vars.Sigma_n**2*Cons.delta*Cons.kappa_m + 15*Vars.Sigma_n**2*Cons.kappa_p - 48*Vars.Sigma_n**2)/Cons.M**4
    Vars.E_SQ_7 = Vars.S_ell**3*(2*Cons.kappa_p + 4*Cons.lambda_p - 20) + Vars.S_ell**2*Vars.Sigma_ell*(2*Cons.delta*Cons.kappa_p + 6*Cons.delta*Cons.lambda_p - 32*Cons.delta + 4*Cons.kappa_m - 6*Cons.lambda_m) + Vars.S_ell*Vars.Sigma_ell**2*(5*Cons.delta*Cons.kappa_m - 6*Cons.delta*Cons.lambda_m - 5*Cons.kappa_p + 6*Cons.lambda_p + Cons.nu*(-2*Cons.kappa_p - 12*Cons.lambda_p + 68) - 12) + Vars.Sigma_ell**3*(-3*Cons.delta*Cons.kappa_p + 2*Cons.delta*Cons.lambda_p + 3*Cons.kappa_m - 2*Cons.lambda_m + Cons.nu*(-2*Cons.delta*Cons.lambda_p + 12*Cons.delta - 6*Cons.kappa_m + 6*Cons.lambda_m))
    Vars.E_SO_3 = (14*Vars.S_ell/3 + 2*Vars.Sigma_ell*Cons.delta)/Cons.M**2
    Vars.E_SO_5 = (Vars.S_ell*(11 - 61*Cons.nu/9) + Vars.Sigma_ell*Cons.delta*(3 - 10*Cons.nu/3))/Cons.M**2
    Vars.E_SO_7 = (Vars.S_ell*(29*Cons.nu**2/12 - 367*Cons.nu/4 + 135/4) + Vars.Sigma_ell*Cons.delta*(5*Cons.nu**2/4 - 39*Cons.nu + 27/4))/Cons.M**2
    Vars.MDot_Alvi_5 = -Cons.M1**3*Vars.chi1_ell*(3*Cons.chi1chi1 + 1)/(4*Cons.M**3) - Cons.M2**3*Vars.chi2_ell*(3*Cons.chi2chi2 + 1)/(4*Cons.M**3)

@njit
def OmegaVec_chiVec_1_5p0(Cons,Vars):
    Omega1_coeff = Vars.v**5/Cons.M
    return Omega1_coeff*(Vars.S2*Vars.v**3*((2.0 - 0.5*Cons.delta)/Cons.M**2 + Vars.v*(-5.0*Vars.S_ell - 3.0*Vars.Sigma_ell*Cons.delta/Cons.M)/Cons.M**4) + Vars.ellHat*(-0.75*Cons.delta + 0.5*Cons.nu + Vars.v**2*(Cons.delta*(0.625*Cons.nu - 0.5625) + Cons.nu*(1.25 - 0.0416666666666667*Cons.nu) + Vars.v*(Vars.v*(Cons.nu*(Cons.nu*(0.229166666666667*Cons.nu - 4.94791666666667) + 2.89583333333333) + 0.96875 + (Cons.delta*(Cons.nu*(5.91666666666667 - 0.489583333333333*Cons.nu) - 0.96875) + (Vars.S_ell*(Vars.S_ell*(Cons.kappa_p*(-0.25*Cons.nu - 0.375) - 0.5*Cons.nu - 0.75) + Vars.Sigma_ell*(Cons.delta*(Cons.kappa_p*(-0.25*Cons.nu - 0.375) - 0.5*Cons.nu - 0.75) + Cons.kappa_m*(0.25*Cons.nu + 0.375))) + Vars.Sigma_ell**2*(Cons.delta*Cons.kappa_m*(0.125*Cons.nu + 0.1875) + Cons.kappa_p*(Cons.nu*(0.25*Cons.nu + 0.25) - 0.1875) + Cons.nu*(0.5*Cons.nu + 0.75)) + (Vars.S_ell*(Vars.S_ell*Cons.delta*(0.375*Cons.kappa_p + 0.75) + Vars.Sigma_ell*Cons.delta*(Cons.delta*(0.375*Cons.kappa_p + 0.75) - 0.375*Cons.kappa_m)) + Vars.Sigma_ell**2*Cons.delta*(-0.1875*Cons.delta*Cons.kappa_m + Cons.kappa_p*(0.1875 - 0.375*Cons.nu) - 0.75*Cons.nu))/Cons.M)/Cons.M**7)/Cons.M) + (Vars.S_ell*(0.833333333333333*Cons.nu + 1.25) + (-1.25*Vars.S_ell*Cons.delta + Vars.Sigma_ell*Cons.delta*(0.5*Cons.nu + 0.75) - 0.75*Vars.Sigma_ell*Cons.delta**2/Cons.M)/Cons.M)/Cons.M**2) + 0.5625) + 0.75) + Vars.lambdaHat*Vars.v**3*(Vars.S2_lambda*(0.5*Cons.delta + 2.0)/Cons.M**2 + Vars.v*(Vars.S1_lambda*(-0.75*Cons.delta - 3.0*Cons.kappa_1 + 2.25) + Vars.S1_lambda*(Cons.delta*(0.75 - 1.5*Cons.kappa_1) + 1.5*Cons.kappa_1 - 0.75)/Cons.nu)/Cons.M**2) + Vars.nHat*Vars.v*(3.0*Vars.chi1_n*Cons.nu + Vars.v**2*((Vars.S1_n*(Cons.delta*(-2.25*Cons.kappa_1 - 0.25) + Cons.kappa_1*(1.5*Cons.nu + 0.75) + 0.5*Cons.nu + 3.25) + Vars.S1_n*(Cons.delta*(1.5 - 0.75*Cons.kappa_1) + 0.75*Cons.kappa_1 - 1.5)/Cons.nu + Vars.S2_n*(1.5*Cons.delta - 1.0*Cons.nu + 3.0))/Cons.M**2 + Vars.v*(-15.0*Vars.S1_n*Vars.S_ell*Cons.kappa_1 + Vars.S1_n*Vars.S_ell*(-7.5*Cons.delta*Cons.kappa_1 + 7.5*Cons.kappa_1)/Cons.nu + 15.0*Vars.S2_n*Vars.S_ell + (-9.0*Vars.S1_n*Vars.Sigma_ell*Cons.delta*Cons.kappa_1 + Vars.S1_n*Vars.Sigma_ell*Cons.delta*(-4.5*Cons.delta*Cons.kappa_1 + 4.5*Cons.kappa_1)/Cons.nu + 9.0*Vars.S2_n*Vars.Sigma_ell*Cons.delta)/Cons.M)/Cons.M**4) + 3.0*Cons.M2**2*Vars.chi2_n/Cons.M**2) - Cons.M2**2*Vars.chiVec2*Vars.v/Cons.M**2)

@njit
def OmegaVec_chiVec_2_5p0(Cons,Vars):
    Omega2_coeff = Vars.v**5/Cons.M
    return Omega2_coeff*(Vars.S1*Vars.v**3*((0.5*Cons.delta + 2.0)/Cons.M**2 + Vars.v*(-5.0*Vars.S_ell + 3.0*Vars.Sigma_ell*Cons.delta/Cons.M)/Cons.M**4) + Vars.ellHat*(0.75*Cons.delta + 0.5*Cons.nu + Vars.v**2*(Cons.delta*(0.5625 - 0.625*Cons.nu) + Cons.nu*(1.25 - 0.0416666666666667*Cons.nu) + Vars.v*(Vars.v*(Cons.nu*(Cons.nu*(0.229166666666667*Cons.nu - 4.94791666666667) + 2.89583333333333) + 0.96875 + (Cons.delta*(Cons.nu*(0.489583333333333*Cons.nu - 5.91666666666667) + 0.96875) + (Vars.S_ell*(Vars.S_ell*(Cons.kappa_p*(-0.25*Cons.nu - 0.375) - 0.5*Cons.nu - 0.75) + Vars.Sigma_ell*(Cons.delta*(Cons.kappa_p*(0.25*Cons.nu + 0.375) + 0.5*Cons.nu + 0.75) + Cons.kappa_m*(-0.25*Cons.nu - 0.375))) + Vars.Sigma_ell**2*(Cons.delta*Cons.kappa_m*(0.125*Cons.nu + 0.1875) + Cons.kappa_p*(Cons.nu*(0.25*Cons.nu + 0.25) - 0.1875) + Cons.nu*(0.5*Cons.nu + 0.75)) + (Vars.S_ell*(Vars.S_ell*Cons.delta*(-0.375*Cons.kappa_p - 0.75) + Vars.Sigma_ell*Cons.delta*(Cons.delta*(0.375*Cons.kappa_p + 0.75) - 0.375*Cons.kappa_m)) + Vars.Sigma_ell**2*Cons.delta*(0.1875*Cons.delta*Cons.kappa_m + Cons.kappa_p*(0.375*Cons.nu - 0.1875) + 0.75*Cons.nu))/Cons.M)/Cons.M**7)/Cons.M) + (Vars.S_ell*(0.833333333333333*Cons.nu + 1.25) + (1.25*Vars.S_ell*Cons.delta + Vars.Sigma_ell*Cons.delta*(-0.5*Cons.nu - 0.75) - 0.75*Vars.Sigma_ell*Cons.delta**2/Cons.M)/Cons.M)/Cons.M**2) + 0.5625) + 0.75) + Vars.lambdaHat*Vars.v**3*(Vars.S1_lambda*(2.0 - 0.5*Cons.delta)/Cons.M**2 + Vars.v*(Vars.S2_lambda*(0.75*Cons.delta - 3.0*Cons.kappa_2 + 2.25) + Vars.S2_lambda*(Cons.delta*(1.5*Cons.kappa_2 - 0.75) + 1.5*Cons.kappa_2 - 0.75)/Cons.nu)/Cons.M**2) + Vars.nHat*Vars.v*(3.0*Vars.chi2_n*Cons.nu + Vars.v**2*((Vars.S1_n*(-1.5*Cons.delta - 1.0*Cons.nu + 3.0) + Vars.S2_n*(Cons.delta*(2.25*Cons.kappa_2 + 0.25) + Cons.kappa_2*(1.5*Cons.nu + 0.75) + 0.5*Cons.nu + 3.25) + Vars.S2_n*(Cons.delta*(0.75*Cons.kappa_2 - 1.5) + 0.75*Cons.kappa_2 - 1.5)/Cons.nu)/Cons.M**2 + Vars.v*(15.0*Vars.S1_n*Vars.S_ell - 15.0*Vars.S2_n*Vars.S_ell*Cons.kappa_2 + Vars.S2_n*Vars.S_ell*(7.5*Cons.delta*Cons.kappa_2 + 7.5*Cons.kappa_2)/Cons.nu + (-9.0*Vars.S1_n*Vars.Sigma_ell*Cons.delta + 9.0*Vars.S2_n*Vars.Sigma_ell*Cons.delta*Cons.kappa_2 + Vars.S2_n*Vars.Sigma_ell*Cons.delta*(-4.5*Cons.delta*Cons.kappa_2 - 4.5*Cons.kappa_2)/Cons.nu)/Cons.M)/Cons.M**4) + 3.0*Cons.M1**2*Vars.chi1_n/Cons.M**2) - Cons.M1**2*Vars.chiVec1*Vars.v/Cons.M**2)

@njit
def OmegaVec_5p0(Cons,Vars):
    gamma_PN_5 = (0.888888888888889*Vars.S_ell*Cons.nu + 3.33333333333333*Vars.S_ell)/Cons.M**2 + 2.0*Vars.Sigma_ell*Cons.delta/Cons.M**3
    gamma_PN_6 = Vars.S_ell**2*(-0.916666666666667*Cons.delta*Cons.kappa_m - 0.916666666666667*Cons.kappa_p + Cons.nu*(-0.166666666666667*Cons.kappa_p - 0.333333333333333) + 1.55555555555556) + Vars.S_ell*Vars.Sigma_ell*(1.66666666666667*Cons.delta + Cons.nu*(-0.166666666666667*Cons.delta*Cons.kappa_p - 0.333333333333333*Cons.delta + 3.83333333333333*Cons.kappa_m)) + Vars.Sigma_ell**2*(Cons.nu**2*(0.166666666666667*Cons.kappa_p + 0.333333333333333) + Cons.nu*(Cons.delta*Cons.kappa_m - Cons.kappa_p - 2.0) + 1.0) + 0.0123456790123457*Cons.nu**3 + 6.36111111111111*Cons.nu**2 - 2.98177812235564*Cons.nu + 1.0
    gamma_PN_0 = 1.00000000000000
    a_ell_4 = Vars.S_n*(5.77777777777778*Cons.nu**2 + 14.75*Cons.nu + 1.5) + Vars.Sigma_n*Cons.delta*(2.83333333333333*Cons.nu**2 + 9.125*Cons.nu + 1.5)
    a_ell_1 = 3.0*Vars.Sigma_ell*Vars.Sigma_n*Cons.nu*(Cons.kappa_p + 2.0)/Cons.M**2 - 3.0*Cons.delta*(Vars.S_ell*Vars.Sigma_n + Vars.S_n*Vars.Sigma_ell)/Cons.M**2 - 0.5*(3.0*Vars.S_n*(2.0*Vars.S_ell*(Cons.kappa_p + 2.0) + Vars.Sigma_ell*(Cons.delta*Cons.kappa_p - Cons.kappa_m)) + 3.0*Vars.Sigma_n*(Vars.S_ell*(Cons.delta*Cons.kappa_p - Cons.kappa_m) + Vars.Sigma_ell*(-Cons.delta*Cons.kappa_m + Cons.kappa_p)))/Cons.M**2
    a_ell_2 = Vars.S_n*(-9.66666666666667*Cons.nu - 10.0) + Vars.Sigma_n*Cons.delta*(-4.5*Cons.nu - 6.0)
    gamma_PN_7 = (Vars.S_ell*(-6.0*Cons.nu**2 - 10.5833333333333*Cons.nu + 5.0) - 2.66666666666667*Vars.Sigma_ell*Cons.delta*Cons.nu**2 + Vars.Sigma_ell*Cons.delta*(3.0 - 10.1666666666667*Cons.nu))/Cons.M**2
    gamma_PN_2 = 1.0 - 0.333333333333333*Cons.nu
    a_ell_3 = Cons.nu*(4.0*Cons.delta*(Vars.S_ell*Vars.Sigma_n + Vars.S_n*Vars.Sigma_ell)/Cons.M**2 + 0.5*(2.0*Vars.S_n*(4.0*Vars.S_ell*(Cons.kappa_p + 2.0) + Vars.Sigma_ell*(2.0*Cons.delta*Cons.kappa_p + 7.0*Cons.kappa_m)) + Vars.Sigma_n*(2.0*Vars.S_ell*(2.0*Cons.delta*Cons.kappa_p + 7.0*Cons.kappa_m) + Vars.Sigma_ell*(5.0*Cons.delta*Cons.kappa_m - 17.0*Cons.kappa_p - 26.0)))/Cons.M**2) - 4.0*Vars.Sigma_ell*Vars.Sigma_n*Cons.nu**2*(Cons.kappa_p + 2.0)/Cons.M**2 + 6.0*Cons.delta*(Vars.S_ell*Vars.Sigma_n + Vars.S_n*Vars.Sigma_ell)/Cons.M**2 + 0.5*(Vars.S_n*(Vars.S_ell*(-9.0*Cons.delta*Cons.kappa_m + 3.0*Cons.kappa_p + 22.0) + Vars.Sigma_ell*(6.0*Cons.delta*Cons.kappa_p - 6.0*Cons.kappa_m)) + 6.0*Vars.Sigma_n*(Vars.S_ell*(Cons.delta*Cons.kappa_p - Cons.kappa_m) + Vars.Sigma_ell*(-Cons.delta*Cons.kappa_m + Cons.kappa_p + 1.0)))/Cons.M**2
    gamma_PN_3 = (1.66666666666667*Vars.S_ell + Vars.Sigma_ell*Cons.delta)/Cons.M**2
    a_ell_0 = 7.0*Vars.S_n + 3.0*Vars.Sigma_n*Cons.delta
    gamma_PN_4 = -5.41666666666667*Cons.nu + 1.0 + (Vars.S_ell**2*(-0.5*Cons.kappa_p - 1.0) + Vars.S_ell*Vars.Sigma_ell*(-0.5*Cons.delta*Cons.kappa_p - Cons.delta + 0.5*Cons.kappa_m) + Vars.Sigma_ell**2*(0.25*Cons.delta*Cons.kappa_m - 0.25*Cons.kappa_p + Cons.nu*(0.5*Cons.kappa_p + 1.0)))/Cons.M**8
    return Vars.ellHat*Vars.v**3/Cons.M + Vars.nHat*Vars.v**6*(a_ell_0 + Vars.v*(a_ell_1 + Vars.v*(a_ell_2 + Vars.v*(a_ell_3 + a_ell_4*Vars.v))))*(gamma_PN_0 + Vars.v**2*(gamma_PN_2 + Vars.v*(gamma_PN_3 + Vars.v*(gamma_PN_4 + Vars.v*(gamma_PN_5 + Vars.v*(gamma_PN_6 + gamma_PN_7*Vars.v))))))/Cons.M**3


@njit(cache=True)
def TaylorT1_5p0(Cons,Vars):
    Flux = Vars.Fcal_coeff*(Cons.Fcal_0 + Vars.v**2*(Cons.Fcal_2 + Vars.v*(Cons.Fcal_3 + Vars.Fcal_SO_3 + Vars.v*(Cons.Fcal_4 + Vars.Fcal_SQ_4 + Vars.v*(Cons.Fcal_5 + Vars.Fcal_SO_5 + Vars.v*(Cons.Fcal_6 + Vars.Fcal_SO_6 + Vars.Fcal_SQ_6 + Cons.Fcal_lnv_6*Vars.logv + Vars.v*(Cons.Fcal_7 + Vars.Fcal_SO_7 + Vars.Fcal_SQ_7 + Vars.v*(Cons.Fcal_8 + Vars.Fcal_SO_8 + Cons.Fcal_lnv_8*Vars.logv + Vars.v*(Cons.Fcal_9 + Cons.Fcal_lnv_9*Vars.logv + Vars.v*(Cons.Fcal_10 + Cons.Fcal_lnv_10*Vars.logv))))))))))
    dEdV = -0.5*Cons.M*Cons.nu*Vars.v*(2.0*Cons.E_0 + Vars.v**2*(4.0*Cons.E_2 + Vars.v*(5.0*Vars.E_SO_3 + Vars.v*(6.0*Cons.E_4 + 6.0*Vars.E_SQ_4 + Vars.v*(7.0*Vars.E_SO_5 + Vars.v*(8.0*Cons.E_6 + 8.0*Vars.E_SQ_6 + Vars.v*(9.0*Vars.E_SO_7 + 9.0*Vars.E_SQ_7 + Vars.v*(10.0*Cons.E_8 + Cons.E_lnv_8*(10.0*Vars.logv + 1.0) + Vars.v**2*(12.0*Cons.E_10 + Cons.E_lnv_10*(12.0*Vars.logv + 1.0))))))))))
    Absorption = Vars.Fcal_coeff*Vars.v**5*(Vars.MDot_Alvi_5 + Cons.MDot_Alvi_8*Vars.v**3)
    dvdt_T1 = (-Absorption - Flux)/dEdV
    dydt=np.zeros(9)
    [dydt[5],dydt[6],dydt[7],dydt[8]] = FrameFromAngularVelocityIntegrand(Vars.R, OmegaVec_5p0(Cons,Vars)[1:])
    dydt[0] = dvdt_T1[0] 
    if(Cons.EvolveSpin1):
        dydt[1], dydt[2]=FrameFromAngularVelocity_2D_Integrand(Vars.rfrak_chi1[0], Vars.rfrak_chi1[1],(mul(mul(inverse(Cons.S_chi1),OmegaVec_chiVec_1_5p0(Cons,Vars)),Cons.S_chi1))[1:])
    else:
        dydt[1], dydt[2] = 0.0, 0.0
    if(Cons.EvolveSpin2):
        dydt[3], dydt[4]=FrameFromAngularVelocity_2D_Integrand(Vars.rfrak_chi2[0], Vars.rfrak_chi2[1],(mul(mul(inverse(Cons.S_chi2),OmegaVec_chiVec_2_5p0(Cons,Vars)),Cons.S_chi2))[1:])
    else:
        dydt[3], dydt[4] = 0.0, 0.0
    return dydt

@njit(cache=True)
def TaylorT4_5p0(Cons,Vars):
    Flux = Vars.Fcal_coeff*(Cons.Fcal_0 + Vars.v**2*(Cons.Fcal_2 + Vars.v*(Cons.Fcal_3 + Vars.Fcal_SO_3 + Vars.v*(Cons.Fcal_4 + Vars.Fcal_SQ_4 + Vars.v*(Cons.Fcal_5 + Vars.Fcal_SO_5 + Vars.v*(Cons.Fcal_6 + Vars.Fcal_SO_6 + Vars.Fcal_SQ_6 + Cons.Fcal_lnv_6*Vars.logv + Vars.v*(Cons.Fcal_7 + Vars.Fcal_SO_7 + Vars.Fcal_SQ_7 + Vars.v*(Cons.Fcal_8 + Vars.Fcal_SO_8 + Cons.Fcal_lnv_8*Vars.logv + Vars.v*(Cons.Fcal_9 + Cons.Fcal_lnv_9*Vars.logv + Vars.v*(Cons.Fcal_10 + Cons.Fcal_lnv_10*Vars.logv))))))))))
    dEdV = -0.5*Cons.M*Cons.nu*Vars.v*(2.0*Cons.E_0 + Vars.v**2*(4.0*Cons.E_2 + Vars.v*(5.0*Vars.E_SO_3 + Vars.v*(6.0*Cons.E_4 + 6.0*Vars.E_SQ_4 + Vars.v*(7.0*Vars.E_SO_5 + Vars.v*(8.0*Cons.E_6 + 8.0*Vars.E_SQ_6 + Vars.v*(9.0*Vars.E_SO_7 + 9.0*Vars.E_SQ_7 + Vars.v*(10.0*Cons.E_8 + Cons.E_lnv_8*(10.0*Vars.logv + 1.0) + Vars.v**2*(12.0*Cons.E_10 + Cons.E_lnv_10*(12.0*Vars.logv + 1.0))))))))))
    Absorption = Vars.Fcal_coeff*Vars.v**5*(Vars.MDot_Alvi_5 + Cons.MDot_Alvi_8*Vars.v**3)
    dvdt_T1 = (-Absorption - Flux)/dEdV
    dydt=np.zeros(9)
    [dydt[5],dydt[6],dydt[7],dydt[8]] = FrameFromAngularVelocityIntegrand(Vars.R, OmegaVec_5p0(Cons,Vars)[1:])
    dydt[0] = dvdt_T4[0]
    if(Cons.EvolveSpin1):
        dydt[1], dydt[2]=FrameFromAngularVelocity_2D_Integrand(Vars.rfrak_chi1[0], Vars.rfrak_chi1[1],(mul(mul(inverse(Cons.S_chi1),OmegaVec_chiVec_1_5p0(Cons,Vars)),Cons.S_chi1))[1:])
    else:
        dydt[1], dydt[2] = 0.0, 0.0
    if(Cons.EvolveSpin2):
        dydt[3], dydt[4]=FrameFromAngularVelocity_2D_Integrand(Vars.rfrak_chi2[0], Vars.rfrak_chi2[1],(mul(mul(inverse(Cons.S_chi2),OmegaVec_chiVec_2_5p0(Cons,Vars)),Cons.S_chi2))[1:])
    else:
        dydt[3], dydt[4] = 0.0, 0.0
    return dydt      

@njit(cache=True)
def TaylorT5_5p0(Cons,Vars):
    Flux = Vars.Fcal_coeff*(Cons.Fcal_0 + Vars.v**2*(Cons.Fcal_2 + Vars.v*(Cons.Fcal_3 + Vars.Fcal_SO_3 + Vars.v*(Cons.Fcal_4 + Vars.Fcal_SQ_4 + Vars.v*(Cons.Fcal_5 + Vars.Fcal_SO_5 + Vars.v*(Cons.Fcal_6 + Vars.Fcal_SO_6 + Vars.Fcal_SQ_6 + Cons.Fcal_lnv_6*Vars.logv + Vars.v*(Cons.Fcal_7 + Vars.Fcal_SO_7 + Vars.Fcal_SQ_7 + Vars.v*(Cons.Fcal_8 + Vars.Fcal_SO_8 + Cons.Fcal_lnv_8*Vars.logv + Vars.v*(Cons.Fcal_9 + Cons.Fcal_lnv_9*Vars.logv + Vars.v*(Cons.Fcal_10 + Cons.Fcal_lnv_10*Vars.logv))))))))))
    dEdV = -0.5*Cons.M*Cons.nu*Vars.v*(2.0*Cons.E_0 + Vars.v**2*(4.0*Cons.E_2 + Vars.v*(5.0*Vars.E_SO_3 + Vars.v*(6.0*Cons.E_4 + 6.0*Vars.E_SQ_4 + Vars.v*(7.0*Vars.E_SO_5 + Vars.v*(8.0*Cons.E_6 + 8.0*Vars.E_SQ_6 + Vars.v*(9.0*Vars.E_SO_7 + 9.0*Vars.E_SQ_7 + Vars.v*(10.0*Cons.E_8 + Cons.E_lnv_8*(10.0*Vars.logv + 1.0) + Vars.v**2*(12.0*Cons.E_10 + Cons.E_lnv_10*(12.0*Vars.logv + 1.0))))))))))
    Absorption = Vars.Fcal_coeff*Vars.v**5*(Vars.MDot_Alvi_5 + Cons.MDot_Alvi_8*Vars.v**3)
    dvdt_T1 = (-Absorption - Flux)/dEdV
    dydt=np.zeros(8)
    [dydt[5],dydt[6],dydt[7],dydt[8]] = FrameFromAngularVelocityIntegrand(Vars.R, OmegaVec_5p0(Cons,Vars)[1:])
    dydt[0] = dvdt_T5[0]
    if(Cons.EvolveSpin1):
        dydt[1], dydt[2]=FrameFromAngularVelocity_2D_Integrand(Vars.rfrak_chi1[0], Vars.rfrak_chi1[1],(mul(mul(inverse(Cons.S_chi1),OmegaVec_chiVec_1_5p0(Cons,Vars)),Cons.S_chi1))[1:])
    else:
        dydt[1], dydt[2] = 0.0, 0.0
    if(Cons.EvolveSpin2):
        dydt[3], dydt[4]=FrameFromAngularVelocity_2D_Integrand(Vars.rfrak_chi2[0], Vars.rfrak_chi2[1],(mul(mul(inverse(Cons.S_chi2),OmegaVec_chiVec_2_5p0(Cons,Vars)),Cons.S_chi2))[1:])
    else:
        dydt[3], dydt[4] = 0.0, 0.0
    return dydt      

@njit(cache=True)
def Recalculate_5p5(Cons,Vars,y):
    Vars.v = np.array([y[0]])
    Vars.rfrak_chi1 = np.array([y[1],y[2]])
    Vars.rfrak_chi2 = np.array([y[3],y[4]])
    Vars.R = y[5]*Cons.wHat + y[6]*Cons.xHat + y[7]*Cons.yHat + y[8]*Cons.zHat
    Vars.nHat = mul(mul(Vars.R,Cons.xHat),inverse(Vars.R))
    Vars.lambdaHat = mul(mul(Vars.R,Cons.yHat),inverse(Vars.R))
    Vars.ellHat = mul(mul(Vars.R,Cons.zHat),inverse(Vars.R))
    Vars.R_S1 = exp(Vars.rfrak_chi1[0]*Cons.xHat + Vars.rfrak_chi1[1]*Cons.yHat)
    Vars.R_S2 = exp(Vars.rfrak_chi2[0]*Cons.xHat + Vars.rfrak_chi2[1]*Cons.yHat)
    Vars.chiVec1 = mul(mul(mul(Cons.S_chi1,Vars.R_S1),Cons.zHat),mul(conjugate(Vars.R_S1),conjugate(Cons.S_chi1)))
    Vars.chiVec2 = mul(mul(mul(Cons.S_chi2,Vars.R_S2),Cons.zHat),mul(conjugate(Vars.R_S2),conjugate(Cons.S_chi2)))
    Vars.chi1_n = np.array([dot(Vars.chiVec1[1:],Vars.nHat[1:])])
    Vars.chi1_lambda = np.array([dot(Vars.chiVec1[1:],Vars.lambdaHat[1:])])
    Vars.chi1_ell = np.array([dot(Vars.chiVec1[1:],Vars.ellHat[1:])])
    Vars.chi2_n = np.array([dot(Vars.chiVec2[1:],Vars.nHat[1:])])
    Vars.chi2_lambda = np.array([dot(Vars.chiVec2[1:],Vars.lambdaHat[1:])])
    Vars.chi2_ell = np.array([dot(Vars.chiVec2[1:],Vars.ellHat[1:])])
    Vars.S = Cons.M1**2*Vars.chiVec1 + Cons.M2**2*Vars.chiVec2
    Vars.S_ell = Cons.M1**2*Vars.chi1_ell + Cons.M2**2*Vars.chi2_ell
    Vars.S_n = Cons.M1**2*Vars.chi1_n + Cons.M2**2*Vars.chi2_n
    Vars.S_lambda = Cons.M1**2*Vars.chi1_lambda + Cons.M2**2*Vars.chi2_lambda
    Vars.Sigma = Cons.M*(-Cons.M1*Vars.chiVec1 + Cons.M2*Vars.chiVec2)
    Vars.Sigma_ell = Cons.M*(-Cons.M1*Vars.chi1_ell + Cons.M2*Vars.chi2_ell)
    Vars.Sigma_n = Cons.M*(-Cons.M1*Vars.chi1_n + Cons.M2*Vars.chi2_n)
    Vars.Sigma_lambda = Cons.M*(-Cons.M1*Vars.chi1_lambda + Cons.M2*Vars.chi2_lambda)
    Vars.SS = np.array([dot(Vars.S[1:],Vars.S[1:])])
    Vars.SigmaSigma = np.array([dot(Vars.Sigma[1:],Vars.Sigma[1:])])
    Vars.SSigma = np.array([dot(Vars.S[1:],Vars.Sigma[1:])])
    Vars.S1 = Cons.M1**2*Vars.chiVec1
    Vars.S1_n = Cons.M1**2*Vars.chi1_n
    Vars.S1_lambda = Cons.M1**2*Vars.chi1_lambda
    Vars.S2 = Cons.M2**2*Vars.chiVec2
    Vars.S2_n = Cons.M2**2*Vars.chi2_n
    Vars.S2_lambda = Cons.M2**2*Vars.chi2_lambda
    Vars.chi_s_ell = Vars.chi1_ell/2 + Vars.chi2_ell/2
    Vars.chi_a_ell = Vars.chi1_ell/2 - Vars.chi2_ell/2
    Vars.logv = log(Vars.v)
    Vars.Fcal_coeff = 32*Cons.nu**2*Vars.v**10/5
    Vars.Fcal_SQ_4 = Cons.chi1chi1*(-89*Cons.delta/192 + 89*Cons.nu/96 - 89/192) - 103*Cons.chi1chi2*Cons.nu/48 + Cons.chi2chi2*(89*Cons.delta/192 + 89*Cons.nu/96 - 89/192) + Vars.chi_a_ell*(Vars.chi_a_ell*(287/96 - 12*Cons.nu) + 287*Vars.chi_s_ell*Cons.delta/48) + Vars.chi_s_ell**2*(Cons.nu/24 + 287/96)
    Vars.Fcal_SQ_6 = Cons.nu*(-0.0714285714285714*Cons.delta*(353*Vars.SSigma + 86*Vars.S_lambda*Vars.Sigma_lambda - 1145*Vars.S_n*Vars.Sigma_n)/Cons.M**4 - 0.00297619047619048*(4236*Vars.SS*Cons.kappa_p + 8472*Vars.SS + 4236*Vars.SSigma*Cons.delta*Cons.kappa_p + 1672*Vars.SSigma*Cons.kappa_m + 1032*Vars.S_lambda**2*(Cons.kappa_p + 2) + 8*Vars.S_lambda*Vars.Sigma_lambda*(129*Cons.delta*Cons.kappa_p + 620*Cons.kappa_m) - 13740*Vars.S_n**2*(Cons.kappa_p + 2) - 4*Vars.S_n*Vars.Sigma_n*(3435*Cons.delta*Cons.kappa_p + 2494*Cons.kappa_m) - 641*Vars.SigmaSigma*Cons.delta*Cons.kappa_m - 2713*Vars.SigmaSigma*Cons.kappa_p - 4146*Vars.SigmaSigma + 982*Vars.Sigma_lambda**2*Cons.delta*Cons.kappa_m - 6944*Vars.Sigma_lambda**2*Cons.kappa_p - 3736*Vars.Sigma_lambda**2 + 941*Vars.Sigma_n**2*Cons.delta*Cons.kappa_m + 15083*Vars.Sigma_n**2*Cons.kappa_p + 44524*Vars.Sigma_n**2)/Cons.M**4) - 0.00595238095238095*Cons.delta*(911*Vars.SSigma + 3643*Vars.S_lambda*Vars.Sigma_lambda - 22399*Vars.S_n*Vars.Sigma_n)/Cons.M**4 + 0.0357142857142857*Cons.nu**2*(Cons.kappa_p + 2)*(353*Vars.SigmaSigma + 86*Vars.Sigma_lambda**2 - 1145*Vars.Sigma_n**2)/Cons.M**4 + 0.00297619047619048*(1477*Vars.SS*Cons.delta*Cons.kappa_m - 1877*Vars.SS*Cons.kappa_p + 502*Vars.SS - 3354*Vars.SSigma*Cons.delta*Cons.kappa_p + 3354*Vars.SSigma*Cons.kappa_m + Vars.S_lambda**2*(1498*Cons.delta*Cons.kappa_m - 4464*Cons.kappa_p - 10811) - 5962*Vars.S_lambda*Vars.Sigma_lambda*(Cons.delta*Cons.kappa_p - Cons.kappa_m) + Vars.S_n**2*(-5929*Cons.delta*Cons.kappa_m + 10095*Cons.kappa_p + 45299) + 16024*Vars.S_n*Vars.Sigma_n*(Cons.delta*Cons.kappa_p - Cons.kappa_m) + 1677*Vars.SigmaSigma*Cons.delta*Cons.kappa_m - 1677*Vars.SigmaSigma*Cons.kappa_p - 168*Vars.SigmaSigma + 2981*Vars.Sigma_lambda**2*Cons.delta*Cons.kappa_m - 2981*Vars.Sigma_lambda**2*Cons.kappa_p + 1313*Vars.Sigma_lambda**2 - 8012*Vars.Sigma_n**2*Cons.delta*Cons.kappa_m + 8012*Vars.Sigma_n**2*Cons.kappa_p + 6387*Vars.Sigma_n**2)/Cons.M**4
    Vars.Fcal_SQ_7 = Vars.S_ell**3*(-16*Cons.kappa_p/3 - 4*Cons.lambda_p + 40/3) + Vars.S_ell**2*Vars.Sigma_ell*(-35*Cons.delta*Cons.kappa_p/6 - 6*Cons.delta*Cons.lambda_p + 73*Cons.delta/3 - 3*Cons.kappa_m/4 + 6*Cons.lambda_m) + Vars.S_ell*Vars.Sigma_ell**2*(-35*Cons.delta*Cons.kappa_m/12 + 6*Cons.delta*Cons.lambda_m + 35*Cons.kappa_p/12 - 6*Cons.lambda_p + Cons.nu*(22*Cons.kappa_p/3 + 12*Cons.lambda_p - 172/3) + 32/3) + Vars.Sigma_ell**3*(67*Cons.delta*Cons.kappa_p/24 - 2*Cons.delta*Cons.lambda_p - Cons.delta/8 - 67*Cons.kappa_m/24 + 2*Cons.lambda_m + Cons.nu*(Cons.delta*Cons.kappa_p/2 + 2*Cons.delta*Cons.lambda_p - 11*Cons.delta + 61*Cons.kappa_m/12 - 6*Cons.lambda_m))
    Vars.Fcal_SO_3 = (-4*Vars.S_ell - 5*Vars.Sigma_ell*Cons.delta/4)/Cons.M**2
    Vars.Fcal_SO_5 = (Vars.S_ell*(272*Cons.nu/9 - 9/2) + Vars.Sigma_ell*Cons.delta*(43*Cons.nu/4 - 13/16))/Cons.M**2
    Vars.Fcal_SO_6 = (-16*Vars.S_ell*pi - 31*Vars.Sigma_ell*Cons.delta*pi/6)/Cons.M**2
    Vars.Fcal_SO_7 = (Vars.S_ell*(-2810*Cons.nu**2/27 + 6172*Cons.nu/189 + 476645/6804) + Vars.Sigma_ell*Cons.delta*(-1501*Cons.nu**2/36 + 1849*Cons.nu/126 + 9535/336))/Cons.M**2
    Vars.Fcal_SO_8 = (Vars.S_ell*pi*(13879*Cons.nu/72 - 3485/96) + Vars.Sigma_ell*Cons.delta*pi*(130583*Cons.nu/2016 - 7163/672))/Cons.M**2
    Vars.E_SQ_4 = -3*Vars.chi_a_ell**2/2 - 3*Vars.chi_s_ell**2/2 - Cons.delta*(Cons.chi2chi2/2 + 3*Vars.chi_a_ell*Vars.chi_s_ell) + Cons.nu*(Cons.chi1chi2 + 6*Vars.chi_a_ell**2) + (Cons.chi1chi1 + Cons.chi2chi2)*(Cons.delta - 2*Cons.nu + 1)/4
    Vars.E_SQ_6 = Cons.nu*(-3*Cons.delta*(Vars.SSigma - 3*Vars.S_n*Vars.Sigma_n)/Cons.M**4 + 0.25*(-6*Vars.SS*Cons.kappa_p - 12*Vars.SS - 6*Vars.SSigma*Cons.delta*Cons.kappa_p - 26*Vars.SSigma*Cons.kappa_m + 24*Vars.S_lambda*Vars.Sigma_lambda*Cons.kappa_m + 18*Vars.S_n**2*(Cons.kappa_p + 2) + 18*Vars.S_n*Vars.Sigma_n*(Cons.delta*Cons.kappa_p + 3*Cons.kappa_m) + 19*Vars.SigmaSigma*Cons.kappa_p - 80*Vars.SigmaSigma - 18*Vars.Sigma_lambda**2*Cons.kappa_p + 76*Vars.Sigma_lambda**2 - 39*Vars.Sigma_n**2*Cons.kappa_p + 168*Vars.Sigma_n**2 + Cons.delta*Cons.kappa_m*(-5*Vars.SigmaSigma + 6*Vars.Sigma_lambda**2 + 9*Vars.Sigma_n**2))/Cons.M**4) + Cons.delta*(15*Vars.SSigma - 13*Vars.S_lambda*Vars.Sigma_lambda - 29*Vars.S_n*Vars.Sigma_n)/Cons.M**4 + 1.5*Cons.nu**2*(Vars.SigmaSigma - 3*Vars.Sigma_n**2)*(Cons.kappa_p + 2)/Cons.M**4 + 0.25*(8*Vars.SS*Cons.delta*Cons.kappa_m - 6*Vars.SS*Cons.kappa_p + 36*Vars.SS - 14*Vars.SSigma*Cons.delta*Cons.kappa_p + 14*Vars.SSigma*Cons.kappa_m + Vars.S_lambda**2*(-6*Cons.delta*Cons.kappa_m + 6*Cons.kappa_p - 28) + 12*Vars.S_lambda*Vars.Sigma_lambda*(Cons.delta*Cons.kappa_p - Cons.kappa_m) - 6*Vars.S_n**2*(3*Cons.delta*Cons.kappa_m - 2*Cons.kappa_p + 10) + 30*Vars.S_n*Vars.Sigma_n*(Cons.delta*Cons.kappa_p - Cons.kappa_m) + 7*Vars.SigmaSigma*Cons.delta*Cons.kappa_m - 7*Vars.SigmaSigma*Cons.kappa_p + 24*Vars.SigmaSigma - 6*Vars.Sigma_lambda**2*Cons.delta*Cons.kappa_m + 6*Vars.Sigma_lambda**2*Cons.kappa_p - 24*Vars.Sigma_lambda**2 - 15*Vars.Sigma_n**2*Cons.delta*Cons.kappa_m + 15*Vars.Sigma_n**2*Cons.kappa_p - 48*Vars.Sigma_n**2)/Cons.M**4
    Vars.E_SQ_7 = Vars.S_ell**3*(2*Cons.kappa_p + 4*Cons.lambda_p - 20) + Vars.S_ell**2*Vars.Sigma_ell*(2*Cons.delta*Cons.kappa_p + 6*Cons.delta*Cons.lambda_p - 32*Cons.delta + 4*Cons.kappa_m - 6*Cons.lambda_m) + Vars.S_ell*Vars.Sigma_ell**2*(5*Cons.delta*Cons.kappa_m - 6*Cons.delta*Cons.lambda_m - 5*Cons.kappa_p + 6*Cons.lambda_p + Cons.nu*(-2*Cons.kappa_p - 12*Cons.lambda_p + 68) - 12) + Vars.Sigma_ell**3*(-3*Cons.delta*Cons.kappa_p + 2*Cons.delta*Cons.lambda_p + 3*Cons.kappa_m - 2*Cons.lambda_m + Cons.nu*(-2*Cons.delta*Cons.lambda_p + 12*Cons.delta - 6*Cons.kappa_m + 6*Cons.lambda_m))
    Vars.E_SO_3 = (14*Vars.S_ell/3 + 2*Vars.Sigma_ell*Cons.delta)/Cons.M**2
    Vars.E_SO_5 = (Vars.S_ell*(11 - 61*Cons.nu/9) + Vars.Sigma_ell*Cons.delta*(3 - 10*Cons.nu/3))/Cons.M**2
    Vars.E_SO_7 = (Vars.S_ell*(29*Cons.nu**2/12 - 367*Cons.nu/4 + 135/4) + Vars.Sigma_ell*Cons.delta*(5*Cons.nu**2/4 - 39*Cons.nu + 27/4))/Cons.M**2
    Vars.MDot_Alvi_5 = -Cons.M1**3*Vars.chi1_ell*(3*Cons.chi1chi1 + 1)/(4*Cons.M**3) - Cons.M2**3*Vars.chi2_ell*(3*Cons.chi2chi2 + 1)/(4*Cons.M**3)

@njit
def OmegaVec_chiVec_1_5p5(Cons,Vars):
    Omega1_coeff = Vars.v**5/Cons.M
    return Omega1_coeff*(Vars.S2*Vars.v**3*((2.0 - 0.5*Cons.delta)/Cons.M**2 + Vars.v*(-5.0*Vars.S_ell - 3.0*Vars.Sigma_ell*Cons.delta/Cons.M)/Cons.M**4) + Vars.ellHat*(-0.75*Cons.delta + 0.5*Cons.nu + Vars.v**2*(Cons.delta*(0.625*Cons.nu - 0.5625) + Cons.nu*(1.25 - 0.0416666666666667*Cons.nu) + Vars.v*(Vars.v*(Cons.nu*(Cons.nu*(0.229166666666667*Cons.nu - 4.94791666666667) + 2.89583333333333) + 0.96875 + (Cons.delta*(Cons.nu*(5.91666666666667 - 0.489583333333333*Cons.nu) - 0.96875) + (Vars.S_ell*(Vars.S_ell*(Cons.kappa_p*(-0.25*Cons.nu - 0.375) - 0.5*Cons.nu - 0.75) + Vars.Sigma_ell*(Cons.delta*(Cons.kappa_p*(-0.25*Cons.nu - 0.375) - 0.5*Cons.nu - 0.75) + Cons.kappa_m*(0.25*Cons.nu + 0.375))) + Vars.Sigma_ell**2*(Cons.delta*Cons.kappa_m*(0.125*Cons.nu + 0.1875) + Cons.kappa_p*(Cons.nu*(0.25*Cons.nu + 0.25) - 0.1875) + Cons.nu*(0.5*Cons.nu + 0.75)) + (Vars.S_ell*(Vars.S_ell*Cons.delta*(0.375*Cons.kappa_p + 0.75) + Vars.Sigma_ell*Cons.delta*(Cons.delta*(0.375*Cons.kappa_p + 0.75) - 0.375*Cons.kappa_m)) + Vars.Sigma_ell**2*Cons.delta*(-0.1875*Cons.delta*Cons.kappa_m + Cons.kappa_p*(0.1875 - 0.375*Cons.nu) - 0.75*Cons.nu))/Cons.M)/Cons.M**7)/Cons.M) + (Vars.S_ell*(0.833333333333333*Cons.nu + 1.25) + (-1.25*Vars.S_ell*Cons.delta + Vars.Sigma_ell*Cons.delta*(0.5*Cons.nu + 0.75) - 0.75*Vars.Sigma_ell*Cons.delta**2/Cons.M)/Cons.M)/Cons.M**2) + 0.5625) + 0.75) + Vars.lambdaHat*Vars.v**3*(Vars.S2_lambda*(0.5*Cons.delta + 2.0)/Cons.M**2 + Vars.v*(Vars.S1_lambda*(-0.75*Cons.delta - 3.0*Cons.kappa_1 + 2.25) + Vars.S1_lambda*(Cons.delta*(0.75 - 1.5*Cons.kappa_1) + 1.5*Cons.kappa_1 - 0.75)/Cons.nu)/Cons.M**2) + Vars.nHat*Vars.v*(3.0*Vars.chi1_n*Cons.nu + Vars.v**2*((Vars.S1_n*(Cons.delta*(-2.25*Cons.kappa_1 - 0.25) + Cons.kappa_1*(1.5*Cons.nu + 0.75) + 0.5*Cons.nu + 3.25) + Vars.S1_n*(Cons.delta*(1.5 - 0.75*Cons.kappa_1) + 0.75*Cons.kappa_1 - 1.5)/Cons.nu + Vars.S2_n*(1.5*Cons.delta - 1.0*Cons.nu + 3.0))/Cons.M**2 + Vars.v*(-15.0*Vars.S1_n*Vars.S_ell*Cons.kappa_1 + Vars.S1_n*Vars.S_ell*(-7.5*Cons.delta*Cons.kappa_1 + 7.5*Cons.kappa_1)/Cons.nu + 15.0*Vars.S2_n*Vars.S_ell + (-9.0*Vars.S1_n*Vars.Sigma_ell*Cons.delta*Cons.kappa_1 + Vars.S1_n*Vars.Sigma_ell*Cons.delta*(-4.5*Cons.delta*Cons.kappa_1 + 4.5*Cons.kappa_1)/Cons.nu + 9.0*Vars.S2_n*Vars.Sigma_ell*Cons.delta)/Cons.M)/Cons.M**4) + 3.0*Cons.M2**2*Vars.chi2_n/Cons.M**2) - Cons.M2**2*Vars.chiVec2*Vars.v/Cons.M**2)

@njit
def OmegaVec_chiVec_2_5p5(Cons,Vars):
    Omega2_coeff = Vars.v**5/Cons.M
    return Omega2_coeff*(Vars.S1*Vars.v**3*((0.5*Cons.delta + 2.0)/Cons.M**2 + Vars.v*(-5.0*Vars.S_ell + 3.0*Vars.Sigma_ell*Cons.delta/Cons.M)/Cons.M**4) + Vars.ellHat*(0.75*Cons.delta + 0.5*Cons.nu + Vars.v**2*(Cons.delta*(0.5625 - 0.625*Cons.nu) + Cons.nu*(1.25 - 0.0416666666666667*Cons.nu) + Vars.v*(Vars.v*(Cons.nu*(Cons.nu*(0.229166666666667*Cons.nu - 4.94791666666667) + 2.89583333333333) + 0.96875 + (Cons.delta*(Cons.nu*(0.489583333333333*Cons.nu - 5.91666666666667) + 0.96875) + (Vars.S_ell*(Vars.S_ell*(Cons.kappa_p*(-0.25*Cons.nu - 0.375) - 0.5*Cons.nu - 0.75) + Vars.Sigma_ell*(Cons.delta*(Cons.kappa_p*(0.25*Cons.nu + 0.375) + 0.5*Cons.nu + 0.75) + Cons.kappa_m*(-0.25*Cons.nu - 0.375))) + Vars.Sigma_ell**2*(Cons.delta*Cons.kappa_m*(0.125*Cons.nu + 0.1875) + Cons.kappa_p*(Cons.nu*(0.25*Cons.nu + 0.25) - 0.1875) + Cons.nu*(0.5*Cons.nu + 0.75)) + (Vars.S_ell*(Vars.S_ell*Cons.delta*(-0.375*Cons.kappa_p - 0.75) + Vars.Sigma_ell*Cons.delta*(Cons.delta*(0.375*Cons.kappa_p + 0.75) - 0.375*Cons.kappa_m)) + Vars.Sigma_ell**2*Cons.delta*(0.1875*Cons.delta*Cons.kappa_m + Cons.kappa_p*(0.375*Cons.nu - 0.1875) + 0.75*Cons.nu))/Cons.M)/Cons.M**7)/Cons.M) + (Vars.S_ell*(0.833333333333333*Cons.nu + 1.25) + (1.25*Vars.S_ell*Cons.delta + Vars.Sigma_ell*Cons.delta*(-0.5*Cons.nu - 0.75) - 0.75*Vars.Sigma_ell*Cons.delta**2/Cons.M)/Cons.M)/Cons.M**2) + 0.5625) + 0.75) + Vars.lambdaHat*Vars.v**3*(Vars.S1_lambda*(2.0 - 0.5*Cons.delta)/Cons.M**2 + Vars.v*(Vars.S2_lambda*(0.75*Cons.delta - 3.0*Cons.kappa_2 + 2.25) + Vars.S2_lambda*(Cons.delta*(1.5*Cons.kappa_2 - 0.75) + 1.5*Cons.kappa_2 - 0.75)/Cons.nu)/Cons.M**2) + Vars.nHat*Vars.v*(3.0*Vars.chi2_n*Cons.nu + Vars.v**2*((Vars.S1_n*(-1.5*Cons.delta - 1.0*Cons.nu + 3.0) + Vars.S2_n*(Cons.delta*(2.25*Cons.kappa_2 + 0.25) + Cons.kappa_2*(1.5*Cons.nu + 0.75) + 0.5*Cons.nu + 3.25) + Vars.S2_n*(Cons.delta*(0.75*Cons.kappa_2 - 1.5) + 0.75*Cons.kappa_2 - 1.5)/Cons.nu)/Cons.M**2 + Vars.v*(15.0*Vars.S1_n*Vars.S_ell - 15.0*Vars.S2_n*Vars.S_ell*Cons.kappa_2 + Vars.S2_n*Vars.S_ell*(7.5*Cons.delta*Cons.kappa_2 + 7.5*Cons.kappa_2)/Cons.nu + (-9.0*Vars.S1_n*Vars.Sigma_ell*Cons.delta + 9.0*Vars.S2_n*Vars.Sigma_ell*Cons.delta*Cons.kappa_2 + Vars.S2_n*Vars.Sigma_ell*Cons.delta*(-4.5*Cons.delta*Cons.kappa_2 - 4.5*Cons.kappa_2)/Cons.nu)/Cons.M)/Cons.M**4) + 3.0*Cons.M1**2*Vars.chi1_n/Cons.M**2) - Cons.M1**2*Vars.chiVec1*Vars.v/Cons.M**2)

@njit
def OmegaVec_5p5(Cons,Vars):
    gamma_PN_5 = (0.888888888888889*Vars.S_ell*Cons.nu + 3.33333333333333*Vars.S_ell)/Cons.M**2 + 2.0*Vars.Sigma_ell*Cons.delta/Cons.M**3
    gamma_PN_6 = Vars.S_ell**2*(-0.916666666666667*Cons.delta*Cons.kappa_m - 0.916666666666667*Cons.kappa_p + Cons.nu*(-0.166666666666667*Cons.kappa_p - 0.333333333333333) + 1.55555555555556) + Vars.S_ell*Vars.Sigma_ell*(1.66666666666667*Cons.delta + Cons.nu*(-0.166666666666667*Cons.delta*Cons.kappa_p - 0.333333333333333*Cons.delta + 3.83333333333333*Cons.kappa_m)) + Vars.Sigma_ell**2*(Cons.nu**2*(0.166666666666667*Cons.kappa_p + 0.333333333333333) + Cons.nu*(Cons.delta*Cons.kappa_m - Cons.kappa_p - 2.0) + 1.0) + 0.0123456790123457*Cons.nu**3 + 6.36111111111111*Cons.nu**2 - 2.98177812235564*Cons.nu + 1.0
    gamma_PN_0 = 1.00000000000000
    a_ell_4 = Vars.S_n*(5.77777777777778*Cons.nu**2 + 14.75*Cons.nu + 1.5) + Vars.Sigma_n*Cons.delta*(2.83333333333333*Cons.nu**2 + 9.125*Cons.nu + 1.5)
    a_ell_1 = 3.0*Vars.Sigma_ell*Vars.Sigma_n*Cons.nu*(Cons.kappa_p + 2.0)/Cons.M**2 - 3.0*Cons.delta*(Vars.S_ell*Vars.Sigma_n + Vars.S_n*Vars.Sigma_ell)/Cons.M**2 - 0.5*(3.0*Vars.S_n*(2.0*Vars.S_ell*(Cons.kappa_p + 2.0) + Vars.Sigma_ell*(Cons.delta*Cons.kappa_p - Cons.kappa_m)) + 3.0*Vars.Sigma_n*(Vars.S_ell*(Cons.delta*Cons.kappa_p - Cons.kappa_m) + Vars.Sigma_ell*(-Cons.delta*Cons.kappa_m + Cons.kappa_p)))/Cons.M**2
    a_ell_2 = Vars.S_n*(-9.66666666666667*Cons.nu - 10.0) + Vars.Sigma_n*Cons.delta*(-4.5*Cons.nu - 6.0)
    gamma_PN_7 = (Vars.S_ell*(-6.0*Cons.nu**2 - 10.5833333333333*Cons.nu + 5.0) - 2.66666666666667*Vars.Sigma_ell*Cons.delta*Cons.nu**2 + Vars.Sigma_ell*Cons.delta*(3.0 - 10.1666666666667*Cons.nu))/Cons.M**2
    gamma_PN_2 = 1.0 - 0.333333333333333*Cons.nu
    a_ell_3 = Cons.nu*(4.0*Cons.delta*(Vars.S_ell*Vars.Sigma_n + Vars.S_n*Vars.Sigma_ell)/Cons.M**2 + 0.5*(2.0*Vars.S_n*(4.0*Vars.S_ell*(Cons.kappa_p + 2.0) + Vars.Sigma_ell*(2.0*Cons.delta*Cons.kappa_p + 7.0*Cons.kappa_m)) + Vars.Sigma_n*(2.0*Vars.S_ell*(2.0*Cons.delta*Cons.kappa_p + 7.0*Cons.kappa_m) + Vars.Sigma_ell*(5.0*Cons.delta*Cons.kappa_m - 17.0*Cons.kappa_p - 26.0)))/Cons.M**2) - 4.0*Vars.Sigma_ell*Vars.Sigma_n*Cons.nu**2*(Cons.kappa_p + 2.0)/Cons.M**2 + 6.0*Cons.delta*(Vars.S_ell*Vars.Sigma_n + Vars.S_n*Vars.Sigma_ell)/Cons.M**2 + 0.5*(Vars.S_n*(Vars.S_ell*(-9.0*Cons.delta*Cons.kappa_m + 3.0*Cons.kappa_p + 22.0) + Vars.Sigma_ell*(6.0*Cons.delta*Cons.kappa_p - 6.0*Cons.kappa_m)) + 6.0*Vars.Sigma_n*(Vars.S_ell*(Cons.delta*Cons.kappa_p - Cons.kappa_m) + Vars.Sigma_ell*(-Cons.delta*Cons.kappa_m + Cons.kappa_p + 1.0)))/Cons.M**2
    gamma_PN_3 = (1.66666666666667*Vars.S_ell + Vars.Sigma_ell*Cons.delta)/Cons.M**2
    a_ell_0 = 7.0*Vars.S_n + 3.0*Vars.Sigma_n*Cons.delta
    gamma_PN_4 = -5.41666666666667*Cons.nu + 1.0 + (Vars.S_ell**2*(-0.5*Cons.kappa_p - 1.0) + Vars.S_ell*Vars.Sigma_ell*(-0.5*Cons.delta*Cons.kappa_p - Cons.delta + 0.5*Cons.kappa_m) + Vars.Sigma_ell**2*(0.25*Cons.delta*Cons.kappa_m - 0.25*Cons.kappa_p + Cons.nu*(0.5*Cons.kappa_p + 1.0)))/Cons.M**8
    return Vars.ellHat*Vars.v**3/Cons.M + Vars.nHat*Vars.v**6*(a_ell_0 + Vars.v*(a_ell_1 + Vars.v*(a_ell_2 + Vars.v*(a_ell_3 + a_ell_4*Vars.v))))*(gamma_PN_0 + Vars.v**2*(gamma_PN_2 + Vars.v*(gamma_PN_3 + Vars.v*(gamma_PN_4 + Vars.v*(gamma_PN_5 + Vars.v*(gamma_PN_6 + gamma_PN_7*Vars.v))))))/Cons.M**3


@njit(cache=True)
def TaylorT1_5p5(Cons,Vars):
    Flux = Vars.Fcal_coeff*(Cons.Fcal_0 + Vars.v**2*(Cons.Fcal_2 + Vars.v*(Cons.Fcal_3 + Vars.Fcal_SO_3 + Vars.v*(Cons.Fcal_4 + Vars.Fcal_SQ_4 + Vars.v*(Cons.Fcal_5 + Vars.Fcal_SO_5 + Vars.v*(Cons.Fcal_6 + Vars.Fcal_SO_6 + Vars.Fcal_SQ_6 + Cons.Fcal_lnv_6*Vars.logv + Vars.v*(Cons.Fcal_7 + Vars.Fcal_SO_7 + Vars.Fcal_SQ_7 + Vars.v*(Cons.Fcal_8 + Vars.Fcal_SO_8 + Cons.Fcal_lnv_8*Vars.logv + Vars.v*(Cons.Fcal_9 + Cons.Fcal_lnv_9*Vars.logv + Vars.v*(Cons.Fcal_10 + Cons.Fcal_lnv_10*Vars.logv + Vars.v*(Cons.Fcal_11 + Cons.Fcal_lnv_11*Vars.logv)))))))))))
    dEdV = -0.5*Cons.M*Cons.nu*Vars.v*(2.0*Cons.E_0 + Vars.v**2*(4.0*Cons.E_2 + Vars.v*(5.0*Vars.E_SO_3 + Vars.v*(6.0*Cons.E_4 + 6.0*Vars.E_SQ_4 + Vars.v*(7.0*Vars.E_SO_5 + Vars.v*(8.0*Cons.E_6 + 8.0*Vars.E_SQ_6 + Vars.v*(9.0*Vars.E_SO_7 + 9.0*Vars.E_SQ_7 + Vars.v*(10.0*Cons.E_8 + Cons.E_lnv_8*(10.0*Vars.logv + 1.0) + Vars.v**2*(12.0*Cons.E_10 + 13.0*Cons.E_11*Vars.v + Cons.E_lnv_10*(12.0*Vars.logv + 1.0))))))))))
    Absorption = Vars.Fcal_coeff*Vars.v**5*(Vars.MDot_Alvi_5 + Cons.MDot_Alvi_8*Vars.v**3)
    dvdt_T1 = (-Absorption - Flux)/dEdV
    dydt=np.zeros(9)
    [dydt[5],dydt[6],dydt[7],dydt[8]] = FrameFromAngularVelocityIntegrand(Vars.R, OmegaVec_5p5(Cons,Vars)[1:])
    dydt[0] = dvdt_T1[0] 
    if(Cons.EvolveSpin1):
        dydt[1], dydt[2]=FrameFromAngularVelocity_2D_Integrand(Vars.rfrak_chi1[0], Vars.rfrak_chi1[1],(mul(mul(inverse(Cons.S_chi1),OmegaVec_chiVec_1_5p5(Cons,Vars)),Cons.S_chi1))[1:])
    else:
        dydt[1], dydt[2] = 0.0, 0.0
    if(Cons.EvolveSpin2):
        dydt[3], dydt[4]=FrameFromAngularVelocity_2D_Integrand(Vars.rfrak_chi2[0], Vars.rfrak_chi2[1],(mul(mul(inverse(Cons.S_chi2),OmegaVec_chiVec_2_5p5(Cons,Vars)),Cons.S_chi2))[1:])
    else:
        dydt[3], dydt[4] = 0.0, 0.0
    return dydt

@njit(cache=True)
def TaylorT4_5p5(Cons,Vars):
    Flux = Vars.Fcal_coeff*(Cons.Fcal_0 + Vars.v**2*(Cons.Fcal_2 + Vars.v*(Cons.Fcal_3 + Vars.Fcal_SO_3 + Vars.v*(Cons.Fcal_4 + Vars.Fcal_SQ_4 + Vars.v*(Cons.Fcal_5 + Vars.Fcal_SO_5 + Vars.v*(Cons.Fcal_6 + Vars.Fcal_SO_6 + Vars.Fcal_SQ_6 + Cons.Fcal_lnv_6*Vars.logv + Vars.v*(Cons.Fcal_7 + Vars.Fcal_SO_7 + Vars.Fcal_SQ_7 + Vars.v*(Cons.Fcal_8 + Vars.Fcal_SO_8 + Cons.Fcal_lnv_8*Vars.logv + Vars.v*(Cons.Fcal_9 + Cons.Fcal_lnv_9*Vars.logv + Vars.v*(Cons.Fcal_10 + Cons.Fcal_lnv_10*Vars.logv + Vars.v*(Cons.Fcal_11 + Cons.Fcal_lnv_11*Vars.logv)))))))))))
    dEdV = -0.5*Cons.M*Cons.nu*Vars.v*(2.0*Cons.E_0 + Vars.v**2*(4.0*Cons.E_2 + Vars.v*(5.0*Vars.E_SO_3 + Vars.v*(6.0*Cons.E_4 + 6.0*Vars.E_SQ_4 + Vars.v*(7.0*Vars.E_SO_5 + Vars.v*(8.0*Cons.E_6 + 8.0*Vars.E_SQ_6 + Vars.v*(9.0*Vars.E_SO_7 + 9.0*Vars.E_SQ_7 + Vars.v*(10.0*Cons.E_8 + Cons.E_lnv_8*(10.0*Vars.logv + 1.0) + Vars.v**2*(12.0*Cons.E_10 + 13.0*Cons.E_11*Vars.v + Cons.E_lnv_10*(12.0*Vars.logv + 1.0))))))))))
    Absorption = Vars.Fcal_coeff*Vars.v**5*(Vars.MDot_Alvi_5 + Cons.MDot_Alvi_8*Vars.v**3)
    dvdt_T1 = (-Absorption - Flux)/dEdV
    dydt=np.zeros(9)
    [dydt[5],dydt[6],dydt[7],dydt[8]] = FrameFromAngularVelocityIntegrand(Vars.R, OmegaVec_5p5(Cons,Vars)[1:])
    dydt[0] = dvdt_T4[0]
    if(Cons.EvolveSpin1):
        dydt[1], dydt[2]=FrameFromAngularVelocity_2D_Integrand(Vars.rfrak_chi1[0], Vars.rfrak_chi1[1],(mul(mul(inverse(Cons.S_chi1),OmegaVec_chiVec_1_5p5(Cons,Vars)),Cons.S_chi1))[1:])
    else:
        dydt[1], dydt[2] = 0.0, 0.0
    if(Cons.EvolveSpin2):
        dydt[3], dydt[4]=FrameFromAngularVelocity_2D_Integrand(Vars.rfrak_chi2[0], Vars.rfrak_chi2[1],(mul(mul(inverse(Cons.S_chi2),OmegaVec_chiVec_2_5p5(Cons,Vars)),Cons.S_chi2))[1:])
    else:
        dydt[3], dydt[4] = 0.0, 0.0
    return dydt      

@njit(cache=True)
def TaylorT5_5p5(Cons,Vars):
    Flux = Vars.Fcal_coeff*(Cons.Fcal_0 + Vars.v**2*(Cons.Fcal_2 + Vars.v*(Cons.Fcal_3 + Vars.Fcal_SO_3 + Vars.v*(Cons.Fcal_4 + Vars.Fcal_SQ_4 + Vars.v*(Cons.Fcal_5 + Vars.Fcal_SO_5 + Vars.v*(Cons.Fcal_6 + Vars.Fcal_SO_6 + Vars.Fcal_SQ_6 + Cons.Fcal_lnv_6*Vars.logv + Vars.v*(Cons.Fcal_7 + Vars.Fcal_SO_7 + Vars.Fcal_SQ_7 + Vars.v*(Cons.Fcal_8 + Vars.Fcal_SO_8 + Cons.Fcal_lnv_8*Vars.logv + Vars.v*(Cons.Fcal_9 + Cons.Fcal_lnv_9*Vars.logv + Vars.v*(Cons.Fcal_10 + Cons.Fcal_lnv_10*Vars.logv + Vars.v*(Cons.Fcal_11 + Cons.Fcal_lnv_11*Vars.logv)))))))))))
    dEdV = -0.5*Cons.M*Cons.nu*Vars.v*(2.0*Cons.E_0 + Vars.v**2*(4.0*Cons.E_2 + Vars.v*(5.0*Vars.E_SO_3 + Vars.v*(6.0*Cons.E_4 + 6.0*Vars.E_SQ_4 + Vars.v*(7.0*Vars.E_SO_5 + Vars.v*(8.0*Cons.E_6 + 8.0*Vars.E_SQ_6 + Vars.v*(9.0*Vars.E_SO_7 + 9.0*Vars.E_SQ_7 + Vars.v*(10.0*Cons.E_8 + Cons.E_lnv_8*(10.0*Vars.logv + 1.0) + Vars.v**2*(12.0*Cons.E_10 + 13.0*Cons.E_11*Vars.v + Cons.E_lnv_10*(12.0*Vars.logv + 1.0))))))))))
    Absorption = Vars.Fcal_coeff*Vars.v**5*(Vars.MDot_Alvi_5 + Cons.MDot_Alvi_8*Vars.v**3)
    dvdt_T1 = (-Absorption - Flux)/dEdV
    dydt=np.zeros(8)
    [dydt[5],dydt[6],dydt[7],dydt[8]] = FrameFromAngularVelocityIntegrand(Vars.R, OmegaVec_5p5(Cons,Vars)[1:])
    dydt[0] = dvdt_T5[0]
    if(Cons.EvolveSpin1):
        dydt[1], dydt[2]=FrameFromAngularVelocity_2D_Integrand(Vars.rfrak_chi1[0], Vars.rfrak_chi1[1],(mul(mul(inverse(Cons.S_chi1),OmegaVec_chiVec_1_5p5(Cons,Vars)),Cons.S_chi1))[1:])
    else:
        dydt[1], dydt[2] = 0.0, 0.0
    if(Cons.EvolveSpin2):
        dydt[3], dydt[4]=FrameFromAngularVelocity_2D_Integrand(Vars.rfrak_chi2[0], Vars.rfrak_chi2[1],(mul(mul(inverse(Cons.S_chi2),OmegaVec_chiVec_2_5p5(Cons,Vars)),Cons.S_chi2))[1:])
    else:
        dydt[3], dydt[4] = 0.0, 0.0
    return dydt      

@njit(cache=True)
def Recalculate_6p0(Cons,Vars,y):
    Vars.v = np.array([y[0]])
    Vars.rfrak_chi1 = np.array([y[1],y[2]])
    Vars.rfrak_chi2 = np.array([y[3],y[4]])
    Vars.R = y[5]*Cons.wHat + y[6]*Cons.xHat + y[7]*Cons.yHat + y[8]*Cons.zHat
    Vars.nHat = mul(mul(Vars.R,Cons.xHat),inverse(Vars.R))
    Vars.lambdaHat = mul(mul(Vars.R,Cons.yHat),inverse(Vars.R))
    Vars.ellHat = mul(mul(Vars.R,Cons.zHat),inverse(Vars.R))
    Vars.R_S1 = exp(Vars.rfrak_chi1[0]*Cons.xHat + Vars.rfrak_chi1[1]*Cons.yHat)
    Vars.R_S2 = exp(Vars.rfrak_chi2[0]*Cons.xHat + Vars.rfrak_chi2[1]*Cons.yHat)
    Vars.chiVec1 = mul(mul(mul(Cons.S_chi1,Vars.R_S1),Cons.zHat),mul(conjugate(Vars.R_S1),conjugate(Cons.S_chi1)))
    Vars.chiVec2 = mul(mul(mul(Cons.S_chi2,Vars.R_S2),Cons.zHat),mul(conjugate(Vars.R_S2),conjugate(Cons.S_chi2)))
    Vars.chi1_n = np.array([dot(Vars.chiVec1[1:],Vars.nHat[1:])])
    Vars.chi1_lambda = np.array([dot(Vars.chiVec1[1:],Vars.lambdaHat[1:])])
    Vars.chi1_ell = np.array([dot(Vars.chiVec1[1:],Vars.ellHat[1:])])
    Vars.chi2_n = np.array([dot(Vars.chiVec2[1:],Vars.nHat[1:])])
    Vars.chi2_lambda = np.array([dot(Vars.chiVec2[1:],Vars.lambdaHat[1:])])
    Vars.chi2_ell = np.array([dot(Vars.chiVec2[1:],Vars.ellHat[1:])])
    Vars.S = Cons.M1**2*Vars.chiVec1 + Cons.M2**2*Vars.chiVec2
    Vars.S_ell = Cons.M1**2*Vars.chi1_ell + Cons.M2**2*Vars.chi2_ell
    Vars.S_n = Cons.M1**2*Vars.chi1_n + Cons.M2**2*Vars.chi2_n
    Vars.S_lambda = Cons.M1**2*Vars.chi1_lambda + Cons.M2**2*Vars.chi2_lambda
    Vars.Sigma = Cons.M*(-Cons.M1*Vars.chiVec1 + Cons.M2*Vars.chiVec2)
    Vars.Sigma_ell = Cons.M*(-Cons.M1*Vars.chi1_ell + Cons.M2*Vars.chi2_ell)
    Vars.Sigma_n = Cons.M*(-Cons.M1*Vars.chi1_n + Cons.M2*Vars.chi2_n)
    Vars.Sigma_lambda = Cons.M*(-Cons.M1*Vars.chi1_lambda + Cons.M2*Vars.chi2_lambda)
    Vars.SS = np.array([dot(Vars.S[1:],Vars.S[1:])])
    Vars.SigmaSigma = np.array([dot(Vars.Sigma[1:],Vars.Sigma[1:])])
    Vars.SSigma = np.array([dot(Vars.S[1:],Vars.Sigma[1:])])
    Vars.S1 = Cons.M1**2*Vars.chiVec1
    Vars.S1_n = Cons.M1**2*Vars.chi1_n
    Vars.S1_lambda = Cons.M1**2*Vars.chi1_lambda
    Vars.S2 = Cons.M2**2*Vars.chiVec2
    Vars.S2_n = Cons.M2**2*Vars.chi2_n
    Vars.S2_lambda = Cons.M2**2*Vars.chi2_lambda
    Vars.chi_s_ell = Vars.chi1_ell/2 + Vars.chi2_ell/2
    Vars.chi_a_ell = Vars.chi1_ell/2 - Vars.chi2_ell/2
    Vars.logv = log(Vars.v)
    Vars.Fcal_coeff = 32*Cons.nu**2*Vars.v**10/5
    Vars.Fcal_SQ_4 = Cons.chi1chi1*(-89*Cons.delta/192 + 89*Cons.nu/96 - 89/192) - 103*Cons.chi1chi2*Cons.nu/48 + Cons.chi2chi2*(89*Cons.delta/192 + 89*Cons.nu/96 - 89/192) + Vars.chi_a_ell*(Vars.chi_a_ell*(287/96 - 12*Cons.nu) + 287*Vars.chi_s_ell*Cons.delta/48) + Vars.chi_s_ell**2*(Cons.nu/24 + 287/96)
    Vars.Fcal_SQ_6 = Cons.nu*(-0.0714285714285714*Cons.delta*(353*Vars.SSigma + 86*Vars.S_lambda*Vars.Sigma_lambda - 1145*Vars.S_n*Vars.Sigma_n)/Cons.M**4 - 0.00297619047619048*(4236*Vars.SS*Cons.kappa_p + 8472*Vars.SS + 4236*Vars.SSigma*Cons.delta*Cons.kappa_p + 1672*Vars.SSigma*Cons.kappa_m + 1032*Vars.S_lambda**2*(Cons.kappa_p + 2) + 8*Vars.S_lambda*Vars.Sigma_lambda*(129*Cons.delta*Cons.kappa_p + 620*Cons.kappa_m) - 13740*Vars.S_n**2*(Cons.kappa_p + 2) - 4*Vars.S_n*Vars.Sigma_n*(3435*Cons.delta*Cons.kappa_p + 2494*Cons.kappa_m) - 641*Vars.SigmaSigma*Cons.delta*Cons.kappa_m - 2713*Vars.SigmaSigma*Cons.kappa_p - 4146*Vars.SigmaSigma + 982*Vars.Sigma_lambda**2*Cons.delta*Cons.kappa_m - 6944*Vars.Sigma_lambda**2*Cons.kappa_p - 3736*Vars.Sigma_lambda**2 + 941*Vars.Sigma_n**2*Cons.delta*Cons.kappa_m + 15083*Vars.Sigma_n**2*Cons.kappa_p + 44524*Vars.Sigma_n**2)/Cons.M**4) - 0.00595238095238095*Cons.delta*(911*Vars.SSigma + 3643*Vars.S_lambda*Vars.Sigma_lambda - 22399*Vars.S_n*Vars.Sigma_n)/Cons.M**4 + 0.0357142857142857*Cons.nu**2*(Cons.kappa_p + 2)*(353*Vars.SigmaSigma + 86*Vars.Sigma_lambda**2 - 1145*Vars.Sigma_n**2)/Cons.M**4 + 0.00297619047619048*(1477*Vars.SS*Cons.delta*Cons.kappa_m - 1877*Vars.SS*Cons.kappa_p + 502*Vars.SS - 3354*Vars.SSigma*Cons.delta*Cons.kappa_p + 3354*Vars.SSigma*Cons.kappa_m + Vars.S_lambda**2*(1498*Cons.delta*Cons.kappa_m - 4464*Cons.kappa_p - 10811) - 5962*Vars.S_lambda*Vars.Sigma_lambda*(Cons.delta*Cons.kappa_p - Cons.kappa_m) + Vars.S_n**2*(-5929*Cons.delta*Cons.kappa_m + 10095*Cons.kappa_p + 45299) + 16024*Vars.S_n*Vars.Sigma_n*(Cons.delta*Cons.kappa_p - Cons.kappa_m) + 1677*Vars.SigmaSigma*Cons.delta*Cons.kappa_m - 1677*Vars.SigmaSigma*Cons.kappa_p - 168*Vars.SigmaSigma + 2981*Vars.Sigma_lambda**2*Cons.delta*Cons.kappa_m - 2981*Vars.Sigma_lambda**2*Cons.kappa_p + 1313*Vars.Sigma_lambda**2 - 8012*Vars.Sigma_n**2*Cons.delta*Cons.kappa_m + 8012*Vars.Sigma_n**2*Cons.kappa_p + 6387*Vars.Sigma_n**2)/Cons.M**4
    Vars.Fcal_SQ_7 = Vars.S_ell**3*(-16*Cons.kappa_p/3 - 4*Cons.lambda_p + 40/3) + Vars.S_ell**2*Vars.Sigma_ell*(-35*Cons.delta*Cons.kappa_p/6 - 6*Cons.delta*Cons.lambda_p + 73*Cons.delta/3 - 3*Cons.kappa_m/4 + 6*Cons.lambda_m) + Vars.S_ell*Vars.Sigma_ell**2*(-35*Cons.delta*Cons.kappa_m/12 + 6*Cons.delta*Cons.lambda_m + 35*Cons.kappa_p/12 - 6*Cons.lambda_p + Cons.nu*(22*Cons.kappa_p/3 + 12*Cons.lambda_p - 172/3) + 32/3) + Vars.Sigma_ell**3*(67*Cons.delta*Cons.kappa_p/24 - 2*Cons.delta*Cons.lambda_p - Cons.delta/8 - 67*Cons.kappa_m/24 + 2*Cons.lambda_m + Cons.nu*(Cons.delta*Cons.kappa_p/2 + 2*Cons.delta*Cons.lambda_p - 11*Cons.delta + 61*Cons.kappa_m/12 - 6*Cons.lambda_m))
    Vars.Fcal_SO_3 = (-4*Vars.S_ell - 5*Vars.Sigma_ell*Cons.delta/4)/Cons.M**2
    Vars.Fcal_SO_5 = (Vars.S_ell*(272*Cons.nu/9 - 9/2) + Vars.Sigma_ell*Cons.delta*(43*Cons.nu/4 - 13/16))/Cons.M**2
    Vars.Fcal_SO_6 = (-16*Vars.S_ell*pi - 31*Vars.Sigma_ell*Cons.delta*pi/6)/Cons.M**2
    Vars.Fcal_SO_7 = (Vars.S_ell*(-2810*Cons.nu**2/27 + 6172*Cons.nu/189 + 476645/6804) + Vars.Sigma_ell*Cons.delta*(-1501*Cons.nu**2/36 + 1849*Cons.nu/126 + 9535/336))/Cons.M**2
    Vars.Fcal_SO_8 = (Vars.S_ell*pi*(13879*Cons.nu/72 - 3485/96) + Vars.Sigma_ell*Cons.delta*pi*(130583*Cons.nu/2016 - 7163/672))/Cons.M**2
    Vars.E_SQ_4 = -3*Vars.chi_a_ell**2/2 - 3*Vars.chi_s_ell**2/2 - Cons.delta*(Cons.chi2chi2/2 + 3*Vars.chi_a_ell*Vars.chi_s_ell) + Cons.nu*(Cons.chi1chi2 + 6*Vars.chi_a_ell**2) + (Cons.chi1chi1 + Cons.chi2chi2)*(Cons.delta - 2*Cons.nu + 1)/4
    Vars.E_SQ_6 = Cons.nu*(-3*Cons.delta*(Vars.SSigma - 3*Vars.S_n*Vars.Sigma_n)/Cons.M**4 + 0.25*(-6*Vars.SS*Cons.kappa_p - 12*Vars.SS - 6*Vars.SSigma*Cons.delta*Cons.kappa_p - 26*Vars.SSigma*Cons.kappa_m + 24*Vars.S_lambda*Vars.Sigma_lambda*Cons.kappa_m + 18*Vars.S_n**2*(Cons.kappa_p + 2) + 18*Vars.S_n*Vars.Sigma_n*(Cons.delta*Cons.kappa_p + 3*Cons.kappa_m) + 19*Vars.SigmaSigma*Cons.kappa_p - 80*Vars.SigmaSigma - 18*Vars.Sigma_lambda**2*Cons.kappa_p + 76*Vars.Sigma_lambda**2 - 39*Vars.Sigma_n**2*Cons.kappa_p + 168*Vars.Sigma_n**2 + Cons.delta*Cons.kappa_m*(-5*Vars.SigmaSigma + 6*Vars.Sigma_lambda**2 + 9*Vars.Sigma_n**2))/Cons.M**4) + Cons.delta*(15*Vars.SSigma - 13*Vars.S_lambda*Vars.Sigma_lambda - 29*Vars.S_n*Vars.Sigma_n)/Cons.M**4 + 1.5*Cons.nu**2*(Vars.SigmaSigma - 3*Vars.Sigma_n**2)*(Cons.kappa_p + 2)/Cons.M**4 + 0.25*(8*Vars.SS*Cons.delta*Cons.kappa_m - 6*Vars.SS*Cons.kappa_p + 36*Vars.SS - 14*Vars.SSigma*Cons.delta*Cons.kappa_p + 14*Vars.SSigma*Cons.kappa_m + Vars.S_lambda**2*(-6*Cons.delta*Cons.kappa_m + 6*Cons.kappa_p - 28) + 12*Vars.S_lambda*Vars.Sigma_lambda*(Cons.delta*Cons.kappa_p - Cons.kappa_m) - 6*Vars.S_n**2*(3*Cons.delta*Cons.kappa_m - 2*Cons.kappa_p + 10) + 30*Vars.S_n*Vars.Sigma_n*(Cons.delta*Cons.kappa_p - Cons.kappa_m) + 7*Vars.SigmaSigma*Cons.delta*Cons.kappa_m - 7*Vars.SigmaSigma*Cons.kappa_p + 24*Vars.SigmaSigma - 6*Vars.Sigma_lambda**2*Cons.delta*Cons.kappa_m + 6*Vars.Sigma_lambda**2*Cons.kappa_p - 24*Vars.Sigma_lambda**2 - 15*Vars.Sigma_n**2*Cons.delta*Cons.kappa_m + 15*Vars.Sigma_n**2*Cons.kappa_p - 48*Vars.Sigma_n**2)/Cons.M**4
    Vars.E_SQ_7 = Vars.S_ell**3*(2*Cons.kappa_p + 4*Cons.lambda_p - 20) + Vars.S_ell**2*Vars.Sigma_ell*(2*Cons.delta*Cons.kappa_p + 6*Cons.delta*Cons.lambda_p - 32*Cons.delta + 4*Cons.kappa_m - 6*Cons.lambda_m) + Vars.S_ell*Vars.Sigma_ell**2*(5*Cons.delta*Cons.kappa_m - 6*Cons.delta*Cons.lambda_m - 5*Cons.kappa_p + 6*Cons.lambda_p + Cons.nu*(-2*Cons.kappa_p - 12*Cons.lambda_p + 68) - 12) + Vars.Sigma_ell**3*(-3*Cons.delta*Cons.kappa_p + 2*Cons.delta*Cons.lambda_p + 3*Cons.kappa_m - 2*Cons.lambda_m + Cons.nu*(-2*Cons.delta*Cons.lambda_p + 12*Cons.delta - 6*Cons.kappa_m + 6*Cons.lambda_m))
    Vars.E_SO_3 = (14*Vars.S_ell/3 + 2*Vars.Sigma_ell*Cons.delta)/Cons.M**2
    Vars.E_SO_5 = (Vars.S_ell*(11 - 61*Cons.nu/9) + Vars.Sigma_ell*Cons.delta*(3 - 10*Cons.nu/3))/Cons.M**2
    Vars.E_SO_7 = (Vars.S_ell*(29*Cons.nu**2/12 - 367*Cons.nu/4 + 135/4) + Vars.Sigma_ell*Cons.delta*(5*Cons.nu**2/4 - 39*Cons.nu + 27/4))/Cons.M**2
    Vars.MDot_Alvi_5 = -Cons.M1**3*Vars.chi1_ell*(3*Cons.chi1chi1 + 1)/(4*Cons.M**3) - Cons.M2**3*Vars.chi2_ell*(3*Cons.chi2chi2 + 1)/(4*Cons.M**3)

@njit
def OmegaVec_chiVec_1_6p0(Cons,Vars):
    Omega1_coeff = Vars.v**5/Cons.M
    return Omega1_coeff*(Vars.S2*Vars.v**3*((2.0 - 0.5*Cons.delta)/Cons.M**2 + Vars.v*(-5.0*Vars.S_ell - 3.0*Vars.Sigma_ell*Cons.delta/Cons.M)/Cons.M**4) + Vars.ellHat*(-0.75*Cons.delta + 0.5*Cons.nu + Vars.v**2*(Cons.delta*(0.625*Cons.nu - 0.5625) + Cons.nu*(1.25 - 0.0416666666666667*Cons.nu) + Vars.v*(Vars.v*(Cons.nu*(Cons.nu*(0.229166666666667*Cons.nu - 4.94791666666667) + 2.89583333333333) + 0.96875 + (Cons.delta*(Cons.nu*(5.91666666666667 - 0.489583333333333*Cons.nu) - 0.96875) + (Vars.S_ell*(Vars.S_ell*(Cons.kappa_p*(-0.25*Cons.nu - 0.375) - 0.5*Cons.nu - 0.75) + Vars.Sigma_ell*(Cons.delta*(Cons.kappa_p*(-0.25*Cons.nu - 0.375) - 0.5*Cons.nu - 0.75) + Cons.kappa_m*(0.25*Cons.nu + 0.375))) + Vars.Sigma_ell**2*(Cons.delta*Cons.kappa_m*(0.125*Cons.nu + 0.1875) + Cons.kappa_p*(Cons.nu*(0.25*Cons.nu + 0.25) - 0.1875) + Cons.nu*(0.5*Cons.nu + 0.75)) + (Vars.S_ell*(Vars.S_ell*Cons.delta*(0.375*Cons.kappa_p + 0.75) + Vars.Sigma_ell*Cons.delta*(Cons.delta*(0.375*Cons.kappa_p + 0.75) - 0.375*Cons.kappa_m)) + Vars.Sigma_ell**2*Cons.delta*(-0.1875*Cons.delta*Cons.kappa_m + Cons.kappa_p*(0.1875 - 0.375*Cons.nu) - 0.75*Cons.nu))/Cons.M)/Cons.M**7)/Cons.M) + (Vars.S_ell*(0.833333333333333*Cons.nu + 1.25) + (-1.25*Vars.S_ell*Cons.delta + Vars.Sigma_ell*Cons.delta*(0.5*Cons.nu + 0.75) - 0.75*Vars.Sigma_ell*Cons.delta**2/Cons.M)/Cons.M)/Cons.M**2) + 0.5625) + 0.75) + Vars.lambdaHat*Vars.v**3*(Vars.S2_lambda*(0.5*Cons.delta + 2.0)/Cons.M**2 + Vars.v*(Vars.S1_lambda*(-0.75*Cons.delta - 3.0*Cons.kappa_1 + 2.25) + Vars.S1_lambda*(Cons.delta*(0.75 - 1.5*Cons.kappa_1) + 1.5*Cons.kappa_1 - 0.75)/Cons.nu)/Cons.M**2) + Vars.nHat*Vars.v*(3.0*Vars.chi1_n*Cons.nu + Vars.v**2*((Vars.S1_n*(Cons.delta*(-2.25*Cons.kappa_1 - 0.25) + Cons.kappa_1*(1.5*Cons.nu + 0.75) + 0.5*Cons.nu + 3.25) + Vars.S1_n*(Cons.delta*(1.5 - 0.75*Cons.kappa_1) + 0.75*Cons.kappa_1 - 1.5)/Cons.nu + Vars.S2_n*(1.5*Cons.delta - 1.0*Cons.nu + 3.0))/Cons.M**2 + Vars.v*(-15.0*Vars.S1_n*Vars.S_ell*Cons.kappa_1 + Vars.S1_n*Vars.S_ell*(-7.5*Cons.delta*Cons.kappa_1 + 7.5*Cons.kappa_1)/Cons.nu + 15.0*Vars.S2_n*Vars.S_ell + (-9.0*Vars.S1_n*Vars.Sigma_ell*Cons.delta*Cons.kappa_1 + Vars.S1_n*Vars.Sigma_ell*Cons.delta*(-4.5*Cons.delta*Cons.kappa_1 + 4.5*Cons.kappa_1)/Cons.nu + 9.0*Vars.S2_n*Vars.Sigma_ell*Cons.delta)/Cons.M)/Cons.M**4) + 3.0*Cons.M2**2*Vars.chi2_n/Cons.M**2) - Cons.M2**2*Vars.chiVec2*Vars.v/Cons.M**2)

@njit
def OmegaVec_chiVec_2_6p0(Cons,Vars):
    Omega2_coeff = Vars.v**5/Cons.M
    return Omega2_coeff*(Vars.S1*Vars.v**3*((0.5*Cons.delta + 2.0)/Cons.M**2 + Vars.v*(-5.0*Vars.S_ell + 3.0*Vars.Sigma_ell*Cons.delta/Cons.M)/Cons.M**4) + Vars.ellHat*(0.75*Cons.delta + 0.5*Cons.nu + Vars.v**2*(Cons.delta*(0.5625 - 0.625*Cons.nu) + Cons.nu*(1.25 - 0.0416666666666667*Cons.nu) + Vars.v*(Vars.v*(Cons.nu*(Cons.nu*(0.229166666666667*Cons.nu - 4.94791666666667) + 2.89583333333333) + 0.96875 + (Cons.delta*(Cons.nu*(0.489583333333333*Cons.nu - 5.91666666666667) + 0.96875) + (Vars.S_ell*(Vars.S_ell*(Cons.kappa_p*(-0.25*Cons.nu - 0.375) - 0.5*Cons.nu - 0.75) + Vars.Sigma_ell*(Cons.delta*(Cons.kappa_p*(0.25*Cons.nu + 0.375) + 0.5*Cons.nu + 0.75) + Cons.kappa_m*(-0.25*Cons.nu - 0.375))) + Vars.Sigma_ell**2*(Cons.delta*Cons.kappa_m*(0.125*Cons.nu + 0.1875) + Cons.kappa_p*(Cons.nu*(0.25*Cons.nu + 0.25) - 0.1875) + Cons.nu*(0.5*Cons.nu + 0.75)) + (Vars.S_ell*(Vars.S_ell*Cons.delta*(-0.375*Cons.kappa_p - 0.75) + Vars.Sigma_ell*Cons.delta*(Cons.delta*(0.375*Cons.kappa_p + 0.75) - 0.375*Cons.kappa_m)) + Vars.Sigma_ell**2*Cons.delta*(0.1875*Cons.delta*Cons.kappa_m + Cons.kappa_p*(0.375*Cons.nu - 0.1875) + 0.75*Cons.nu))/Cons.M)/Cons.M**7)/Cons.M) + (Vars.S_ell*(0.833333333333333*Cons.nu + 1.25) + (1.25*Vars.S_ell*Cons.delta + Vars.Sigma_ell*Cons.delta*(-0.5*Cons.nu - 0.75) - 0.75*Vars.Sigma_ell*Cons.delta**2/Cons.M)/Cons.M)/Cons.M**2) + 0.5625) + 0.75) + Vars.lambdaHat*Vars.v**3*(Vars.S1_lambda*(2.0 - 0.5*Cons.delta)/Cons.M**2 + Vars.v*(Vars.S2_lambda*(0.75*Cons.delta - 3.0*Cons.kappa_2 + 2.25) + Vars.S2_lambda*(Cons.delta*(1.5*Cons.kappa_2 - 0.75) + 1.5*Cons.kappa_2 - 0.75)/Cons.nu)/Cons.M**2) + Vars.nHat*Vars.v*(3.0*Vars.chi2_n*Cons.nu + Vars.v**2*((Vars.S1_n*(-1.5*Cons.delta - 1.0*Cons.nu + 3.0) + Vars.S2_n*(Cons.delta*(2.25*Cons.kappa_2 + 0.25) + Cons.kappa_2*(1.5*Cons.nu + 0.75) + 0.5*Cons.nu + 3.25) + Vars.S2_n*(Cons.delta*(0.75*Cons.kappa_2 - 1.5) + 0.75*Cons.kappa_2 - 1.5)/Cons.nu)/Cons.M**2 + Vars.v*(15.0*Vars.S1_n*Vars.S_ell - 15.0*Vars.S2_n*Vars.S_ell*Cons.kappa_2 + Vars.S2_n*Vars.S_ell*(7.5*Cons.delta*Cons.kappa_2 + 7.5*Cons.kappa_2)/Cons.nu + (-9.0*Vars.S1_n*Vars.Sigma_ell*Cons.delta + 9.0*Vars.S2_n*Vars.Sigma_ell*Cons.delta*Cons.kappa_2 + Vars.S2_n*Vars.Sigma_ell*Cons.delta*(-4.5*Cons.delta*Cons.kappa_2 - 4.5*Cons.kappa_2)/Cons.nu)/Cons.M)/Cons.M**4) + 3.0*Cons.M1**2*Vars.chi1_n/Cons.M**2) - Cons.M1**2*Vars.chiVec1*Vars.v/Cons.M**2)

@njit
def OmegaVec_6p0(Cons,Vars):
    gamma_PN_5 = (0.888888888888889*Vars.S_ell*Cons.nu + 3.33333333333333*Vars.S_ell)/Cons.M**2 + 2.0*Vars.Sigma_ell*Cons.delta/Cons.M**3
    gamma_PN_6 = Vars.S_ell**2*(-0.916666666666667*Cons.delta*Cons.kappa_m - 0.916666666666667*Cons.kappa_p + Cons.nu*(-0.166666666666667*Cons.kappa_p - 0.333333333333333) + 1.55555555555556) + Vars.S_ell*Vars.Sigma_ell*(1.66666666666667*Cons.delta + Cons.nu*(-0.166666666666667*Cons.delta*Cons.kappa_p - 0.333333333333333*Cons.delta + 3.83333333333333*Cons.kappa_m)) + Vars.Sigma_ell**2*(Cons.nu**2*(0.166666666666667*Cons.kappa_p + 0.333333333333333) + Cons.nu*(Cons.delta*Cons.kappa_m - Cons.kappa_p - 2.0) + 1.0) + 0.0123456790123457*Cons.nu**3 + 6.36111111111111*Cons.nu**2 - 2.98177812235564*Cons.nu + 1.0
    gamma_PN_0 = 1.00000000000000
    a_ell_4 = Vars.S_n*(5.77777777777778*Cons.nu**2 + 14.75*Cons.nu + 1.5) + Vars.Sigma_n*Cons.delta*(2.83333333333333*Cons.nu**2 + 9.125*Cons.nu + 1.5)
    a_ell_1 = 3.0*Vars.Sigma_ell*Vars.Sigma_n*Cons.nu*(Cons.kappa_p + 2.0)/Cons.M**2 - 3.0*Cons.delta*(Vars.S_ell*Vars.Sigma_n + Vars.S_n*Vars.Sigma_ell)/Cons.M**2 - 0.5*(3.0*Vars.S_n*(2.0*Vars.S_ell*(Cons.kappa_p + 2.0) + Vars.Sigma_ell*(Cons.delta*Cons.kappa_p - Cons.kappa_m)) + 3.0*Vars.Sigma_n*(Vars.S_ell*(Cons.delta*Cons.kappa_p - Cons.kappa_m) + Vars.Sigma_ell*(-Cons.delta*Cons.kappa_m + Cons.kappa_p)))/Cons.M**2
    a_ell_2 = Vars.S_n*(-9.66666666666667*Cons.nu - 10.0) + Vars.Sigma_n*Cons.delta*(-4.5*Cons.nu - 6.0)
    gamma_PN_7 = (Vars.S_ell*(-6.0*Cons.nu**2 - 10.5833333333333*Cons.nu + 5.0) - 2.66666666666667*Vars.Sigma_ell*Cons.delta*Cons.nu**2 + Vars.Sigma_ell*Cons.delta*(3.0 - 10.1666666666667*Cons.nu))/Cons.M**2
    gamma_PN_2 = 1.0 - 0.333333333333333*Cons.nu
    a_ell_3 = Cons.nu*(4.0*Cons.delta*(Vars.S_ell*Vars.Sigma_n + Vars.S_n*Vars.Sigma_ell)/Cons.M**2 + 0.5*(2.0*Vars.S_n*(4.0*Vars.S_ell*(Cons.kappa_p + 2.0) + Vars.Sigma_ell*(2.0*Cons.delta*Cons.kappa_p + 7.0*Cons.kappa_m)) + Vars.Sigma_n*(2.0*Vars.S_ell*(2.0*Cons.delta*Cons.kappa_p + 7.0*Cons.kappa_m) + Vars.Sigma_ell*(5.0*Cons.delta*Cons.kappa_m - 17.0*Cons.kappa_p - 26.0)))/Cons.M**2) - 4.0*Vars.Sigma_ell*Vars.Sigma_n*Cons.nu**2*(Cons.kappa_p + 2.0)/Cons.M**2 + 6.0*Cons.delta*(Vars.S_ell*Vars.Sigma_n + Vars.S_n*Vars.Sigma_ell)/Cons.M**2 + 0.5*(Vars.S_n*(Vars.S_ell*(-9.0*Cons.delta*Cons.kappa_m + 3.0*Cons.kappa_p + 22.0) + Vars.Sigma_ell*(6.0*Cons.delta*Cons.kappa_p - 6.0*Cons.kappa_m)) + 6.0*Vars.Sigma_n*(Vars.S_ell*(Cons.delta*Cons.kappa_p - Cons.kappa_m) + Vars.Sigma_ell*(-Cons.delta*Cons.kappa_m + Cons.kappa_p + 1.0)))/Cons.M**2
    gamma_PN_3 = (1.66666666666667*Vars.S_ell + Vars.Sigma_ell*Cons.delta)/Cons.M**2
    a_ell_0 = 7.0*Vars.S_n + 3.0*Vars.Sigma_n*Cons.delta
    gamma_PN_4 = -5.41666666666667*Cons.nu + 1.0 + (Vars.S_ell**2*(-0.5*Cons.kappa_p - 1.0) + Vars.S_ell*Vars.Sigma_ell*(-0.5*Cons.delta*Cons.kappa_p - Cons.delta + 0.5*Cons.kappa_m) + Vars.Sigma_ell**2*(0.25*Cons.delta*Cons.kappa_m - 0.25*Cons.kappa_p + Cons.nu*(0.5*Cons.kappa_p + 1.0)))/Cons.M**8
    return Vars.ellHat*Vars.v**3/Cons.M + Vars.nHat*Vars.v**6*(a_ell_0 + Vars.v*(a_ell_1 + Vars.v*(a_ell_2 + Vars.v*(a_ell_3 + a_ell_4*Vars.v))))*(gamma_PN_0 + Vars.v**2*(gamma_PN_2 + Vars.v*(gamma_PN_3 + Vars.v*(gamma_PN_4 + Vars.v*(gamma_PN_5 + Vars.v*(gamma_PN_6 + gamma_PN_7*Vars.v))))))/Cons.M**3


@njit(cache=True)
def TaylorT1_6p0(Cons,Vars):
    Flux = Vars.Fcal_coeff*(Cons.Fcal_0 + Vars.v**2*(Cons.Fcal_2 + Vars.v*(Cons.Fcal_3 + Vars.Fcal_SO_3 + Vars.v*(Cons.Fcal_4 + Vars.Fcal_SQ_4 + Vars.v*(Cons.Fcal_5 + Vars.Fcal_SO_5 + Vars.v*(Cons.Fcal_6 + Vars.Fcal_SO_6 + Vars.Fcal_SQ_6 + Cons.Fcal_lnv_6*Vars.logv + Vars.v*(Cons.Fcal_7 + Vars.Fcal_SO_7 + Vars.Fcal_SQ_7 + Vars.v*(Cons.Fcal_8 + Vars.Fcal_SO_8 + Cons.Fcal_lnv_8*Vars.logv + Vars.v*(Cons.Fcal_9 + Cons.Fcal_lnv_9*Vars.logv + Vars.v*(Cons.Fcal_10 + Cons.Fcal_lnv_10*Vars.logv + Vars.v*(Cons.Fcal_11 + Cons.Fcal_lnv_11*Vars.logv + Vars.v*(Cons.Fcal_12 + Cons.Fcal_lnv2_12*Vars.logv**2 + Cons.Fcal_lnv_12*Vars.logv))))))))))))
    dEdV = -0.5*Cons.M*Cons.nu*Vars.v*(2.0*Cons.E_0 + Vars.v**2*(4.0*Cons.E_2 + Vars.v*(5.0*Vars.E_SO_3 + Vars.v*(6.0*Cons.E_4 + 6.0*Vars.E_SQ_4 + Vars.v*(7.0*Vars.E_SO_5 + Vars.v*(8.0*Cons.E_6 + 8.0*Vars.E_SQ_6 + Vars.v*(9.0*Vars.E_SO_7 + 9.0*Vars.E_SQ_7 + Vars.v*(10.0*Cons.E_8 + Cons.E_lnv_8*(10.0*Vars.logv + 1.0) + Vars.v**2*(12.0*Cons.E_10 + Cons.E_lnv_10*(12.0*Vars.logv + 1.0) + Vars.v*(13.0*Cons.E_11 + Vars.v*(14.0*Cons.E_12 + Cons.E_lnv_12*(14.0*Vars.logv + 1.0))))))))))))
    Absorption = Vars.Fcal_coeff*Vars.v**5*(Vars.MDot_Alvi_5 + Cons.MDot_Alvi_8*Vars.v**3)
    dvdt_T1 = (-Absorption - Flux)/dEdV
    dydt=np.zeros(9)
    [dydt[5],dydt[6],dydt[7],dydt[8]] = FrameFromAngularVelocityIntegrand(Vars.R, OmegaVec_6p0(Cons,Vars)[1:])
    dydt[0] = dvdt_T1[0] 
    if(Cons.EvolveSpin1):
        dydt[1], dydt[2]=FrameFromAngularVelocity_2D_Integrand(Vars.rfrak_chi1[0], Vars.rfrak_chi1[1],(mul(mul(inverse(Cons.S_chi1),OmegaVec_chiVec_1_6p0(Cons,Vars)),Cons.S_chi1))[1:])
    else:
        dydt[1], dydt[2] = 0.0, 0.0
    if(Cons.EvolveSpin2):
        dydt[3], dydt[4]=FrameFromAngularVelocity_2D_Integrand(Vars.rfrak_chi2[0], Vars.rfrak_chi2[1],(mul(mul(inverse(Cons.S_chi2),OmegaVec_chiVec_2_6p0(Cons,Vars)),Cons.S_chi2))[1:])
    else:
        dydt[3], dydt[4] = 0.0, 0.0
    return dydt

@njit(cache=True)
def TaylorT4_6p0(Cons,Vars):
    Flux = Vars.Fcal_coeff*(Cons.Fcal_0 + Vars.v**2*(Cons.Fcal_2 + Vars.v*(Cons.Fcal_3 + Vars.Fcal_SO_3 + Vars.v*(Cons.Fcal_4 + Vars.Fcal_SQ_4 + Vars.v*(Cons.Fcal_5 + Vars.Fcal_SO_5 + Vars.v*(Cons.Fcal_6 + Vars.Fcal_SO_6 + Vars.Fcal_SQ_6 + Cons.Fcal_lnv_6*Vars.logv + Vars.v*(Cons.Fcal_7 + Vars.Fcal_SO_7 + Vars.Fcal_SQ_7 + Vars.v*(Cons.Fcal_8 + Vars.Fcal_SO_8 + Cons.Fcal_lnv_8*Vars.logv + Vars.v*(Cons.Fcal_9 + Cons.Fcal_lnv_9*Vars.logv + Vars.v*(Cons.Fcal_10 + Cons.Fcal_lnv_10*Vars.logv + Vars.v*(Cons.Fcal_11 + Cons.Fcal_lnv_11*Vars.logv + Vars.v*(Cons.Fcal_12 + Cons.Fcal_lnv2_12*Vars.logv**2 + Cons.Fcal_lnv_12*Vars.logv))))))))))))
    dEdV = -0.5*Cons.M*Cons.nu*Vars.v*(2.0*Cons.E_0 + Vars.v**2*(4.0*Cons.E_2 + Vars.v*(5.0*Vars.E_SO_3 + Vars.v*(6.0*Cons.E_4 + 6.0*Vars.E_SQ_4 + Vars.v*(7.0*Vars.E_SO_5 + Vars.v*(8.0*Cons.E_6 + 8.0*Vars.E_SQ_6 + Vars.v*(9.0*Vars.E_SO_7 + 9.0*Vars.E_SQ_7 + Vars.v*(10.0*Cons.E_8 + Cons.E_lnv_8*(10.0*Vars.logv + 1.0) + Vars.v**2*(12.0*Cons.E_10 + Cons.E_lnv_10*(12.0*Vars.logv + 1.0) + Vars.v*(13.0*Cons.E_11 + Vars.v*(14.0*Cons.E_12 + Cons.E_lnv_12*(14.0*Vars.logv + 1.0))))))))))))
    Absorption = Vars.Fcal_coeff*Vars.v**5*(Vars.MDot_Alvi_5 + Cons.MDot_Alvi_8*Vars.v**3)
    dvdt_T1 = (-Absorption - Flux)/dEdV
    dydt=np.zeros(9)
    [dydt[5],dydt[6],dydt[7],dydt[8]] = FrameFromAngularVelocityIntegrand(Vars.R, OmegaVec_6p0(Cons,Vars)[1:])
    dydt[0] = dvdt_T4[0]
    if(Cons.EvolveSpin1):
        dydt[1], dydt[2]=FrameFromAngularVelocity_2D_Integrand(Vars.rfrak_chi1[0], Vars.rfrak_chi1[1],(mul(mul(inverse(Cons.S_chi1),OmegaVec_chiVec_1_6p0(Cons,Vars)),Cons.S_chi1))[1:])
    else:
        dydt[1], dydt[2] = 0.0, 0.0
    if(Cons.EvolveSpin2):
        dydt[3], dydt[4]=FrameFromAngularVelocity_2D_Integrand(Vars.rfrak_chi2[0], Vars.rfrak_chi2[1],(mul(mul(inverse(Cons.S_chi2),OmegaVec_chiVec_2_6p0(Cons,Vars)),Cons.S_chi2))[1:])
    else:
        dydt[3], dydt[4] = 0.0, 0.0
    return dydt      

@njit(cache=True)
def TaylorT5_6p0(Cons,Vars):
    Flux = Vars.Fcal_coeff*(Cons.Fcal_0 + Vars.v**2*(Cons.Fcal_2 + Vars.v*(Cons.Fcal_3 + Vars.Fcal_SO_3 + Vars.v*(Cons.Fcal_4 + Vars.Fcal_SQ_4 + Vars.v*(Cons.Fcal_5 + Vars.Fcal_SO_5 + Vars.v*(Cons.Fcal_6 + Vars.Fcal_SO_6 + Vars.Fcal_SQ_6 + Cons.Fcal_lnv_6*Vars.logv + Vars.v*(Cons.Fcal_7 + Vars.Fcal_SO_7 + Vars.Fcal_SQ_7 + Vars.v*(Cons.Fcal_8 + Vars.Fcal_SO_8 + Cons.Fcal_lnv_8*Vars.logv + Vars.v*(Cons.Fcal_9 + Cons.Fcal_lnv_9*Vars.logv + Vars.v*(Cons.Fcal_10 + Cons.Fcal_lnv_10*Vars.logv + Vars.v*(Cons.Fcal_11 + Cons.Fcal_lnv_11*Vars.logv + Vars.v*(Cons.Fcal_12 + Cons.Fcal_lnv2_12*Vars.logv**2 + Cons.Fcal_lnv_12*Vars.logv))))))))))))
    dEdV = -0.5*Cons.M*Cons.nu*Vars.v*(2.0*Cons.E_0 + Vars.v**2*(4.0*Cons.E_2 + Vars.v*(5.0*Vars.E_SO_3 + Vars.v*(6.0*Cons.E_4 + 6.0*Vars.E_SQ_4 + Vars.v*(7.0*Vars.E_SO_5 + Vars.v*(8.0*Cons.E_6 + 8.0*Vars.E_SQ_6 + Vars.v*(9.0*Vars.E_SO_7 + 9.0*Vars.E_SQ_7 + Vars.v*(10.0*Cons.E_8 + Cons.E_lnv_8*(10.0*Vars.logv + 1.0) + Vars.v**2*(12.0*Cons.E_10 + Cons.E_lnv_10*(12.0*Vars.logv + 1.0) + Vars.v*(13.0*Cons.E_11 + Vars.v*(14.0*Cons.E_12 + Cons.E_lnv_12*(14.0*Vars.logv + 1.0))))))))))))
    Absorption = Vars.Fcal_coeff*Vars.v**5*(Vars.MDot_Alvi_5 + Cons.MDot_Alvi_8*Vars.v**3)
    dvdt_T1 = (-Absorption - Flux)/dEdV
    dydt=np.zeros(8)
    [dydt[5],dydt[6],dydt[7],dydt[8]] = FrameFromAngularVelocityIntegrand(Vars.R, OmegaVec_6p0(Cons,Vars)[1:])
    dydt[0] = dvdt_T5[0]
    if(Cons.EvolveSpin1):
        dydt[1], dydt[2]=FrameFromAngularVelocity_2D_Integrand(Vars.rfrak_chi1[0], Vars.rfrak_chi1[1],(mul(mul(inverse(Cons.S_chi1),OmegaVec_chiVec_1_6p0(Cons,Vars)),Cons.S_chi1))[1:])
    else:
        dydt[1], dydt[2] = 0.0, 0.0
    if(Cons.EvolveSpin2):
        dydt[3], dydt[4]=FrameFromAngularVelocity_2D_Integrand(Vars.rfrak_chi2[0], Vars.rfrak_chi2[1],(mul(mul(inverse(Cons.S_chi2),OmegaVec_chiVec_2_6p0(Cons,Vars)),Cons.S_chi2))[1:])
    else:
        dydt[3], dydt[4] = 0.0, 0.0
    return dydt      

class PNEv:
    def Integrand(t,y):
        PNEv.Recalculate.get(2*PNEv.PNEvolutionOrder)(PNEv.Cons,PNEv.Vars,y)
        dydt=PNEv.Taylor.get(PNEv.TaylorTn+20*PNEv.PNEvolutionOrder)(PNEv.Cons,PNEv.Vars)
        if PNEv.Vars.v>=1.0 and PNEv.NotBackward:
        #    print("Beyond domain of PN validity, this is a good way to terminate.")
            PNEv.terminal1=False
        if y[0]>=PNEv.v_end or y[0]<=PNEv.v_start:
        #    print("Integration finished.")
            PNEv.terminal1=False
        #if dydt[0]<1.0e-12 and PNEv.NotBackward:
        #    print("v is decreasing, which is not an uncommon way to stop.")
        #    PNEv.terminal2=False
        return dydt
        
    def Evolution(wHat_i, xHat_i, yHat_i, zHat_i, M1_i, M2_i, v_i, S_chi1_i, S_chi2_i, R_i,
        omega_start=None, omega_end=None, t_PNStart=False, t_PNEnd=False, PNEvolutionOrder=3.5, TaylorTn=1, StepsPerOrbit=32, dt=None, tol=1e-8, MinStep=1e-7): 
        # Initialization of constants
        PNEv.terminal1=True
        PNEv.terminal2=True
        PNEv.NotBackward=True
        PNEv.PNEvolutionOrder=PNEvolutionOrder
        PNEv.TaylorTn=TaylorTn
        PNEv.Recalculate={            0:Recalculate_0,
            1:Recalculate_0p50,
            2:Recalculate_1p0,
            3:Recalculate_1p5,
            4:Recalculate_2p0,
            5:Recalculate_2p5,
            6:Recalculate_3p0,
            7:Recalculate_3p5,
            8:Recalculate_4p0,
            9:Recalculate_4p5,
            10:Recalculate_5p0,
            11:Recalculate_5p5,
            12:Recalculate_6p0}
        PNEv.Taylor={
            1:TaylorT1_0,
            11:TaylorT1_0p50,
            21:TaylorT1_1p0,
            31:TaylorT1_1p5,
            41:TaylorT1_2p0,
            51:TaylorT1_2p5,
            61:TaylorT1_3p0,
            71:TaylorT1_3p5,
            81:TaylorT1_4p0,
            91:TaylorT1_4p5,
            101:TaylorT1_5p0,
            111:TaylorT1_5p5,
            121:TaylorT1_6p0,
            4:TaylorT4_0,
            14:TaylorT4_0p50,
            24:TaylorT4_1p0,
            34:TaylorT4_1p5,
            44:TaylorT4_2p0,
            54:TaylorT4_2p5,
            64:TaylorT4_3p0,
            74:TaylorT4_3p5,
            84:TaylorT4_4p0,
            94:TaylorT4_4p5,
            104:TaylorT4_5p0,
            114:TaylorT4_5p5,
            124:TaylorT4_6p0,
            5:TaylorT5_0,
            15:TaylorT5_0p50,
            25:TaylorT5_1p0,
            35:TaylorT5_1p5,
            45:TaylorT5_2p0,
            55:TaylorT5_2p5,
            65:TaylorT5_3p0,
            75:TaylorT5_3p5,
            85:TaylorT5_4p0,
            95:TaylorT5_4p5,
            105:TaylorT5_5p0,
            115:TaylorT5_5p5,
            125:TaylorT5_6p0}
        if omega_start is None:
            omega_start=0.0
        if omega_end is None:
            omega_end=1.0/(M1_i+M2_i)
        PNEv.v_start=(omega_start*(M1_i+M2_i))**(1/3)
        PNEv.v_end=(omega_end*(M1_i+M2_i))**(1/3)
        z=np.array([0.0])
        kappa = 1.0
        PNEv.Cons=Cons(z,z,z,z,z,z,z,z,z,z,z,z,z,z,z,z,z,z,z,z,z,z,z,z,z,z,z,z,z,z,z,z,z,z,z,z,z,z,z,z,z,z,z,z,z,z,z,z,z,z,z,z,z,True,True)
        PNEv.Vars=Vars(z,z,z,z,z,z,z,z,z,z,z,z,z,z,z,z,z,z,z,z,z,z,z,z,z,z,z,z,z,z,z,z,z,z,z,z,z,z,z,z,z,z,z,z,z,z,z,z,z,z,z,z,z)
        Initialization(PNEv.Cons, wHat_i, xHat_i, yHat_i, zHat_i, M1_i, M2_i, v_i, S_chi1_i, S_chi2_i, kappa, kappa, kappa, kappa)
    
        def terminate(t,y):
            return 1.0*PNEv.terminal1*PNEv.terminal2
        terminate.terminal=True
        TMerger=5.0/(256.0*PNEv.Cons.nu*v_i**8)
        TEnd=TMerger
        if t_PNEnd:
            TEnd=t_PNEnd
        time=[0.0]
        if dt is None:
            while time[-1]<TEnd and 2*PNEv.Cons.M*(256*PNEv.Cons.nu*(TMerger-time[-1])/5)**(3/8)/StepsPerOrbit>MinStep:
                time.append(time[-1]+(2*PNEv.Cons.M*(256*PNEv.Cons.nu*(TMerger-time[-1])/5)**(3/8)/StepsPerOrbit)[0])
        else:
            while time[-1]<TEnd:
                time.append(time[-1]+dt)
        time=np.delete(time, -1)
       
        # Integrate
        try:
            yy=solve_ivp(PNEv.Integrand, [time[0],time[-1]], [v_i,0.0,
                0.0,0.0,0.0,R_i[0],R_i[1],R_i[2],R_i[3]], method='DOP853',
                t_eval=time, dense_output=True, events=terminate, rtol=tol, atol=tol)
        except:
            PNEv.terminal1=True
            PNEv.terminal2=True
            yy=solve_ivp(PNEv.Integrand, [time[0],time[-1]], [v_i,0.0,
                0.0,0.0,0.0,R_i[0],R_i[1],R_i[2],R_i[3]], method='DOP853',
                events=terminate, rtol=tol, atol=tol)
            PNEv.terminal1=True
            PNEv.terminal2=True
            time=time[time<yy.t[-1]]
            yy=solve_ivp(PNEv.Integrand, [time[0],time[-1]], [v_i,0.0,
                0.0,0.0,0.0,R_i[0],R_i[1],R_i[2],R_i[3]], method='DOP853',
                t_eval=time, dense_output=True, events=terminate, rtol=tol, atol=tol)
        if omega_start!=0.0 or t_PNStart:
            PNEv.terminal1=True
            PNEv.terminal2=True
            PNEv.NotBackward=False
            time=[0.0]
            TStart=-3*TMerger
            if t_PNStart:
                TStart=t_PNStart
            if omega_start!=0.0:
                TStart=TMerger-10.0/(256*PNEv.Cons.nu)*(np.pi/(omega_start*PNEv.Cons.M))**(8/3)
            if dt is None:
                while time[-1]>TStart:
                    time.append(time[-1]-(2*PNEv.Cons.M*(256*PNEv.Cons.nu*(TMerger-time[-1])/5)**(3/8)/StepsPerOrbit)[0])
            else:
                while time[-1]>TStart:
                    time.append(time[-1]-dt)
            time=np.array(time)
            try:
                yyForward=solve_ivp(PNEv.Integrand, [time[0],time[-1]], [v_i,0.0,
                    0.0,0.0,0.0,R_i[0],R_i[1],R_i[2],R_i[3]], method='DOP853',
                    t_eval=time, dense_output=True, events=terminate, rtol=tol, atol=tol)
            except:
                PNEv.terminal1=True
                PNEv.terminal2=True
                yyForward=solve_ivp(PNEv.Integrand, [time[0],time[-1]], [v_i,0.0,
                    0.0,0.0,0.0,R_i[0],R_i[1],R_i[2],R_i[3]], method='DOP853',
                    events=terminate, rtol=tol, atol=tol)
                PNEv.terminal1=True
                PNEv.terminal2=True
                time=time[time>yyForward.t[-1]]
                yyForward=solve_ivp(PNEv.Integrand, [time[0],time[-1]], [v_i,0.0,
                    0.0,0.0,0.0,R_i[0],R_i[1],R_i[2],R_i[3]], method='DOP853',
                    t_eval=time, dense_output=True, events=terminate, rtol=tol, atol=tol)
            yy.t=np.append(yyForward.t[1:][::-1],yy.t)
            data=np.empty((9,len(yy.t)))
            for i in range(9):
                data[i]=np.append(yyForward.y[i][1:][::-1],yy.y[i])
            yy.y=data
             
        return yy
