import scri
from scipy.fft import fft, ifft, fftfreq
from scipy.ndimage import gaussian_filter
from scipy.signal import hilbert
from PYPostNewtonian.Code import PostNewtonian
import PNBMS
import numpy as np
import math
import quaternion
import h5py
import time
import matplotlib
matplotlib.use('AGG')
import matplotlib.pyplot as plt
from scipy.optimize import least_squares
from scipy.optimize import basinhopping, minimize, Bounds
from scipy.integrate import simpson
from scipy.interpolate import InterpolatedUnivariateSpline as Spline
from scipy.interpolate import CubicSpline
from scipy.integrate import solve_ivp
import copy
from memory_profiler import profile

class SplineArray:
    def __init__(self, x, y):
        self.complex = np.iscomplexobj(y)
        if self.complex:
            y = y.view(dtype=float)
        self.splines = [Spline(x, y[:, i]) for i in range(y.shape[1])]
    def __call__(self, xprime):
        yprime = np.concatenate([spl(xprime)[:, np.newaxis] for spl in self.splines], axis=1)
        if self.complex:
            yprime = yprime.view(dtype=complex)
        return yprime

def nOrbits_to_length(nOrbits, t_end, Omega_array, time_array):
    tol=1e-3
    number=0
    a=time_array[0]
    b=t_end
    t0=(a+b)/2
    orbits=simpson(Omega_array[(time_array>t0)&(time_array<t_end)],
        time_array[(time_array>t0)&(time_array<t_end)])/2/np.pi
    while abs(orbits-nOrbits)>tol and number<=100:
        orbits=simpson(Omega_array[(time_array>t0)&(time_array<t_end)],
        time_array[(time_array>t0)&(time_array<t_end)])/2/np.pi
        if orbits<nOrbits:
            b=t0
        else:
            a=t0
        t0=(a+b)/2
        number+=1
    if number>100:
        message=("The waveform before endtime only has {0} orbits, which is smaller than required {1}.")
        raise ValueError(message.format(orbits, norbits))
    return t_end-t0

def Physical_to_Parameterize(x):
    """
    Reparameterize physical parameters x=[q, M, chi1_x, chi1_y, chi1_z, chi2_x, chi2_y, chi2_z, t_delta,
    logR_delta_x, logR_delta_y, logR_delta_z] for optimization. The last four components of x are optional.
    """
    if len(x)!=12 and len(x)!=8 and np.linalg.norm(np.cross(chi1_i,chi2_i))<1e-8:
        message=("Incorrect format for x={0}. x should be [q, M, chi1_x, chi1_y, chi1_z, chi2_x, chi2_y, chi2_z, t_delta,"
            +"logR_delta_x, logR_delta_y, logR_delta_z], where the last four components are optional.")
        raise ValueError(message.format(x))
    x_2=np.linalg.norm(x[2:5])/np.linalg.norm(chi1_i)
    x_3=np.linalg.norm(x[5:8])/np.linalg.norm(chi2_i)
    rotation=quaternion.optimal_alignment_in_Euclidean_metric(x[2:5],x_2*chi1_i)
    axis=quaternion.from_float_array(np.append(0,x[2:5])).normalized()
    chi1_hat=x[2:5]/np.linalg.norm(x[2:5])
    chi2_i_proj=x_3*(rotation*quaternion.quaternion(0,chi2_i[0],chi2_i[1],chi2_i[2])*rotation.conjugate()).vec
    chi2_i_proj=chi2_i_proj-np.dot(chi2_i_proj,chi1_hat)*chi1_hat
    chi2_0_proj=x[5:8]-np.dot(x[5:8],chi1_hat)*chi1_hat
    angle=np.arccos(np.dot(chi2_i_proj,chi2_0_proj)/(np.linalg.norm(chi2_i_proj)*np.linalg.norm(chi2_0_proj)))
    if np.dot(chi2_i_proj,chi2_0_proj)/(np.linalg.norm(chi2_i_proj)*np.linalg.norm(chi2_0_proj))>1:
        angle=0.0
    elif np.dot(chi2_i_proj,chi2_0_proj)/(np.linalg.norm(chi2_i_proj)*np.linalg.norm(chi2_0_proj))<-1:
        angle=np.pi
    sign=np.sign(np.dot(chi1_hat,np.cross(chi2_i_proj,chi2_0_proj)))
    rotation=np.exp(sign*axis*angle/2)*rotation
    if len(x)==12:
        x_47=(quaternion.as_float_array(np.log(rotation))-np.append(0.0, omega_mean*x[8]/2))[1:]
    else:
        x_47=(quaternion.as_float_array(np.log(rotation)))[1:]
    x_47[x_47>np.pi]=x_47[x_47>np.pi]-np.pi
    chi2_0=x_3*(rotation*quaternion.quaternion(0,chi2_i[0],chi2_i[1],chi2_i[2])*rotation.conjugate()).vec
    sign=np.sign(np.dot(np.cross(chi2_0,chi1_hat),np.cross(chi2_0,x[5:8])))
    x_7=sign*np.arccos(np.dot(chi2_0,x[5:8])/(np.linalg.norm(x[5:8])*np.linalg.norm(chi2_0)))
    if np.dot(chi2_0,x[5:8])/(np.linalg.norm(x[5:8])*np.linalg.norm(chi2_0))>1:
        x_7=0.0
    x[0]=x[0]/q_0
    x[1]=x[1]/Mc_0
    x[2]=x_2
    x[3]=x_3
    x[4:7]=x_47
    x[7]=x_7
    if len(x)==12:
        x[9:]=x[9:]-omega_mean*x[8]/2
    return x
        
def Parameterize_to_Physical(x):
    """
    Output optimization parameters as physical quantities x=[q, M, chi1_x, chi1_y, chi1_z, chi2_x, chi2_y, chi2_z, t_delta,
    logR_delta_x, logR_delta_y, logR_delta_z]. The last four components of x are optional.
    """
    if len(x)!=12 and len(x)!=8:
        message=("Incorrect format for x={0}")
        raise ValueError(message.format(x))
    x[0]=x[0]*q_0
    x[1]=x[1]*Mc_0
    if len(x)==12:
        phase=quaternion.quaternion(0.0, omega_mean[0]*x[8]/2, omega_mean[1]*x[8]/2, omega_mean[2]*x[8]/2)
        rotation=np.exp(xHat*x[4]+yHat*x[5]+zHat*x[6]+phase)
    else:
        rotation=np.exp(xHat*x[4]+yHat*x[5]+zHat*x[6])
    chi1_0=x[2]*(rotation*quaternion.quaternion(0,chi1_i[0],chi1_i[1],chi1_i[2])*rotation.conjugate()).vec
    chi2_0=x[3]*rotation*quaternion.quaternion(0,chi2_i[0],chi2_i[1],chi2_i[2])*rotation.conjugate()
    axis=quaternion.from_float_array(np.append(0,np.cross(chi2_0.vec,chi1_0))).normalized()
    angle=np.exp(axis*x[7]/2)
    chi2_0=(angle*chi2_0*angle.conjugate()).vec
    x[2:5]=chi1_0
    x[5:8]=chi2_0
    if len(x)==12:
        x[9:]=x[9:]+omega_mean*x[8]/2
    return x

def InitialT(x):
    return simpson((omega_NR_mag_matching-omega_PN_mag_spline(matchingt+x))**2, matchingt)/NormalizationOmega

def InitialR(theta):
    R_temp=R_delta*np.exp(theta/2*quaternion.quaternion(0.0, omega_NR_hat[0,0], omega_NR_hat[0,1],omega_NR_hat[0,2]))
    W_temp=scri.rotate_physical_system(W_NR_matching_in.copy(), R_temp)
    cost=(np.angle(W_temp.data[int(len(matchingt)/2),4])\
        -np.angle(PNData_spline(mint.x+W_temp.t[int(len(matchingt)/2)])[0,4]))**2
    return cost

def Optimize4D(x):
    phase=quaternion.quaternion(0.0, omega_mean[0]*x[0]/2, omega_mean[1]*x[0]/2, omega_mean[2]*x[0]/2)
    R_delta=np.exp(quaternion.quaternion(0.0,x[1],x[2],x[3])+phase)
    W_temp=scri.rotate_physical_system(W_NR_matching_in.copy(), R_delta)
    temp=PNData_spline(matchingt+x[0])
    h1h2=np.sum(simpson(abs(temp-W_temp.data)**2.0, matchingt, axis=0))
    cost=0.5*h1h2/Normalization
    return 10.0+np.log10(cost)

def Mismatch(W1,W2,t1,t2):
    """
    Calculate the mismatch of W1 and W2 between t1 and t2.
    """
    W2_spline=SplineArray(W2.t, W2.data)
    matchingtt=W1.t[(W1.t>=t1)&(W1.t<=t2)]
    h1h2=np.sum(simpson(W2_spline(matchingtt)*np.conjugate(W1.data[(W1.t>=t1)&(W1.t<=t2),:]), matchingtt, axis=0))
    h1h1=np.sum(simpson(W2_spline(matchingtt)*np.conjugate(W2_spline(matchingtt)), matchingtt, axis=0))
    h2h2=np.sum(simpson(W1.data[(W1.t>=t1)&(W1.t<=t2),:]*np.conjugate(W1.data[(W1.t>=t1)&(W1.t<=t2),:]), matchingtt, axis=0))
    return 1-h1h2/np.sqrt(h1h1*h2h2)

def SquaredError(W1,W2,t1,t2,mode=None):
    """
    Calculate the residue of W1 and W2 between t1 and t2.
    """
    W2_spline=SplineArray(W2.t, W2.data)
    matchingtt=W1.t[(W1.t>=t1)&(W1.t<=t2)]
    h1h2=np.sum(simpson(abs(W2_spline(matchingtt)-W1.data[(W1.t>=t1)&(W1.t<=t2),:])**2.0, matchingtt, axis=0))
    h1h1=np.sum(simpson(abs(W1.data[(W1.t>=t1)&(W1.t<=t2),:])**2.0, matchingtt, axis=0))
    if type(mode) != type(None):
        h1h2=np.sum(simpson(abs(W2_spline(matchingtt)-W1.data[(W1.t>=t1)&(W1.t<=t2),:])[:,mode]**2.0, matchingtt, axis=0))
    return 0.5*h1h2/h1h1

def SquaredErrorNorm(W1,W2,t1,t2,mode=None):
    """
    Calculate the residue of W1 and W2 between t1 and t2.
    """
    W2_spline=SplineArray(W2.t, W2.data)
    matchingtt=W1.t[(W1.t>=t1)&(W1.t<=t2)]
    h1h2=np.sum(simpson(abs(np.abs(W2_spline(matchingtt))-np.abs(W1.data[(W1.t>=t1)&(W1.t<=t2),:]))**2.0, matchingtt, axis=0))
    h1h1=np.sum(simpson(abs(W1.data[(W1.t>=t1)&(W1.t<=t2),:])**2.0, matchingtt, axis=0))
    if mode != None:
        h1h2=np.sum(simpson(abs(np.abs(W2_spline(matchingtt))-np.abs(W1.data[(W1.t>=t1)&(W1.t<=t2),:]))[:,mode]**2.0, matchingtt, axis=0))
    return 0.5*h1h2/h1h1

def SquaredErrorScalar(t_NR,t_PN,f_NR,f_PN,t1,t2):
    fPN_spline=Spline(t_PN, f_PN)
    matchingtt=t_NR[(t_NR>=t1)&(t_NR<=t2)]
    f1f2=np.sum(simpson(abs(fPN_spline(matchingtt)-f_NR[(t_NR>=t1)&(t_NR<=t2)])**2.0, matchingtt, axis=0))
    f1f1=np.sum(simpson(abs(f_NR[(t_NR>=t1)&(t_NR<=t2)])**2.0, matchingtt, axis=0))
    return 0.5*f1f2/f1f1

def Align(x):
    """
    Generate PN waveform and align it with NR waveform.
    x=[q, chi1_x, chi1_y, chi1_z, chi2_x, chi2_y, chi2_z]
    """
    clock1=time.time()
    global mint, minima, R_delta, W_PN, W_PN_corot, t_end0, t_pre, t_end, omega_0, omega_PN,\
        omega_PN_spline, omega_PN_mag, omega_PN_mag_spline, PNData_spline, matchingt,\
        omega_NR_mag_matching, omega_mean, iter_num, W_NR_matching_in, omega_NR_hat,\
        Normalization, lowbound, upbound, scale, spin1, spin2,t_start, chiA, chiB, chi1_PN_spline, chi2_PN_spline,q,PNNorm_spline,cost1,cost2,cost3,cost4,cost5
    t_start=t_starts
    iter_num+=1
    
    Phys=Parameterize_to_Physical(np.copy(x))
    print(("Call # {4}, generating PN with parameters q={0}, M={7}, omega_0={1}, chi1_0={2}, chi2_0={3},"
        +"t_PNstart={5}, t_PNend={6}.").format(Phys[0], omega_00,Phys[2:5], Phys[5:8],iter_num, t_PNStart, t_PNEnd, Phys[1]))
    W_PN_corot=PostNewtonian.PNWaveform(Phys[0], omega_00.copy()*Phys[1], Phys[2:5], Phys[5:8], frame_0, t_start.copy()/Phys[1], t_PNStart, t_PNEnd)
    if PNIter==0:
        ZeroModes=[2,8,16,26,38,52,68] # Memory modes
        W_PN_corot.data[:,ZeroModes]=0.0*W_PN_corot.data[:,ZeroModes] # Not cosider memory effect since NR dosen't have corrrect memory.
    W_PN=scri.to_inertial_frame(W_PN_corot.copy())
    W_PN.t=W_PN.t*Phys[1]
    # Set up the matching region data for PN, and get the corresponding angular velocity and frame
    omega_PN=W_PN.angular_velocity()
    omega_PN_spline=SplineArray(W_PN.t, omega_PN)
    omega_PN_mag=gaussian_filter(np.linalg.norm(omega_PN, axis=1),100)
    omega_PN_mag_spline=Spline(W_PN.t, omega_PN_mag)
    PNData_spline=SplineArray(W_PN.t, W_PN.data)
    PNNorm_spline=Spline(W_PN.t,gaussian_filter(np.linalg.norm(W_PN.data,axis=1),100))

    # Get initial guess
    Initial1=np.array([0.0,0.0,0.0,0.0])
    Initial2=np.array([0.0,0.0,0.0,0.0])
    if len(Phys)==12:
        Initial1[0]=Phys[8]
        Initial1[1:]=Phys[9:]-omega_mean*Phys[8]/2
        R_delta=np.exp(quaternion.quaternion(0.0,Phys[9],Phys[10],Phys[11]))
        R_delta=R_delta*np.exp(np.pi/2*quaternion.quaternion(0.0, omega_NR_hat[0,0], omega_NR_hat[0,1],omega_NR_hat[0,2]))
        Initial2[0]=Phys[8]
        Initial2[1:]=quaternion.as_float_array(np.log(R_delta))[1:]-omega_mean*Phys[8]/2
    else:
        # Get initial guess of time alignment by matching angular velocity
        mint=least_squares(InitialT, 0.0, bounds=[-10*np.pi/omega_0,10*np.pi/omega_0],ftol=1e-14, xtol=1e-14, gtol=1e-14)
        # Get initial guess of frame alignment
        R_delta = quaternion.optimal_alignment_in_Euclidean_metric(omega_NR[(W_NR.t>=t_start)&(W_NR.t<=t_end0)],\
            omega_PN_spline(matchingt+mint.x), matchingt)
        minf=least_squares(InitialR, 0.0, bounds=[-np.pi,np.pi],ftol=1e-14, xtol=1e-14, gtol=1e-14)
        # Pi degeneracy
        phase=quaternion.quaternion(0.0, omega_mean[0]*mint.x/2, omega_mean[1]*mint.x/2, omega_mean[2]*mint.x/2)
        R_delta1=R_delta*np.exp(minf.x/2*quaternion.quaternion(0.0, omega_NR_hat[0,0], omega_NR_hat[0,1],omega_NR_hat[0,2]))
        R_delta2=R_delta1*np.exp(np.pi/2*quaternion.quaternion(0.0, omega_NR_hat[0,0], omega_NR_hat[0,1],omega_NR_hat[0,2]))
        logR_delta1=quaternion.as_float_array(np.log(R_delta1))-np.append(0.0, omega_mean*mint.x/2)
        logR_delta2=quaternion.as_float_array(np.log(R_delta2))-np.append(0.0, omega_mean*mint.x/2)
        Initial1[0]=mint.x[0]
        Initial2[0]=mint.x[0]
        Initial1[1:]=logR_delta1[0][1:]
        Initial2[1:]=logR_delta2[0][1:]

    # Alignment of time and frame
    scale=[np.pi/omega_0/2,np.pi/4.0,np.pi/4.0,np.pi/4.0]
    lowbound1=Initial1-scale
    upbound1=Initial1+scale
    lowbound2=Initial2-scale
    upbound2=Initial2+scale
    minima=least_squares(Optimize4D, Initial1,bounds=(lowbound1,upbound1),ftol=1e-16, xtol=1e-16, gtol=1e-14,x_scale='jac',max_nfev=2000)
    minima2=least_squares(Optimize4D, Initial2,bounds=(lowbound2,upbound2),ftol=1e-16, xtol=1e-16, gtol=1e-14,x_scale='jac',max_nfev=2000)
    cost1=minima.fun
    cost2=minima2.fun
    if cost1<=cost2:
        cost=cost1
    else:
        cost=cost2
        minima.x=minima2.x
        minima.fun=minima2.fun
    print(minima)
    return cost

def Optimize11D(x):
    """
    Generate PN waveform and align it with NR waveform.
    x=[q, chi1_x, chi1_y, chi1_z, chi2_x, chi2_y, chi2_z]
    """
    clock1=time.time()
    global mint, minima, R_delta, W_PN, W_PN_corot, t_end0, t_pre, t_end, omega_0, omega_PN,\
        omega_PN_spline, omega_PN_mag, omega_PN_mag_spline, PNData_spline, matchingt,\
        omega_NR_mag_matching, omega_mean, iter_num, W_NR_matching_in, omega_NR_hat,\
        Normalization, lowbound, upbound, scale, spin1, spin2,t_start, chiA, chiB, chi1_PN_spline, chi2_PN_spline,cost1,cost2,cost3,cost4,cost5,PNNorm_spline
    t_start=t_starts
    iter_num+=1
    
    Phys=Parameterize_to_Physical(np.copy(x))
    W_PN_corot=PostNewtonian.PNWaveform(Phys[0], omega_00.copy()*Phys[1], Phys[2:5], Phys[5:8], frame_0, t_start.copy()/Phys[1], t_PNStart, t_PNEnd)
    if PNIter==0:
        ZeroModes=[2,8,16,26,38,52,68] # Memory modes
        W_PN_corot.data[:,ZeroModes]=0.0*W_PN_corot.data[:,ZeroModes] # Not cosider memory effect since NR dosen't have corrrect memory.
    W_PN=scri.to_inertial_frame(W_PN_corot.copy())
    W_PN.t=W_PN.t*Phys[1]

    # Set up the matching region data for PN, and get the corresponding angular velocity and frame
    omega_PN=W_PN.angular_velocity()
    omega_PN_spline=SplineArray(W_PN.t, omega_PN)
    omega_PN_mag=gaussian_filter(np.linalg.norm(omega_PN, axis=1),100)
    omega_PN_mag_spline=Spline(W_PN.t, omega_PN_mag)
    PNData_spline=SplineArray(W_PN.t, W_PN.data)
    PNNorm_spline=Spline(W_PN.t,gaussian_filter(np.linalg.norm(W_PN.data,axis=1),100))

    cost=Optimize4D(x[8:])

    return cost

def Optimize1D(x):
    return Optimize11D(PNParas+x*direction)

#@profile
def Hybridize(t_end00, data_dir, cce_dir, out_dir, length, nOrbits, debug=0, OptimizePNParas=0, truncate=None):
    """
    Align and hybridize given NR waveform with PN waveform.
    """
    clock0=time.time()
    global mint, minima, W_NR, W_PN, W_PN_corot, t_starts, t_end0, t_pre, t_end, omega_0, omega_PN,\
        omega_NR, omega_PN_spline, omega_PN_mag, omega_PN_mag_spline, PNData_spline, matchingt,\
        omega_NR_mag_matching, omega_mean, iter_num, W_NR_matching_in, R_delta, omega_NR_hat,\
        Normalization, xHat, yHat, zHat, t_PNStart, t_PNEnd, lowbound, upbound, scale, frame_0,\
        spin1,spin2,chiA,chiB,NormalizationChi,chi1_i,chi2_i,q_0,NormalizationOmega,omega_NR_matching,t_delta,Mc_0,cc,coeff2,coeff3,coeff4,coeff5,WaveformNorm_NR_matching,PNParas,direction,PNIter,omega_00,PhyParas,length_global,Output
    print('Iteration #',PNIter,'between BMS and PN parameters fixing.')
    t_end0=t_end00
    length=length_global
    t_start=t_end0-length
# Get NR waveform
    abd=scri.SpEC.create_abd_from_h5("CCE",h=cce_dir+'/Strain.h5',Psi4=cce_dir+'/Psi4.h5',Psi3=cce_dir+'/Psi3.h5',Psi2=cce_dir+'/Psi2.h5',Psi1=cce_dir+'/Psi1.h5',Psi0=cce_dir+'/Psi0.h5')
    t0=-abd.t[np.argmax(np.linalg.norm(abd.sigma.bar,axis=1))]
    abd.t=abd.t+t0
    if truncate!=None:
        abd1 = scri.asymptotic_bondi_data.AsymptoticBondiData(time = abd.t[(abd.t>=truncate[0])&(abd.t<truncate[1])],ell_max = 8,multiplication_truncator = max)
        abd1.sigma=abd.sigma[(abd.t>=truncate[0])&(abd.t<truncate[1])]
        abd1.psi0=abd.psi0[(abd.t>=truncate[0])&(abd.t<truncate[1])]
        abd1.psi1=abd.psi1[(abd.t>=truncate[0])&(abd.t<truncate[1])]
        abd1.psi2=abd.psi2[(abd.t>=truncate[0])&(abd.t<truncate[1])]
        abd1.psi3=abd.psi3[(abd.t>=truncate[0])&(abd.t<truncate[1])]
        abd1.psi4=abd.psi4[(abd.t>=truncate[0])&(abd.t<truncate[1])]
        abd=abd1
    if nOrbits!=None:   
        W_temp=scri.WaveformModes()
        W_temp.t=abd.t
        W_temp.data=2*abd.sigma.bar
        W_temp.data=np.copy(W_temp.data[:,4:])
        W_temp.ells=2,8
        omega_NR=W_temp.angular_velocity()
        omega_NR_mag = np.linalg.norm(omega_NR, axis=1)
        length_global=nOrbits_to_length(nOrbits,t_end0,omega_NR_mag,W_temp.t)
        length=length_global
        t_start=t_end0-length
    
    W_NR=scri.WaveformModes()
    if PNIter==0:
        abd1,trans=abd.map_to_superrest_frame(t_0=t_start+length/2)
        W_NR.t=abd1.t
        W_NR.data=2*abd1.sigma.bar
        W_NR.data=np.copy(W_NR.data[:,4:])
        W_NR.ells=2,8
        W_NR_corot=scri.to_corotating_frame(W_NR.copy())
        ZeroModes=[2,8,16,26,38,52,68]
        W_NR_corot.data[:,ZeroModes]=0.0*W_NR_corot.data[:,ZeroModes]
        W_NR=scri.to_inertial_frame(W_NR_corot.copy())
        W_NR_corot=scri.to_corotating_frame(W_NR.copy())
    else:
        W_PN_corot=PostNewtonian.PNWaveform(PhyParas[0], omega_00.copy()*PhyParas[1], PhyParas[2:5], PhyParas[5:8], frame_0, t_start.copy()/PhyParas[1])
        W_PN=scri.to_inertial_frame(W_PN_corot.copy())
        W_PN.t=W_PN.t*PhyParas[1]
        W_PN.data=np.append(0.0*W_PN.data[:,0:4],np.copy(W_PN.data),axis=1)
        W_PN.ells=0,8
        W_PN.dataType=scri.h
        W_PN_PsiM_corot=PostNewtonian.PNWaveform(PhyParas[0], omega_00.copy()*PhyParas[1], PhyParas[2:5], PhyParas[5:8], frame_0, t_start.copy()/PhyParas[1], datatype="Psi_M")
        W_PN_PsiM=scri.to_inertial_frame(W_PN_PsiM_corot.copy())
        W_PN_PsiM.t=W_PN_PsiM.t*PhyParas[1]
        tp1, W_NR, tp2, idx=PNBMS.PN_BMS_w_time_phase(abd,W_PN,W_PN_PsiM,t_start,t_start+length,None)
        W_NR_corot=scri.to_corotating_frame(W_NR.copy())
    print("Map to superrest frame used ",time.time()-clock0)
    
# Get the initial angular velocity in matching region
    omega_NR=W_NR.angular_velocity()
    omega_NR_mag = np.linalg.norm(omega_NR, axis=1)
    omega_NR_mag = gaussian_filter(omega_NR_mag, sigma=100)
    if W_NR.t[0]>=t_start-10 or W_NR.t[-1]<=t_start+10:
        message=("t_start {0} should be much larger than the start time of NR"
            +" waveform {1} and much smaller than the ending time {2}.")
        raise ValueError(message.format(t_start, W_NR.t[0], W_NR.t[-1]))
    omega_0=omega_NR_mag[(W_NR.t>=t_start-10)&(W_NR.t<=t_start+10)]
    omega_0=np.mean(omega_0)
    if PNIter==0:
        omega_00=omega_0
        if nOrbits!=None:  
            length_global=nOrbits_to_length(nOrbits,t_end0,omega_NR_mag,W_NR.t)
    if W_NR.t[0]>=t_start-2.0*np.pi/omega_0 or W_NR.t[-1]<=t_start+2.0*np.pi/omega_0:
        message=("t_start {0} should be at least {1} larger than the start time of NR waveform {2}"
            +" and smaller than the ending time of NR waveform {3}.")
        raise ValueError(message.format(t_start, 2.0*np.pi/omega_0, W_NR.t[0], W_NR.t[-1]))
    t_start=t_end0-length
    t_pre=t_start-1000
    t_end=t_start+length+1000
    t_PNStart=-5000
    t_PNEnd=t_end-t_start
    matchingt=W_NR.t[(W_NR.t>=t_start)&(W_NR.t<=t_end0)]
    omega_NR_mag_matching=omega_NR_mag[(W_NR.t>=t_start)&(W_NR.t<=t_end0)]
    W_NR_matching_in=W_NR.interpolate(matchingt)
    omega_NR_matching = omega_NR[(W_NR.t>=t_start) & (W_NR.t<=t_end0)]
    omega_mean = np.mean(omega_NR_matching, axis=0)
    omega_NR_hat = omega_NR_matching / np.linalg.norm(omega_NR_matching, axis=1)[:, np.newaxis]
    Normalization=np.sum(simpson(abs(W_NR_matching_in.data)**2.0, matchingt, axis=0))
    NormalizationOmega=simpson(np.linalg.norm(omega_NR_matching,axis=1), matchingt)

# Get initial guess of PN Parameters
    xHat=quaternion.quaternion(0.0,1.0,0.0,0.0)
    yHat=quaternion.quaternion(0.0,0.0,1.0,0.0)
    zHat=quaternion.quaternion(0.0,0.0,0.0,1.0)
    with h5py.File(data_dir+'/Horizons.h5', 'r') as f:
        tA=f['AhA.dir/CoordCenterInertial.dat'][:,0]+t0
        xA=f['AhA.dir/CoordCenterInertial.dat'][:,1:]
        mA=f['AhA.dir/ArealMass.dat'][:,1]
        chiA=f['AhA.dir/chiInertial.dat'][:,1:]
        xB=f['AhB.dir/CoordCenterInertial.dat'][:,1:]
        mB=f['AhB.dir/ArealMass.dat'][:,1]
        chiB=f['AhB.dir/chiInertial.dat'][:,1:]
    i_1=abs(tA-t_start).argmin()
    m1=mA[i_1]
    m2=mB[i_1]
    q_0=m1/m2
    if m1<m2:
        q_0=1.0/q_0
    Mc_0=1.0
    d=xA-xB
    nHat=np.empty((len(d),4))
    for i in range(len(d)):
        nHat[i]=np.append(0,d[i])
    nHatArray=nHat
    nHat=quaternion.from_float_array(nHat)
    for i in range(len(d)):
        nHat[i]=nHat[i].normalized()
    dnHatdt=CubicSpline(tA, nHatArray).derivative()(tA)
    lambdaHat=quaternion.from_float_array(dnHatdt)
    for i in range(len(d)):
        lambdaHat[i]=lambdaHat[i].normalized()
    Ra=np.sqrt(-nHat[i_1]*xHat)
    beta=math.atan2(np.dot((Ra*zHat*Ra.inverse()).vec,lambdaHat[i_1].vec),\
        np.dot((Ra*yHat*Ra.inverse()).vec,lambdaHat[i_1].vec))
    frame_0=Ra*np.exp((beta)/2*xHat)
    print(frame_0)
    chi1_i=chiA[i_1]
    chi2_i=chiB[i_1]
    chi_p=max(np.sqrt(np.dot(chi1_i,chi1_i)-(np.dot(chi1_i,omega_NR[W_NR.t>t_start][0])/np.linalg.norm(omega_NR[W_NR.t>t_start][0]))**2),np.sqrt(np.dot(chi2_i,chi2_i)-(np.dot(chi2_i,omega_NR[W_NR.t>t_start][0])/np.linalg.norm(omega_NR[W_NR.t>t_start][0]))**2)*1/q_0*(4/q_0+3)/(4+3/q_0))
    print("chi_p = ",chi_p,np.sqrt(np.dot(chi1_i,chi1_i)))
    chiA_spline=SplineArray(tA, chiA)
    chiB_spline=SplineArray(tA, chiB)
    chiA=chiA_spline(matchingt)
    chiB=chiB_spline(matchingt)
    NormalizationChi=simpson(np.linalg.norm(chiA,axis=1)**2.0, matchingt)+simpson(np.linalg.norm(chiB,axis=1)**2.0, matchingt)
    chiA=quaternion.from_float_array(np.column_stack((0.0*matchingt,chiA)))
    chiB=quaternion.from_float_array(np.column_stack((0.0*matchingt,chiB)))

# Align Waveforms
    t_starts=t_start
    if PNIter==0:
        PNParas=np.array([1.0,1.0,1.0,1.0,0.0,0.0,0.0,0.0])
        PhyParas=Parameterize_to_Physical(np.copy(PNParas))
    else:
        PNParas=Physical_to_Parameterize(np.copy(PhyParas))
    iter_num=0
    if OptimizePNParas:
        COST=Align(PNParas)
        logR_delta=np.append(minima.x[0],minima.x[1:]+omega_mean*minima.x[0]/2)#quaternion.as_float_array(np.log(R_delta))
        if len(PhyParas)==12:
            PhyParas[8:]=logR_delta
        else:
            PhyParas=np.append(PhyParas,logR_delta)
        PNParas=Physical_to_Parameterize(np.copy(PhyParas))
        lowbound12D=PNParas-[0.2,0.1,0.2,0.2,np.pi/4,np.pi/4,np.pi/4,np.pi/4,np.pi/omega_0/2.0,np.pi/4,np.pi/4,np.pi/4]
        upbound12D=PNParas+[0.2,0.1,0.2,0.2,np.pi/4,np.pi/4,np.pi/4,np.pi/4,np.pi/omega_0/2.0,np.pi/4,np.pi/4,np.pi/4]
        minima12D=minimize(Optimize11D,PNParas,method='Nelder-Mead',bounds=Bounds(lowbound12D,upbound12D),options={'return_all':True,'xatol': 3e-16, 'fatol': 1e-15,'maxfev':2000,'adaptive':True})
        #minima12D=least_squares(Optimize11D, minima12D.x,bounds=(lowbound12D,upbound12D),ftol=3e-16, xtol=3e-16, gtol=1e-15,x_scale='jac',max_nfev=5000)
        #print(minima12D)
        PNParas=minima12D.x

    t_PNStart, t_PNEnd=-80000, -500-t_start
    Align(PNParas)
    t_delta=minima.x[0]
    print("Time shift=", t_delta)
    phase=quaternion.quaternion(0.0, omega_mean[0]*minima.x[0]/2,omega_mean[1]*minima.x[0]/2, omega_mean[2]*minima.x[0]/2)
    R_delta=np.exp(quaternion.quaternion(0.0,minima.x[1],minima.x[2],minima.x[3])+phase)
    print("R_delta=",R_delta)
    W_PN.t=W_PN.t-t_delta
    W_NR=scri.rotate_physical_system(W_NR, R_delta)
    f1=np.linalg.norm(W_NR.angular_velocity(), axis=1)
    f2=np.linalg.norm(W_PN.angular_velocity(), axis=1)
    modes = np.arange(len(W_NR.data[0,:]))
    modes = np.delete(modes,2)
    print("SquaredError over matching window: ",SquaredError(W_NR,W_PN,t_start,t_start+length), " ", SquaredError(W_NR,W_PN,t_start,t_start+length,mode=modes), " ",SquaredErrorNorm(W_NR,W_PN,t_start,t_start+length)," ",SquaredErrorScalar(W_NR.t,W_PN.t,f1,f2,t_start,t_start+length))
    Output=SquaredError(W_NR,W_PN,t_start,t_start+length)
    """
    t1=t_end0-nOrbits_to_length(nOrbits/2,t_end0,omega_NR_mag,W_NR.t)-nOrbits_to_length(25,t_end0,omega_NR_mag,W_NR.t)
    length=nOrbits_to_length(10,t1,omega_NR_mag,W_NR.t)#################################################################################
    Output=np.append(Output,np.array([SquaredError(W_NR,W_PN,t1-1.0*length,t1),SquaredError(W_NR,W_PN,t1-4.0*length,t1)]))
    print("SquaredError over longer region: ","t=40 ",t1-4.0*length," ",SquaredError(W_NR,W_PN,t1-4.0*length,t1), " ",SquaredError(W_NR,W_PN,t1-4.0*length,t1,mode=modes), " ", SquaredErrorNorm(W_NR,W_PN,t1-4.0*length,t1)," ",SquaredErrorScalar(W_NR.t,W_PN.t,f1,f2,t1-4.0*length,t1))
    print("SquaredError over longer region: ","t=30 ",t1-3.0*length," ",SquaredError(W_NR,W_PN,t1-3.0*length,t1), " ",SquaredError(W_NR,W_PN,t1-3.0*length,t1,mode=modes), " ", SquaredErrorNorm(W_NR,W_PN,t1-3.0*length,t1)," ",SquaredErrorScalar(W_NR.t,W_PN.t,f1,f2,t1-3.0*length,t1))
    print("SquaredError over longer region: ","t=20 ",t1-2.0*length," ",SquaredError(W_NR,W_PN,t1-2.0*length,t1), " ",SquaredError(W_NR,W_PN,t1-2.0*length,t1,mode=modes), " ", SquaredErrorNorm(W_NR,W_PN,t1-2.0*length,t1)," ",SquaredErrorScalar(W_NR.t,W_PN.t,f1,f2,t1-2.0*length,t1))
    print("SquaredError over longer region: ","t=10 ",t1-1.0*length," ",SquaredError(W_NR,W_PN,t1-1.0*length,t1), " ",SquaredError(W_NR,W_PN,t1-1.0*length,t1,mode=modes), " ",SquaredErrorNorm(W_NR,W_PN,t1-1.0*length,t1)," ",SquaredErrorScalar(W_NR.t,W_PN.t,f1,f2,t1-1.0*length,t1))
    print("SquaredError of (2,2) mode: ",SquaredError(W_NR,W_PN,t_start,t_start+length,mode=4))
    print("SquaredError (2,1) mode: ",SquaredError(W_NR,W_PN,t_start,t_start+length,mode=3))
    print("SquaredError (2,0) mode: ",SquaredError(W_NR,W_PN,t_start,t_start+length,mode=2))
    for i in range(len(W_NR.data[0,:])):
        print(i,SquaredError(W_NR,W_PN,t_start,t_start+length,mode=i)," ",SquaredError(W_NR,W_PN,t1-1.0*length,t1,mode=i)," ",SquaredError(W_NR,W_PN,t1-2.0*length,t1,mode=i)," ",SquaredError(W_NR,W_PN,t1-3.0*length,t1,mode=i)," ",SquaredError(W_NR,W_PN,t1-4.0*length,t1,mode=i))
    length=length_global###################################################################################
    """
    PhyParas=Parameterize_to_Physical(np.copy(PNParas))
    logR_delta=np.append(minima.x[0],minima.x[1:]+omega_mean*minima.x[0]/2)
    if len(PhyParas)==12:
        PhyParas[8:]=logR_delta
    else:
        PhyParas=np.append(PhyParas,logR_delta)
    print("Physical PN parameters are: ",PhyParas)

    if debug:
        plt.subplot(211)
        plt.plot(W_NR.t, omega_NR_mag, label='Angular velocity')
        plt.plot(W_NR.t, omega_PN_mag_spline(W_NR.t+t_delta), label='Angular velocity')
        plt.axvline(t_start, linestyle='dotted')
        plt.axvline(t_end0, linestyle='dotted')
        plt.legend(['NR', 'PN'], loc="upper right")
        plt.ylabel("Omega Magnititude")
        plt.xlabel("Time")
        #plt.xlim((-95000,-2000))
        plt.ylim((0.005,0.006))
        plt.subplot(212)
        plt.plot(W_NR.t, np.abs(omega_NR_mag-omega_PN_mag_spline(W_NR.t+t_delta)))
        plt.yscale('log')
        plt.savefig(out_dir+"/hybridCheckOmega"+data_dir[-8:-5],dpi=1000)
        plt.clf()

    # Hybridize waveform
    PNData_spline=SplineArray(W_PN.t, W_PN.data)
    tTemp=np.array(np.append(W_PN.t[W_PN.t<t_start], W_NR.t[W_NR.t>=t_start]))
    dataTemp=np.empty((len(tTemp), len(W_NR.LM)), dtype=complex)
    Smooth=np.resize(scri.utilities.transition_function(matchingt, t_start, t_end0), (len(W_NR.LM), len(matchingt))).T
    # Hybridize data
    matching_data=(1-Smooth)*PNData_spline(matchingt)+Smooth*W_NR.data[(W_NR.t>=t_start)&(W_NR.t<=t_end0),:]
    dataTemp[tTemp<t_start,:]=W_PN.data[W_PN.t<t_start,:]
    dataTemp[(tTemp>=t_start)&(tTemp<=t_end0),:]=matching_data
    dataTemp[tTemp>t_end0,:]=W_NR.data[W_NR.t>t_end0,:]
    # Delete indices that cause tiny time step
    minstep=min(min(np.diff(W_NR.t[(W_NR.t>t_pre)&(W_NR.t<t_end)])),min(np.diff(W_PN.t[(W_PN.t>t_pre)&(W_PN.t<t_end)])))
    BadIndices=np.nonzero(np.append(np.diff(tTemp)<minstep,0)&(tTemp>t_pre)&(tTemp<t_end))
    while len(BadIndices[0])>0:
        tTemp=np.delete(tTemp,BadIndices)
        dataTemp=np.delete(dataTemp,BadIndices,axis=0)
        BadIndices=np.nonzero(np.append(np.diff(tTemp)<minstep,0)&(tTemp>t_pre)&(tTemp<t_end))
    # Construct Hybrid waveform
    W_H=scri.WaveformModes()
    W_H.t=tTemp
    W_H.data=dataTemp
    ell_min, ell_max = min(W_NR.LM[:, 0]), max(W_NR.LM[:, 0])
    W_H.ells = ell_min, ell_max

    # Plot results   
    t1=-80000
    t2=0
    fig, (ax1, ax2, ax3) = plt.subplots(3,1)
    ax1.plot(W_NR.t, W_NR.data[:,4].real-W_NR.data[:,4].imag, label='NR', linewidth=1)
    ax1.plot(W_PN.t, W_PN.data[:,4].real-W_PN.data[:,4].imag, label='PN', linewidth=1)
    #ax1.plot(W_NR.t, np.abs(W_NR.data[:,4]), label='NR', linewidth=1)################3#########
    #ax1.plot(W_PN.t, np.abs(W_PN.data[:,4]), label='PN', linewidth=1)####################
    ax1.plot(W_H.t, W_H.data[:,4].real-W_H.data[:,4].imag, ls='--', label='Hybrid', linewidth=0.5)
    ax1.set_xlim((t1,t2))
    ax1.set_ylim((-0.15,0.15))
    ax1.set_ylabel("(2,2) mode")
    ax1.legend(['NR', 'PN'], loc="upper right")
    ax1.axvline(t_start, linestyle='dotted')
    ax1.axvline(t_end0, linestyle='dotted')
    ax2.plot(W_NR.t, W_NR.data[:,3].real-W_NR.data[:,3].imag, label='NR', linewidth=1)
    ax2.plot(W_PN.t, W_PN.data[:,3].real-W_PN.data[:,3].imag, label='PN', linewidth=1)
    ax2.plot(W_H.t, W_H.data[:,3].real-W_H.data[:,3].imag, ls='--', label='Hybrid', linewidth=0.5)
    ax2.set_xlim((t1,t2))
    ax2.set_ylim((-0.03,0.03))
    ax2.set_ylabel("(2,1) mode")
    ax2.axvline(t_start, linestyle='dotted')
    ax2.axvline(t_end0, linestyle='dotted')
    ax3.plot(W_NR.t, W_NR.data[:,2].real-W_NR.data[:,2].imag, label='NR', linewidth=1)
    ax3.plot(W_PN.t, W_PN.data[:,2].real-W_PN.data[:,2].imag, label='PN', linewidth=1)
    ax3.plot(W_H.t, W_H.data[:,2].real-W_H.data[:,2].imag, ls='--', label='Hybrid', linewidth=0.5)
    ax3.set_xlim((t1,t2))
    ax3.set_ylim((-0.01,0.02))
    ax3.set_ylabel("(2,0) mode")
    ax3.set_xlabel("Time")
    ax3.axvline(t_start, linestyle='dotted')
    ax3.axvline(t_end0, linestyle='dotted')
    fig.savefig(out_dir+"/hybridCheckResults1"+data_dir[-8:-5]+str(t_start)[:3],dpi=1000)
    fig.clf()
    global W_NR1    
    W_NR1=W_NR.copy()
    W_NR=W_NR.interpolate(W_PN.t)
    fig, (ax1, ax2, ax3) = plt.subplots(3,1)
    ax1.plot(W_PN.t, np.abs(W_NR.data[:,4]-W_PN.data[:,4]), linewidth=1)
    ax1.set_yscale('log')
    ax1.set_xlim((t1,t2))
    ax1.set_ylabel("(2,2) mode Error")
    ax1.axvline(t_start, linestyle='dotted')
    ax1.axvline(t_end0, linestyle='dotted')
    ax2.plot(W_PN.t, np.abs(W_NR.data[:,3]-W_PN.data[:,3]), linewidth=1)
    ax2.set_yscale('log')
    ax2.set_xlim((t1,t2))
    ax2.set_ylabel("(2,1) mode Error")
    ax2.axvline(t_start, linestyle='dotted')
    ax2.axvline(t_end0, linestyle='dotted')
    ax3.plot(W_PN.t, np.abs(W_NR.data[:,2]-W_PN.data[:,2]), linewidth=1)
    ax3.set_yscale('log')
    ax3.set_xlim((t1,t2))
    ax3.set_ylabel("(2,0) mode Error")
    ax3.set_xlabel("Time")
    ax3.axvline(t_start, linestyle='dotted')
    ax3.axvline(t_end0, linestyle='dotted')
    """
    plt.subplot(414)
    plt.plot(W_PN.t, np.linalg.norm((W_NR.data.real-W_NR.data.imag)\
        -(W_PN.data.real-W_PN.data.imag),axis=1), linewidth=1)
    plt.yscale('log')
    plt.xlim((t1,t2))#plt.xlim((t_start-25*np.pi/omega_0, t_end0+25*np.pi/omega_0))
    plt.ylabel("Norm Error")
    plt.xlabel("Time")
    plt.axvline(t_start, linestyle='dotted')
    plt.axvline(t_end0, linestyle='dotted')
    """
    fig.savefig(out_dir+"/hybridCheckResultsError1"+data_dir[-8:-5]+str(t_start)[:3],dpi=1000)
    fig.clf()
    np.savetxt(out_dir+"/hybridCheckResultsErrort"+data_dir[-8:-5]+str(t_start)[:3]+'.txt', W_PN.t, delimiter=',')
    np.savetxt(out_dir+"/hybridCheckResultsError22"+data_dir[-8:-5]+str(t_start)[:3]+'.txt', np.abs(W_NR.data[:,4]-W_PN.data[:,4]), delimiter=',')
    np.savetxt(out_dir+"/hybridCheckResultsError21"+data_dir[-8:-5]+str(t_start)[:3]+'.txt', np.abs(W_NR.data[:,3]-W_PN.data[:,3]), delimiter=',')
    np.savetxt(out_dir+"/hybridCheckResultsError20"+data_dir[-8:-5]+str(t_start)[:3]+'.txt', np.abs(W_NR.data[:,2]-W_PN.data[:,2]), delimiter=',')

    # Output results 
    """
    outname=out_dir+'/hybridPN.h5'
    scri.SpEC.write_to_h5(W_PN, outname, file_write_mode='w')
    outname=out_dir+'/hybridHybrid'+str(t_start)+'.h5'
    scri.SpEC.write_to_h5(W_H, outname, file_write_mode='w')
    """
    print("All done, total time:",time.time()-clock0)

# Run the code
import os
import argparse
parser = argparse.ArgumentParser()
parser.add_argument('--t',type=float, default=-7000.0,help='Start time of matching window')
parser.add_argument('--SimDir', default='/panfs/ds09/sxs/dzsun/SimAnnex/Public/HybTest/015/Lev3',help='Path in which to find the extropolated waveform data')
parser.add_argument('--CCEDir', default='/home/dzsun/CCEAnnex/Public/HybTest/015_CCE/Lev3/CCE',help='Path in which to find the CCE waveform data')
parser.add_argument('--OutDir', default='/home/dzsun',help='Path in which to output results')
parser.add_argument('--length',type=float, default=5000.0,help='Length of matching window')
parser.add_argument('--nOrbits',type=float, default=None,help='Length of matching window in orbits, will disable "length" option if not None')
parser.add_argument('--truncate',nargs=2,type=float, default=None,help='--truncate t1 t2. If specified, it will truncate the abd object and keep only data between t1 and t2')
args = vars(parser.parse_args())

global mismatch, W_NR1, PNIter, omega_00,PhyParas,length_global, Output
PNIter=0
length_global=args['length']
truncate=args['truncate']
OptArg=1
maxiter=3
while PNIter<=maxiter:
    print("PNIter=: ",PNIter)
    Hybridize(args['t'],args['SimDir'],args['CCEDir'],args['OutDir'], args['length'], args['nOrbits'], debug=0, OptimizePNParas=OptArg, truncate=truncate)
    PNIter=PNIter+1
    print("t_start = ",args['t']-length_global)
    print(Output)
"""
Some garbage

chi1temp=quaternion.as_float_array(R_delta*chiA*R_delta.conjugate())[:,1:]
    chi2temp=quaternion.as_float_array(R_delta*chiB*R_delta.conjugate())[:,1:]
    cost3=simpson(-(omega_NR_mag_matching-omega_PN_mag_spline(matchingt+x[0])),matchingt)/NormalizationOmega
    WaveformNorm_PN=PNNorm_spline(matchingt+x[0])
    cost5=simpson((WaveformNorm_NR_matching-WaveformNorm_PN)**2.0,matchingt)/Normalization
    cost4=(np.mean(omega_PN_mag_spline(matchingt+x[0])[-10:])-np.mean(omega_PN_mag_spline(matchingt+x[0])[0:10])-np.mean(omega_NR_mag_matching[-10:])+np.mean(omega_NR_mag_matching[0:10]))

def Plot_Cost_Funciton():
    xx = np.linspace(0.98,1.02, 200)
    yy = np.linspace(0.99,1.01, 200)
    X, Y = np.meshgrid(xx, yy)
    Z=np.empty((len(xx),len(yy)))
    for i in range(len(xx)):
        for j in range(len(yy)):
            Z[j,i] = Optimize11D([xx[i],yy[j],1.0,1.0,0.0,0.0,0.0,0.0,-1.18811841,  0.00257369, -0.01022434, -0.03249287])
    fig=plt.pcolor(X,Y,np.log(Z),cmap='rainbow')
    plt.xlabel("q/q_optimum")
    plt.ylabel("Mc/Mc_optimum")
    plt.colorbar(fig)
    plt.savefig("/home/dzsun/hybridO_009qMc",dpi=1000)
    plt.clf()
    
    cc=0
    x0=np.array([1.5,1,0.1,0.2,0.3,0.2,0.4,0.1,0,0,0,0])
    x1=np.array([1.4972630070632238,0.9997536456747469,0.12225649,0.21137109,0.27426194,0.20982456,0.37112871,0.14216262,-4.96687163e+00,  4.10469907e-03,  1.95525101e-03,  2.04323446e-03])
    dx=(x0-x1)/200
    print(dx)
    dx=dx*[0,1,0,0,0,0,0,0,0,0,0,0]
    Z=np.empty(300)
    for i in range(-50,250):
        x=x1+i*dx
        Z[i+50]=abs((Optimize11D(x+dx)-Optimize11D(x-dx))/(2*np.linalg.norm(dx)))
        print(Z[i+50],i)
    plt.plot(np.array(range(-50,250))/200,Z)
    plt.xlabel("lambda")
    plt.ylabel("derivative of cost1 on profile along M axis")
    plt.yscale('log')
    #plt.ylim(1e-11,1e-5)
    plt.title('Derivative profile along M axis')
    plt.savefig("/home/dzsun/hybridO_PNDevM",dpi=1000)
    plt.clf()
    
    t1=t_start-10000#-25000
    t2=t_start+10000#-5000
    t3=t2#-80000#t_start+10000
    dt=(cost3-0.5*cost4/(t_end-t_start)*(0.5*t_start-t1+0.5*t_end)/omega_0)*(0.5*t_start-t1+0.5*t_end)
    matchingt=W_NR.t[(W_NR.t>=t1)&(W_NR.t<=t2)]
    h1h2=np.linalg.norm(np.real(simpson(PNData_spline(matchingt)*np.conjugate(PNData_spline(matchingt+dt)), matchingt, axis=0)))
    h1h1=np.linalg.norm(np.real(simpson(PNData_spline(matchingt)*np.conjugate(PNData_spline(matchingt)), matchingt, axis=0)))
    h2h2=np.linalg.norm(np.real(simpson(PNData_spline(matchingt+dt)*np.conjugate(PNData_spline(matchingt+dt)), matchingt, axis=0)))
    print(h1h2,h1h1,h2h2,1-h1h2/np.sqrt(h1h1*h2h2))
"""
