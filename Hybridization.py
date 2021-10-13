import scri
import numpy as np
import quaternion
import h5py
import time
from scipy.optimize import least_squares
from scipy.interpolate import InterpolatedUnivariateSpline as Spline

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

def Hybridize(t_start, data_dir, out_dir):
    """
    Align and hybridize given NR waveform with PN waveform, the matching region starts at t_start, and last 3 orbits
    """

    clock0=time.time()
# Get NR waveform
    NRFileName=data_dir+'/rhOverM_Asymptotic_GeometricUnits.h5/Extrapolated_N2.dir'
    W_NR=scri.SpEC.read_from_h5(NRFileName)
    #W_NR.to_coprecessing_frame()####################################################################################3
    W_NR.t=W_NR.t-W_NR.max_norm_time()
    W_NR.data=-W_NR.data
    W_NR_corot=scri.to_corotating_frame(W_NR.copy())

# Get PN waveform
    PNFileName=data_dir+'/rhOverM_Inertial_PN.h5'
    W_PN=scri.SpEC.read_from_h5(PNFileName)
    #W_PN.to_coprecessing_frame()####################################################################################
    W_PN_corot=scri.to_corotating_frame(W_PN.copy())

# Get the initial angular velocity in matching region
    omega_NR=W_NR.copy().angular_velocity()
    omega_NR_mag = np.linalg.norm(omega_NR, axis=1)
    omega_0=omega_NR_mag[(W_NR.t>=t_start-10)&(W_NR.t<=t_start+10)]
    omega_0=np.mean(omega_0)

# Set up the matching region data for PN, and get the corresponding angular velocity and frame
    t_pre=t_start-10*np.pi/omega_0
    t_end=t_start+10*np.pi/omega_0
    t_end0=t_start+6*np.pi/omega_0
    omega_PN=W_PN.copy().angular_velocity()
    omega_PN_spline=SplineArray(W_PN.t, omega_PN)
    omega_PN_mag=np.linalg.norm(omega_PN, axis=1)
    omega_PN_mag_spline=Spline(W_PN.t, omega_PN_mag)
    PNData_spline=SplineArray(W_PN.t, W_PN.data)

# Get initial guess of time alignment by matching angular velocity
    matchingt=W_NR.t[(W_NR.t>=t_start)&(W_NR.t<=t_end0)]
    temp=np.append(np.diff(W_NR.t),0.0)
    dt=temp[(W_NR.t>=t_start)&(W_NR.t<=t_end0)]
    omega_NR_mag_matching=omega_NR_mag[(W_NR.t>=t_start)&(W_NR.t<=t_end0)]
    def InitialT(x):
        print(x)
        return np.sum((omega_NR_mag_matching-omega_PN_mag_spline(matchingt+x))**2*dt)/(omega_0**2.0)
    mint=least_squares(InitialT, 0.0, bounds=[-10*np.pi/omega_0,10*np.pi/omega_0])
    print(mint)
    print("Initial guess of t:", mint.x)

# Get initial guess of frame alignment
    R_delta = quaternion.optimal_alignment_in_Euclidean_metric(omega_NR[(W_NR.t>=t_start)&(W_NR.t<=t_end0)], omega_PN_spline(matchingt+mint.x), matchingt)
    W_NR_matching_in=W_NR.copy().interpolate(matchingt)
    omega_NR_matching = omega_NR[(W_NR.t>=t_start) & (W_NR.t<=t_end0)]
    omega_NR_hat = omega_NR_matching / np.linalg.norm(omega_NR_matching, axis=1)[:, np.newaxis]
    def InitialR(theta):
        print(theta)
        R_temp=R_delta*np.exp(theta/2*quaternion.quaternion(0.0, omega_NR_hat[0,0], omega_NR_hat[0,1],omega_NR_hat[0,2]))
        W_temp=scri.rotate_decomposition_basis(W_NR_matching_in.copy(), R_temp)
        temp=(np.angle(W_temp.data[int(len(matchingt)/2),4])-np.angle(PNData_spline(mint.x+W_temp.t[int(len(matchingt)/2)])[0,4]))**2
        print(temp)
        return temp
    minf=least_squares(InitialR, 0.0, bounds=[-np.pi,np.pi])
    print(minf)
    # Get rid of pi degeneracy
    R_delta=R_delta*np.exp(minf.x/2*quaternion.quaternion(0.0, omega_NR_hat[0,0], omega_NR_hat[0,1],omega_NR_hat[0,2]))
    R_delta2=R_delta*np.exp(np.pi/2*quaternion.quaternion(0.0, omega_NR_hat[0,0], omega_NR_hat[0,1],omega_NR_hat[0,2]))
    W_temp1=scri.rotate_decomposition_basis(W_NR_matching_in.copy(), R_delta)
    W_temp2=scri.rotate_decomposition_basis(W_NR_matching_in.copy(), R_delta2)
    temp1=np.linalg.norm(PNData_spline(matchingt+mint.x)-W_temp1.data[(W_temp1.t>=t_start)&(W_temp1.t<=t_end0)],axis=1)**2.0*dt
    temp2=np.linalg.norm(PNData_spline(matchingt+mint.x)-W_temp2.data[(W_temp2.t>=t_start)&(W_temp2.t<=t_end0)],axis=1)**2.0*dt
    temp1=np.real(sum(temp1))
    temp2=np.real(sum(temp2))
    if temp2<temp1:
        R_delta=R_delta2
    print("Initial guess of R_delta:",R_delta)
    logR_delta=quaternion.as_float_array(np.log(R_delta))
    print(logR_delta)

# Alignment of time and frame
    clock1=time.time()
    def Optimize4D(x):
        print(x)
        R_delta=np.exp(quaternion.quaternion(0.0,x[0],x[1],x[2]))
        W_temp=scri.rotate_decomposition_basis(W_NR_matching_in.copy(), R_delta)
        def Optimize1D(x):
            print(x)
            temp=np.linalg.norm(PNData_spline(matchingt+x)-W_temp.data[(W_temp.t>=t_start)&(W_temp.t<=t_end0)],axis=1)**2.0*dt
            temp1=np.linalg.norm(PNData_spline(matchingt+x),axis=1)*np.linalg.norm(W_temp.data[(W_temp.t>=t_start)&(W_temp.t<=t_end0)],axis=1) # Normalization factor
            temp1=sum(temp1)
            temp=np.real(sum(temp))/temp1
            print(temp)
            return temp
        optimizeT=least_squares(Optimize1D, mint.x, bounds=(mint.x-np.pi/omega_0,mint.x+np.pi/omega_0))
        return optimizeT.fun, optimizeT.x
    def _Optimize4D(x):
        return Optimize4D(x)[0]
    minima=least_squares(_Optimize4D, [logR_delta[0][1],logR_delta[0][2],logR_delta[0][3]], bounds=([-np.pi,-np.pi,-np.pi], [np.pi,np.pi,np.pi]))
    print(minima)
    t_delta=Optimize4D(minima.x)[1]
    print("Time shift=", t_delta)
    R_delta=np.exp(quaternion.quaternion(0.0,minima.x[0],minima.x[1],minima.x[2]))
    print("R_delta=",R_delta)
    print("Optimization time used:",time.time()-clock1)
    W_PN.t=W_PN.t-t_delta
    W_NR=scri.rotate_decomposition_basis(W_NR, R_delta)

# Hybridize waveform
    PNData_spline=SplineArray(W_PN.t, W_PN.data)
    W_H=scri.WaveformModes()
    W_H.t=np.append(W_PN.t[W_PN.t<t_start], W_NR.t[W_NR.t>=t_start])
    W_H.data=1j*np.empty((len(W_H.t), len(W_NR.LM)))
    ell_min, ell_max = min(W_NR.LM[:, 0]), max(W_NR.LM[:, 0])
    W_H.ells = ell_min, ell_max
    N=len(matchingt)
    xx=np.arange(N)/N
    # Hybridize data
    for i_m in range(W_NR.LM.shape[0]):
        matching_data=(1-scri.utilities.transition_function(xx,0.0,1.0,0.0,1.0))*PNData_spline(matchingt)[:,i_m]+scri.utilities.transition_function(xx,0.0,1.0,0.0,1.0)*W_NR.data[(W_NR.t>=t_start)&(W_NR.t<=t_end0),i_m]
        W_H.data[:,i_m]=np.append(np.append(W_PN.data[W_PN.t<t_start,i_m],matching_data),W_NR.data[W_NR.t>t_end0,i_m])
    print("finished")

# Output results
    outname=out_dir+'/hybridHybrid'+str(t_start)+'.h5'
    scri.SpEC.write_to_h5(W_H, outname, file_write_mode='w')
    outname=out_dir+'/hybridNR'+str(t_start)+'.h5'
    scri.SpEC.write_to_h5(W_NR, outname, file_write_mode='w')
    outname=out_dir+'/hybridPN'+str(t_start)+'.h5'
    scri.SpEC.write_to_h5(W_PN, outname, file_write_mode='w')
    print("Total time:",time.time()-clock0)

def Run():
    import os
    for i in [-30000]:
        Hybridize(i,'/home/dzsun/SimAnnex/Public/HybTest/006/Lev3','/home/dzsun')
        os.rename('/home/dzsun/rhOverM_hybridNR'+str(i)+'.h5','/home/dzsun/hybridNR'+str(i)+'.h5')
        os.rename('/home/dzsun/rhOverM_hybridPN'+str(i)+'.h5','/home/dzsun/hybridPN'+str(i)+'.h5')
        os.rename('/home/dzsun/UnknownDataType_hybridHybrid'+str(i)+'.h5','/home/dzsun/hybridHybrid'+str(i)+'.h5')
