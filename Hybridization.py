import scri
import numpy as np
import quaternion
import h5py
import time
import matplotlib.pyplot as plt
from scipy.optimize import least_squares
from scipy.optimize import minimize
from scipy.integrate import simpson
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

def Hybridize(t_start, data_dir, out_dir, debug=0):
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
    omega_NR=W_NR.angular_velocity()
    omega_NR_mag = np.linalg.norm(omega_NR, axis=1)
    if W_NR.t[0]>=t_start-10 or W_NR.t[-1]<=t_start+10:
        message=("t_start {0} should be much larger than the start time of NR waveform {1} and much smaller than the ending time {2}.")
        raise ValueError(message.format(t_start, W_NR.t[0], W_NR.t[-1]))
    omega_0=omega_NR_mag[(W_NR.t>=t_start-10)&(W_NR.t<=t_start+10)]
    omega_0=np.mean(omega_0)

# Set up the matching region data for PN, and get the corresponding angular velocity and frame
    if W_NR.t[0]>=t_start-10*np.pi/omega_0 or W_NR.t[-1]<=t_start+10*np.pi/omega_0\
        or W_PN.t[0]>=t_start-10*np.pi/omega_0 or W_PN.t[-1]<=t_start+10*np.pi/omega_0:
        message=("t_start {0} should be at least {1} larger than the start time of NR waveform {2} and start of PN waveform {3},"
                  +" and smaller than the ending time of NR waveform {4} and PN waveform {5}.")
        raise ValueError(message.format(t_start, 10*np.pi/omega_0, W_NR.t[0], W_NR.t[-1], W_PN.t[0], W_PN.t[-1]))
    t_pre=t_start-10*np.pi/omega_0
    t_end=t_start+10*np.pi/omega_0
    t_end0=t_start+6*np.pi/omega_0
    omega_PN=W_PN.angular_velocity()
    omega_PN_spline=SplineArray(W_PN.t, omega_PN)
    omega_PN_mag=np.linalg.norm(omega_PN, axis=1)
    omega_PN_mag_spline=Spline(W_PN.t, omega_PN_mag)
    PNData_spline=SplineArray(W_PN.t, W_PN.data)
    PN22Data_spline=Spline(W_PN.t, W_PN.data[:,4])
    matchingt=W_NR.t[(W_NR.t>=t_start)&(W_NR.t<=t_end0)]
    omega_NR_mag_matching=omega_NR_mag[(W_NR.t>=t_start)&(W_NR.t<=t_end0)]
    if debug:
        plt.plot(matchingt, omega_NR_mag_matching, label='Angular velocity')
        plt.savefig(out_dir+"/hybridCheckOmega")
        plt.clf()

# Get initial guess of time alignment by matching angular velocity
    def InitialT(x):
        print(x)
        return simpson((omega_NR_mag_matching-omega_PN_mag_spline(matchingt+x))**2, x=matchingt)/(omega_0**2.0)
    mint=least_squares(InitialT, 0.0, bounds=[-10*np.pi/omega_0,10*np.pi/omega_0])
    print(mint)
    print("Initial guess of t:", mint.x)

# Get initial guess of frame alignment
    R_delta = quaternion.optimal_alignment_in_Euclidean_metric(omega_NR[(W_NR.t>=t_start)&(W_NR.t<=t_end0)],\
        omega_PN_spline(matchingt+mint.x), matchingt)
    W_NR_matching_in=W_NR.interpolate(matchingt)
    omega_NR_matching = omega_NR[(W_NR.t>=t_start) & (W_NR.t<=t_end0)]
    omega_NR_hat = omega_NR_matching / np.linalg.norm(omega_NR_matching, axis=1)[:, np.newaxis]
    def InitialR(theta):
        print(theta)
        R_temp=R_delta*np.exp(theta/2*quaternion.quaternion(0.0, omega_NR_hat[0,0], omega_NR_hat[0,1],omega_NR_hat[0,2]))
        W_temp=scri.rotate_physical_system(W_NR_matching_in.copy(), R_temp)
        cost=(np.angle(W_temp.data[int(len(matchingt)/2),4])-np.angle(PNData_spline(mint.x+W_temp.t[int(len(matchingt)/2)])[0,4]))**2
        print(cost)
        return cost
    minf=least_squares(InitialR, 0.0, bounds=[-np.pi,np.pi])
    print(minf)
    # Get rid of pi degeneracy
    R_delta=R_delta*np.exp(minf.x/2*quaternion.quaternion(0.0, omega_NR_hat[0,0], omega_NR_hat[0,1],omega_NR_hat[0,2]))
    R_delta2=R_delta*np.exp(np.pi/2*quaternion.quaternion(0.0, omega_NR_hat[0,0], omega_NR_hat[0,1],omega_NR_hat[0,2]))
    W_temp1=scri.rotate_physical_system(W_NR_matching_in.copy(), R_delta)
    W_temp2=scri.rotate_physical_system(W_NR_matching_in.copy(), R_delta2)
    cost1=np.linalg.norm(PNData_spline(matchingt+mint.x)-W_temp1.data,axis=1)**2.0
    cost2=np.linalg.norm(PNData_spline(matchingt+mint.x)-W_temp2.data,axis=1)**2.0
    cost1=np.real(simpson(cost1, matchingt))
    cost2=np.real(simpson(cost2, matchingt))
    if cost2<cost1:
        R_delta=R_delta2
    print("Initial guess of R_delta:",R_delta)
    logR_delta=quaternion.as_float_array(np.log(R_delta))

# Alignment of time and frame
    clock1=time.time()
    Normalization=simpson(np.linalg.norm(W_NR_matching_in.data,axis=1)**2.0, matchingt)
    def Optimize4D(x):
        print(x)
        R_delta=np.exp(quaternion.quaternion(0.0,x[1],x[2],x[3]))
        W_temp=scri.rotate_physical_system(W_NR_matching_in.copy(), R_delta)
        cost=np.linalg.norm(PNData_spline(matchingt+x[0])-W_temp.data,axis=1)**2.0
        cost=np.real(simpson(cost, matchingt))/Normalization
        print(cost)
        return cost
    lowbound=np.append(mint.x-np.pi/omega_0/2.0, logR_delta[0][1:]-np.pi/8)
    upbound=np.append(mint.x+np.pi/omega_0/2.0,logR_delta[0][1:]+np.pi/8)
    scale=[np.pi/omega_0,np.pi/4.0,np.pi/4.0,np.pi/4.0]
    minima=minimize(Optimize4D, np.append(mint.x, logR_delta[0][1:]), bounds=((lowbound[0], upbound[0]),\
        (lowbound[1], upbound[1]), (lowbound[2], upbound[2]), (lowbound[3], upbound[3])), method='BFGS', tol=1e-8)
    if min((minima.x-lowbound)/scale)<1e-2 or min((upbound-minima.x)/scale)<1e-2:
        message=("Minima {0} near bounds {1}, {2}.")
        raise ValueError(message.format(Minima.x, lowbound, upbound))
    print(minima)
    t_delta=minima.x[0]
    print("Time shift=", t_delta)
    R_delta=np.exp(quaternion.quaternion(0.0,minima.x[1],minima.x[2],minima.x[3]))
    print("R_delta=",R_delta)
    print("Optimization time used:",time.time()-clock1)
    W_PN.t=W_PN.t-t_delta
    W_NR=scri.rotate_physical_system(W_NR, R_delta)
    if debug:
        from mpl_toolkits.mplot3d import Axes3D
        xx = np.linspace(t_delta-10, t_delta+10, 40)
        yy = np.linspace(minima.x[3]-0.05, minima.x[3]+0.05, 40)
        X, Y = np.meshgrid(xx, yy)
        Z=np.empty((len(xx),len(yy)))
        for i in range(len(xx)):
            for j in range(len(yy)):
                Z[i,j] = Optimize4D([xx[i], minima.x[1], minima.x[2], yy[j]])
        fig=plt.figure()
        ax3 = Axes3D(fig)
        ax3.view_init(elev=90,azim=0)
        ax3.plot_surface(X,Y,np.exp(Z),cmap='rainbow')
        plt.xlabel("Time")
        plt.ylabel("One component of quarternion")
        plt.savefig(out_dir+"/hybridCheck4DOptimization")
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
    if debug:
        plt.plot(W_NR.t, W_NR.data[:,4].real-W_NR.data[:,4].imag, label='NR')
        plt.plot(W_PN.t, W_PN.data[:,4].real-W_PN.data[:,4].imag, label='PN')
        plt.plot(W_H.t, W_H.data[:,4].real-W_H.data[:,4].imag, ls='--', label='Hybrid')
        plt.xlim((t_start-5000, t_start+2000))
        plt.legend(['NR', 'PN', 'Hybrid'])
        plt.savefig(out_dir+"/hybridCheckResults")
        plt.clf()

# Output results
    outname=out_dir+'/hybridHybrid'+str(t_start)+'.h5'
    scri.SpEC.write_to_h5(W_H, outname, file_write_mode='w')
    outname=out_dir+'/hybridNR'+str(t_start)+'.h5'
    scri.SpEC.write_to_h5(W_NR, outname, file_write_mode='w')
    outname=out_dir+'/hybridPN'+str(t_start)+'.h5'
    scri.SpEC.write_to_h5(W_PN, outname, file_write_mode='w')
    print("Finished, total time:",time.time()-clock0)

def Run():
    import os
    for i in [-30000]:
        Hybridize(i,'/home/dzsun/SimAnnex/Public/HybTest/006/Lev3','/home/dzsun', debug=0)
        os.rename('/home/dzsun/rhOverM_hybridNR'+str(i)+'.h5','/home/dzsun/hybridNR'+str(i)+'.h5')
        os.rename('/home/dzsun/rhOverM_hybridPN'+str(i)+'.h5','/home/dzsun/hybridPN'+str(i)+'.h5')
        os.rename('/home/dzsun/UnknownDataType_hybridHybrid'+str(i)+'.h5','/home/dzsun/hybridHybrid'+str(i)+'.h5')
