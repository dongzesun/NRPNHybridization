import scri
import PNEvolution
import PNWaveformModes
import numpy as np
import math
import quaternion
import h5py
import time
import matplotlib.pyplot as plt
from scipy.optimize import least_squares
from scipy.integrate import simpson
from scipy.interpolate import InterpolatedUnivariateSpline as Spline
from scipy.interpolate import CubicSpline
from scipy.integrate import solve_ivp

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

def InitialT(x):
    return simpson((omega_NR_mag_matching-omega_PN_mag_spline(matchingt+x))**2, matchingt)/(omega_0**2.0)

def InitialR(theta):
    R_temp=R_delta*np.exp(theta/2*quaternion.quaternion(0.0, omega_NR_hat[0,0], omega_NR_hat[0,1],omega_NR_hat[0,2]))
    W_temp=scri.rotate_physical_system(W_NR_matching_in.copy(), R_temp)
    cost=(np.angle(W_temp.data[int(len(matchingt)/2),4])\
        -np.angle(PNData_spline(mint.x+W_temp.t[int(len(matchingt)/2)])[0,4]))**2
    return cost

def Optimize5D(x):
    phase=quaternion.quaternion(0.0, omega_mean[0]*x[0]/2, omega_mean[1]*x[0]/2, omega_mean[2]*x[0]/2)
    R_delta=np.exp(quaternion.quaternion(0.0,x[1],x[2],x[3])+phase)
    W_temp=scri.rotate_physical_system(W_NR_matching_in.copy(), R_delta)
    temp=PNData_spline(matchingt+x[0])
    temp[:,2]=temp[:,2]+x[4]
    cost=np.linalg.norm(temp-W_temp.data,axis=1)**2.0
    cost=np.real(simpson(cost, matchingt))/Normalization
    return cost

def Align(x):
    """
    Generate PN waveform and align it with NR waveform.
    x=[delta, chi1_x, chi1_y, chi1_z, chi2_x, chi2_y, chi2_z]
    """
    clock1=time.time()
    global mint, minima, R_delta, W_PN, W_PN_corot, t_end0, t_pre, t_end, omega_0, omega_PN,\
        omega_PN_spline, omega_PN_mag, omega_PN_mag_spline, PNData_spline, matchingt,\
        omega_NR_mag_matching, omega_mean, iter_num, W_NR_matching_in, omega_NR_hat, Normalization
    t_start=t_starts
    iter_num+=1
    [delta, chi1_x, chi1_y, chi1_z, chi2_x, chi2_y, chi2_z]=x
    m1=(1+delta)/2.0
    m2=(1-delta)/2.0
    v_i=omega_0**(1/3)
    chi1_i=[chi1_x, chi1_y, chi1_z]
    chi2_i=[chi2_x, chi2_y, chi2_z]
    chi1Mag=quaternion.quaternion(0,chi1_i[0],chi1_i[1],chi1_i[2]).abs()
    chi2Mag=quaternion.quaternion(0,chi2_i[0],chi2_i[1],chi2_i[2]).abs()
    S_chi1_i=quaternion.quaternion(0.0,0.0,0.0,0.0)
    S_chi2_i=S_chi1_i
    if chi1Mag>1e-12:
        S_chi1_i=np.sqrt(chi1Mag)*np.sqrt(\
            -quaternion.quaternion(0,chi1_i[0],chi1_i[1],chi1_i[2]).normalized()*zHat).normalized()
    if chi2Mag>1e-12:
        S_chi2_i=np.sqrt(chi2Mag)*np.sqrt(\
            -quaternion.quaternion(0,chi2_i[0],chi2_i[1],chi2_i[2]).normalized()*zHat).normalized()
    print(("Call # {8}, generating PN with parameters m1={0}, m2={1}, v_i={2}, S_chi1_i={3}, S_chi2_i={4},"
        +"rfrak_frame_i=[{5}, {6}, {7}].").format(m1, m2, v_i,S_chi1_i, S_chi2_i,\
        rfrak_frame_i[0], rfrak_frame_i[1], rfrak_frame_i[2], iter_num))
    PN=PNEvolution.TaylorTn_4p0PN_Q.TaylorTn_4p0PN_Q(xHat, yHat, zHat, m1, m2, v_i,S_chi1_i, S_chi2_i,\
        0.0, 0.0, 0.0, 0.0,rfrak_frame_i[0], rfrak_frame_i[1], rfrak_frame_i[2])
    W_PN_corot=scri.WaveformModes()
    W_PN_corot.t=PN.t-PN.t[-1]
    frame=np.empty(len(PN.t), dtype=quaternion.quaternion)
    for i in range (len(PN.t)):
        frame[i]=np.exp(quaternion.quaternion(0.0,PN.y[5,i],PN.y[6,i],PN.y[7,i]))
    W_PN_corot.frame=frame
    W_PN_corot.data=PNWaveformModes.WaveformModes_3p5PN.WaveformModes_3p5PN(xHat, yHat, zHat,\
        m1, m2, v_i,S_chi1_i, S_chi2_i,0.0, 0.0, 0.0, 0.0,\
        rfrak_frame_i[0], rfrak_frame_i[1], rfrak_frame_i[2],PN.y)
    ell_min, ell_max = 2, 8
    W_PN_corot.ells = ell_min, ell_max
    W_PN=scri.to_inertial_frame(W_PN_corot.copy())
    
# Set up the matching region data for PN, and get the corresponding angular velocity and frame
    if W_PN.t[0]>=t_start-10*np.pi/omega_0 or W_PN.t[-1]<=t_start+10*np.pi/omega_0:
        message=("t_start {0} should be at least {1} larger than the start time of PN waveform {2}"
            +" and smaller than the ending time of PN waveform {3}.")
        raise ValueError(message.format(t_start, 10*np.pi/omega_0, W_PN.t[0], W_PN.t[-1]))
    omega_PN=W_PN.angular_velocity()
    omega_PN_spline=SplineArray(W_PN.t, omega_PN)
    omega_PN_mag=np.linalg.norm(omega_PN, axis=1)
    omega_PN_mag_spline=Spline(W_PN.t, omega_PN_mag)
    PNData_spline=SplineArray(W_PN.t, W_PN.data)

# Get initial guess of time alignment by matching angular velocity
    mint=least_squares(InitialT, 0.0, bounds=[-10*np.pi/omega_0,10*np.pi/omega_0])

# Get initial guess of frame alignment
    R_delta = quaternion.optimal_alignment_in_Euclidean_metric(omega_NR[(W_NR.t>=t_start)&(W_NR.t<=t_end0)],\
        omega_PN_spline(matchingt+mint.x), matchingt)
    minf=least_squares(InitialR, 0.0, bounds=[-np.pi,np.pi])
    # Pi degeneracy
    R_delta1=R_delta*np.exp(minf.x/2*quaternion.quaternion(0.0, omega_NR_hat[0,0], omega_NR_hat[0,1],omega_NR_hat[0,2]))
    R_delta2=R_delta1*np.exp(np.pi/2*quaternion.quaternion(0.0, omega_NR_hat[0,0], omega_NR_hat[0,1],omega_NR_hat[0,2]))
    W_temp1=scri.rotate_physical_system(W_NR_matching_in.copy(), R_delta1)
    W_temp2=scri.rotate_physical_system(W_NR_matching_in.copy(), R_delta2)
    temp=PNData_spline(matchingt+mint.x)
    temp[:,2]=temp[:,2]+0.009
    cost1=np.linalg.norm(temp-W_temp1.data,axis=1)**2.0
    cost1=np.real(simpson(cost1, matchingt))/Normalization
    cost2=np.linalg.norm(temp-W_temp2.data,axis=1)**2.0
    cost2=np.real(simpson(cost2, matchingt))/Normalization
    if cost1<=cost2:
        R_delta=R_delta1
    else:
        R_delta=R_delta2
    logR_delta=quaternion.as_float_array(np.log(R_delta))-np.append(0.0, omega_mean*mint.x/2)

# Alignment of time and frame
    lowbound=np.append(np.append(mint.x-np.pi/omega_0/2.0, logR_delta[0][1:]-np.pi/4),0.0)
    upbound=np.append(np.append(mint.x+np.pi/omega_0/2.0, logR_delta[0][1:]+np.pi/4),0.015)
    scale=[np.pi/omega_0,np.pi/4.0,np.pi/4.0,np.pi/4.0,0.01]
    minima=least_squares(Optimize5D, np.append(np.append(mint.x, logR_delta[0][1:]),0.009),\
        bounds=(lowbound,upbound),x_scale='jac')
    if min((minima.x-lowbound)/scale)<1e-2 or min((upbound-minima.x)/scale)<1e-2:
        message=("Minima {0} near bounds {1}, {2}.")
        raise ValueError(message.format(minima.x, lowbound, upbound))
    if minima.success==0:
        raise ValueError("Optimized unsuccessful.")
    print(minima)
    print("{0}th call of cost function used time:".format(iter_num),time.time()-clock1)
    return minima.cost

def Hybridize(t_start, data_dir, out_dir, debug=0, OptimizePNParas=1):
    """
    Align and hybridize given NR waveform with PN waveform.
    """
    clock0=time.time()
    global mint, minima, W_NR, W_PN, W_PN_corot, t_starts, t_end0, t_pre, t_end, omega_0, omega_PN,\
        omega_NR, omega_PN_spline, omega_PN_mag, omega_PN_mag_spline, PNData_spline, matchingt,\
        omega_NR_mag_matching, omega_mean, iter_num, W_NR_matching_in, R_delta, omega_NR_hat,\
        Normalization, xHat, yHat, zHat, rfrak_frame_i
# Get NR waveform
    NRFileName=data_dir+'/rhOverM_Asymptotic_GeometricUnits_CoM.h5/Extrapolated_N2.dir'
    W_NR=scri.SpEC.read_from_h5(NRFileName)
    t0=-W_NR.max_norm_time()
    W_NR.t=W_NR.t+t0
    W_NR.data=-W_NR.data
    W_NR_corot=scri.to_corotating_frame(W_NR.copy())

# Get the initial angular velocity in matching region
    omega_NR=W_NR.angular_velocity()
    omega_NR_mag = np.linalg.norm(omega_NR, axis=1)
    if W_NR.t[0]>=t_start-10 or W_NR.t[-1]<=t_start+10:
        message=("t_start {0} should be much larger than the start time of NR"
            +" waveform {1} and much smaller than the ending time {2}.")
        raise ValueError(message.format(t_start, W_NR.t[0], W_NR.t[-1]))
    omega_0=omega_NR_mag[(W_NR.t>=t_start-10)&(W_NR.t<=t_start+10)]
    omega_0=np.mean(omega_0)
    if W_NR.t[0]>=t_start-10*np.pi/omega_0 or W_NR.t[-1]<=t_start+10*np.pi/omega_0:
        message=("t_start {0} should be at least {1} larger than the start time of NR waveform {2}"
            +" and smaller than the ending time of NR waveform {3}.")
        raise ValueError(message.format(t_start, 10*np.pi/omega_0, W_NR.t[0], W_NR.t[-1]))
    elif omega_NR_mag[W_NR.t[-1]-W_NR.t>10*np.pi/omega_0][-1]<=1.11*omega_0\
        or omega_NR_mag[W_NR.t-W_NR.t[0]>10*np.pi/omega_0][0]>=0.99*omega_0:
        message=("The range of angular velocity magnitude should at least contain (0.99,1.11)*omega_0,"
            +" where omega_0 is the angular velocity magnitude at t_start. Current omega_0 = {0},"
            +" while angular velocity magnitude range = ({1}, {2}).")
        raise ValueError(message.format(omega_0, omega_NR_mag[W_NR.t-W_NR.t[0]>10*np.pi/omega_0][0],\
            omega_NR_mag[W_NR.t[-1]-W_NR.t>10*np.pi/omega_0][-1]))
    t_end0=W_NR.t[(omega_NR_mag>1.1*omega_0)&(W_NR.t-W_NR.t[0]>10*np.pi/omega_0)\
        &(W_NR.t[-1]-W_NR.t>10*np.pi/omega_0)][0]
    t_pre=W_NR.t[(omega_NR_mag<0.99*omega_0)&(W_NR.t-W_NR.t[0]>10*np.pi/omega_0)\
        &(W_NR.t[-1]-W_NR.t>10*np.pi/omega_0)][-1]
    t_end=W_NR.t[(omega_NR_mag>1.11*omega_0)&(W_NR.t-W_NR.t[0]>10*np.pi/omega_0)\
        &(W_NR.t[-1]-W_NR.t>10*np.pi/omega_0)][0]
    matchingt=W_NR.t[(W_NR.t>=t_start)&(W_NR.t<=t_end0)]
    omega_NR_mag_matching=omega_NR_mag[(W_NR.t>=t_start)&(W_NR.t<=t_end0)]
    W_NR_matching_in=W_NR.interpolate(matchingt)
    omega_NR_matching = omega_NR[(W_NR.t>=t_start) & (W_NR.t<=t_end0)]
    omega_mean = np.mean(omega_NR_matching, axis=0)
    omega_NR_hat = omega_NR_matching / np.linalg.norm(omega_NR_matching, axis=1)[:, np.newaxis]
    Normalization=1e-4*simpson(np.linalg.norm(W_NR_matching_in.data,axis=1)**2.0, matchingt)

# Get PN Parameters
    xHat=quaternion.quaternion(0.0,1.0,0.0,0.0)
    yHat=quaternion.quaternion(0.0,0.0,1.0,0.0)
    zHat=quaternion.quaternion(0.0,0.0,0.0,1.0)
    with h5py.File(data_dir+'/Horizons.h5', 'r') as f:
        tA=f['AhA.dir/CoordCenterInertial.dat'][:,0]+t0
        xA=f['AhA.dir/CoordCenterInertial.dat'][:,1:]
        mA=f['AhA.dir/ArealMass.dat'][:,1]
        chiA=f['AhA.dir/chiInertial.dat'][:,1:]
        tB=f['AhB.dir/CoordCenterInertial.dat'][:,0]+t0
        xB=f['AhB.dir/CoordCenterInertial.dat'][:,1:]
        mB=f['AhB.dir/ArealMass.dat'][:,1]
        chiB=f['AhB.dir/chiInertial.dat'][:,1:]
    i_1=abs(tA-t_start).argmin()
    m1=mA[i_1]
    m2=mB[i_1]
    delta=(m1-m2)/(m1+m2)
    v_i=omega_0**(1/3)
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
    R_frame_i=Ra*np.exp((beta)/2*xHat)
    rfrak_frame_i=np.log(R_frame_i).vec
    
# Align Waveforms
    t_starts=t_start
    PNParas=np.append(np.append(delta, chiA[i_1]), chiB[i_1])
    iter_num=0
    if OptimizePNParas:
        PNPara=least_squares(Align, np.append(np.append(delta, chiA[i_1]), chiB[i_1]),\
            bounds=(np.append(np.append(delta-0.1, chiA[i_1]-0.2), chiB[i_1]-0.2),\
            np.append(np.append(delta+0.1, chiA[i_1]+0.2), chiB[i_1]+0.2)), x_scale='jac')
        print(PNPara)
        PNParas=PNPara.x
    Align(PNParas)
    t_delta=minima.x[0]
    print("Time shift=", t_delta)
    phase=quaternion.quaternion(0.0, omega_mean[0]*minima.x[0]/2,\
        omega_mean[1]*minima.x[0]/2, omega_mean[2]*minima.x[0]/2)
    R_delta=np.exp(quaternion.quaternion(0.0,minima.x[1],minima.x[2],minima.x[3])+phase)
    print("R_delta=",R_delta)
    W_PN.t=W_PN.t-t_delta
    W_NR=scri.rotate_physical_system(W_NR, R_delta)
    W_NR.data[:,2]=W_NR.data[:,2]-minima.x[4]
    if debug:
        plt.plot(matchingt, omega_NR_mag_matching, label='Angular velocity')
        plt.plot(matchingt, omega_PN_mag_spline(matchingt+mint.x), label='Angular velocity')
        plt.legend(['NR', 'PN'], loc="upper right")
        plt.ylabel("Omega Magnititude")
        plt.xlabel("Time")
        plt.savefig(out_dir+"/hybridCheckOmega")
        plt.clf()
        xx = np.linspace(minima.x[0]-30*np.pi/omega_0, minima.x[0]+30*np.pi/omega_0, 40)
        yy = np.linspace(minima.x[3]-2.5, minima.x[3]+2.5, 40)
        X, Y = np.meshgrid(xx, yy)
        Z=np.empty((len(xx),len(yy)))
        for i in range(len(xx)):
            for j in range(len(yy)):
                Z[j,i] = Optimize4D([xx[i], minima.x[1], minima.x[2], yy[j]])
        fig=plt.pcolor(X,Y,Z,cmap='rainbow')
        plt.xlabel("Time")
        plt.ylabel("One component of log of quarternion")
        plt.colorbar(fig)
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
    plt.subplot(311)
    plt.plot(W_NR.t, W_NR.data[:,4].real-W_NR.data[:,4].imag, label='NR', linewidth=1)
    plt.plot(W_PN.t, W_PN.data[:,4].real-W_PN.data[:,4].imag, label='PN', linewidth=1)
    plt.plot(W_H.t, W_H.data[:,4].real-W_H.data[:,4].imag, ls='--', label='Hybrid', linewidth=1)
    plt.xlim((t_start-25*np.pi/omega_0, t_end0+25*np.pi/omega_0))
    plt.ylim((-0.15,0.15))
    plt.ylabel("(2,2) mode")
    plt.legend(['NR', 'PN', 'Hybrid'], loc="upper right")
    plt.axvline(t_start, linestyle='dotted')
    plt.axvline(t_end0, linestyle='dotted')
    plt.subplot(312)
    plt.plot(W_NR.t, W_NR.data[:,3].real-W_NR.data[:,3].imag, label='NR', linewidth=1)
    plt.plot(W_PN.t, W_PN.data[:,3].real-W_PN.data[:,3].imag, label='PN', linewidth=1)
    plt.plot(W_H.t, W_H.data[:,3].real-W_H.data[:,3].imag, ls='--', label='Hybrid', linewidth=1)
    plt.xlim((t_start-25*np.pi/omega_0, t_end0+25*np.pi/omega_0))
    plt.ylim((-0.03,0.03))
    plt.ylabel("(2,1) mode")
    plt.axvline(t_start, linestyle='dotted')
    plt.axvline(t_end0, linestyle='dotted')
    plt.subplot(313)
    plt.plot(W_NR.t, W_NR.data[:,2].real-W_NR.data[:,2].imag, label='NR', linewidth=1)
    plt.plot(W_PN.t, W_PN.data[:,2].real-W_PN.data[:,2].imag, label='PN', linewidth=1)
    plt.plot(W_H.t, W_H.data[:,2].real-W_H.data[:,2].imag, ls='--', label='Hybrid', linewidth=1)
    plt.xlim((t_start-25*np.pi/omega_0, t_end0+25*np.pi/omega_0))
    plt.ylim((-0.015,0.01))
    plt.ylabel("(2,2) mode")
    plt.xlabel("Time")
    plt.axvline(t_start, linestyle='dotted')
    plt.axvline(t_end0, linestyle='dotted')
    plt.savefig(out_dir+"/hybridCheckResults",dpi=1000)
    plt.clf()

# Output results
    outname=out_dir+'/hybridHybrid'+str(t_start)+'.h5'
    scri.SpEC.write_to_h5(W_H, outname, file_write_mode='w')
    outname=out_dir+'/hybridNR'+str(t_start)+'.h5'
    scri.SpEC.write_to_h5(W_NR, outname, file_write_mode='w')
    outname=out_dir+'/hybridPN'+str(t_start)+'.h5'
    scri.SpEC.write_to_h5(W_PN, outname, file_write_mode='w')
    print("All done, total time:",time.time()-clock0)

# Run the code
import os
for i in [-30000]:
    Hybridize(i,'/home/dzsun/SimAnnex/Public/HybTest/006/Lev3','/home/dzsun', debug=0, OptimizePNParas=1)
    #Hybridize(i,'/home/dzsun/SimAnnex/Public/NonSpinningSurrogate/BBH_SKS_d17.5_q2_sA_0_0_0_sB_0_0_0/Lev4',\
    #    '/home/dzsun',debug=0)
    os.rename('/home/dzsun/rhOverM_hybridNR'+str(i)+'.h5','/home/dzsun/hybridNR'+str(i)+'.h5')
    os.rename('/home/dzsun/rhOverM_hybridPN'+str(i)+'.h5','/home/dzsun/hybridPN'+str(i)+'.h5')
    os.rename('/home/dzsun/UnknownDataType_hybridHybrid'+str(i)+'.h5','/home/dzsun/hybridHybrid'+str(i)+'.h5')
