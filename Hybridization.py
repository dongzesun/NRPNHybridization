import scri
from PYPostNewtonian.Code import PostNewtonian
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
import pyswarms as ps
from pyswarms.utils.plotters import (plot_cost_history, plot_contour, plot_surface)

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
    cost=(np.angle(W_temp.data[int(len(matchingt)/2),4])
        -np.angle(PNData_spline(mint.x+W_temp.t[int(len(matchingt)/2)])[0,4]))**2
    return cost

def Optimize4D(x):
    """
    Align 3D rotation and time shift.
    """
    phase=quaternion.quaternion(0.0, omega_mean[0]*x[0]/2, omega_mean[1]*x[0]/2, omega_mean[2]*x[0]/2)
    R_delta=np.exp(quaternion.quaternion(0.0,x[1],x[2],x[3])+phase)
    W_temp=scri.rotate_physical_system(W_NR_matching_in.copy(), R_delta)
    temp=PNData_spline(matchingt+x[0])
    cost=np.linalg.norm(temp-W_temp.data,axis=1)**2.0
    cost=np.real(simpson(cost, matchingt))/Normalization
    return cost

def OptimizePNParas(x):
    """
    Optimize PN parameters.
    """
    global mint, minima, R_delta, W_PN, W_PN_corot, t_end0, t_pre, t_end, omega_0, omega_PN,\
        omega_PN_spline, omega_PN_mag, omega_PN_mag_spline, PNData_spline, matchingt,\
        omega_NR_mag_matching, omega_mean, iter_num, W_NR_matching_in, omega_NR_hat,\
        Normalization, lowbound, upbound, scale, spin1, spin2, t_start, chiA, chiB, chi1_PN_spline, chi2_PN_spline
    clock1=time.time()
    t_start=t_starts
    iter_num+=1
    q=q_0*x[0]
    rotation=np.exp(xHat*x[3]/2+yHat*x[4]/2+zHat*x[5]/2)
    axis=quaternion.from_float_array(np.append(0,np.cross(chi1_i,chi2_i)))
    angle=np.exp(axis*x[6]/2)
    chi1_0=x[1]*(rotation*quaternion.quaternion(0,chi1_i[0],chi1_i[1],chi1_i[2])*rotation.conjugate()).vec
    chi2_0=x[2]*rotation*quaternion.quaternion(0,chi2_i[0],chi2_i[1],chi2_i[2])*rotation.conjugate()
    chi2_0=(angle*chi2_0*angle.conjugate()).vec
    chi1Mag=quaternion.quaternion(0,chi1_0[0],chi1_0[1],chi1_0[2]).abs()
    chi2Mag=quaternion.quaternion(0,chi2_0[0],chi2_0[1],chi2_0[2]).abs()
    S_chi1_0=quaternion.quaternion(0.0,0.0,0.0,0.0)
    S_chi2_0=S_chi1_0
    if chi1Mag>1e-12:
        S_chi1_0=np.sqrt(chi1Mag)*np.sqrt(
            -quaternion.quaternion(0,chi1_0[0],chi1_0[1],chi1_0[2]).normalized()*zHat).normalized()
    if chi2Mag>1e-12:
        S_chi2_0=np.sqrt(chi2Mag)*np.sqrt(
            -quaternion.quaternion(0,chi2_0[0],chi2_0[1],chi2_0[2]).normalized()*zHat).normalized()
    print(("Call # {4}, generating PN with parameters q={0}, omega_0={1}, chi1_0={2}, chi2_0={3},"
        +"t_PNstart={5}, t_PNend={6}.").format(q, omega_0,chi1_0, chi2_0,iter_num, t_PNStart, t_PNEnd))
    W_PN_corot, spin1, spin2=PostNewtonian.PNWaveform(q, omega_0, chi1_0, chi2_0, frame_0, t_start, t_PNStart, t_PNEnd)
    W_PN_corot.data[:,2]=0.0*W_PN_corot.data[:,2] # Not cosider memory effect since NR dosen't have corrrect memory.
    W_PN=scri.to_inertial_frame(W_PN_corot.copy())
    PNData_spline=SplineArray(W_PN.t, W_PN.data)
    chi1_PN_spline=SplineArray(W_PN.t, quaternion.as_float_array(spin1)[:,1:])
    chi2_PN_spline=SplineArray(W_PN.t, quaternion.as_float_array(spin2)[:,1:])
    omega_PN=W_PN.angular_velocity()
    omega_PN_spline=SplineArray(W_PN.t, omega_PN)
    omega_PN_mag=np.linalg.norm(omega_PN, axis=1)
    omega_PN_mag_spline=Spline(W_PN.t, omega_PN_mag)
    PNData_spline=SplineArray(W_PN.t, W_PN.data)

    # Get initial guess of time alignment by matching angular velocity
    mint=least_squares(InitialT, 0.0, bounds=[-10*np.pi/omega_0,10*np.pi/omega_0])

    # Get initial guess of frame alignment
    R_delta = quaternion.optimal_alignment_in_Euclidean_metric(omega_NR[(W_NR.t>=t_start)&(W_NR.t<=t_end0)],
        omega_PN_spline(matchingt+mint.x), matchingt)
    minf=least_squares(InitialR, 0.0, bounds=[-np.pi,np.pi])
    # Pi degeneracy
    phase=quaternion.quaternion(0.0, omega_mean[0]*mint.x/2, omega_mean[1]*mint.x/2, omega_mean[2]*mint.x/2)
    R_delta1=R_delta*np.exp(minf.x/2*quaternion.quaternion(0.0, omega_NR_hat[0,0], omega_NR_hat[0,1],omega_NR_hat[0,2]))
    R_delta2=R_delta1*np.exp(np.pi/2*quaternion.quaternion(0.0, omega_NR_hat[0,0], omega_NR_hat[0,1],omega_NR_hat[0,2]))
    W_temp1=scri.rotate_physical_system(W_NR_matching_in.copy(), R_delta1)
    W_temp2=scri.rotate_physical_system(W_NR_matching_in.copy(), R_delta2)
    temp=PNData_spline(matchingt+mint.x)
    cost1=np.linalg.norm(temp-W_temp1.data,axis=1)**2.0
    cost1=np.real(simpson(cost1, matchingt))/Normalization
    cost2=np.linalg.norm(temp-W_temp2.data,axis=1)**2.0
    cost2=np.real(simpson(cost2, matchingt))/Normalization
    if cost1<=cost2:
        R_delta=R_delta1
    else:
        R_delta=R_delta2
    print('Initial guess: ', mint.x, R_delta)
    W_temp=scri.rotate_physical_system(W_NR_matching_in.copy(), R_delta)
    temp=PNData_spline(matchingt+mint.x[0])
    cost1=np.linalg.norm(temp-W_temp.data,axis=1)**2.0
    cost1=np.real(simpson(cost1, matchingt))/Normalization
    chi1temp=quaternion.as_float_array(R_delta*chiA*R_delta.conjugate())[:,1:]
    chi2temp=quaternion.as_float_array(R_delta*chiB*R_delta.conjugate())[:,1:]
    chi1PN=chi1_PN_spline(matchingt+mint.x)
    chi2PN=chi2_PN_spline(matchingt+mint.x)
    cost2=simpson(np.linalg.norm(chi1temp-chi1PN,axis=1)**2.0+np.linalg.norm(chi2temp-chi2PN,axis=1)**2.0,matchingt)
    cost2=cost2/NormalizationChi
    print("{0}th call of cost function used time:".format(iter_num),time.time()-clock1)
    return cost1#+cost2

def AlignSwarms(x):
    y=np.empty(len(x))
    for i in range(len(x)):
        y[i]=OptimizePNParas(x[i])
    return y
     
def Align(x):
    """
    Generate PN waveform and align it with NR waveform.
    """
    clock1=time.time()
    global mint, minima, R_delta, W_PN, W_PN_corot, t_end0, t_pre, t_end, omega_0, omega_PN,\
        omega_PN_spline, omega_PN_mag, omega_PN_mag_spline, PNData_spline, matchingt,\
        omega_NR_mag_matching, omega_mean, iter_num, W_NR_matching_in, omega_NR_hat,\
        Normalization, lowbound, upbound, scale, spin1, spin2
    t_start=t_starts
    iter_num+=1
    q=q_0*x[0]
    rotation=np.exp(xHat*x[3]/2+yHat*x[4]/2+zHat*x[5]/2)
    axis=quaternion.from_float_array(np.append(0,np.cross(chi1_i,chi2_i)))
    angle=np.exp(axis*x[6]/2)
    chi1_0=x[1]*(rotation*quaternion.quaternion(0,chi1_i[0],chi1_i[1],chi1_i[2])*rotation.conjugate()).vec
    chi2_0=x[2]*rotation*quaternion.quaternion(0,chi2_i[0],chi2_i[1],chi2_i[2])*rotation.conjugate()
    chi2_0=(angle*chi2_0*angle.conjugate()).vec
    chi1Mag=quaternion.quaternion(0,chi1_0[0],chi1_0[1],chi1_0[2]).abs()
    chi2Mag=quaternion.quaternion(0,chi2_0[0],chi2_0[1],chi2_0[2]).abs()
    S_chi1_0=quaternion.quaternion(0.0,0.0,0.0,0.0)
    S_chi2_0=S_chi1_0
    if chi1Mag>1e-12:
        S_chi1_0=np.sqrt(chi1Mag)*np.sqrt(
            -quaternion.quaternion(0,chi1_0[0],chi1_0[1],chi1_0[2]).normalized()*zHat).normalized()
    if chi2Mag>1e-12:
        S_chi2_0=np.sqrt(chi2Mag)*np.sqrt(
            -quaternion.quaternion(0,chi2_0[0],chi2_0[1],chi2_0[2]).normalized()*zHat).normalized()
    print(("Call # {4}, generating PN with parameters q={0}, omega_0={1}, chi1_0={2}, chi2_0={3},"
        +"t_PNstart={5}, t_PNend={6}.").format(q, omega_0,chi1_0, chi2_0,iter_num, t_PNStart, t_PNEnd))
    W_PN_corot, spin1, spin2=PostNewtonian.PNWaveform(q, omega_0, chi1_0, chi2_0, frame_0, t_start, t_PNStart, t_PNEnd)
    W_PN_corot.data[:,2]=0.0*W_PN_corot.data[:,2] # Not cosider memory effect since NR dosen't have corrrect memory.
    W_PN=scri.to_inertial_frame(W_PN_corot.copy())
    
    # Set up the matching region data for PN, and get the corresponding angular velocity and frame
    if W_PN.t[0]>=t_start-10*np.pi/omega_0 or W_PN.t[-1]<=t_start+10*np.pi/omega_0:
        message=("t_start {0} should be at least {1} larger than the start time of PN waveform:{2},"
            +" and smaller than the ending time of PN waveform:{3}.")
        raise ValueError(message.format(t_start, 10*np.pi/omega_0, W_PN.t[0], W_PN.t[-1]))
    omega_PN=W_PN.angular_velocity()
    omega_PN_spline=SplineArray(W_PN.t, omega_PN)
    omega_PN_mag=np.linalg.norm(omega_PN, axis=1)
    omega_PN_mag_spline=Spline(W_PN.t, omega_PN_mag)
    PNData_spline=SplineArray(W_PN.t, W_PN.data)

    # Get initial guess of time alignment by matching angular velocity
    mint=least_squares(InitialT, 0.0, bounds=[-10*np.pi/omega_0,10*np.pi/omega_0])

    # Get initial guess of frame alignment
    R_delta = quaternion.optimal_alignment_in_Euclidean_metric(omega_NR[(W_NR.t>=t_start)&(W_NR.t<=t_end0)],
        omega_PN_spline(matchingt+mint.x), matchingt)
    minf=least_squares(InitialR, 0.0, bounds=[-np.pi,np.pi])
    # Pi degeneracy
    phase=quaternion.quaternion(0.0, omega_mean[0]*mint.x/2, omega_mean[1]*mint.x/2, omega_mean[2]*mint.x/2)
    R_delta1=R_delta*np.exp(minf.x/2*quaternion.quaternion(0.0, omega_NR_hat[0,0], omega_NR_hat[0,1],omega_NR_hat[0,2]))
    R_delta2=R_delta1*np.exp(np.pi/2*quaternion.quaternion(0.0, omega_NR_hat[0,0], omega_NR_hat[0,1],omega_NR_hat[0,2]))
    W_temp1=scri.rotate_physical_system(W_NR_matching_in.copy(), R_delta1)
    W_temp2=scri.rotate_physical_system(W_NR_matching_in.copy(), R_delta2)
    temp=PNData_spline(matchingt+mint.x)
    cost1=np.linalg.norm(temp-W_temp1.data,axis=1)**2.0
    cost1=np.real(simpson(cost1, matchingt))/Normalization
    cost2=np.linalg.norm(temp-W_temp2.data,axis=1)**2.0
    cost2=np.real(simpson(cost2, matchingt))/Normalization
    if cost1<=cost2:
        R_delta=R_delta1
    else:
        R_delta=R_delta2
    print("initial guess", mint.x, R_delta)
    logR_delta=quaternion.as_float_array(np.log(R_delta))-np.append(0.0, omega_mean*mint.x/2)

    # Alignment of time and frame
    lowbound=np.append(mint.x-np.pi/omega_0/2.0, logR_delta[0][1:]-np.pi/4)
    upbound=np.append(mint.x+np.pi/omega_0/2.0, logR_delta[0][1:]+np.pi/4)
    scale=[np.pi/omega_0,np.pi/4.0,np.pi/4.0,np.pi/4.0]
    minima=least_squares(Optimize4D, np.append(mint.x, logR_delta[0][1:]),
        bounds=(lowbound,upbound),x_scale='jac')
    print(minima)
    print("{0}th call of cost function used time:".format(iter_num),time.time()-clock1)
    return minima.fun

def Hybridize(t_start, data_dir, out_dir, debug=0, OptimizePNParas=0):
    """
    Align and hybridize given NR waveform with PN waveform.
    """
    clock0=time.time()
    global mint, minima, W_NR, W_PN, W_PN_corot, t_starts, t_end0, t_pre, t_end, omega_0, omega_PN,\
        omega_NR, omega_PN_spline, omega_PN_mag, omega_PN_mag_spline, PNData_spline, matchingt,\
        omega_NR_mag_matching, omega_mean, iter_num, W_NR_matching_in, R_delta, omega_NR_hat,\
        Normalization, xHat, yHat, zHat, t_PNStart, t_PNEnd, lowbound, upbound, scale, frame_0,\
        spin1,spin2,chiA,chiB,NormalizationChi,chi1_i,chi2_i,q_0
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
        raise ValueError(message.format(omega_0, omega_NR_mag[W_NR.t-W_NR.t[0]>10*np.pi/omega_0][0],
            omega_NR_mag[W_NR.t[-1]-W_NR.t>10*np.pi/omega_0][-1]))
    t_end0=-4000#W_NR.t[(omega_NR_mag>2.1*omega_0)&(W_NR.t-W_NR.t[0]>10*np.pi/omega_0)\
        #&(W_NR.t[-1]-W_NR.t>10*np.pi/omega_0)][0]
    t_pre=W_NR.t[(omega_NR_mag<0.99*omega_0)&(W_NR.t-W_NR.t[0]>10*np.pi/omega_0)\
        &(W_NR.t[-1]-W_NR.t>10*np.pi/omega_0)][-1]
    t_end=-5000#W_NR.t[(omega_NR_mag>2.11*omega_0)&(W_NR.t-W_NR.t[0]>10*np.pi/omega_0)\
        #&(W_NR.t[-1]-W_NR.t>10*np.pi/omega_0)][0]
    t_PNStart=W_NR.t[(omega_NR_mag<0.9*omega_0)&(W_NR.t[-1]-W_NR.t>10*np.pi/omega_0)][-1]-t_start
    t_PNEnd=-1000-t_start#W_NR.t[(omega_NR_mag>1.4*omega_0)&(W_NR.t-W_NR.t[0]>10*np.pi/omega_0)][0]-t_start
    matchingt=W_NR.t[(W_NR.t>=t_start)&(W_NR.t<=t_end0)]
    omega_NR_mag_matching=omega_NR_mag[(W_NR.t>=t_start)&(W_NR.t<=t_end0)]
    W_NR_matching_in=W_NR.interpolate(matchingt)
    omega_NR_matching = omega_NR[(W_NR.t>=t_start) & (W_NR.t<=t_end0)]
    omega_mean = np.mean(omega_NR_matching, axis=0)
    omega_NR_hat = omega_NR_matching / np.linalg.norm(omega_NR_matching, axis=1)[:, np.newaxis]
    Normalization=0.1*simpson(np.linalg.norm(W_NR_matching_in.data,axis=1)**2.0, matchingt)

# Get PN Parameters
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
        q=1/q
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
    beta=math.atan2(np.dot((Ra*zHat*Ra.inverse()).vec,lambdaHat[i_1].vec),
        np.dot((Ra*yHat*Ra.inverse()).vec,lambdaHat[i_1].vec))
    frame_0=Ra*np.exp((beta)/2*xHat)
    chi1_i=chiA[i_1]
    chi2_i=chiB[i_1]
    chiA_spline=SplineArray(tA, chiA)
    chiB_spline=SplineArray(tA, chiB)
    chiA=chiA_spline(matchingt)
    chiB=chiB_spline(matchingt)
    NormalizationChi=0.1*(simpson(np.linalg.norm(chiA,axis=1)**2.0, matchingt)+simpson(np.linalg.norm(chiB,axis=1)**2.0, matchingt))
    chiA=quaternion.from_float_array(np.column_stack((0.0*matchingt,chiA)))
    chiB=quaternion.from_float_array(np.column_stack((0.0*matchingt,chiB)))

# Align Waveforms
    t_starts=t_start
    PNParas=[1.0,1.0,1.0,0.0,0.0,0.0,0.0]
    iter_num=0
    if OptimizePNParas:
        xmax=[1.02,1.02,1.02,np.pi/10,np.pi/10,np.pi/10,np.pi/20]
        xmin=[0.98,0.98,0.98,-np.pi/10,-np.pi/10,-np.pi/10,-np.pi/20]
        bounds=(xmin,xmax)
        options = {'c1': 0.5, 'c2': 0.3, 'w':0.9}
        optimizer = ps.single.GlobalBestPSO(n_particles=10, dimensions=7, options=options,bounds=bounds)
        cost, pos = optimizer.optimize(AlignSwarms, iters=100, n_processes=10)
        print(cost,pos)
        plot_cost_history(cost_history=optimizer.cost_history)
        plt.savefig(out_dir+"/hybridOptimizeHistory{0}".format(t_start))
        plt.clf()
        PNPara=least_squares(AlignChi, pos, bounds=(xmin,xmax), x_scale='jac')
        print(PNPara)
        PNParas=PNPara.x
    t_PNStart, t_PNEnd=-90000, -1000-t_start#t_end-t_start###########################################
    Align(PNParas)

    if min((minima.x-lowbound)/scale)<1e-2 or min((upbound-minima.x)/scale)<1e-2:
        message=("Minima {0} near bounds {1}, {2}.")
        raise ValueError(message.format(minima.x, lowbound, upbound))
    if minima.success==0:
        raise ValueError("Optimized unsuccessful.")

    t_delta=minima.x[0]
    print("Time shift=", t_delta)
    phase=quaternion.quaternion(0.0, omega_mean[0]*minima.x[0]/2,\
        omega_mean[1]*minima.x[0]/2, omega_mean[2]*minima.x[0]/2)
    R_delta=np.exp(quaternion.quaternion(0.0,minima.x[1],minima.x[2],minima.x[3])+phase)
    print("R_delta=",R_delta)
    W_PN.t=W_PN.t-t_delta
    W_NR=scri.rotate_physical_system(W_NR, R_delta)
    if debug:
        spin1=R_delta.conjugate()*spin1*R_delta
        spin2=R_delta.conjugate()*spin2*R_delta
        plt.plot(matchingt,quaternion.as_float_array(chiA)[:,1])
        plt.plot(matchingt,quaternion.as_float_array(chiB)[:,1])
        plt.plot(W_PN.t,quaternion.as_float_array(spin1)[:,1],linestyle='dotted')
        plt.plot(W_PN.t,quaternion.as_float_array(spin2)[:,1],linestyle='dotted')
        plt.xlabel("Time")
        plt.ylabel("x component of chi")
        plt.xlim((-95000,-15000))
        plt.legend(['NR chi1', 'NR chi2', 'PN chi1', 'PN chi2'], loc="upper right")
        plt.savefig(out_dir+"/hybridCheckChi")
        plt.clf()

        plt.plot(matchingt, omega_NR_mag_matching, label='Angular velocity')
        plt.plot(matchingt, omega_PN_mag_spline(matchingt+mint.x), label='Angular velocity')
        plt.legend(['NR', 'PN'], loc="upper right")
        plt.ylabel("Omega Magnititude")
        plt.xlabel("Time")
        plt.savefig(out_dir+"/hybridCheckOmega")
        plt.clf()
        xx = np.linspace(minima.x[0]-5*np.pi/omega_0, minima.x[0]+5*np.pi/omega_0, 40)
        yy = np.linspace(minima.x[3]-2.0, minima.x[3]+2.0, 40)
        X, Y = np.meshgrid(xx, yy)
        Z=np.empty((len(xx),len(yy)))
        for i in range(len(xx)):
            for j in range(len(yy)):
                Z[j,i] = Optimize4D([xx[i], minima.x[1], minima.x[2], yy[j]])
        fig=plt.pcolor(X,Y,np.log(Z),cmap='rainbow')
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

    """#Calculate mismatch
    matchingt=W_NR.t[(W_NR.t>=-95000)&(W_NR.t<=-5000)]
    h1h2=np.linalg.norm(np.real(simpson(PNData_spline(matchingt)*np.conjugate(W_NR.data[(W_NR.t>=-95000)&(W_NR.t<=-5000),:]), matchingt, axis=0)))**2.0
    h1h1=np.linalg.norm(np.real(simpson(PNData_spline(matchingt)*np.conjugate(PNData_spline(matchingt)), matchingt, axis=0)))**2.0
    h2h2=np.linalg.norm(np.real(simpson(W_NR.data[(W_NR.t>=-95000)&(W_NR.t<=-5000),:]*np.conjugate(W_NR.data[(W_NR.t>=-95000)&(W_NR.t<=-5000),:]), matchingt, axis=0)))**2.0
    mismatch.append(1-h1h2/np.sqrt(h1h1*h2h2))
    print(h1h2,h1h1,h2h2,1-h1h2/np.sqrt(h1h1*h2h2))
    matchingt=W_NR.t[(W_NR.t>=-25000)&(W_NR.t<=-5000)]
    h1h2=np.linalg.norm(np.real(simpson(PNData_spline(matchingt)*np.conjugate(W_NR.data[(W_NR.t>=-25000)&(W_NR.t<=-5000),:]), matchingt, axis=0)))**2.0
    h1h1=np.linalg.norm(np.real(simpson(PNData_spline(matchingt)*np.conjugate(PNData_spline(matchingt)), matchingt, axis=0)))**2.0
    h2h2=np.linalg.norm(np.real(simpson(W_NR.data[(W_NR.t>=-25000)&(W_NR.t<=-5000),:]*np.conjugate(W_NR.data[(W_NR.t>=-25000)&(W_NR.t<=-5000),:]), matchingt, axis=0)))**2.0
    print(h1h2,h1h1,h2h2,1-h1h2/np.sqrt(h1h1*h2h2))
    """

# Plot results############################################################################    
    plt.subplot(311)
    plt.plot(W_NR.t, W_NR.data[:,4].real-W_NR.data[:,4].imag, label='NR', linewidth=1)
    plt.plot(W_PN.t, W_PN.data[:,4].real-W_PN.data[:,4].imag, label='PN', linewidth=1)
    plt.plot(W_H.t, W_H.data[:,4].real-W_H.data[:,4].imag, ls='--', label='Hybrid', linewidth=1)
    plt.xlim((t_start-25*np.pi/omega_0, t_end0+25*np.pi/omega_0))
    #plt.xlim((-105000,5000))
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
    #plt.xlim((-105000,5000))
    plt.ylim((-0.03,0.03))
    plt.ylabel("(2,1) mode")
    plt.axvline(t_start, linestyle='dotted')
    plt.axvline(t_end0, linestyle='dotted')
    plt.subplot(313)
    plt.plot(W_NR.t, W_NR.data[:,2].real-W_NR.data[:,2].imag, label='NR', linewidth=1)
    plt.plot(W_PN.t, W_PN.data[:,2].real-W_PN.data[:,2].imag, label='PN', linewidth=1)
    plt.plot(W_H.t, W_H.data[:,2].real-W_H.data[:,2].imag, ls='--', label='Hybrid', linewidth=1)
    plt.xlim((t_start-25*np.pi/omega_0, t_end0+25*np.pi/omega_0))
    #plt.xlim((-105000,5000))
    plt.ylim((-0.01,0.01))
    plt.ylabel("(2,0) mode")
    plt.xlabel("Time")
    plt.axvline(t_start, linestyle='dotted')
    plt.axvline(t_end0, linestyle='dotted')
    plt.savefig(out_dir+"/hybridCheckResults",dpi=1000)
    plt.clf()

    global W_NR1    
    W_NR1=W_NR.copy()#########################################################
    W_NR=W_NR.interpolate(W_PN.t)
    plt.subplot(311)
    plt.plot(W_PN.t, np.abs((W_NR.data[:,4].real-W_NR.data[:,4].imag)\
        -(W_PN.data[:,4].real-W_PN.data[:,4].imag)), linewidth=1)
    plt.yscale('log')
    plt.xlim((t_start-25*np.pi/omega_0, t_end0+25*np.pi/omega_0))
    plt.ylabel("(2,2) mode Error")
    plt.axvline(t_start, linestyle='dotted')
    plt.axvline(t_end0, linestyle='dotted')
    plt.subplot(312)
    plt.plot(W_PN.t, np.abs((W_NR.data[:,3].real-W_NR.data[:,3].imag)\
        -(W_PN.data[:,3].real-W_PN.data[:,3].imag)), linewidth=1)
    plt.yscale('log')
    plt.xlim((t_start-25*np.pi/omega_0, t_end0+25*np.pi/omega_0))
    plt.ylabel("(2,1) mode Error")
    plt.axvline(t_start, linestyle='dotted')
    plt.axvline(t_end0, linestyle='dotted')
    plt.subplot(313)
    plt.plot(W_NR.t, np.abs((W_NR.data[:,2].real-W_NR.data[:,2].imag)\
        -(W_PN.data[:,2].real-W_PN.data[:,2].imag)), linewidth=1)
    plt.yscale('log')
    plt.xlim((t_start-25*np.pi/omega_0, t_end0+25*np.pi/omega_0))
    plt.ylabel("(2,0) mode Error")
    plt.xlabel("Time")
    plt.axvline(t_start, linestyle='dotted')
    plt.axvline(t_end0, linestyle='dotted')
    plt.savefig(out_dir+"/hybridCheckResultsError",dpi=1000)
    plt.clf()

    plt.subplot(211)
    plt.plot(W_PN.t[W_PN.t>t0],np.linalg.norm(W_NR.data[W_PN.t>t0],axis=1))
    plt.plot(W_PN.t[W_PN.t>t0],np.linalg.norm(W_PN.data[W_PN.t>t0],axis=1))
    plt.ylabel("Waveform norm")
    plt.legend(['NR','PN'],loc="upper right")
    plt.subplot(212)
    plt.plot(W_PN.t[W_PN.t>t0],np.abs(np.linalg.norm(W_PN.data,axis=1)-np.linalg.norm(W_NR.data,axis=1))[W_PN.t>t0])
    plt.ylabel("Error of waveform norm")
    plt.xlabel("Time")
    plt.yscale('log')
    plt.axvline(t_start, linestyle='dotted')
    plt.axvline(t_end0, linestyle='dotted')
    plt.savefig(out_dir+"/hybridCheckWaveformNormError",dpi=1000)
    plt.clf()

# Output results 
    #outname=out_dir+'/hybridPN.h5'
    #scri.SpEC.write_to_h5(W_PN, outname, file_write_mode='w')
    outname=out_dir+'/hybridHybrid'+str(t_start)+'.h5'
    scri.SpEC.write_to_h5(W_H, outname, file_write_mode='w')
    print("All done, total time:",time.time()-clock0)

# Run the code
import os
global mismatch, W_NR1
mismatch=[]

for i in [-7000,-10000,-15000,-25000,-45000,-85000]:
    Hybridize(i,'data_dir','out_dir', debug=0, OptimizePNParas=1)
    print(mismatch)
