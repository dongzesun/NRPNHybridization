import scri
import numpy as np
import quaternion
import h5py
import time
from scipy.signal import savgol_filter
from scipy.optimize import least_squares

def time_indices(T, t1, t2):
    """
    Given an array of times, return the indices between t1 and t2
    """

    Ind = range(len(T))
    Size = len(Ind)
    i = 1
    while i < Size:
        if T[Ind[i]] <= t1 and T[Ind[i+1]] > t1:
            Ind = np.delete(Ind, range(i))
            Size = len(Ind)
            i = 1
        if T[Ind[i]] >= t2 and T[Ind[i-1]] < t2:
            Ind = np.delete(Ind, range(i+1, Size))
            Size = len(Ind)
        i += 1
    return Ind

def MatchingRegion(W, t1, t2):
    """
    Interpolate waveform_base object W between t1 and t2, and calculate angular velocity
    """

    W_matching=W.copy().interpolate(np.arange(t1, t2, 1.0))
    # Calculate angular velocity
    omega=W_matching.copy().angular_velocity()
    omega=np.sqrt(omega[:,0]**2+omega[:,1]**2+omega[:,2]**2)
    omega=savgol_filter(omega,31,1)
    omega_prime=np.diff(omega)
    omega_prime=savgol_filter(omega_prime,31,1)
    return W_matching, omega, omega_prime

def Hybridize(t_start, data_dir, out_dir):
    """
    Align and hybridize given NR waveform with PN waveform, the matching region starts at t_start, and last 3 orbits
    """

    clock0=time.time()
# Get NR waveform
    NRFileName=data_dir+'/rhOverM_Asymptotic_GeometricUnits.h5/Extrapolated_N2.dir'
    W_NR=scri.SpEC.read_from_h5(NRFileName)
    W_NR.t=W_NR.t-W_NR.max_norm_time()
    W_NR.data=-W_NR.data
    W_NR_corot=scri.to_corotating_frame(W_NR.copy())
# Get PN waveform
    PNFileName=data_dir+'/rhOverM_Inertial_PN.h5'
    W_PN=scri.SpEC.read_from_h5(PNFileName)
    W_PN_corot=scri.to_corotating_frame(W_PN.copy())

# Get the initial angular velocity in matching region
    temp1, omega_NR, temp2=MatchingRegion(W_NR, t_start-1000, t_start+1000) # Here 1000 is an arbitrary large number
    omega_0=omega_NR[int(len(omega_NR)/2)]

# Set up the matching region data for NR, and get the corresponding angular velocity and frame
    t_pre=t_start-10*np.pi/omega_0
    t_end=t_start+10*np.pi/omega_0
    W_NR_matching_in, omega_NR, temp1=MatchingRegion(W_NR, t_pre, t_end)
    W_NR_matching_corot, temp1, temp2=MatchingRegion(W_NR_corot, t_pre, t_end)
    W_PN_matching_in, omega_PN, omega_PN_prime=MatchingRegion(W_PN, t_pre, t_end)
    W_PN_matching_corot, temp1, temp2=MatchingRegion(W_PN_corot, t_pre, t_end)
    print("After interpolate:",time.time()-clock0)

# Get initial guess of time alignment by matching angular velocity
    t_start_indice=time_indices(W_NR_matching_in.t, t_start, t_start+0.01)[0]
    def minix(x):
        print(x)
        Indices_PN=time_indices(W_PN_matching_in.t, t_start+x, t_start+x+6*np.pi/omega_0)
        dt=t_start+x-W_PN_matching_in.t[Indices_PN[0]]
        omega_PN_matching=omega_PN[Indices_PN]+omega_PN_prime[Indices_PN]*dt
        return np.sum((omega_NR[t_start_indice:t_start_indice+len(Indices_PN)]-omega_PN_matching)**2)
    mint=least_squares(minix, 0.0, bounds=[-10*np.pi/omega_0,10*np.pi/omega_0], gtol=2.23e-16)
    print(mint)
    t_delta=-mint.x
    print("Initial guess of t:",t_delta)
    clock1=time.time()

# Get initial guess of frame alignment
    Indices_PN=time_indices(W_PN_matching_in.t, t_start+mint.x, t_start+mint.x+6*np.pi/omega_0)
    R_delta=quaternion.optimal_alignment_in_chordal_metric(W_PN_matching_corot.frame[Indices_PN].conjugate(), W_NR_matching_corot.frame[t_start_indice:t_start_indice+len(Indices_PN)].conjugate()).conjugate()
    print("Initial guess of R_delta:",R_delta)
    R_delta=quaternion.as_float_array(R_delta)

# Alignment of time and frame
    def Optimize4D(x):
        print(x)
        temp=0j;
        R_delta=np.exp(quaternion.quaternion(0.0,x[1],x[2],x[3]))
        Indices_PN=time_indices(W_PN_matching_in.t, t_start+x[0], t_start+x[0]+6*np.pi/omega_0)
        W_temp=W_NR_matching_in.copy()
        W_temp=scri.rotate_decomposition_basis(W_temp, R_delta)
        temp=np.linalg.norm(W_PN_matching_in.data[Indices_PN]-W_temp.data[t_start_indice:t_start_indice+len(Indices_PN)],axis=1)**2.0
        temp=np.real(sum(temp))
        return temp, R_delta
    def _Optimize4D(x):
        return Optimize4D(x)[0]
    mini=least_squares(_Optimize4D, [-t_delta,R_delta[1],R_delta[2],R_delta[3]], bounds=([-t_delta-np.pi/omega_0,-1.0,-1.0,-1.0],[-t_delta+np.pi/omega_0,1.0,1.0,1.0]))
    print("Optimization time used:",time.time()-clock1)
    print(mini)
    print("Time shift=", -mini.x[0])
    R_delta=Optimize4D(mini.x)[1]
    print("R_delta=",R_delta)
    W_PN.t=W_PN.t-mini.x[0]
    W_NR=scri.rotate_decomposition_basis(W_NR, R_delta)

# Hybridize waveform
    t_end0=t_start+6*np.pi/omega_0
    W_matching_NR, temp1, temp2=MatchingRegion(W_NR, t_pre, t_end)
    W_matching_PN, temp1, temp2=MatchingRegion(W_PN, t_pre, t_end)
    t_start_indicePN=time_indices(W_PN.t, t_start, t_start+0.01)[0]
    t_start_indiceNR=time_indices(W_NR.t, t_end0, t_end0+0.01)[1]
    t_start_indiceMatching=time_indices(W_matching_NR.t, t_start, t_start+0.01)[0]
    t_end_indiceMatching=time_indices(W_matching_NR.t, t_end0, t_end0+0.01)[1]
    Ind=np.arange(t_start_indiceMatching,t_end_indiceMatching)
    W_H=scri.WaveformModes()
    W_H.t=np.append(np.append(W_PN.t[0:t_start_indicePN], W_matching_NR.t[Ind]), W_NR.t[t_start_indiceNR:-1])
    W_H.data=1j*np.empty((len(W_H.t), len(W_NR.LM)))
    ell_min, ell_max = min(W_NR.LM[:, 0]), max(W_NR.LM[:, 0])
    W_H.ells = ell_min, ell_max
    N=len(Ind)
    xx=np.arange(N)/N
    # Hybridize data
    for i_m in range(W_NR.LM.shape[0]):
        matching_data=(1-scri.utilities.transition_function(xx,0.0,1.0,0.0,1.0))*W_matching_PN.data[Ind,i_m]+scri.utilities.transition_function(xx,0.0,1.0,0.0,1.0)*W_matching_NR.data[Ind,i_m]
        W_H.data[:,i_m]=np.append(np.append(W_PN.data[0:t_start_indicePN,i_m],matching_data),W_NR.data[t_start_indiceNR:-1,i_m])
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
    for i in [-30000]:
           Hybridize(i,'/home/dzsun/SimAnnex/Public/HybTest/006/Lev3','/home/dzsun') 
