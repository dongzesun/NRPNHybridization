import scri
import numpy as np
import quaternion

def monotonic_indices(T, MinTimeStep=1.0e-3):
    """
    Given an array of times, return the indices that make the array strictly monotonic.
    """
    from numpy import delete

    Ind = range(len(T))
    Size = len(Ind)
    i = 1
    while i < Size:
        if T[Ind[i]] <= T[Ind[i - 1]] + MinTimeStep:
            j = 0
            while T[Ind[j]] + MinTimeStep < T[Ind[i]]:
                j += 1
            # erase data from j (inclusive) to i (exclusive)
            Ind = delete(Ind, range(j, i))
            Size = len(Ind)
            i = j - 1
        i += 1
    return Ind

def ReadFromNRAR(FileName) :
    """
    Read data from an H5 file in NRAR format.
    """
    import re
    import h5py
    #from Quaternions import Quaternion
    import numpy

    YlmRegex = re.compile(r"""Y_l(?P<L>[0-9]+)_m(?P<M>[-+0-9]+)\.dat""")
# Initialize the Waveform object
    W = scri.waveform_base.WaveformBase()
# Record the filename being read in
    W._append_history(str("*this = GWFrames.ReadFromNRAR(FileName='{0}')\n".format(FileName)))
    try :
        FileName, RootGroup = FileName.rsplit('.h5', 1)
        FileName += '.h5'
    except ValueError :
        RootGroup = '' # FileName is just a file, not a group in a file
    try :
        f_h5 = h5py.File(FileName, 'r')
    except IOError :
        print("ReadFromNRAR could not open the file '{0}'\n\n".format(FileName))
        raise
    if(RootGroup) :
        f = f_h5[RootGroup]
    else :
        f = f_h5
    try :
        try :
            for r in f['Frame'] :
                print(Quaternion(r))
                W.frame=Quaternion(r)
        except KeyError :
            pass # There was no frame data
# Get the descriptive items
        try :
            W.frameType=int(f.attrs['FrameType'])
            W.dataType=int(f.attrs['DataType'])
            W.r_is_scaled_out=bool(f.attrs['RIsScaledOut'])
            W.m_is_scaled_out=bool(f.attrs['MIsScaledOut'])
        except KeyError :
            pass
# Get the names of all the datasets in the h5 file, and check for matches
        YLMdata = [DataSet for DataSet in list(f) for m in [YlmRegex.search(DataSet)] if m]
        if(len(YLMdata)==0) :
            raise ValueError("Couldn't understand dataset names in '{0}'.".format(FileName))
# Sort the dataset names by increasing ell, then increasing m
        YLMdata = sorted(YLMdata, key=lambda DataSet : [int(YlmRegex.search(DataSet).group('L')), int(YlmRegex.search(DataSet).group('M'))])
# List just the ell and m numbers
        LM = sorted([[int(m.group('L')), int(m.group('M'))] for DataSet in YLMdata for m in [YlmRegex.search(DataSet)] if m])
# Get the time data (assuming all are equal)
        NModes = len(LM)
        Wdata = f[YLMdata[0]]
        NTimes = Wdata.shape[0]
        T = Wdata[:,0]
# Set up storage
        Re = numpy.empty((NTimes, NModes))
        Im = numpy.empty((NTimes, NModes))
        m = 0
# Loop through, getting each mode
        for DataSet in YLMdata :
            if( not (f[DataSet].shape[0]==NTimes) ) :
                raise ValueError("The number of time steps in this dataset should be {0}; ".format(NTimes) +
                                 "it is {0} in '{1}'.".format(f[DataSet].shape[0], DataSet))
            Re[:,m] = f[DataSet][:,1]
            Im[:,m] = f[DataSet][:,2]
            m += 1
# Make sure time is monotonic and set the data
        Indices = monotonic_indices(T)
        BadIndices = numpy.setdiff1d(range(len(T)), Indices)
        W.NTimes=NTimes
        W.t=T[Indices]
        W.LM=np.array(LM)
        W.data=numpy.delete(Re, BadIndices, 0)+1j*numpy.delete(Im, BadIndices, 0)
        #Ws = scri.WaveformModes(W)
    except KeyError :
        print("This H5 file appears to have not stored all the required information.\n\n")
        raise # Re-raise the exception after adding our information
    finally : # Use `finally` to make sure this happens:
        f_h5.close()
    return W

def OutputToNRAR(W, FileName, FileWriteMode='w') :
    """
    Output the Waveform in NRAR format.

    Note that the FileName is prepended with some descriptive
    information involving the data type and the frame type, such as
    'rhOverM_Corotating_' or 'rMPsi4_Aligned_'.

    """
    from h5py import File
    Group = None
    if('.h5' in FileName and not FileName.endswith('.h5')) :
        FileName,Group = FileName.split('.h5')
        FileName += '.h5'
# Open the file for output
    try :
        F = File(FileName, FileWriteMode)
    except IOError : # If that did not work...
        print("OutputToNRAR was unable to open the file '{0}'.\n\n".format(FileName))
        raise # re-raise the exception after the informative message above
    try :
# If we are writing to a group within the file, create it
        if(Group) :
            G = F.create_group(Group)
        else :
            G = F
# Now write data
        for i_m in range(W.LM.shape[0]) :
            ell,m = W.LM[i_m]
            Data_m = G.create_dataset("Y_l{0}_m{1}.dat".format(ell, m), data=[[t, d.real, d.imag] for t,d in zip(W.t,W.data[:,i_m])],
                                      compression="gzip", shuffle=True)
            Data_m.attrs['ell'] = ell
            Data_m.attrs['m'] = m
    finally : # Use `finally` to make sure this happens:
# Close the file and we are done
        F.close()

def time_indices(T, t1, t2):
    """
    Given an array of times, return the indices between t1 and t2
    """
    from numpy import delete

    Ind = range(len(T))
    Size = len(Ind)
    i = 1
    while i < Size:
        if T[Ind[i]] <= t1 and T[Ind[i+1]] > t1:
            Ind = delete(Ind, range(i))
            Size = len(Ind)
            i = 1
        if T[Ind[i]] >= t2 and T[Ind[i-1]] < t2:
            Ind = delete(Ind, range(i+1, Size))
            Size = len(Ind)
        i += 1
    return Ind

def smooth(x):
    import numpy as np

    for i in range(len(x)):
        if (x[i]==0):
            x[i]=0.001
        elif (x[i]==1):
            x[i]=0.999
    return 1.0/(1.0+np.exp(1.0/(x-1.0) + 1.0/x))

def MatchingRegion(W, t1, t2):
    """
    Interpolate waveform_base object W between t1 and t2, and calculate angular velocity
    """
    import h5py
    from scipy.interpolate import interp1d
    import numpy as np
    import quaternion
    import os

    Indices=time_indices(W.t, t1, t2)
    W_matching=scri.waveform_base.WaveformBase()
    W_matching.t=np.arange(t1, t2, 1.0)
    W_matching.LM=W.LM
    W_matching.data=np.empty((len(W_matching.t), len(W_matching.LM)), dtype=complex)
    # Interpolate waveform data
    for i_m in range(W.LM.shape[0]):
        tempfunc = interp1d(W.t[Indices],W.data[Indices,i_m])
        W_matching.data[:,i_m]=tempfunc(W_matching.t)
    # Interpolate frame
    if len(W.frame)==len(W.t):
        temp = quaternion.as_float_array(W.frame)
        tempfunc0 = interp1d(W.t[Indices],temp[Indices,0])
        tempfunc1 = interp1d(W.t[Indices],temp[Indices,1])
        tempfunc2 = interp1d(W.t[Indices],temp[Indices,2])
        tempfunc3 = interp1d(W.t[Indices],temp[Indices,3])
        temp0=tempfunc0(W_matching.t)
        temp1=tempfunc1(W_matching.t)
        temp2=tempfunc2(W_matching.t)
        temp3=tempfunc3(W_matching.t)
        temp=W.frame[range(len(W_matching.t))]
        for i in range(len(W_matching.t)):
            temp[i]=quaternion.quaternion(temp0[i], temp1[i], temp2[i], temp3[i])
        W_matching.frame=temp
    outname='./W_matching.h5' # Output the matching data, in case it is changed when calculating angular velocity or rotating frame
    OutputToNRAR(W_matching, outname, FileWriteMode='w')
    # Calculate angular velocity
    temp,omega=scri.corotating_frame(W_matching,return_omega=True)
    W_matching=ReadFromNRAR(outname)
    os.remove(outname)
    omega=np.sqrt(omega[:,0]**2+omega[:,1]**2+omega[:,2]**2)
    omega_prime=np.diff(omega)
    return W_matching, omega, omega_prime

def Hybridize(t_start, data_dir, out_dir):
    """
    Align and hybridize given NR waveform with PN waveform, the matching region starts at t_start, and last 3 orbits
    """
    import h5py
    from scipy.optimize import basinhopping
    from scipy.optimize import least_squares
    import numpy as np
    import quaternion
    import time
    import os

    clock0=time.time()
# Get NR waveform
    NRFileName=data_dir+'/rhOverM_Asymptotic_GeometricUnits.h5/Extrapolated_N2.dir'
    W_NR=ReadFromNRAR(NRFileName)
    W_NR.t=W_NR.t-W_NR.max_norm_time()
    W_NR.data=-W_NR.data
    W_NR_corot=scri.to_corotating_frame(W_NR)
    W_NR=ReadFromNRAR(NRFileName)
    W_NR.t=W_NR.t-W_NR.max_norm_time()
    W_NR.data=-W_NR.data
# Get PN waveform
    PNFileName=data_dir+'/rhOverM_Inertial_PN.h5'
    W_PN=ReadFromNRAR(PNFileName)
    W_PN_corot=scri.to_corotating_frame(W_PN)
    W_PN=ReadFromNRAR(PNFileName)

# Get the initial angular velocity in matching region
    temp1, omega_NR, temp2=MatchingRegion(W_NR, t_start-1000, t_start)
    omega_0=omega_NR[-1]

# Set up the matching region data for NR, and get the corresponding angular velocity and frame
    t_pre=t_start-10*3.1416/omega_0
    t_end=t_start+10*3.1416/omega_0
    W_NR_matching_in, omega_NR, temp1=MatchingRegion(W_NR, t_pre, t_end)
    W_NR_matching_in.LM=W_NR.LM
    matching_outname='./W_NR_matching_in.h5'
    OutputToNRAR(W_NR_matching_in, matching_outname, FileWriteMode='w')
    W_NR_matching_corot, temp1, temp2=MatchingRegion(W_NR_corot, t_pre, t_end)
    W_PN_matching_in, omega_PN, omega_PN_prime=MatchingRegion(W_PN, t_pre, t_end)
    W_PN_matching_corot, temp1, temp2=MatchingRegion(W_PN_corot, t_pre, t_end)
    print("After interpolate:",time.time()-clock0)

# Get initial guess of time alignment by matching angular velocity
    t_start_indice=time_indices(W_NR_matching_corot.t, t_start, t_start+0.01)[0]
    def minix(x):
        print(x)
        Indices_PN=time_indices(W_PN_matching_corot.t, t_start+x, t_start+x+6*3.1416/omega_0)
        dt=t_start+x-W_PN_matching_corot.t[Indices_PN[0]]
        omega_PN_matching=omega_PN[Indices_PN]+omega_PN_prime[Indices_PN]*dt
        return 1e6*np.sum((omega_NR[t_start_indice:t_start_indice+len(Indices_PN)]-omega_PN_matching)**2)
    mint=basinhopping(minix, 0.0, niter=200, T=0.001, stepsize=40.0)
    mint=least_squares(minix, mint.x, bounds=(mint.x-3.1416/omega_0,mint.x+3.1416/omega_0),ftol=2.23e-16,xtol=2.23e-16,gtol=2.23e-16)
    print(mint)
    t_delta=-mint.x
    print("Initial guess of t:",t_delta)
    clock1=time.time()

# Alignment of time and frame
    def Optimize4D(x):
        print(x)
        temp=1j;
        R_delta=quaternion.quaternion(1.0,0.0,0.0,0.0)
        if x[1]**2-x[2]**2-x[3]**2>1 or x[0]+t_delta>3.1416/omega_0 or x[0]+t_delta<-3.1416/omega_0:
            temp=10000*(x[0]+t_delta)**2
        else:
            R_delta1=quaternion.quaternion(np.sqrt(abs(1-x[1]**2-x[2]**2-x[3]**2)),x[1],x[2],x[3])
            R_delta2=quaternion.quaternion(-np.sqrt(abs(1-x[1]**2-x[2]**2-x[3]**2)),x[1],x[2],x[3])
            R_delta3=-quaternion.quaternion(np.sqrt(abs(1-x[1]**2-x[2]**2-x[3]**2)),x[1],x[2],x[3])
            R_delta4=-quaternion.quaternion(-np.sqrt(abs(1-x[1]**2-x[2]**2-x[3]**2)),x[1],x[2],x[3])
            Indices_PN=time_indices(W_PN_matching_corot.t, t_start+x[0], t_start+x[0]+6*3.1416/omega_0)
            W_temp1=ReadFromNRAR(matching_outname)
            W_temp2=ReadFromNRAR(matching_outname)
            W_temp3=ReadFromNRAR(matching_outname)
            W_temp4=ReadFromNRAR(matching_outname)
            W_temp1=scri.rotate_decomposition_basis(W_temp1, R_delta1)
            W_temp2=scri.rotate_decomposition_basis(W_temp2, R_delta2)
            W_temp3=scri.rotate_decomposition_basis(W_temp3, R_delta3)
            W_temp4=scri.rotate_decomposition_basis(W_temp4, R_delta4)
            temp1=0j
            temp2=0j
            temp3=0j
            temp4=0j
            for i_m in range(W_NR.LM.shape[0]):
                temp11=W_PN_matching_in.data[Indices_PN][i_m]-W_temp1.data[t_start_indice:t_start_indice+len(Indices_PN)][i_m]
                temp22=W_PN_matching_in.data[Indices_PN][i_m]-W_temp2.data[t_start_indice:t_start_indice+len(Indices_PN)][i_m]
                temp33=W_PN_matching_in.data[Indices_PN][i_m]-W_temp3.data[t_start_indice:t_start_indice+len(Indices_PN)][i_m]
                temp44=W_PN_matching_in.data[Indices_PN][i_m]-W_temp4.data[t_start_indice:t_start_indice+len(Indices_PN)][i_m]
                temp11=temp11*temp11.conjugate()
                temp22=temp22*temp22.conjugate()
                temp33=temp33*temp33.conjugate()
                temp44=temp44*temp44.conjugate()
                temp1+=temp11
                temp2+=temp22
                temp3+=temp33
                temp4+=temp44
            temp1=np.real(sum(temp1))
            temp2=np.real(sum(temp2))
            temp3=np.real(sum(temp3))
            temp4=np.real(sum(temp4))
            temp=min(temp1,temp2,temp3,temp4)
            if temp1==temp:
                R_delta=R_delta1
            elif temp2==temp:
                R_delta=R_delta2
            elif temp3==temp:
                R_delta=R_delta3
            else:
                R_delta=R_delta4
            print(temp)
        return temp, R_delta
    def _Optimize4D(x):
        return Optimize4D(x)[0]
    mini=least_squares(_Optimize4D, [-t_delta,0.0,0.0,0.0], bounds=([-t_delta-3.1416/omega_0,-1.0,-1.0,-1.0],[-t_delta+3.1416/omega_0,1.0,1.0,1.0]))
    print("Optimization time used:",time.time()-clock1)
    print(mini)
    print("Time shift=", -mini.x[0])
    R_delta=Optimize4D(mini.x)[1]
    print("R_delta=",R_delta)
    W_PN.t=W_PN.t-mini.x[0]
    W_NR=scri.rotate_decomposition_basis(W_NR, R_delta)
    os.remove(matching_outname)

# Hybridize waveform
    t_end0=t_start+6*3.1416/omega_0
    W_matching_NR, temp1, temp2=MatchingRegion(W_NR, t_pre, t_end)
    W_matching_PN, temp1, temp2=MatchingRegion(W_PN, t_pre, t_end)
    t_start_indicePN=time_indices(W_PN.t, t_start, t_start+0.01)[0]
    t_start_indiceNR=time_indices(W_NR.t, t_end0, t_end0+0.01)[1]
    t_start_indiceMatching=time_indices(W_matching_NR.t, t_start, t_start+0.01)[0]
    t_end_indiceMatching=time_indices(W_matching_NR.t, t_end0, t_end0+0.01)[1]
    Ind=np.arange(t_start_indiceMatching,t_end_indiceMatching)
    W_H=scri.waveform_base.WaveformBase()
    W_H.t=np.append(np.append(W_PN.t[0:t_start_indicePN], W_matching_NR.t[Ind]), W_NR.t[t_start_indiceNR:-1])
    W_H.LM=W_PN.LM
    W_H.data=1j*np.empty((len(W_H.t), len(W_H.LM)))
    N=len(Ind)
    xx=np.arange(N)/N
    # Hybridize data
    for i_m in range(W_H.LM.shape[0]):
        matching_data=(1-smooth(xx))*W_matching_PN.data[Ind,i_m]+smooth(xx)*W_matching_NR.data[Ind,i_m]
        W_H.data[:,i_m]=np.append(np.append(W_PN.data[0:t_start_indicePN,i_m],matching_data),W_NR.data[t_start_indiceNR:-1,i_m])
    print("finished")

# Output results
    outname=out_dir+'/hybridHybrid'+str(t_start)+'.h5'
    OutputToNRAR(W_H, outname, FileWriteMode='w')
    outname=out_dir+'/hybridNR'+str(t_start)+'.h5'
    OutputToNRAR(W_NR, outname, FileWriteMode='w')
    outname=out_dir+'/hybridPN'+str(t_start)+'.h5'
    OutputToNRAR(W_PN, outname, FileWriteMode='w')
    print("Total time:",time.time()-clock0)

def Run():
    for i in [-32000]:
           Hybridize(i,'/home/dzsun/SimAnnex/Public/HybTest/006/Lev3','/home/dzsun') 
