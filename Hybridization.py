"""NR/PN hybridization driver.

This module keeps the historical numerical algorithms intact while exposing
the command-line workflow through an import-safe ``main`` function.
"""

import argparse
import hashlib
import json
import os

import scri
from scipy.ndimage import gaussian_filter
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
from scipy.optimize import approx_fprime
from scipy.integrate import simpson
from scipy.interpolate import InterpolatedUnivariateSpline as Spline
from scipy.interpolate import CubicSpline
import copy
from Provenance import (
    BMSProvenanceRecorder,
    abd_fingerprint,
    optimizer_summary,
    transform_summary,
    waveform_fingerprint,
)
#from memory_profiler import profile

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


class hyb_quantites:
    def __init__(self, t_end, length):
        self.PNIter = 0
        self.cost=[]
        self.length = length
        self.t_start = t_end - length
        self.t_end = t_end
        self.t_PNStart = -5000.0
        self.t_PNEnd = length + 1000.0
        self.omega_NR_mag = None
        self.omega_i = None
        self.matchingt = None
        self.omega_NR_mag_matching = None
        self.W_NR_matching_in = None
        self.omega_NR_matching = None
        self.omega_mean = None
        self.omega_NR_hat = None
        self.Normalization = None
        self.NormalizationOmega = None
        self.omega_PN_spline = None
        self.omega_PN_mag_spline = None
        self.PNData_spline = None

    def get_window_NR(self, W_NR):
        omega_NR = W_NR.angular_velocity()
        self.omega_NR_mag = np.linalg.norm(omega_NR, axis=1)
        omega_NR_mag = gaussian_filter(self.omega_NR_mag, sigma=100)

        if W_NR.t[0]>=self.t_start-10 or W_NR.t[-1]<=self.t_start+10:
            message = ("t_start {0} should be much larger than the start time of NR"
                +" waveform {1} and much smaller than the ending time {2}.")
            raise ValueError(message.format(self.omega_it_start, W_NR.t[0], W_NR.t[-1]))
        omega_0 = omega_NR_mag[(W_NR.t>=self.t_start-10)&(W_NR.t<=self.t_start+10)]
        self.omega_i = np.mean(omega_0)

        if W_NR.t[0]>=self.t_start-2.0*np.pi/self.omega_i or W_NR.t[-1]<=self.t_start+2.0*np.pi/self.omega_i:
            message = ("t_start {0} should be at least {1} larger than the start time of NR waveform {2}"
                +" and smaller than the ending time of NR waveform {3}.")
            raise ValueError(message.format(self.t_start, 2.0*np.pi/self.omega_i, W_NR.t[0], W_NR.t[-1]))

        self.matchingt = W_NR.t[(W_NR.t>=self.t_start)&(W_NR.t<=self.t_end)]
        self.omega_NR_mag_matching = omega_NR_mag[(W_NR.t>=self.t_start)&(W_NR.t<=self.t_end)]
        self.W_NR_matching_in = W_NR.interpolate(self.matchingt)
        self.omega_NR_matching = omega_NR[(W_NR.t>=self.t_start) & (W_NR.t<=self.t_end)]
        self.omega_mean = np.mean(self.omega_NR_matching, axis=0)
        self.omega_NR_hat = self.omega_NR_matching / np.linalg.norm(self.omega_NR_matching, axis=1)[:, np.newaxis]
        self.Normalization = np.sum(simpson(abs(self.W_NR_matching_in.data)**2.0, x=self.matchingt, axis=0))
        self.NormalizationOmega = simpson(np.linalg.norm(self.omega_NR_matching, axis=1), x=self.matchingt)

    def get_window_PN(self, W_PN):
        omega_PN = W_PN.angular_velocity()
        self.omega_PN_spline = SplineArray(W_PN.t, omega_PN)
        omega_PN_mag = gaussian_filter(np.linalg.norm(omega_PN, axis=1),100)
        self.omega_PN_mag_spline = Spline(W_PN.t, omega_PN_mag)
        self.PNData_spline = SplineArray(W_PN.t, W_PN.data)


class PNParameters:
    def __init__(self, data_dir, hyb, t0):
        self.xHat = quaternion.quaternion(0.0,1.0,0.0,0.0)
        self.yHat = quaternion.quaternion(0.0,0.0,1.0,0.0)
        self.zHat = quaternion.quaternion(0.0,0.0,0.0,1.0)

        with h5py.File(data_dir+'/Horizons.h5', 'r') as f:
            tA = f['AhA.dir/CoordCenterInertial.dat'][:,0] + t0
            xA = f['AhA.dir/CoordCenterInertial.dat'][:,1:]
            mA = f['AhA.dir/ChristodoulouMass.dat'][:,1]
            chiA = f['AhA.dir/chiInertial.dat'][:,1:]
            xB = f['AhB.dir/CoordCenterInertial.dat'][:,1:]
            mB = f['AhB.dir/ChristodoulouMass.dat'][:,1]
            chiB = f['AhB.dir/chiInertial.dat'][:,1:]
        self.i = abs(tA - hyb.t_start).argmin()
        self.tA = tA

        m1 = mA[self.i]
        m2 = mB[self.i]
        self.q_i = m1/m2
        if m1 < m2:
            self.q_i = 1.0/self.q_i

        self.M_i = m1 + m2

        # calculate frame_i
        d = xA - xB
        nHat = np.empty((len(d),4))
        for i in range(len(d)):
            nHat[i] = np.append(0,d[i])
        nHatArray = nHat
        nHat = quaternion.from_float_array(nHat)
        for i in range(len(d)):
            nHat[i] = nHat[i].normalized()
        dnHatdt = CubicSpline(tA, nHatArray).derivative()(tA)
        lambdaHat = quaternion.from_float_array(dnHatdt)
        for i in range(len(d)):
            lambdaHat[i] = lambdaHat[i].normalized()
        Ra = np.sqrt(-nHat[self.i]*self.xHat)
        beta = math.atan2(np.dot((Ra*self.zHat*Ra.inverse()).vec, lambdaHat[self.i].vec),
            np.dot((Ra*self.yHat*Ra.inverse()).vec, lambdaHat[self.i].vec))
        self.frame_i = Ra*np.exp((beta)/2*self.xHat)
        self.frame_i = quaternion.as_float_array(self.frame_i)

        self.chiA = chiA
        self.chiB = chiB
        self.chi1_i = chiA[self.i]
        self.chi2_i = chiB[self.i]

        self.OptParas = np.array([1.0,1.0,1.0,1.0,0.0,0.0,0.0,0.0])
        self.PhyParas = np.array([1.0,1.0,1.0,1.0,0.0,0.0,0.0,0.0])


    def rotate(self, trans):
        self.frame_i = quaternion.from_float_array(self.frame_i)
        self.frame_i = self.frame_i*quaternion.from_float_array(trans["transformations"]["frame_rotation"])
        self.frame_i = quaternion.as_float_array(self.frame_i)

        chiA = quaternion.from_float_array(np.column_stack((0.0*self.tA, self.chiA)))
        chiB = quaternion.from_float_array(np.column_stack((0.0*self.tA, self.chiB)))
        chiA = quaternion.from_float_array(trans["transformations"]["frame_rotation"])*chiA*quaternion.from_float_array(trans["transformations"]["frame_rotation"]).inverse()
        chiB = quaternion.from_float_array(trans["transformations"]["frame_rotation"])*chiB*quaternion.from_float_array(trans["transformations"]["frame_rotation"]).inverse()
        self.chi1_i = chiA[self.i].vec
        self.chi2_i = chiB[self.i].vec

    def Physical_to_Parameterize(self, hyb):
        """
        Reparameterize physical parameters x=[q, M, chi1_x, chi1_y, chi1_z, chi2_x, chi2_y, chi2_z, t_delta,
        logR_delta_x, logR_delta_y, logR_delta_z] for optimization. The last four components of x are optional.
        """
        x = np.copy(self.PhyParas)
        if len(x)!=12 and len(x)!=8 and np.linalg.norm(np.cross(self.chi1_i, self.chi2_i))<1e-8:
            message = ("Incorrect format for x={0}. x should be [q, M, chi1_x, chi1_y, chi1_z, chi2_x, chi2_y, chi2_z, t_delta,"
                +"logR_delta_x, logR_delta_y, logR_delta_z], where the last four components are optional.")
            raise ValueError(message.format(x))

        x_2 = np.linalg.norm(x[2:5])/np.linalg.norm(self.chi1_i)
        x_3 = np.linalg.norm(x[5:8])/np.linalg.norm(self.chi2_i)

        rotation = quaternion.optimal_alignment_in_Euclidean_metric(x[2:5],x_2*self.chi1_i)
        axis = quaternion.from_float_array(np.append(0,x[2:5])).normalized()
        chi1_hat = x[2:5]/np.linalg.norm(x[2:5])
        chi2_i_proj = x_3*(rotation*quaternion.quaternion(0,self.chi2_i[0],self.chi2_i[1],self.chi2_i[2])*rotation.conjugate()).vec
        chi2_i_proj = chi2_i_proj - np.dot(chi2_i_proj, chi1_hat)*chi1_hat
        chi2_0_proj = x[5:8] - np.dot(x[5:8], chi1_hat)*chi1_hat
        angle = np.arccos(np.dot(chi2_i_proj, chi2_0_proj)/(np.linalg.norm(chi2_i_proj)*np.linalg.norm(chi2_0_proj)))
        if np.dot(chi2_i_proj, chi2_0_proj)/(np.linalg.norm(chi2_i_proj)*np.linalg.norm(chi2_0_proj)) > 1:
            angle = 0.0
        elif np.dot(chi2_i_proj, chi2_0_proj)/(np.linalg.norm(chi2_i_proj)*np.linalg.norm(chi2_0_proj)) < -1:
            angle = np.pi
        sign = np.sign(np.dot(chi1_hat,np.cross(chi2_i_proj, chi2_0_proj)))
        rotation = np.exp(sign*axis*angle/2)*rotation

        x_47 = (quaternion.as_float_array(np.log(rotation)))[1:]
        x_47[x_47>np.pi] = x_47[x_47>np.pi] - np.pi

        chi2_0 = x_3*(rotation*quaternion.quaternion(0,self.chi2_i[0],self.chi2_i[1],self.chi2_i[2])*rotation.conjugate()).vec
        sign = np.sign(np.dot(np.cross(chi2_0, chi1_hat), np.cross(chi2_0, x[5:8])))
        x_7 = sign*np.arccos(np.dot(chi2_0, x[5:8])/(np.linalg.norm(x[5:8])*np.linalg.norm(chi2_0)))
        if np.dot(chi2_0,x[5:8])/(np.linalg.norm(x[5:8])*np.linalg.norm(chi2_0)) > 1:
            x_7 = 0.0

        x[0] = x[0]/self.q_i
        x[1] = x[1]/self.M_i
        x[2] = x_2
        x[3] = x_3
        x[4:7] = x_47
        x[7] = x_7
        if len(x) == 12:
            x[9:] = x[9:] - hyb.omega_mean*x[8]/2
        self.OptParas = x


    def Parameterize_to_Physical(self, hyb):
        """
        Output optimization parameters as physical quantities x=[q, M, chi1_x, chi1_y, chi1_z, chi2_x, chi2_y, chi2_z, t_delta,
        logR_delta_x, logR_delta_y, logR_delta_z]. The last four components of x are optional.
        """
        x = np.copy(self.OptParas)
        if len(x)!=12 and len(x)!=8:
            message = ("Incorrect format for x={0}")
            raise ValueError(message.format(x))

        x[0] = x[0]*self.q_i
        x[1] = x[1]*self.M_i

        chi1_0 = quaternion.quaternion(0,self.chi1_i[0],self.chi1_i[1],self.chi1_i[2])
        chi2_0 = quaternion.quaternion(0,self.chi2_i[0],self.chi2_i[1],self.chi2_i[2])
        rotation = np.exp(self.xHat*x[4] + self.yHat*x[5] + self.zHat*x[6])
        chi1_0 = x[2]*(rotation*chi1_0*rotation.conjugate()).vec
        chi2_0 = x[3]*rotation*chi2_0*rotation.conjugate()
        axis = quaternion.from_float_array(np.append(0, np.cross(chi2_0.vec, chi1_0))).normalized()
        angle = np.exp(axis*x[7]/2)
        chi2_0 = (angle*chi2_0*angle.conjugate()).vec
        x[2:5] = chi1_0
        x[5:8] = chi2_0

        if len(x) == 12:
            x[9:] = x[9:] + hyb.omega_mean*x[8]/2
        self.PhyParas = x


def get_extrapolated_NR(data_dir, nOrbits, t_end, length):
    W = scri.SpEC.read_from_h5(data_dir+'/rhOverM_Asymptotic_GeometricUnits_CoM.h5/Extrapolated_N4.dir')
    t0 = -W.max_norm_time()
    W.t = W.t+t0
    W_corot = scri.to_corotating_frame(W.copy())
    ZeroModes = [2,8,16,26,38,52,68]
    W_corot.data[:,ZeroModes] = 0.0*W_corot.data[:,ZeroModes]
    W = scri.to_inertial_frame(W_corot.copy())
    if nOrbits != None:
        omega = W.angular_velocity()
        omega_mag = np.linalg.norm(omega, axis=1)
        length = nOrbits_to_length(nOrbits, t_end, omega_mag,W.t)
    return W, t0, length


def get_abd(cce_dir, truncate, abd_format="RPDMB", provenance_recorder=None, outer_iteration=None):
    """
    Create AsymptoticBondiData from CCE data and apply truncation if required.
    """
    abd = scri.SpEC.create_abd_from_h5(
        abd_format,
        h = cce_dir + '/Strain.h5',
        Psi4 = cce_dir + '/Psi4.h5',
        Psi3 = cce_dir + '/Psi3.h5',
        Psi2 = cce_dir + '/Psi2.h5',
        Psi1 = cce_dir + '/Psi1.h5',
        Psi0 = cce_dir + '/Psi0.h5')
    if provenance_recorder:
        input_files = []
        for filename in ("Strain.h5", "Psi4.h5", "Psi3.h5", "Psi2.h5", "Psi1.h5", "Psi0.h5"):
            path = os.path.realpath(os.path.join(cce_dir, filename))
            stat = os.stat(path)
            input_files.append({"path": path, "size_bytes": int(stat.st_size), "mtime_ns": int(stat.st_mtime_ns)})
        provenance_recorder.record(
            "load_raw_cce_bondi_data", outer_iteration=outer_iteration,
            abd_format=abd_format, input_files=input_files, output_abd=abd_fingerprint(abd),
        )
    t0 = -abd.t[np.argmax(np.linalg.norm(abd.sigma.bar, axis=1))]
    before_time_shift = abd_fingerprint(abd) if provenance_recorder else None
    abd.t = abd.t + t0
    if provenance_recorder:
        provenance_recorder.record(
            "shift_abd_time_origin", outer_iteration=outer_iteration, time_shift=float(t0),
            definition="negative time of maximum strain norm", input_abd=before_time_shift,
            output_abd=abd_fingerprint(abd),
        )
    interpolation_input = abd_fingerprint(abd) if provenance_recorder else None
    interpolation_grid = np.arange(abd.t[0],abd.t[-1],1.0)
    abd = abd.interpolate(interpolation_grid)
    if provenance_recorder:
        provenance_recorder.record(
            "interpolate_abd", outer_iteration=outer_iteration, input_abd=interpolation_input,
            output_abd=abd_fingerprint(abd), grid_start=float(interpolation_grid[0]),
            grid_stop_exclusive=float(interpolation_grid[-1] + 1.0), grid_step=1.0,
            grid_size=int(len(interpolation_grid)),
        )

    if truncate != None:
        abd_prime = scri.asymptotic_bondi_data.AsymptoticBondiData(
            time = abd.t[(abd.t>=truncate[0])&(abd.t<truncate[1])],
            ell_max = 8,
            multiplication_truncator = max
        )
        abd_prime.sigma = abd.sigma[(abd.t>=truncate[0])&(abd.t<truncate[1])]
        abd_prime.psi0 = abd.psi0[(abd.t>=truncate[0])&(abd.t<truncate[1])]
        abd_prime.psi1 = abd.psi1[(abd.t>=truncate[0])&(abd.t<truncate[1])]
        abd_prime.psi2 = abd.psi2[(abd.t>=truncate[0])&(abd.t<truncate[1])]
        abd_prime.psi3 = abd.psi3[(abd.t>=truncate[0])&(abd.t<truncate[1])]
        abd_prime.psi4 = abd.psi4[(abd.t>=truncate[0])&(abd.t<truncate[1])]
        abd = abd_prime
        if provenance_recorder:
            provenance_recorder.record(
                "truncate_abd", outer_iteration=outer_iteration,
                interval=[float(truncate[0]), float(truncate[1])], output_abd=abd_fingerprint(abd),
            )

    return abd, t0


def abd_to_WM(abd, lmin = 2):
    W = scri.WaveformModes()
    W.t = abd.t
    W.data = 2*abd.sigma.bar
    if len(W.data[0]) > 77 and lmin == 2:
        W.data = np.copy(W.data[:,4:])
        W.ells = 2,8
    else:
        W.ells = 0,8
    W.frameType = scri.Inertial
    W.dataType = scri.h
    return W


def nOrbits_to_length(nOrbits, t_end, Omega_array, time_array):
    """
    Calculate the waveform length in M, given number of orbits.
    """
    tol = 1e-2
    number = 0
    a = time_array[0]
    b = t_end
    t0 = (a + b)/2
    orbits = simpson(Omega_array[(time_array>t0)&(time_array<t_end)],
        x=time_array[(time_array>t0)&(time_array<t_end)])/2/np.pi
    while abs(orbits-nOrbits) > tol and number <= 100:
        orbits = simpson(Omega_array[(time_array>t0)&(time_array<t_end)],
        x=time_array[(time_array>t0)&(time_array<t_end)])/2/np.pi
        if orbits < nOrbits:
            b = t0
        else:
            a = t0
        t0 = (a + b)/2
        number += 1
    if number > 100:
        message = ("The waveform before endtime only has {0} orbits, which is smaller than required {1}.")
        raise ValueError(message.format(orbits, nOrbits))
    return t_end - t0


def get_length_from_abd(abd, nOrbits, t_end):
    W = abd_to_WM(abd)
    omega = W.angular_velocity()
    omega_mag = np.linalg.norm(omega, axis=1)
    length = nOrbits_to_length(nOrbits, t_end, omega_mag, W.t)
    return length


#@profile
def fix_BMS(
    abd,
    hyb,
    PN,
    bms_iterations=4,
    bms_align_nprocs=None,
    bms_bruteforce_dt=1000,
    bms_bruteforce_phi=20,
    pn_evolution_order=6.0,
    pn_waveform_mode_order=6.0,
    provenance_recorder=None,
):
    if hyb.PNIter == 0:
        temp = time.time()
        abd_before = abd_fingerprint(abd) if provenance_recorder else None
        abd_prime, trans, abd_err = abd.map_to_superrest_frame(t_0=hyb.t_start+hyb.length/2)
        if provenance_recorder:
            provenance_recorder.record(
                "initial_superrest_map", objects={"local_transform": trans, "cumulative_transform": trans},
                outer_iteration=int(hyb.PNIter), map_center_time=float(hyb.t_start + hyb.length / 2),
                local_transform=transform_summary(trans), cumulative_transform=transform_summary(trans),
                abd_error=np.asarray(abd_err), abd_before=abd_before, abd_after=abd_fingerprint(abd_prime),
                selected_inner_index=0,
            )
        W_NR = abd_to_WM(abd_prime)
        waveform_before_zeroing = waveform_fingerprint(W_NR) if provenance_recorder else None
        W_NR_corot = scri.to_corotating_frame(W_NR.copy())
        ZeroModes = [2,8,16,26,38,52,68]
        W_NR_corot.data[:,ZeroModes] = 0.0*W_NR_corot.data[:,ZeroModes]
        W_NR = scri.to_inertial_frame(W_NR_corot.copy())
        if provenance_recorder:
            provenance_recorder.record(
                "zero_corotating_m0_modes", outer_iteration=int(hyb.PNIter),
                zero_mode_indices=ZeroModes, waveform_before=waveform_before_zeroing,
                waveform_after=waveform_fingerprint(W_NR),
            )
        waveform_before_rotation = waveform_fingerprint(W_NR) if provenance_recorder else None
        W_NR = scri.rotate_physical_system(W_NR, trans.frame_rotation.conjugate())
        if provenance_recorder:
            provenance_recorder.record(
                "apply_conjugate_selected_bms_rotation", outer_iteration=int(hyb.PNIter),
                rotation_wxyz=quaternion.as_float_array(trans.frame_rotation.conjugate()),
                source_transform="selected local initial superrest transform",
                waveform_before=waveform_before_rotation, waveform_after=waveform_fingerprint(W_NR),
            )
        print('BMS frame mapping time used:',time.time()-temp)
    else:
        Phys = PN.PhyParas
        W_PN = PostNewtonian.PNWaveform(
            Phys[0], np.copy(hyb.omega_i)*Phys[1], Phys[2:5], Phys[5:8], PN.frame_i, np.copy(hyb.t_start)/Phys[1],
            t_PNStart=hyb.t_PNStart, t_PNEnd=hyb.t_PNEnd,
            PNEvolutionOrder=pn_evolution_order, PNWaveformModeOrder=pn_waveform_mode_order,
        )
        W_PN.t = W_PN.t*Phys[1]
        W_PN.data = W_PN.data*Phys[1]#DQ
        #W_PN.data = np.append(0.0*W_PN.data[:,0:4], np.copy(W_PN.data), axis=1)
        W_PN.ells = 2,8
        W_PN.dataType = scri.h

        W_PN_PsiM = PostNewtonian.PNWaveform(
            Phys[0], np.copy(hyb.omega_i)*Phys[1], Phys[2:5], Phys[5:8], PN.frame_i, np.copy(hyb.t_start)/Phys[1],
            t_PNStart=hyb.t_PNStart, t_PNEnd=hyb.t_PNEnd, datatype="Psi_M",
            PNEvolutionOrder=pn_evolution_order, PNWaveformModeOrder=pn_waveform_mode_order,
        )
        W_PN_PsiM.t = W_PN_PsiM.t*Phys[1]
        W_PN_PsiM.data = W_PN_PsiM.data*Phys[1]#DQ

        if provenance_recorder:
            provenance_recorder.record(
                "construct_pn_bms_target", outer_iteration=int(hyb.PNIter),
                physical_parameters=np.asarray(Phys), omega_i=float(hyb.omega_i),
                frame_i_wxyz=quaternion.as_float_array(PN.frame_i),
                pn_evolution_order=float(pn_evolution_order),
                pn_waveform_mode_order=float(pn_waveform_mode_order),
                strain_target=waveform_fingerprint(W_PN), moreschi_target=waveform_fingerprint(W_PN_PsiM),
            )

        temp = time.time()
        tp1, W_NR, trans, idx = PNBMS.PN_BMS_w_time_phase(
            abd,
            W_PN,
            W_PN_PsiM,
            hyb.t_start,
            hyb.t_start+hyb.length,
            None,
            N=bms_iterations,
            nprocs=bms_align_nprocs,
            n_brute_force_δt=bms_bruteforce_dt,
            n_brute_force_δϕ=bms_bruteforce_phi,
            provenance_recorder=provenance_recorder,
            outer_iteration=int(hyb.PNIter),
        )
        waveform_before_rotation = waveform_fingerprint(W_NR) if provenance_recorder else None
        W_NR = scri.rotate_physical_system(W_NR, trans.frame_rotation.conjugate())
        W_NR.frame = np.array([])
        if provenance_recorder:
            provenance_recorder.record(
                "apply_conjugate_selected_bms_rotation", outer_iteration=int(hyb.PNIter),
                rotation_wxyz=quaternion.as_float_array(trans.frame_rotation.conjugate()),
                source_transform="selected local inner transform",
                waveform_before=waveform_before_rotation, waveform_after=waveform_fingerprint(W_NR),
            )
        print('BMS frame mapping time used:',time.time()-temp)
        print('BMS mapping error is: ', tp1)

    return W_NR, trans


def fix_BMS_NRNR(abd, abd2, hyb, PN, bms_iterations=4, bms_align_nprocs=None, bms_bruteforce_dt=1000, bms_bruteforce_phi=20):
    abd_prime, trans, abd_err = abd.map_to_superrest_frame(t_0=hyb.t_start+hyb.length/2)
    W_NR = abd_to_WM(abd_prime)

    Psi_M = PNBMS.MT_to_WM(abd_prime.supermomentum('Moreschi'), False, dataType = scri.psi2)
    tp1, W_NR2, tp2, idx = PNBMS.PN_BMS_w_time_phase(
        abd2,
        W_NR,
        Psi_M,
        hyb.t_start,
        hyb.t_start+hyb.length,
        None,
        N=bms_iterations,
        nprocs=bms_align_nprocs,
        n_brute_force_δt=bms_bruteforce_dt,
        n_brute_force_δϕ=bms_bruteforce_phi,
    )

    return W_NR, W_NR2


def InitialT(x, hyb):
    return simpson((hyb.omega_NR_mag_matching - hyb.omega_PN_mag_spline(hyb.matchingt + x))**2, x=hyb.matchingt)/hyb.NormalizationOmega


def InitialR(theta, hyb, t_delta, R_delta):
    R_temp = R_delta*np.exp(theta/2*quaternion.quaternion(0.0, hyb.omega_NR_hat[0,0], hyb.omega_NR_hat[0,1], hyb.omega_NR_hat[0,2]))
    W_temp = scri.rotate_physical_system(hyb.W_NR_matching_in.copy(), R_temp)
    cost = (np.angle(W_temp.data[int(len(hyb.matchingt)/2),4]) - np.angle(hyb.PNData_spline(t_delta + W_temp.t[int(len(hyb.matchingt)/2)])[0,4]))**2
    return cost


def Optimize4D(x, hyb):
    phase = quaternion.quaternion(0.0, hyb.omega_mean[0]*x[0]/2, hyb.omega_mean[1]*x[0]/2, hyb.omega_mean[2]*x[0]/2)
    R_delta = np.exp(quaternion.quaternion(0.0,x[1],x[2],x[3]) + phase)
    W_temp = scri.rotate_physical_system(hyb.W_NR_matching_in.copy(), R_delta)
    temp = hyb.PNData_spline(hyb.matchingt + x[0])
    cost = temp - W_temp.data
    return np.array([cost.real, cost.imag]).flatten('F')/np.sqrt(np.sum(abs(temp)**2.0))


def SquaredError(W1, W2, t1, t2, mode=None):
    """
    Calculate the residue of W1 and W2 between t1 and t2.
    """
    W2_spline = SplineArray(W2.t, W2.data)
    matchingt = W1.t[(W1.t>=t1)&(W1.t<=t2)]
    h1h2 = np.sum(simpson(abs(W2_spline(matchingt) - W1.data[(W1.t>=t1)&(W1.t<=t2),:])**2.0, x=matchingt, axis=0))
    h1h1 = np.sum(simpson(abs(W1.data[(W1.t>=t1)&(W1.t<=t2),:])**2.0, x=matchingt, axis=0))
    if type(mode) != type(None):
        h1h2 = np.sum(simpson(abs(W2_spline(matchingt) - W1.data[(W1.t>=t1)&(W1.t<=t2),:])[:,mode]**2.0, x=matchingt, axis=0))
    return 0.5*h1h2/h1h1


def ModeSquaredErrors(W1, W2, t1, t2):
    W2_spline = SplineArray(W2.t, W2.data)
    mask = (W1.t>=t1)&(W1.t<=t2)
    matchingt = W1.t[mask]
    diff = W2_spline(matchingt) - W1.data[mask, :]
    h1h1 = np.sum(simpson(abs(W1.data[mask, :])**2.0, x=matchingt, axis=0))
    return list(0.5*simpson(abs(diff)**2.0, x=matchingt, axis=0)/h1h1)


def StandardError(minima, PN):
    """
    Calculate standard error from least_squares optimization.
    """
    J = minima.jac
    cov = np.linalg.inv(J.T.dot(J))
    var = np.sqrt(np.diagonal(cov))
    var = np.sqrt(np.sum(minima.fun**2)/(minima.fun.size/77.0 - minima.x.size))*var

    if np.linalg.norm(PN.chi1_i)<1e-4 or np.linalg.norm(PN.chi2_i)<1e-4:
        var[7] = 1e-5
        if np.linalg.norm(PN.chi1_i)<1e-4:
            var[2] = 1e-5
        if np.linalg.norm(PN.chi2_i)<1e-4:
            var[3] = 1e-5

    return var


def Align(PN, hyb, W_NR2=None, pn_evolution_order=6.0, pn_waveform_mode_order=6.0):
    """
    Generate PN waveform and align it with NR waveform.
    """
    Phys = PN.PhyParas

    print(("Generating PN with parameters q={0}, M={4}, omega_0={1}, chi1_0={2}, chi2_0={3}.").format(
        Phys[0], hyb.omega_i, Phys[2:5], Phys[5:8], Phys[1]))
    W_PN_corot, chi1, chi2 = PostNewtonian.PNWaveform(
        Phys[0], np.copy(hyb.omega_i)*Phys[1], Phys[2:5], Phys[5:8], PN.frame_i, np.copy(hyb.t_start)/Phys[1],
        t_PNStart=hyb.t_PNStart, t_PNEnd=hyb.t_PNEnd, frametype="corotating", return_chi=True,
        PNEvolutionOrder=pn_evolution_order, PNWaveformModeOrder=pn_waveform_mode_order,
    )

    if hyb.PNIter == 0:
        ZeroModes = [2,8,16,26,38,52,68] # Memory modes
        W_PN_corot.data[:,ZeroModes] = 0.0*W_PN_corot.data[:,ZeroModes] # Not consider memory effect since NR dosen't have corrrect memory.
    W_PN = scri.to_inertial_frame(W_PN_corot.copy())
    W_PN.t = W_PN.t*Phys[1]
    W_PN.data = W_PN.data*Phys[1]#DQ

    if W_NR2 != None:
        W_PN = W_NR2

    # Set up the matching region data for PN, and get the corresponding angular velocity and frame
    hyb.get_window_PN(W_PN)

    # Get initial guess
    Initial1 = np.array([0.0,0.0,0.0,0.0])
    Initial2 = np.array([0.0,0.0,0.0,0.0])
    if len(Phys) == 12:
        Initial1[0] = Phys[8]
        Initial1[1:] = Phys[9:] - hyb.omega_mean*Phys[8]/2

        R_delta = np.exp(quaternion.quaternion(0.0,Phys[9],Phys[10],Phys[11]))
        R_delta = R_delta*np.exp(np.pi/2*quaternion.quaternion(0.0, hyb.omega_NR_hat[0,0], hyb.omega_NR_hat[0,1],hyb.omega_NR_hat[0,2]))
        Initial2[0] = Phys[8]
        Initial2[1:] = quaternion.as_float_array(np.log(R_delta))[1:] - hyb.omega_mean*Phys[8]/2
    else:
        # Get initial guess of time alignment by matching angular velocity
        mint = least_squares(InitialT, 0.0, bounds=[-10*np.pi/hyb.omega_i, 10*np.pi/hyb.omega_i], ftol=1e-14, xtol=1e-14, gtol=1e-14, args=([hyb]))

        # Get initial guess of frame alignment
        R_delta = quaternion.optimal_alignment_in_Euclidean_metric(
            hyb.omega_NR_matching, hyb.omega_PN_spline(hyb.matchingt + mint.x), hyb.matchingt
        )
        minf = least_squares(InitialR, 0.0, bounds=[-np.pi,np.pi], ftol=1e-14, xtol=1e-14, gtol=1e-14, args=(hyb, mint.x, R_delta))

        # Pi degeneracy
        phase = quaternion.quaternion(0.0, hyb.omega_mean[0]*mint.x/2, hyb.omega_mean[1]*mint.x/2, hyb.omega_mean[2]*mint.x/2)
        R_delta1 = R_delta*np.exp(minf.x/2*quaternion.quaternion(0.0, hyb.omega_NR_hat[0,0], hyb.omega_NR_hat[0,1], hyb.omega_NR_hat[0,2]))
        R_delta2 = R_delta1*np.exp(np.pi/2*quaternion.quaternion(0.0, hyb.omega_NR_hat[0,0], hyb.omega_NR_hat[0,1], hyb.omega_NR_hat[0,2]))
        logR_delta1 = quaternion.as_float_array(np.log(R_delta1)) - np.append(0.0, hyb.omega_mean*mint.x/2)
        logR_delta2 = quaternion.as_float_array(np.log(R_delta2)) - np.append(0.0, hyb.omega_mean*mint.x/2)
        Initial1[0] = mint.x[0]
        Initial2[0] = mint.x[0]
        Initial1[1:] = logR_delta1[0][1:]
        Initial2[1:] = logR_delta2[0][1:]

    # Alignment of time and frame
    scale = [np.pi/hyb.omega_i/2, np.pi, np.pi, np.pi]
    lowbound1 = Initial1 - scale
    upbound1 = Initial1 + scale
    lowbound2 = Initial2 - scale
    upbound2 = Initial2 + scale
    minima1 = least_squares(Optimize4D, Initial1, bounds=(lowbound1,upbound1), ftol=3e-15, xtol=3e-15, gtol=1e-8, x_scale='jac', args=([hyb]))
    minima2 = least_squares(Optimize4D, Initial2, bounds=(lowbound2,upbound2), ftol=3e-15, xtol=3e-15, gtol=1e-8, x_scale='jac', args=([hyb]))
    if minima1.cost > minima2.cost:
        minima1.x = minima2.x
        minima1.cost = minima2.cost
    return minima1, W_PN, chi1, chi2


def Optimize12D(x, PN, hyb, pn_evolution_order=6.0, pn_waveform_mode_order=6.0):
    """
    Generate PN waveform and align it with NR waveform.
    """
    PN.OptParas = x
    PN.Parameterize_to_Physical(hyb)
    Phys = PN.PhyParas
    W_PN_corot = PostNewtonian.PNWaveform(
        Phys[0], np.copy(hyb.omega_i)*Phys[1], Phys[2:5], Phys[5:8], PN.frame_i, np.copy(hyb.t_start)/Phys[1],
        t_PNStart=hyb.t_PNStart, t_PNEnd=hyb.t_PNEnd, frametype="corotating",
        PNEvolutionOrder=pn_evolution_order, PNWaveformModeOrder=pn_waveform_mode_order,
    )

    if hyb.PNIter == 0:
        ZeroModes = [2,8,16,26,38,52,68] # Memory modes
        W_PN_corot.data[:,ZeroModes] = 0.0*W_PN_corot.data[:,ZeroModes] # Not cosider memory effect since NR dosen't have corrrect memory.
    W_PN = scri.to_inertial_frame(W_PN_corot.copy())
    W_PN.t = W_PN.t*Phys[1]
    W_PN.data = W_PN.data*Phys[1]

    hyb.get_window_PN(W_PN)

    return Optimize4D(x[8:], hyb)


def Stitch(W_PN, W_NR, hyb):
    """
    Stitch PN and NR waveforms.
    """
    PNData_spline = SplineArray(W_PN.t, W_PN.data)
    tTemp = np.array(np.append(W_PN.t[W_PN.t<hyb.t_start], W_NR.t[W_NR.t>=hyb.t_start]))
    dataTemp = np.empty((len(tTemp), len(W_NR.LM)), dtype=complex)
    Smooth = np.resize(scri.utilities.transition_function(hyb.matchingt, hyb.t_start, hyb.t_end), (len(W_NR.LM), len(hyb.matchingt))).T

    # Stitch data
    matching_data = (1 - Smooth)*PNData_spline(hyb.matchingt) + Smooth*W_NR.data[(W_NR.t>=hyb.t_start)&(W_NR.t<=hyb.t_end),:]
    dataTemp[tTemp<hyb.t_start,:] = W_PN.data[W_PN.t<hyb.t_start,:]
    dataTemp[(tTemp>=hyb.t_start)&(tTemp<=hyb.t_end),:] = matching_data
    dataTemp[tTemp>hyb.t_end,:] = W_NR.data[W_NR.t>hyb.t_end,:]

    # Delete indices that cause tiny time step
    minstep = min(min(np.diff(W_NR.t[(W_NR.t>hyb.t_start-10)&(W_NR.t<hyb.t_end+10)])),
                min(np.diff(W_PN.t[(W_PN.t>hyb.t_start-10)&(W_PN.t<hyb.t_end+10)])))
    BadIndices = np.nonzero(np.append(np.diff(tTemp)<minstep,0)&(tTemp>hyb.t_start-10)&(tTemp<hyb.t_end+10))
    while len(BadIndices[0]) > 0:
        tTemp = np.delete(tTemp, BadIndices)
        dataTemp = np.delete(dataTemp, BadIndices, axis=0)
        BadIndices = np.nonzero(np.append(np.diff(tTemp)<minstep,0)&(tTemp>hyb.t_start-10)&(tTemp<hyb.t_end+10))

    # Construct Hybrid waveform
    W_H = scri.WaveformModes()
    W_H.t = tTemp
    W_H.data = dataTemp
    W_H.ells = min(W_NR.LM[:, 0]), max(W_NR.LM[:, 0])
    W_H.dataType = scri.h
    W_H.frameType = scri.Inertial
    return W_H


def Output(out_name, W_NR, W_PN, W_H, trans, minima12D, PN, hyb, nOrbits, TestWindowError=False, Uncertainty=True):
    outname = 'hybridNR.h5'
    scri.SpEC.write_to_h5(W_NR, outname, file_write_mode='w')
    outname = 'hybridPN.h5'
    scri.SpEC.write_to_h5(W_PN, outname, file_write_mode='w')
    outname = 'hybrid.h5'
    scri.SpEC.write_to_h5(W_H, outname, file_write_mode='w')

    ErrorMatchingWindow = SquaredError(W_NR, W_PN, hyb.t_start, hyb.t_start + hyb.length)
    hyb.cost.append(ErrorMatchingWindow)

    checkpoint = np.array([hyb.PNIter, hyb.cost, hyb.omega_i, hyb.t_start, hyb.length, hyb.t_PNStart, hyb.t_PNEnd, PN.frame_i], dtype="object")

    ErrorTestWindow = []
    if TestWindowError:
        t1 = hyb.t_end - nOrbits_to_length(25+nOrbits/2, hyb.t_end, hyb.omega_NR_mag, W_NR.t)
        length = nOrbits_to_length(10, t1, hyb.omega_NR_mag, W_NR.t)
        ErrorTestWindow = SquaredError(W_NR, W_PN, t1-length, t1)
        print('Test Window Error = ', ErrorTestWindow)

    ModeError = ModeSquaredErrors(W_NR, W_PN, hyb.t_start, hyb.t_start + hyb.length)
    print(ModeError)

    change = []
    if Uncertainty:
        #var = StandardError(minima12D, PN)
        change.append(PN.PhyParas[:2]*(1 - 1/PN.OptParas[:2]))
        change.append(np.linalg.norm(PN.PhyParas[2:5])*(1 - 1/PN.OptParas[2]))
        change.append(np.linalg.norm(PN.PhyParas[5:8])*(1 - 1/PN.OptParas[3]))
        change.append(PN.OptParas[7])
        #change.append(PN.PhyParas[:2]/PN.OptParas[:2]*var[:2])
        #change.append(np.linalg.norm(PN.PhyParas[2:5])/PN.OptParas[2]*var[2])
        #change.append(np.linalg.norm(PN.PhyParas[5:8])/PN.OptParas[3]*var[3])
        #change.append(var[7])
        change = np.array(change, dtype="object")

    np.savez(
        out_name,
        checkpoint = checkpoint,
        ErrorMatchingWindow = ErrorMatchingWindow,
        ErrorTestWindow = ErrorTestWindow,
        BMSTrans = trans,
        OptParas = PN.OptParas,
        PhyParas = PN.PhyParas,
        ModeError = ModeError,
        change = change
    )


#@profile
def Hybridize(
    WaveformType,
    t_end,
    sim_dir,
    cce_dir,
    out_name,
    length,
    nOrbits,
    hyb,
    PN,
    debug=0,
    OptimizePNParas=0,
    truncate=None,
    cce_dir2=None,
    abd_format="RPDMB",
    bms_iterations=4,
    bms_align_nprocs=None,
    bms_bruteforce_dt=1000,
    bms_bruteforce_phi=20,
    pn_order=6.0,
    provenance_recorder=None,
    fixed_intrinsic_parameters=None,
    fixed_omega_i=None,
    pn_evolution_order=None,
    pn_waveform_mode_order=None,
):
    """
    Align and hybridize given NR waveform with PN waveform.
    """

    pn_evolution_order = pn_order if pn_evolution_order is None else pn_evolution_order
    pn_waveform_mode_order = pn_order if pn_waveform_mode_order is None else pn_waveform_mode_order

    # Get NR waveform
    if WaveformType == 'cce':
        abd, t0 = get_abd(cce_dir, truncate, abd_format, provenance_recorder, int(hyb.PNIter))
        if cce_dir2 != None:
            abd2, t02 = get_abd(cce_dir2, truncate, abd_format)
            OptimizePNParas = 0
    else:
        W_NR, t0, length = get_extrapolated_NR(data_dir, nOrbits, t_end, length)


    # BMS tranformations
    trans = []
    if WaveformType == 'cce':
        if cce_dir2 == None:
            W_NR, trans = fix_BMS(
                abd,
                hyb,
                PN,
                bms_iterations=bms_iterations,
                bms_align_nprocs=bms_align_nprocs,
                bms_bruteforce_dt=bms_bruteforce_dt,
                bms_bruteforce_phi=bms_bruteforce_phi,
                pn_evolution_order=pn_evolution_order,
                pn_waveform_mode_order=pn_waveform_mode_order,
                provenance_recorder=provenance_recorder,
            )
            #PN.rotate(trans)
        else:
            W_NR, W_NR2 = fix_BMS_NRNR(
                abd,
                abd2,
                hyb,
                PN,
                bms_iterations=bms_iterations,
                bms_align_nprocs=bms_align_nprocs,
                bms_bruteforce_dt=bms_bruteforce_dt,
                bms_bruteforce_phi=bms_bruteforce_phi,
            )


    # Get matching window
    hyb.get_window_NR(W_NR)
    if fixed_omega_i is not None:
        hyb.omega_i = float(fixed_omega_i)


    # Optimize parameters
    if fixed_intrinsic_parameters is not None:
        if len(PN.PhyParas) == 12:
            PN.PhyParas[:8] = fixed_intrinsic_parameters
        else:
            PN.PhyParas = np.asarray(fixed_intrinsic_parameters, dtype=float).copy()
        PN.Physical_to_Parameterize(hyb)
        if provenance_recorder:
            provenance_recorder.record(
                "enforce_fixed_intrinsic_parameters",
                outer_iteration=int(hyb.PNIter),
                fixed_physical_parameters=np.asarray(PN.PhyParas[:8]),
                fixed_omega_i=float(hyb.omega_i),
            )
    elif hyb.PNIter == 0:
        PN.Parameterize_to_Physical(hyb)
    else:
        PN.Physical_to_Parameterize(hyb)

    minima12D=[]
    if OptimizePNParas:
        if provenance_recorder:
            provenance_recorder.record(
                "begin_12d_optimization", outer_iteration=int(hyb.PNIter),
                initial_opt_parameters=np.asarray(PN.OptParas),
                initial_physical_parameters=np.asarray(PN.PhyParas),
            )
        minima, W_PN, chiA, chiB = Align(PN, hyb, pn_evolution_order=pn_evolution_order, pn_waveform_mode_order=pn_waveform_mode_order)
        logR_delta = np.append(minima.x[0], minima.x[1:] + hyb.omega_mean*minima.x[0]/2)
        if len(PN.PhyParas) == 12:
            PN.PhyParas[8:] = logR_delta
        else:
            PN.PhyParas = np.append(PN.PhyParas, logR_delta)
        PN.Physical_to_Parameterize(hyb)

        scale = np.array([0.05,0.02,0.1,0.1,np.pi*2,np.pi*2,np.pi,np.pi/4,np.pi/hyb.omega_i/2.0,np.pi/4,np.pi/4,np.pi/4])
        lowbound12D = PN.OptParas - scale
        upbound12D = PN.OptParas + scale

        temp = time.time()
        minima12D = least_squares(
            Optimize12D,
            PN.OptParas,
            bounds=(lowbound12D, upbound12D),
            ftol=3e-15,
            xtol=3e-15,
            gtol=1e-8,
            x_scale='jac',
            args=(PN, hyb, pn_evolution_order, pn_waveform_mode_order),
        )
        print('Optimization time used:',time.time()-temp)
        minima12D.jac = approx_fprime(minima12D.x, Optimize12D, np.full_like(minima12D.x, 1.49e-8), PN, hyb, pn_evolution_order, pn_waveform_mode_order)
        if minima12D.success == False:
            print("12-D Optimization doesn't converge.")

        PN.OptParas = minima12D.x
        PN.Parameterize_to_Physical(hyb)
        if provenance_recorder:
            provenance_recorder.record(
                "complete_12d_optimization", outer_iteration=int(hyb.PNIter),
                optimizer=optimizer_summary(minima12D), optimized_parameters=np.asarray(PN.OptParas),
                physical_parameters=np.asarray(PN.PhyParas),
            )


    # Get aligned NR and PN waveforms
    hyb.t_PNStart, hyb.t_PNEnd = -90000, 1000 - hyb.t_start
    if cce_dir2 != None:
        minima, W_PN, chiA, chiB = Align(PN, hyb, W_NR2, pn_evolution_order=pn_evolution_order, pn_waveform_mode_order=pn_waveform_mode_order)
    else:
        minima, W_PN, chiA, chiB = Align(PN, hyb, pn_evolution_order=pn_evolution_order, pn_waveform_mode_order=pn_waveform_mode_order)

    t_delta = minima.x[0]
    logR_delta = np.append(0.0, minima.x[1:] + hyb.omega_mean*minima.x[0]/2)
    if len(PN.PhyParas) == 12:
        PN.PhyParas[8] = t_delta
        PN.PhyParas[9:] = logR_delta[1:]
    else:
        PN.PhyParas = np.append(PN.PhyParas, logR_delta)
    PN.Physical_to_Parameterize(hyb)
    R_delta = np.exp(quaternion.from_float_array(logR_delta))

    pn_before_time_shift = waveform_fingerprint(W_PN) if provenance_recorder else None
    W_PN.t = W_PN.t - t_delta
    nr_before_rotation = waveform_fingerprint(W_NR) if provenance_recorder else None
    W_NR = scri.rotate_physical_system(W_NR, R_delta)
    chiA = R_delta*chiA*R_delta.conjugate()
    chiB = R_delta*chiB*R_delta.conjugate()

    if provenance_recorder:
        provenance_recorder.record(
            "apply_final_residual_alignment", outer_iteration=int(hyb.PNIter),
            align_result=optimizer_summary(minima), pn_time_shift=float(-t_delta),
            pn_time_action="W_PN.t = W_PN.t - t_delta; NR time is unchanged",
            residual_log_rotation_wxyz=quaternion.as_float_array(quaternion.from_float_array(logR_delta)),
            residual_rotation_wxyz=quaternion.as_float_array(R_delta),
            pn_before=pn_before_time_shift, pn_after=waveform_fingerprint(W_PN),
            nr_before=nr_before_rotation, nr_after=waveform_fingerprint(W_NR),
        )

    print("SquaredError over matching window: ",SquaredError(W_NR, W_PN, hyb.t_start, hyb.t_start + hyb.length))
    if cce_dir2 == None:
        print(minima12D)


    # Stitch PN and NR waveforms
    W_H = Stitch(W_PN, W_NR, hyb)


    # Plot results
    if debug:
        t1=-80000
        t2=0
        fig, (ax1, ax2, ax3) = plt.subplots(3,1)
        ax1.plot(W_NR.t, W_NR.data[:,4].real-W_NR.data[:,4].imag, label='NR', linewidth=1)
        ax1.plot(W_PN.t, W_PN.data[:,4].real-W_PN.data[:,4].imag, label='PN', linewidth=1)
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
        fig.savefig( "hybridCheckResults")
        fig.clf()

    if cce_dir2 == None:
        Output(out_name, W_NR, W_PN, W_H, trans, minima12D, PN, hyb, nOrbits)
        if provenance_recorder:
            provenance_recorder.record(
                "serialize_iteration_outputs", outer_iteration=int(hyb.PNIter),
                output_npz=os.path.realpath(out_name), nr_waveform_file=os.path.realpath("hybridNR.h5"),
                pn_waveform_file=os.path.realpath("hybridPN.h5"), hybrid_waveform_file=os.path.realpath("hybrid.h5"),
                nr_waveform=waveform_fingerprint(W_NR), pn_waveform=waveform_fingerprint(W_PN),
                hybrid_waveform=waveform_fingerprint(W_H), physical_parameters=np.asarray(PN.PhyParas),
                optimized_parameters=np.asarray(PN.OptParas),
            )
        return W_NR, W_PN, W_H, minima12D, trans
    else:
        return W_NR, W_PN


# Run the code

def build_arg_parser():
    parser = argparse.ArgumentParser()
    parser.add_argument('--WaveformType', default='cce', help='cce for CCE waveform, and extrapolated for extrapolated waveform')
    parser.add_argument('--t', type=float, default=-7000.0, help='End time of matching window')
    parser.add_argument('--SimDir', default='SimData/Lev3', help='Path in which to find the extropolated waveform data')
    parser.add_argument('--CCEDir', default='CCEData/Lev3/CCE', help='Path in which to find the CCE waveform data')
    parser.add_argument('--OutName', default='Output.npz', help='Path in which to output results')
    parser.add_argument('--length', type=float, default=5000.0, help='Length of matching window')
    parser.add_argument('--nOrbits', type=float, default=None, help='Length of matching window in orbits, will disable "length" option if not None')
    parser.add_argument('--truncate', nargs=2, type=float, default=None, help='--truncate t1 t2. If specified, it will truncate the abd object and keep only data between t1 and t2')
    parser.add_argument('--CCEDir2', default=None, help='Path in which to find the second CCE waveform data')
    parser.add_argument('--ABDFormat', default='RPDMB', help='Format label passed to scri.SpEC.create_abd_from_h5')
    parser.add_argument('--BMSIterations', type=int, default=4, help='Number of sequential PN-BMS refinement iterations')
    parser.add_argument('--BMSAlignNProcs', type=int, default=None, help='Worker count for the BMS 2D brute-force alignment; defaults to SLURM_CPUS_PER_TASK')
    parser.add_argument('--BMSBruteForceDt', type=int, default=1000, help='Number of time-shift samples in the BMS 2D brute-force alignment')
    parser.add_argument('--BMSBruteForcePhi', type=int, default=20, help='Number of phase samples in the BMS 2D brute-force alignment')
    parser.add_argument('--PNOrder', type=float, default=6.0, help='PN truncation order used for both orbital evolution and waveform modes')
    parser.add_argument('--PNEvolutionOrder', type=float, default=None, help='Override evolution order, up to 22 with Schwarzschild supplementation')
    parser.add_argument('--PNWaveformModeOrder', type=float, default=None, help='Override waveform-mode order independently (0 through 6)')
    parser.add_argument('--RecordBMSProvenance', action='store_true', help='Enable default-off diagnostic BMS provenance recording')
    parser.add_argument('--BMSProvenanceDir', default=None, help='Diagnostic output directory used only with --RecordBMSProvenance')
    parser.add_argument('--BMSProvenanceLabel', default='bms_framefix_diagnostic', help='Run label stored in the diagnostic provenance ledger')
    parser.add_argument('--FixedParameterConfig', default=None, help='Optional JSON config that fixes intrinsic parameters and disables 12-D optimization')
    return parser


def run(args):
    clock_start = time.time()
    WaveformType = args.WaveformType
    t_end = args.t
    data_dir = args.SimDir
    cce_dir = args.CCEDir
    cce_dir2 = args.CCEDir2
    out_name = args.OutName
    length = np.array(args.length)
    nOrbits = args.nOrbits
    truncate = args.truncate
    abd_format = args.ABDFormat
    bms_iterations = args.BMSIterations
    bms_align_nprocs = args.BMSAlignNProcs
    bms_bruteforce_dt = args.BMSBruteForceDt
    bms_bruteforce_phi = args.BMSBruteForcePhi
    pn_order = args.PNOrder
    pn_evolution_order = getattr(args, "PNEvolutionOrder", None)
    pn_waveform_mode_order = getattr(args, "PNWaveformModeOrder", None)
    provenance_recorder = None
    fixed_parameter_config = None
    fixed_intrinsic_parameters = None
    fixed_omega_i = None
    OptArg = 1
    if args.FixedParameterConfig is not None:
        fixed_parameter_path = os.path.realpath(args.FixedParameterConfig)
        with open(fixed_parameter_path, "rb") as stream:
            fixed_parameter_bytes = stream.read()
        fixed_parameter_config = json.loads(fixed_parameter_bytes)
        fixed_intrinsic_parameters = np.asarray([
            fixed_parameter_config["q"], fixed_parameter_config["M"],
            *fixed_parameter_config["chi1"], *fixed_parameter_config["chi2"],
        ], dtype=float)
        fixed_omega_i = float(fixed_parameter_config["initial_orbital_frequency"])
        OptArg = 0
    if args.RecordBMSProvenance:
        if args.BMSProvenanceDir is None:
            raise ValueError("--BMSProvenanceDir is required with --RecordBMSProvenance")
        provenance_recorder = BMSProvenanceRecorder(args.BMSProvenanceDir, args.BMSProvenanceLabel)
        provenance_recorder.record(
            "run_configuration", command_line_arguments=vars(args), working_directory=os.getcwd(),
            source_file=os.path.realpath(__file__),
        )
        if fixed_parameter_config is not None:
            provenance_recorder.record(
                "fixed_parameter_policy",
                fixed_parameter_file=fixed_parameter_path,
                fixed_parameter_sha256=hashlib.sha256(fixed_parameter_bytes).hexdigest(),
                parameter_status=fixed_parameter_config["parameter_status"],
                fixed_intrinsic_parameters=fixed_intrinsic_parameters,
                fixed_omega_i=fixed_omega_i,
                fixed_frame_quaternion_wxyz=fixed_parameter_config["frame_quaternion_wxyz"],
                allowed_variable_policy=fixed_parameter_config["allowed_variable_policy"],
            )
    maxiter = 30

    if WaveformType == 'cce':
        abd, t0 = get_abd(cce_dir, truncate, abd_format, provenance_recorder, "initial_length_probe")
        if nOrbits != None:
            length = get_length_from_abd(abd, nOrbits, t_end)
    else:
        W_NR, t0, length = get_extrapolated_NR(data_dir, nOrbits, t_end, length)
        maxiter = 0
    if cce_dir2 != None:
        maxiter = 0

    hyb = hyb_quantites(t_end, length)
    PN = PNParameters(data_dir, hyb, t0)
    if fixed_parameter_config is not None:
        PN.frame_i = np.asarray(fixed_parameter_config["frame_quaternion_wxyz"], dtype=float)
        PN.PhyParas = fixed_intrinsic_parameters.copy()
    if os.path.exists(out_name):
        Checkpoint = np.load(out_name, allow_pickle=True)
        hyb.PNIter = Checkpoint['checkpoint'][0]
        hyb.cost = Checkpoint['checkpoint'][1]
        hyb.omega_i = Checkpoint['checkpoint'][2]
        hyb.t_start = Checkpoint['checkpoint'][3]
        hyb.length = Checkpoint['checkpoint'][4]
        hyb.t_PNStart = Checkpoint['checkpoint'][5]
        hyb.t_PNEnd = Checkpoint['checkpoint'][6]
        PN.frame_i = Checkpoint['checkpoint'][7]
        PN.OptParas = Checkpoint['OptParas']
        PN.PhyParas = Checkpoint['PhyParas']

    if hyb.PNIter>maxiter:
        hyb.PNIter = 10
    while hyb.PNIter<=maxiter:
        print("PNIter=: ", hyb.PNIter)
        if provenance_recorder:
            provenance_recorder.record(
                "begin_outer_iteration", outer_iteration=int(hyb.PNIter),
                physical_parameters=np.asarray(PN.PhyParas), optimized_parameters=np.asarray(PN.OptParas),
                matching_window=[float(hyb.t_start), float(hyb.t_start + hyb.length)],
            )
        if fixed_intrinsic_parameters is not None:
            if len(PN.PhyParas) == 12:
                PN.PhyParas[:8] = fixed_intrinsic_parameters
            else:
                PN.PhyParas = fixed_intrinsic_parameters.copy()
            PN.frame_i = np.asarray(fixed_parameter_config["frame_quaternion_wxyz"], dtype=float)
            hyb.omega_i = fixed_omega_i
        W_NR, W_PN, W_H, minima12D, BMSTrans = Hybridize(
            WaveformType,
            t_end,
            data_dir,
            cce_dir,
            out_name,
            length,
            nOrbits,
            hyb,
            PN,
            debug=0,
            OptimizePNParas=OptArg,
            truncate=truncate,
            cce_dir2=cce_dir2,
            abd_format=abd_format,
            bms_iterations=bms_iterations,
            bms_align_nprocs=bms_align_nprocs,
            bms_bruteforce_dt=bms_bruteforce_dt,
            bms_bruteforce_phi=bms_bruteforce_phi,
            pn_order=pn_order,
            pn_evolution_order=pn_evolution_order,
            pn_waveform_mode_order=pn_waveform_mode_order,
            provenance_recorder=provenance_recorder,
            fixed_intrinsic_parameters=fixed_intrinsic_parameters,
            fixed_omega_i=fixed_omega_i,
        )

        if hyb.PNIter >= 1 and abs(hyb.cost[-1]-hyb.cost[-2])/hyb.cost[-1]<1e-2:
            hyb.PNIter = maxiter + 1
        else:
            hyb.PNIter += 1

    # Output results
    Output(out_name, W_NR, W_PN, W_H, BMSTrans, minima12D, PN, hyb, nOrbits)
    if provenance_recorder:
        provenance_recorder.record(
            "run_complete", final_nr_waveform=waveform_fingerprint(W_NR),
            final_pn_waveform=waveform_fingerprint(W_PN), final_hybrid_waveform=waveform_fingerprint(W_H),
            final_physical_parameters=np.asarray(PN.PhyParas), final_optimized_parameters=np.asarray(PN.OptParas),
        )
    print("All done, total time:", time.time()-clock_start)


def main(argv=None):
    parser = build_arg_parser()
    run(parser.parse_args(argv))


if __name__ == '__main__':
    main()
