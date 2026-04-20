#%%

# ----------------------------------------------------------------------------------------------------------- #
#                                          IMPORT MODULES                                                     #
# ----------------------------------------------------------------------------------------------------------- #

import sys
sys.path.append('/home/fra/Politecnico Di Torino Studenti Dropbox/Francesco Bellelli/12. LINING/03. Codes/')

import shutil
from pathlib import Path

import io
import re
import h5py
import numpy as np
import pyvista as pv
from PIL import Image
from cycler import cycler
import matplotlib.pyplot as plt
from scipy.io import loadmat, savemat
from fra_toolkit.ade_tdibc import tdim_bc
from scipy.interpolate import interp1d, griddata
from fra_toolkit.spectral_process import get_single_sided_fft
from fra_toolkit.impedance_vector_fit import multipole_to_impedance, impedance_to_multipole
from fra_toolkit.read_external_data import read_impedance_lbm, read_facesheet, remove_zero_signal_part
from fra_toolkit.multi_point_method import mpm_impedance, compute_impedance_multipt, compute_impedance_pairwise, pressure_fft, complex_pressure_amplitude

# ----------------------------------------------------------------------------------------------------------- #
#                                           PLOT SETTINGS                                                     #
# ----------------------------------------------------------------------------------------------------------- #

# Line width, marker size
plt.rcParams['lines.linewidth'] = 2
plt.rcParams['lines.markersize'] = 10

# Font settings
plt.rcParams['font.size'] = 24
plt.rcParams['font.family'] = 'serif'
plt.rcParams['text.usetex'] = True

# Axes and grid
plt.rcParams['axes.labelsize'] = 24
plt.rcParams['axes.titlesize'] = 24
plt.rcParams['axes.grid'] = False

# Legend
plt.rcParams['legend.fontsize'] = 24
plt.rcParams['legend.frameon'] = False

# Figure size (optional default)
plt.rcParams['figure.figsize'] = (12, 10)

# Color cycle
cols = ['b', 'r', 'g', 'm', 'c', 'y', 'k']
plt.rcParams['axes.prop_cycle'] = cycler(color=cols, linestyle=['-', '--', '-.', ':', '-', '--', '-.'])

# ----------------------------------------------------------------------------------------------------------- #
#                                           CALCULATIONS                                                      #
# ----------------------------------------------------------------------------------------------------------- #

def read_probe(f_exc, p_a, x_target):

    dx = 3 * nu * cs / (c0 * (tau - 0.5))
    dt = dx**2 * (tau - 0.5) / (3 * nu)

    it = 0
    p = []
    rho = []
    u = []

    for it in range(0, 630000, save_iter):
        try:
            mesh = pv.read(f'/media/fra/Volume/liner_test/sol_new2/sol/fields_{f_exc}Hz_{p_a}Pa_{it:05d}.vtk')

            dx = mesh.spacing[0]
            Lx = mesh.dimensions[0]
            x = np.arange(0,dx*Lx,dx)

            idx_x = np.argmin(np.abs(x - x_target))

            p.append(mesh['density'][idx_x] * (287.05 * T) - 101325)
            rho.append(mesh['density'][idx_x])
            u.append(mesh['velocity'][idx_x, 0])

            it += 1
        except:
            break

    t = dt*np.arange(0, it, save_iter)

    return t, np.array(p), np.array(u), np.array(rho)


def write_file_probe(t, p, u, rho, x, i):

    # Sample data: A list of lists or tuples containing your numerical values
    p += 101325*np.ones_like(p)
    data = np.column_stack((t, p, u, rho))

    file_path = f'probe_{i}.txt'

    with open(file_path, "w") as f:
        # 1. Write the metadata headers
        f.write("# Time history output for all variables, bands and points\n")
        f.write(f"# Number of frames = {len(data)}\n")
        f.write("# Number of variables = 3\n")
        f.write("# Number of bands = 1\n")
        f.write("# Number of points = 1\n")
        f.write("# Description of columns: \n")
        f.write("# Column, Point_#, (X (m),Y (m),Z (m)), Variable_#, Variable name, Band_#, Freq_range\n")
        f.write("# 1, time(sec)\n")
        f.write(f"# 2, 0, ({x},7.87443e-05,7.87443e-05), 0, ExtractSignal Pressure, 0, Single frequency band\n")
        f.write(f"# 3, 0, ({x},7.87443e-05,7.87443e-05), 1, ExtractSignal X-Velocity, 0, Single frequency band\n")
        f.write(f"# 4, 0, ({x},7.87443e-05,7.87443e-05), 2, ExtractSignal Density, 0, Single frequency band\n")

        # 2. Write the numerical data
        for row in data:
            # Convert each number to a string and join with a space
            line = " ".join(map(str, row))
            f.write(line + "\n")

    return


def read_impedance_optydb(SPL):
    # Define the array of frequencies that have been simulated
    f = np.array([800, 1000, 1400, 2000])

    # Initialize resistance (R) and reactance (X)
    R = np.zeros(len(f))
    X = np.zeros(len(f))

    # Read R and X at each probe
    for j in range(len(f)):
        for i in range(1,12):
            data = np.genfromtxt(f'/run/user/1000/gvfs/sftp:host=hpc-legionlogin.polito.it,user=fbellelli/home/fbellelli/test_kundt/f{f[j]}_SPL{SPL}/probe_{i}_Z.txt',skip_header=1)

            R[j] += data[1]
            X[j] += data[2]

    # Average over probes
    R /= 11
    X /= 11

    # Compute impedance
    Zs0 = R + 1j*X

    return Zs0


def read_impedance_single(SPL, f):
    # Initialize R and X
    R = 0
    X = 0

    # Read R and X at each probe
    for i in range(1,12):
        data = np.genfromtxt(f'/media/fra/Expansion/NIT/Sharp_SingleT_{f}Hz_{SPL}dB/probe_{i}_Z.txt',skip_header=1)

        R += data[1]
        X += data[2]

    # Average over probes
    R /= 11
    X /= 11

    # Compute impedance
    Zs0 = R + 1j*X

    return Zs0

# Parameters
c0 = 343.283
rho = 1.20376
nu = 1.5067e-5
tau = 0.5005
cs = 1 / np.sqrt(3)
save_iter = 500

# Derivate quantities
T = 101325 / (287.05 * rho)
mu = nu * rho

# Define frequency grid and compute angular frequency
f_hz = np.logspace(np.log10(500), np.log10(2500), 100)
omega = 2 * np.pi * f_hz

# Define SPL, acoustic pressure, and excitation frequency arrays
SPL = [130, 145]
f_exc = [800, 1000, 1400, 2000]
p_a = [89, 503]

#%% Compute impedance with multi point method (MPM)
fft_flag = False
k_ref = np.array([[18.92, 18.74, 19, 21], [19.4, 19.6, 19, 21]])

for i in range(len(SPL)):
    for j in range(len(f_exc)):
        print(f'SPL = {SPL[i]} dB - f = {f_exc[j]} Hz...')

        # Compute wavelength
        lam = c0 / f_exc[j]

        if (SPL[i] == 130) & (f_exc[j] == 2000):
            k_ref2 = 0.53
        elif (SPL[i] == 145) & (f_exc[j] == 800):
            k_ref2 = 0.306
        else:
            k_ref2 = 0.25

        x_ref = k_ref[i, j] * lam
        x = x_ref - lam * np.array([0, k_ref2])

        # TDIM BC simulations (quasi 1D)
        t, p1, _, _ = read_probe(f_exc[j], p_a[i], x[0])
        _, p2, _, _ = read_probe(f_exc[j], p_a[i], x[1])

        dt = t[1] - t[0]

        p_time = np.array([p1, p2])

        Z = compute_impedance_multipt(
            p_time,
            x,
            dt,
            f_exc[j],
            c0,
        )

        plt.figure(i)
        plt.subplot(211)
        plt.plot(f_exc[j],np.real(Z),'or')

        plt.subplot(212)
        plt.plot(f_exc[j],np.imag(Z),'or')

    Zs0 = read_impedance_lbm(SPL[i], f_hz)

    plt.figure(i)
    plt.subplot(211)
    plt.plot(f_hz,np.real(Zs0),'-k')
    plt.ylabel('$R$ [-]')

    plt.subplot(212)
    plt.plot(f_hz,np.imag(Zs0),'-k')
    plt.ylabel('$X$ [-]')
    plt.xlabel('$f$ [Hz]')
    plt.suptitle(f'SPL = {SPL[i]} dB')
    plt.tight_layout()
    #plt.savefig(f'plots/impedance_{SPL[i]}dB.png', dpi=300, bbox_inches='tight')

    plt.show()

#%% Input/Output validation on the cavity surface
for i in range(len(SPL)):
    for j in range(len(f_exc)):

        # Read simulated TDIBC data
        data = np.genfromtxt(f'sol_new/tdibc_{f_exc[j]}Hz_{p_a[i]}Pa.txt', invalid_raise=False)
        t = data[:, 0]
        p = data[:, 1]
        v = data[:, 2]

        f_p = interp1d(t, p)
        f_v = interp1d(t, v)

        t = np.linspace(t[0], t[-1], len(t))
        p = f_p(t)
        v = f_v(t)

        # Remove the first acoustic cycle
        mask = (t*f_exc[j] >= 10) & (t*f_exc[j] <= 20)

        t = t[mask]
        t -= t[0]

        p = p[mask]
        v = v[mask]
        
        # Flip the sign except for 1000 Hz to match reference data phase (different initialization)
        if f_exc[j] != 1000:
            p = -p
            v = -v

        p -= np.mean(p)
        v -= np.mean(v)

        # Read reference data
        t_ref, v_ref, p_ref = read_facesheet(f'/media/fra/Volume/liner_test/01_a_posteriori_validation_NIT/nit_{f_exc[j]}Hz_{SPL[i]}dB', remove_zero_part=True)

        p_ref -= np.mean(p_ref)
        v_ref -= np.mean(v_ref)

        # Compute FFT of the time series
        f_ref, p_fft_ref = get_single_sided_fft(p_ref, t_ref)
        _, v_fft_ref = get_single_sided_fft(v_ref, t_ref)

        f, p_fft = get_single_sided_fft(p, t)
        _, v_fft = get_single_sided_fft(v, t)

        # Plot results
        plt.figure()
        plt.subplot(221)
        plt.plot(t*f_exc[j],p)
        try:
            plt.plot(t_ref*f_exc[j],p_ref)
        except:
            print()
        plt.ylabel('$p^\prime$ [Pa]')

        plt.subplot(223)
        plt.plot(t*f_exc[j],v)
        try:
            plt.plot(t_ref*f_exc[j],v_ref)
        except:
            print()
        plt.xlabel('$t/t_{cycle}$ [-]')
        plt.ylabel('$v^\prime$ [m/s]')
        plt.suptitle(f'SPL = {SPL[i]} dB - f = {f_exc[j]} Hz')

        plt.subplot(222)
        plt.plot(f/f_exc[j],20*np.log10(np.abs(p_fft)/2e-5))
        try:
            plt.plot(f_ref/f_exc[j],20*np.log10(np.abs(p_fft_ref)/2e-5))
        except:
            print()
        plt.ylabel('$\\tilde{p}^\prime$ [dB]')
        plt.xlim(0.5, 1.5)
        plt.ylim(0.8*SPL[i], 1.1*SPL[i])

        plt.subplot(224)
        plt.plot(f/f_exc[j],np.abs(v_fft))
        try:
            plt.plot(f_ref/f_exc[j],np.abs(v_fft_ref))
        except:
            print()
        plt.xlabel('$f/f_{exc}$ [-]')
        plt.ylabel('$\\tilde{v}^\prime$ [m/s]')
        plt.xlim(0.5, 1.5)
        plt.figlegend(['TDIBC','LBM'], loc='upper center', ncol=2, bbox_to_anchor=(0.55, 0.965))
        plt.tight_layout()
        # plt.savefig(f'plots/pressure_and_velocity_{f_exc[j]}Hz_{SPL[i]}dB.png', dpi=300, bbox_inches='tight')

        plt.show()

#%% Compare impedance computed with optydb_kundt

# Define the array of frequencies that have been simulated
f = np.array([800, 1000, 1400, 2000])

# Initialize resistance (R) and reactance (X)
R = np.zeros(len(f))
X = np.zeros(len(f))

# Read R and X at each probe
for j in range(len(f)):
    for i in range(1,12):
        data = np.genfromtxt(f'/run/user/1000/gvfs/sftp:host=hpc-legionlogin.polito.it,user=fbellelli/home/fbellelli/test_kundt/test5/f{f[j]}_SPL{130}/probe_{i}_Z.txt',skip_header=1)

        R[j] += data[1]
        X[j] += data[2]

# Average over probes
R /= 11
X /= 11

plt.figure(1)
plt.subplot(211)
plt.plot(f,R,'om')

plt.subplot(212)
plt.plot(f,R,'om')

# Initialize resistance (R) and reactance (X)
R = np.zeros(len(f))
X = np.zeros(len(f))

# Read R and X at each probe
for j in range(len(f)):
    for i in range(1,12):
        data = np.genfromtxt(f'/run/user/1000/gvfs/sftp:host=hpc-legionlogin.polito.it,user=fbellelli/home/fbellelli/test_kundt/test6/f{f[j]}_SPL{130}/probe_{i}_Z.txt',skip_header=1)

        R[j] += data[1]
        X[j] += data[2]

# Average over probes
R /= 11
X /= 11

plt.figure(1)
plt.subplot(211)
plt.plot(f,R,'xg')

plt.subplot(212)
plt.plot(f,X,'xg')

plt.show()

#%% Write probes for optydb_kundt impedance calculation
for i in range(len(SPL)):
    for j in range(len(f_exc)):
        print(f'SPL = {SPL[i]} dB - f = {f_exc[j]} Hz...')

        # Compute wavelength
        lam = c0 / f_exc[j]

        x_liner = 13 * lam
        x_ref = 11 * lam
        x = x_ref - lam * np.linspace(0, 0.6, 11)

        # TDIM BC simulations (quasi 1D)
        t, p1, u1, rho1 = read_probe(f_exc[j], p_a[i], x[0])
        _, p2, u2, rho2 = read_probe(f_exc[j], p_a[i], x[1])
        _, p3, u3, rho3 = read_probe(f_exc[j], p_a[i], x[2])
        _, p4, u4, rho4 = read_probe(f_exc[j], p_a[i], x[3])
        _, p5, u5, rho5 = read_probe(f_exc[j], p_a[i], x[4])
        _, p6, u6, rho6 = read_probe(f_exc[j], p_a[i], x[5])
        _, p7, u7, rho7 = read_probe(f_exc[j], p_a[i], x[6])
        _, p8, u8, rho8 = read_probe(f_exc[j], p_a[i], x[7])
        _, p9, u9, rho9 = read_probe(f_exc[j], p_a[i], x[8])
        _, p10, u10, rho10 = read_probe(f_exc[j], p_a[i], x[9])
        _, p11, u11, rho11 = read_probe(f_exc[j], p_a[i], x[10])

        write_file_probe(t, p1, u1, rho1, 0, 1)
        write_file_probe(t, p2, u2, rho2, 0, 2)
        write_file_probe(t, p3, u3, rho3, 0, 3)
        write_file_probe(t, p4, u4, rho4, 0, 4)
        write_file_probe(t, p5, u5, rho5, 0, 5)
        write_file_probe(t, p6, u6, rho6, 0, 6)
        write_file_probe(t, p7, u7, rho7, 0, 7)
        write_file_probe(t, p8, u8, rho8, 0, 8)
        write_file_probe(t, p9, u9, rho9, 0, 9)
        write_file_probe(t, p10, u10, rho10, 0, 10)
        write_file_probe(t, p11, u11, rho11, 0, 11)

        source_dir = Path('.')
        target_dir = source_dir / f'probes_f_{f_exc[j]}Hz_SPL_{SPL[i]}dB'

        target_dir.mkdir(exist_ok=True)

        for file_path in source_dir.glob('probe_*.txt'):
            shutil.move(str(file_path), str(target_dir / file_path.name))

#%% Find optimal probe location for resistance calculation
SPL = [145]
f_exc = [800]
p_a = [89 if SPL[0] == 130 else 503]

for i in range(len(SPL)):
    for j in range(len(f_exc)):
        print(f'SPL = {SPL[i]} dB - f = {f_exc[j]} Hz...')

        # Compute wavelength
        lam = c0 / f_exc[j]

        for k_ref in np.arange(0.3,0.31,0.002):
            print(f'k_ref = {k_ref}...')

            x_ref = 19.4 * lam
            x = x_ref - lam * np.array([0, k_ref])

            # TDIM BC simulations (quasi 1D)
            t, p1, _, _ = read_probe(f_exc[j], p_a[i], x[0])
            _, p2, _, _ = read_probe(f_exc[j], p_a[i], x[1])

            dt = t[1] - t[0]

            p_time = np.array([p1, p2])

            Z = compute_impedance_multipt(
                p_time,
                x,
                dt,
                f_exc[j],
                c0,
            )

            Zs_single = read_impedance_single(SPL[i], f_exc[j])

            plt.figure(i)
            plt.subplot(211)
            plt.plot(k_ref,np.real(Z)/np.real(Zs_single) - 1,'ok')

            plt.subplot(212)
            plt.plot(k_ref,np.imag(Z)/np.imag(Zs_single) - 1,'ok')
        
        plt.show()