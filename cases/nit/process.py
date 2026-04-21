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

def read_probe(f_exc, p_a, x_target, foldername='/media/fra/Volume/liner_test/02_tdibc_validation_NIT/uniform'):

    dx = 3 * nu * cs / (c0 * (tau - 0.5))
    dt = dx**2 * (tau - 0.5) / (3 * nu)

    it = 0
    p = []
    rho = []
    u = []

    for it in range(0, 630000, save_iter):
        try:
            mesh = pv.read(f'{foldername}/fields_{f_exc}Hz_{p_a}Pa_{it:05d}.vtk')

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


def read_tdibc(foldername,f_exc, p_a, start_cycles=1, ac_cycles=10):

    # Read simulated TDIBC data
    data = np.genfromtxt(f'{foldername}/tdibc_{f_exc}Hz_{p_a}Pa.txt', invalid_raise=False)
    t = data[:, 0]
    p = data[:, 1]
    v = data[:, 2]

    # Interpolate on uniform time grid
    f_p = interp1d(t, p)
    f_v = interp1d(t, v)

    t = np.linspace(t[0], t[-1], len(t))
    p = f_p(t)
    v = f_v(t)

    # Remove the starting acoustic cycles
    mask = (t*f_exc >= start_cycles) & (t*f_exc <= start_cycles + ac_cycles)

    t = t[mask]
    t -= t[0]

    p = p[mask]
    v = v[mask]
    
    # Flip the sign except for 1000 Hz to match reference data phase (different initialization)
    if f_exc != 1000:
        p = -p
        v = -v

    p -= np.mean(p)
    v -= np.mean(v)

    return t, p, v


def read_pf_sim(SPL, f_exc, ac_cycles=10, start_cycles=5.25):
    # Load data from HDF5/MAT file
    with h5py.File(f'/media/fra/Volume/liner_test/01_a_posteriori_validation_NIT/nit_{f_exc}Hz_{SPL}dB_duct.mat', 'r') as f:
        v = np.array(f['fields/x_velocity'], dtype=float)
        p = np.array(f['fields/static_pressure'], dtype=float)
        cc = np.array(f['centroids_m'], dtype=float)

    # Define time array
    dt = 2.111e-05
    t = dt*np.arange(np.shape(p)[-1])

    # Define spatial arrays
    x = cc[:,0]
    y = cc[:,1]

    # Remove x < 0
    mask = x < 0
    x = x[mask]
    y = y[mask]
    p = p[mask, :]
    v = v[mask, :]

    # Remove the starting acoustic cycles
    mask = (t*f_exc >= start_cycles) & (t*f_exc <= start_cycles + ac_cycles)

    t = t[mask]
    t -= t[0]

    p = p[:,mask]
    v = v[:,mask]

    # Remove mean
    v -= np.mean(v, axis=1, keepdims=True)
    p -= np.mean(p, axis=1, keepdims=True)

    # Compute phase average
    len_cycle = len(t[t*f_exc <= 1])
    cycles = int(len(t) / len_cycle)
    p = p[:, 0:cycles * len_cycle]
    v = v[:, 0:cycles * len_cycle]
    t = t[0:len_cycle]

    p = np.mean(p.reshape(len(x), cycles, len_cycle), axis=1)
    v = np.mean(v.reshape(len(x), cycles, len_cycle), axis=1)

    # Resample on uniform 2D grid
    dx = 5e-5
    array_x = np.arange(np.min(x), np.max(x), dx)
    array_y = np.arange(np.min(y), np.max(y), dx)
    grid_x, grid_y = np.mgrid[min(x):max(x):dx, min(y):max(y):dx]

    p_grid = np.zeros((len(array_x), len(array_y), len(t)))
    v_grid = np.zeros((len(array_x), len(array_y), len(t)))
    for i in range(len(t)):
        p_grid[:, :, i] = griddata((x, y), p[:, i], (grid_x, grid_y), method='linear', fill_value=0.0)
        v_grid[:, :, i] = griddata((x, y), v[:, i], (grid_x, grid_y), method='linear', fill_value=0.0)

    return t, p_grid, v_grid, array_x, array_y


def read_lbc_sim(foldername, p_a, f_exc, ac_cycles=10, start_cycles=1, damp_cycles=2, x_cut=0.04):
    # Pre compute dx, dt, and number of time steps
    dx = 3 * nu * cs / (c0 * (tau - 0.5))
    dt = dx**2 * (tau - 0.5) / (3 * nu)
    Nt = int(np.ceil((ac_cycles + start_cycles + damp_cycles) / (f_exc*dt)))

    # Remove start and damping acoustic cycles
    Nt_start = int(np.ceil((start_cycles / f_exc) / dt))
    Nt_start -= np.mod(Nt_start, save_iter)

    Nt_end = int(np.ceil(((start_cycles + ac_cycles) / f_exc) / dt))
    Nt_end -= np.mod(Nt_end, save_iter)

    # Define time array
    Nt = len(range(Nt_start, Nt_end, save_iter))
    t = dt * save_iter * np.arange(0, Nt)

    # Define spatial arrays
    mesh = pv.read(f'{foldername}/fields_{f_exc}Hz_{p_a}Pa_00000.vtk')

    dx = mesh.spacing[0]
    Lx = mesh.dimensions[0]
    Ly = mesh.dimensions[1]

    x = np.arange(0,dx*Lx,dx)
    y = np.arange(0,dx*Ly,dx)

    # Cut along x
    mask = (np.max(x) - x  < x_cut)
    x = x[mask]

    x -= np.max(x)
    y -= np.mean(y)

    Lx = len(x)

    # Read simulated TDIBC data from vtk
    ii = 0
    p = np.zeros((Ly, Lx, Nt))
    v = np.zeros((Ly, Lx, Nt))

    for it in range(Nt_start, Nt_end, save_iter):
        mesh = pv.read(f'{foldername}/fields_{f_exc}Hz_{p_a}Pa_{it:05d}.vtk')

        p_temp = (mesh['density'] * (287.05 * T)).reshape(mesh.dimensions[1], mesh.dimensions[0])
        v_temp = (mesh['velocity'][:, 0]).reshape(mesh.dimensions[1], mesh.dimensions[0])

        p[:, :, ii] = p_temp[:, mask]
        v[:, :, ii] = v_temp[:, mask]

        ii += 1

    # Remove mean
    v -= np.mean(v, axis=2, keepdims=True)
    p -= np.mean(p, axis=2, keepdims=True)

    # Compute phase average
    len_cycle = len(t[t*f_exc <= 1])
    cycles = int(len(t) / len_cycle)
    p = p[:, :, 0:cycles * len_cycle]
    v = v[:, :, 0:cycles * len_cycle]
    t = t[0:len_cycle]

    p = np.mean(p.reshape(Ly, Lx, cycles, len_cycle), axis=2)
    v = np.mean(v.reshape(Ly, Lx, cycles, len_cycle), axis=2)

    return t, p, v, x, y


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
ac_cycles_pf = [6, 7, 10, 14]
start_cycles_pf = 5.25 / f_exc[0] * np.array(f_exc, dtype=float)

#%%

# ----------------------------------------------------------------------------------------------------------- #
#                                               FLOW VIZ                                                      #
# ----------------------------------------------------------------------------------------------------------- #

for i in range(len(SPL)):
    for j in range(len(f_exc)):
        print(f'SPL = {SPL[i]} dB - f = {f_exc[j]} Hz...')

        print('Reading pf...')
        t_pf, p_pf, v_pf, x_pf, y_pf = read_pf_sim(SPL[i], f_exc[j], ac_cycles=ac_cycles_pf[j], start_cycles=start_cycles_pf[j])

        print('Reading lbc...')
        t_lc, p_lc, v_lc, x_lc, y_lc = read_lbc_sim('/media/fra/Volume/liner_test/02_tdibc_validation_NIT/orifices', p_a[i], f_exc[j], ac_cycles=10, start_cycles=1, damp_cycles=2, x_cut=0.04)

        POA = 8 * np.pi * (1.17 / 2)**2 / 9.906**2
        vscale = 9.906/(2 * 1.17)
        v_lc /= (vscale * POA)

        max_vel = np.round(np.max([np.max(np.abs(v_pf)), np.max(np.abs(v_lc))]))

        idx_inflow_pf = np.unravel_index(np.argmin(v_pf),v_pf.shape)[-1]
        idx_outflow_pf = np.unravel_index(np.argmax(v_pf),v_pf.shape)[-1]

        idx_inflow_lc = np.unravel_index(np.argmin(v_lc),v_lc.shape)[-1]
        idx_outflow_lc = np.unravel_index(np.argmax(v_lc),v_lc.shape)[-1]

        print('Plots...')
        plt.figure()
        plt.subplot(221)
        im = plt.pcolormesh(x_pf/1.17e-3,y_pf/1.17e-3,v_pf[:, :, idx_inflow_pf].T,vmin=-max_vel,vmax=max_vel,cmap='bwr')
        plt.title('Inflow - Ref LBM')
        plt.ylabel('$y/d$ [-]')
        plt.xlim(-10,0)
        plt.ylim(-4, 4)

        plt.subplot(222)
        plt.pcolormesh(x_lc/1.17e-3,y_lc/1.17e-3,v_lc[:, :, idx_inflow_lc],vmin=-max_vel,vmax=max_vel,cmap='bwr')
        plt.title('Inflow - TDIM BC')
        plt.xlim(-10,0)
        plt.ylim(-4, 4)

        plt.subplot(223)
        plt.pcolormesh(x_pf/1.17e-3,y_pf/1.17e-3,v_pf[:, :, idx_outflow_pf].T,vmin=-max_vel,vmax=max_vel,cmap='bwr')
        plt.title('Outflow - Ref LBM')
        plt.xlabel('$x/d$ [-]')
        plt.ylabel('$y/d$ [-]')
        plt.xlim(-10,0)
        plt.ylim(-4, 4)

        plt.subplot(224)
        plt.pcolormesh(x_lc/1.17e-3,y_lc/1.17e-3,v_lc[:, :, idx_outflow_lc],vmin=-max_vel,vmax=max_vel,cmap='bwr')
        plt.title('Outflow - TDIM BC')
        plt.xlabel('$x/d$ [-]')
        plt.xlim(-10,0)
        plt.ylim(-4, 4)

        plt.suptitle(f'SPL = {SPL[i]} dB - f = {f_exc[j]} Hz')
        plt.tight_layout()
        plt.colorbar(im, ax=plt.gcf().get_axes(), orientation='vertical', fraction=0.05, pad=0.05, label='$v^\prime$ [m/s]')
        #plt.savefig(f'plots/inflow_outflow_{f_exc[j]}Hz_{SPL[i]}dB_contour.png', dpi=300, bbox_inches='tight')

        plt.figure()
        plt.plot(y_pf/1.17e-3,v_pf[-1, :, idx_inflow_pf],'-b')
        plt.plot(y_pf/1.17e-3,v_pf[-1, :, idx_outflow_pf],'-r')
        plt.plot(y_lc/1.17e-3,v_lc[:, -1, idx_inflow_lc],'--b')
        plt.plot(y_lc/1.17e-3,v_lc[:, -1, idx_outflow_lc],'--r')
        plt.xlabel('$y/d$ [-]')
        plt.ylabel('$v^\prime$ [m/s]')
        plt.title(f'SPL = {SPL[i]} dB - f = {f_exc[j]} Hz')
        plt.xlim(-4, 4)
        plt.ylim(-max_vel-0.5, max_vel+0.5)
        plt.legend(['Inflow - Ref LBM','Outflow - Ref LBM','Inflow - TDIM BC','Outflow - TDIM BC'], ncols=1, loc='lower center')
        plt.tight_layout()
        #plt.savefig(f'plots/inflow_outflow_{f_exc[j]}Hz_{SPL[i]}dB_line.png', dpi=300, bbox_inches='tight')

        plt.show()

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
        print(f'SPL = {SPL[i]} dB - f = {f_exc[j]} Hz...')

        t_unif, p_unif, v_unif = read_tdibc('/media/fra/Volume/liner_test/02_tdibc_validation_NIT/uniform', f_exc[j], p_a[i])
        t_orif, p_orif, v_orif = read_tdibc('/media/fra/Volume/liner_test/02_tdibc_validation_NIT/orifices', f_exc[j], p_a[i])

        # Read reference data
        t_ref, v_ref, p_ref = read_facesheet(f'/media/fra/Volume/liner_test/01_a_posteriori_validation_NIT/nit_{f_exc[j]}Hz_{SPL[i]}dB', remove_zero_part=True)

        p_ref -= np.mean(p_ref)
        v_ref -= np.mean(v_ref)

        # Compute FFT of the time series
        f_ref, p_fft_ref = get_single_sided_fft(p_ref, t_ref)
        _, v_fft_ref = get_single_sided_fft(v_ref, t_ref)

        f_unif, p_fft_unif = get_single_sided_fft(p_unif, t_unif)
        _, v_fft_unif = get_single_sided_fft(v_unif, t_unif)

        f_orif, p_fft_orif = get_single_sided_fft(p_orif, t_orif)
        _, v_fft_orif = get_single_sided_fft(v_orif, t_orif)

        # Plot results
        plt.figure()
        plt.subplot(221)
        plt.plot(t_unif*f_exc[j],p_unif)
        try:
            plt.plot(t_ref*f_exc[j],p_ref)
        except:
            print()
        plt.plot(t_orif*f_exc[j],p_orif)
        plt.ylabel('$p^\prime$ [Pa]')

        plt.subplot(223)
        plt.plot(t_unif*f_exc[j],v_unif)
        try:
            plt.plot(t_ref*f_exc[j],v_ref)
        except:
            print()
        plt.plot(t_orif*f_exc[j],v_orif)
        plt.xlabel('$t/t_{cycle}$ [-]')
        plt.ylabel('$v^\prime$ [m/s]')
        plt.suptitle(f'SPL = {SPL[i]} dB - f = {f_exc[j]} Hz')

        plt.subplot(222)
        plt.plot(f_unif/f_exc[j],20*np.log10(np.abs(p_fft_unif)/2e-5))
        try:
            plt.plot(f_ref/f_exc[j],20*np.log10(np.abs(p_fft_ref)/2e-5))
        except:
            print()
        plt.plot(f_orif/f_exc[j],20*np.log10(np.abs(p_fft_orif)/2e-5))
        plt.ylabel('$\\tilde{p}^\prime$ [dB]')
        plt.xlim(0.5, 1.5)
        plt.ylim(0.8*SPL[i], 1.1*SPL[i])

        plt.subplot(224)
        plt.plot(f_unif/f_exc[j],np.abs(v_fft_unif))
        try:
            plt.plot(f_ref/f_exc[j],np.abs(v_fft_ref))
        except:
            print()
        plt.plot(f_orif/f_exc[j],np.abs(v_fft_orif))
        plt.xlabel('$f/f_{exc}$ [-]')
        plt.ylabel('$\\tilde{v}^\prime$ [m/s]')
        plt.xlim(0.5, 1.5)
        plt.figlegend(['TDIBC (uniform)','LBM','TDIBC (orifices)'], loc='upper center', ncol=3, bbox_to_anchor=(0.55, 0.965))
        plt.tight_layout()
        #plt.savefig(f'plots/pressure_and_velocity_{f_exc[j]}Hz_{SPL[i]}dB_orifices.png', dpi=300, bbox_inches='tight')

        plt.show()


#%%

# ----------------------------------------------------------------------------------------------------------- #
#                                              TRASH CAN                                                      #
# ----------------------------------------------------------------------------------------------------------- #

# Write probes for optydb_kundt impedance calculation

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