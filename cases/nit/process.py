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
from scipy.optimize import curve_fit
from scipy.io import loadmat, savemat
from scipy.integrate import trapezoid
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
#                                         INPUT PARAMETERS                                                    #
# ----------------------------------------------------------------------------------------------------------- #

# Parameters
c0 = 340
rho = 1.184
nu = 1.5e-5

N_o = 2
d = 1.17e-3
L = 9.906e-3

# Derivate quantities
T = 101325 / (287.05 * rho)
mu = nu * rho
POA = N_o * d / L

# Define frequency grid and compute angular frequency
f_hz = np.logspace(np.log10(500), np.log10(2500), 100)
omega = 2 * np.pi * f_hz

# SPL and excitation frequency
SPL = [130, 145]
f_exc = [800, 1000, 1400, 2000]

# Path to folders
foldername_tdibc = '/run/user/1000/gvfs/sftp:host=hpc-legionlogin.polito.it,user=fbellelli/mnt/beegfs/fbellelli/lbc/cases/nit'
foldername_explicit = '/run/user/1000/gvfs/sftp:host=hpc-legionlogin.polito.it,user=fbellelli/mnt/beegfs/fbellelli/lbc/cases/nit_explicit'

# ----------------------------------------------------------------------------------------------------------- #
#                                         ANCILLARY FUNCTIONS                                                 #
# ----------------------------------------------------------------------------------------------------------- #

def plot_mesh(ax, x, y, z, title, vmin, vmax):
    X, Y = np.meshgrid(x / d, y / d)
    im = ax.pcolormesh(X, Y, z, cmap='bwr', vmin=vmin, vmax=vmax)
    ax.set_xlim([-10, 0])
    ax.set_ylim([0, 10])
    ax.set_title(title)
    ax.set_xlabel(r'$x/d$ [-]')
    ax.set_ylabel(r'$y/d$ [-]')
    return im

#%%

# ----------------------------------------------------------------------------------------------------------- #
#                                           IMPEDANCE CALCULATION                                             #
# ----------------------------------------------------------------------------------------------------------- #

Z = {}

for i in range(len(SPL)):
    p_a = 89 if SPL[i] == 130 else 503
    Z[str(SPL[i])] = np.zeros(len(f_exc), dtype=complex)
    for j in range(len(f_exc)):
        # Read explicit simulation
        data_explicit = np.genfromtxt(f'{foldername_explicit}/explicit/tdibc_{f_exc[j]}Hz_{p_a}Pa.txt', invalid_raise=False, skip_header=0)
        t_explicit = data_explicit[:, 0]
        p_explicit = data_explicit[:, 1]
        v_explicit = data_explicit[:, 2]      

        f_p = interp1d(t_explicit, p_explicit)
        f_v = interp1d(t_explicit, v_explicit)

        t = np.linspace(t_explicit[0], t_explicit[-1], len(t_explicit))
        p = f_p(t)
        v = f_v(t)

        f, p_fft = get_single_sided_fft(p, t)
        _, v_fft = get_single_sided_fft(v, t)

        Z[str(SPL[i])][j] = p_fft[np.argmin(np.abs(f-f_exc[j]))] / v_fft[np.argmin(np.abs(f-f_exc[j]))]

for i in range(len(SPL)):
    fR = interp1d(f_exc, np.real(Z[str(SPL[i])]), fill_value='extrapolate')
    fX = interp1d(f_exc, np.imag(Z[str(SPL[i])]), fill_value='extrapolate')

    Z_interp = fR(f_hz) + 1j * fX(f_hz)

    Yinf, lambdas, A, alpha, beta, B, C = impedance_to_multipole(omega, 1/Z_interp, 2, 2)
    Y = multipole_to_impedance(omega, Yinf, lambdas, A, alpha, beta, B, C)

    plt.figure()
    plt.subplot(211)
    plt.plot(f_hz, rho * c0 * np.real(1 / Z_interp), '-b')
    plt.plot(f_hz, rho * c0 * np.real(Y), '--r')
    plt.plot(f_exc, rho * c0 * np.real(1 / Z[str(SPL[i])]), 'ob')
    plt.ylabel('$\\Re(\\tilde{Y})\\,\\rho_0 c_0$ [-]')

    plt.subplot(212)
    plt.plot(f_hz, rho * c0 * np.imag(1 / Z_interp), '-b')
    plt.plot(f_hz, rho * c0 * np.imag(Y), '--r')
    plt.plot(f_exc, rho * c0 * np.imag(1 / Z[str(SPL[i])]), 'sb')
    plt.ylabel('$\\Im(\\tilde{Y})\\,\\rho_0 c_0$ [-]')
    plt.xlabel('$f$ [Hz]')

    plt.figlegend(['Measured','Fit ($N_R=2$, $N_{CC}=2$)'], ncols=2, loc='upper center', bbox_to_anchor=(0.55, 0.965))
    plt.suptitle(f'SPL = {SPL[i]} dB')
    plt.tight_layout()
    plt.savefig(f'plots/admittance_fit_{SPL[i]}dB.png', dpi=300, bbox_inches='tight')

    plt.show()

#%%

# ----------------------------------------------------------------------------------------------------------- #
#                                         A POSTERIORI VALIDATION                                             #
# ----------------------------------------------------------------------------------------------------------- #

for i in range(len(SPL)):
    p_a = 89 if SPL[i] == 130 else 503
    filename_impedance = f'{foldername_tdibc}/impedance_{p_a}Pa_notscaled_test.txt'
    for j in range(len(f_exc)):
        # Read explicit simulation
        data_explicit = np.genfromtxt(f'{foldername_explicit}/explicit/tdibc_{f_exc[j]}Hz_{p_a}Pa.txt', invalid_raise=False, skip_header=0)
        t_explicit = data_explicit[:, 0]
        p_explicit = data_explicit[:, 1]
        v_explicit = data_explicit[:, 2]      

        f_p = interp1d(t_explicit, p_explicit)
        f_v = interp1d(t_explicit, v_explicit)

        t = np.linspace(t_explicit[0], t_explicit[-1], len(t_explicit))
        p = f_p(t)
        v = f_v(t)

        f, p_fft = get_single_sided_fft(p, t)
        _, v_fft = get_single_sided_fft(v, t)

        # Read admittance vector fit parameters
        loaded_arrays = []

        with open(filename_impedance, 'r') as fp:
            Yinf = np.array([float(fp.readline().strip())])
            
            for line in fp:
                if line.strip():
                    arr = np.fromstring(line, sep=' ')
                    loaded_arrays.append(arr)

        A, lambdas, B, C, alpha, beta = loaded_arrays

        # Compute velocity as a function of pressure and admittance with ADE
        t_bc, p_bc, v_bc = tdim_bc(t, lambdas, A, alpha, beta, B, C, Yinf, f_p, vacFlag=False, EEFlag=True)

        f_bc, p_fft_bc = get_single_sided_fft(p_bc, t_bc)
        _, v_fft_bc = get_single_sided_fft(v_bc, t_bc)

        # Plot results
        plt.figure()
        plt.subplot(321)
        plt.plot(t_bc*f_exc[j],p_bc, '-b')
        plt.plot(t*f_exc[j],p, '--r')
        plt.ylabel('$p^\prime$ [Pa]')
        plt.xlabel('$t/t_{cycle}$ [-]')
        plt.xlim(0, 10)

        plt.subplot(323)
        plt.plot(f_bc/f_exc[j],20*np.log10(np.abs(p_fft_bc)/2e-5), '-b')
        plt.plot(f/f_exc[j],20*np.log10(np.abs(p_fft)/2e-5), '--r')
        plt.xlim(0, 2)
        plt.ylim(0.7*SPL[i],1.1*SPL[i])
        plt.ylabel('$\\tilde{p}^\prime$ [dB]')

        plt.subplot(325)
        plt.plot(f_bc/f_exc[j],np.angle(p_fft_bc), '-b')
        plt.plot(f/f_exc[j],np.angle(p_fft), '--r')
        plt.xlim(0, 2)
        plt.ylabel('$\\angle{p}^\prime$ [dB]')
        plt.xlabel('$f/f_{exc}$ [-]')

        plt.subplot(322)
        plt.plot(t_bc*f_exc[j],v_bc, '-b')
        plt.plot(t*f_exc[j],v, '--r')
        plt.ylabel('$v^\prime$ [m/s]')
        plt.xlabel('$t/t_{cycle}$ [-]')
        plt.xlim(0, 10)

        plt.subplot(324)
        plt.plot(f_bc/f_exc[j],np.abs(v_fft_bc), '-b')
        plt.plot(f/f_exc[j],np.abs(v_fft), '--r')
        plt.xlim(0, 2)
        plt.ylabel('$\\tilde{v}^\prime$ [m/s]')

        plt.subplot(326)
        plt.plot(f_bc/f_exc[j],np.angle(v_fft_bc), '-b')
        plt.plot(f/f_exc[j],np.angle(v_fft), '--r')
        plt.xlim(0, 2)
        plt.ylabel('$\\angle{v}^\prime$ [m/s]')
        plt.xlabel('$f/f_{exc}$ [-]')

        plt.figlegend(['ADE','Sample'], ncols=2, loc='upper center', bbox_to_anchor=(0.525, 0.965))
        plt.suptitle(f'SPL = {SPL[i]} dB - f = {f_exc[j]} Hz')
        plt.tight_layout()
        plt.savefig(f'plots/a_posteriori_validation_{f_exc[j]}Hz_{SPL[i]}dB.png', dpi=300, bbox_inches='tight')

        plt.show()
    
#%%

# ----------------------------------------------------------------------------------------------------------- #
#                                               SURFACE                                                       #
# ----------------------------------------------------------------------------------------------------------- #

for i in range(len(SPL)):
    p_a = 89 if SPL[i] == 130 else 503
    for j in range(len(f_exc)):
        # Read TDIBC uniform simulation
        data_uniform = np.genfromtxt(f'{foldername_tdibc}/uniform/tdibc_{f_exc[j]}Hz_{p_a}Pa.txt', invalid_raise=False, skip_header=0)
        t_uniform = data_uniform[:, 0]
        p_uniform = data_uniform[:, 1]
        v_uniform = data_uniform[:, 2]

        f_p = interp1d(t_uniform, p_uniform)
        f_v = interp1d(t_uniform, v_uniform)
        t = np.linspace(t_uniform[0], t_uniform[-1], len(t_uniform))
        p = f_p(t)
        v = f_v(t)

        f_uniform, p_fft_uniform = get_single_sided_fft(p, t)
        _, v_fft_uniform = get_single_sided_fft(v, t)

        # Read TDIBC orifices simulation
        data_orifices = np.genfromtxt(f'{foldername_tdibc}/orifices/tdibc_{f_exc[j]}Hz_{p_a}Pa.txt', invalid_raise=False, skip_header=0)
        t_orifices = data_orifices[:, 0]
        p_orifices = data_orifices[:, 1]
        v_orifices = data_orifices[:, 2]

        f_p = interp1d(t_orifices, p_orifices)
        f_v = interp1d(t_orifices, v_orifices)
        t = np.linspace(t_orifices[0], t_orifices[-1], len(t_orifices))
        p = f_p(t)
        v = f_v(t)

        f_orifices, p_fft_orifices = get_single_sided_fft(p, t)
        _, v_fft_orifices = get_single_sided_fft(v, t)

        # Read explicit simulation
        data_explicit = np.genfromtxt(f'{foldername_explicit}/explicit/tdibc_{f_exc[j]}Hz_{p_a}Pa.txt', invalid_raise=False, skip_header=0)
        t_explicit = data_explicit[:, 0]
        p_explicit = data_explicit[:, 1]
        v_explicit = data_explicit[:, 2]      

        f_p = interp1d(t_explicit, p_explicit)
        f_v = interp1d(t_explicit, v_explicit)
        t = np.linspace(t_explicit[0], t_explicit[-1], len(t_explicit))
        p = f_p(t)
        v = f_v(t)

        f_explicit, p_fft_explicit = get_single_sided_fft(p, t)
        _, v_fft_explicit = get_single_sided_fft(v, t)

        # Plot results
        plt.figure()
        plt.subplot(321)
        plt.plot(t_uniform*f_exc[j],p_uniform, '-.k')
        plt.plot(t_orifices*f_exc[j],p_orifices, '-b')
        plt.plot(t_explicit*f_exc[j],p_explicit, '--r')
        plt.ylabel('$p^\prime$ [Pa]')
        plt.xlabel('$t/t_{cycle}$ [-]')
        plt.xlim(0, 10)

        plt.subplot(323)
        plt.plot(f_uniform/f_exc[j],20*np.log10(np.abs(p_fft_uniform)/2e-5), '-.k')
        plt.plot(f_orifices/f_exc[j],20*np.log10(np.abs(p_fft_orifices)/2e-5), '-b')
        plt.plot(f_explicit/f_exc[j],20*np.log10(np.abs(p_fft_explicit)/2e-5), '--r')
        plt.xlim(0, 2)
        plt.ylim(0.7*SPL[i],1.1*SPL[i])
        plt.ylabel('$\\tilde{p}^\prime$ [dB]')

        plt.subplot(325)
        plt.plot(f_uniform/f_exc[j],np.angle(p_fft_uniform), '-.k')
        plt.plot(f_orifices/f_exc[j],np.angle(p_fft_orifices), '-b')
        plt.plot(f_explicit/f_exc[j],np.angle(p_fft_explicit), '--r')
        plt.xlim(0, 2)
        plt.ylabel('$\\angle{p}^\prime$ [dB]')
        plt.xlabel('$f/f_{exc}$ [-]')

        plt.subplot(322)
        plt.plot(t_uniform*f_exc[j],v_uniform, '-.k')
        plt.plot(t_orifices*f_exc[j],v_orifices, '-b')
        plt.plot(t_explicit*f_exc[j],v_explicit, '--r')
        plt.ylabel('$v^\prime$ [m/s]')
        plt.xlabel('$t/t_{cycle}$ [-]')
        plt.xlim(0, 10)

        plt.subplot(324)
        plt.plot(f_uniform/f_exc[j],np.abs(v_fft_uniform), '-.k')
        plt.plot(f_orifices/f_exc[j],np.abs(v_fft_orifices), '-b')
        plt.plot(f_explicit/f_exc[j],np.abs(v_fft_explicit), '--r')
        plt.xlim(0, 2)
        plt.ylabel('$\\tilde{v}^\prime$ [m/s]')

        plt.subplot(326)
        plt.plot(f_uniform/f_exc[j],np.angle(v_fft_uniform), '-.k')
        plt.plot(f_orifices/f_exc[j],np.angle(v_fft_orifices), '-b')
        plt.plot(f_explicit/f_exc[j],np.angle(v_fft_explicit), '--r')
        plt.xlim(0, 2)
        plt.ylabel('$\\angle{v}^\prime$ [m/s]')
        plt.xlabel('$f/f_{exc}$ [-]')

        plt.figlegend(['TDIBC Uniform','TDIBC Orifices','Explicit'], ncols=3, loc='upper center', bbox_to_anchor=(0.525, 0.965))
        plt.suptitle(f'SPL = {SPL[i]} dB - f = {f_exc[j]} Hz')
        plt.tight_layout()
        plt.savefig(f'plots/pressure_and_velocity_{f_exc[j]}Hz_{SPL[i]}dB.png', dpi=300, bbox_inches='tight')

        plt.show()

#%%

# ----------------------------------------------------------------------------------------------------------- #
#                                               FLOW VIZ                                                      #
# ----------------------------------------------------------------------------------------------------------- #

for i in range(len(SPL)):
    p_a = 89 if SPL[i] == 130 else 503
    for j in range(len(f_exc)):
        # Load inflow and outflow fields
        explicit = loadmat(f'{foldername_explicit}/v_in_out_{f_exc[j]}Hz_{SPL[i]}dB_explicit.mat', squeeze_me=True)
        tdibc = loadmat(f'{foldername_tdibc}/v_in_out_{f_exc[j]}Hz_{SPL[i]}dB_orifices.mat', squeeze_me=True)

        # Compute wavelength
        lam = c0 / f_exc[j]

        # TDIBC Processing
        tdibc['res'] = 10
        tdibc['dx'] = d / tdibc['res']

        tdibc['y'] = tdibc['dx'] * np.arange(tdibc['v_in'].shape[0]) + 5.2650e-04
        tdibc['x'] = tdibc['dx'] * np.arange(tdibc['v_in'].shape[1])
        tdibc['x'] = tdibc['x'] - np.max(tdibc['x'])

        tdibc['mask'] = np.where((tdibc['x'] / d >= -10) & (tdibc['x'] / lam <= 12))[0]

        tdibc['v_in'] = tdibc['v_in'][:, tdibc['mask']]
        tdibc['v_out'] = tdibc['v_out'][:, tdibc['mask']]
        tdibc['x'] = tdibc['x'][tdibc['mask']] - np.max(tdibc['x'][tdibc['mask']])

        val_tdibc_in = trapezoid(tdibc['v_in'][:, -1], tdibc['y']) / L
        val_tdibc_out = trapezoid(tdibc['v_out'][:, -1], tdibc['y']) / L

        # Explicit Processing
        explicit['res'] = 20
        explicit['dx'] = d / explicit['res']

        explicit['y'] = explicit['dx'] * np.arange(explicit['v_in'].shape[0])
        explicit['x'] = explicit['dx'] * np.arange(explicit['v_in'].shape[1])

        explicit['mask'] = np.where((explicit['x'] / d >= -10) & (explicit['x'] / lam <= 12))[0]

        explicit['v_in'] = explicit['v_in'][:, explicit['mask']]
        explicit['v_out'] = explicit['v_out'][:, explicit['mask']]
        explicit['x'] = explicit['x'][explicit['mask']] - np.max(explicit['x'][explicit['mask']])

        explicit_v_in_mean = explicit['v_in'][:, -1]
        explicit_v_out_mean = explicit['v_out'][:, -1]

        val_explicit_in = trapezoid(explicit_v_in_mean, explicit['y']) / L
        val_explicit_out = trapezoid(explicit_v_out_mean, explicit['y']) / L

        # Plot results
        vmin = np.nanmin([np.nanmin(explicit['v_in']), np.nanmin(tdibc['v_in'])])
        vmax = np.nanmax([np.nanmax(explicit['v_out']), np.nanmax(tdibc['v_out'])])

        fig = plt.figure()
        ax1 = fig.add_subplot(231)
        ax1.plot(tdibc['y'] / d, tdibc['v_out'][:, -1], '-b')
        ax1.plot(explicit['y'] / d, explicit_v_out_mean, '-r')
        ax1.plot(explicit['y'] / d, explicit_v_in_mean, '--r')
        ax1.plot(tdibc['y'] / d, tdibc['v_in'][:, -1], '--b')
        ax1.set_xlabel(r'$y/d$ [-]')
        ax1.set_ylabel(r"$v^\prime$ [m/s]")
        ax1.set_title(f'SPL = {SPL[i]} dB, f = {f_exc[j]} Hz')

        ax4 = fig.add_subplot(234)
        ax4.barh(1, val_tdibc_in, color='b', linestyle='--', linewidth=2, alpha=0.5, edgecolor='b')
        ax4.barh(2, val_explicit_in, color='r', linestyle='--', linewidth=2, alpha=0.5, edgecolor='r')
        ax4.barh(1, val_tdibc_out, color='b', linestyle='-', linewidth=2, alpha=0.75, edgecolor='b')
        ax4.barh(2, val_explicit_out, color='r', linestyle='-', linewidth=2, alpha=0.75, edgecolor='r')
        ax4.set_ylim([0.5, 2.5])
        ax4.set_yticks([1, 2])
        ax4.set_yticklabels(['TDIBC', 'Explicit'])
        ax4.set_xlabel(r'$<v_{in/out}>$ [m/s]')

        ax2 = fig.add_subplot(232)
        plot_mesh(ax2, tdibc['x'], tdibc['y'], tdibc['v_in'], 'TDIBC - Inflow', vmin, vmax)

        ax3 = fig.add_subplot(233)
        plot_mesh(ax3, tdibc['x'], tdibc['y'], tdibc['v_out'], 'TDIBC - Outflow', vmin, vmax)

        ax5 = fig.add_subplot(235)
        plot_mesh(ax5, explicit['x'], explicit['y'], explicit['v_in'], 'Explicit - Inflow', vmin, vmax)

        ax6 = fig.add_subplot(236)
        pcm = plot_mesh(ax6, explicit['x'], explicit['y'], explicit['v_out'], 'Explicit - Outflow', vmin, vmax)

        cax = fig.add_axes([0.42, -0.01, 0.55, 0.035])
        fig.colorbar(pcm, cax=cax, orientation='horizontal', label='$v^\prime$ [m/s]')

        plt.tight_layout()
        plt.savefig(f'plots/inflow_outflow_{f_exc[j]}Hz_{SPL[i]}dB.png', dpi=300, bbox_inches='tight')

        plt.show()