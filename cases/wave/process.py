#%%
import io
import numpy as np
import pyvista as pv
from PIL import Image
import matplotlib.pyplot as plt

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
plt.rcParams['axes.grid'] = True
plt.rcParams['grid.linestyle'] = '--'
plt.rcParams['grid.alpha'] = 0.7

# Legend
plt.rcParams['legend.fontsize'] = 24
plt.rcParams['legend.frameon'] = False

# Figure size (optional default)
plt.rcParams['figure.figsize'] = (10, 8)

def compute_A_and_phi(tau, m, A, Nx):

    Nt = int(np.ceil(0.223 / (tau - 0.5) * (Nx / m)**2))
    t = np.arange(Nt)

    rho_cos = np.zeros(Nt)
    rho_sin = np.zeros(Nt)

    for it in range(Nt):
        try:
            # Read numerical solution
            mesh = pv.read(f"sol/fields_tau_{tau:.6f}_m_{m}_A_{A:.6f}_{it:05d}.vtk")
            rho_fluct = mesh['density'][0:Nx] - 1.0

            rho_cos[it] = 2 / Nx * np.sum(np.cos(k * x) * rho_fluct)
            rho_sin[it] = 2 / Nx * np.sum(np.sin(k * x) * rho_fluct)
        except:
            continue

    A_t = np.sqrt(rho_cos**2 + rho_sin**2)
    phi_t = np.arctan2(rho_sin, rho_cos)

    return t, A_t, phi_t


def fit_acoustic_mode(t, A_t, phi_t, k, cs=1/np.sqrt(3), discard_periods=0, amp_threshold=0):

    t = np.asarray(t)
    A_t = np.asarray(A_t)
    phi_t = np.asarray(phi_t)

    # --- unwrap phase ---
    phi_unwrapped = np.unwrap(phi_t)

    # --- estimate acoustic period ---
    omega_est = cs * k
    T_est = 2 * np.pi / omega_est

    # --- discard initial transient ---
    t_min = discard_periods * T_est

    # --- amplitude threshold ---
    A0 = A_t[0]
    valid_amp = A_t > amp_threshold * A0

    # --- combine masks ---
    mask = (t >= t_min) & valid_amp

    t_fit = t[mask]
    A_fit = A_t[mask]
    phi_fit = phi_unwrapped[mask]

    if len(t_fit) < 10:
        raise ValueError("Not enough data points after filtering")

    # --- fit log(A) = log(A0) - Gamma t ---
    logA = np.log(A_fit)
    coeffs_amp = np.polyfit(t_fit, logA, 1)
    Gamma_num = -coeffs_amp[0]

    # --- fit phi = omega t + phi0 ---
    coeffs_phase = np.polyfit(t_fit, phi_fit, 1)
    omega_num = coeffs_phase[0]

    # --- compute phase speed ---
    c_num = omega_num / k

    return Gamma_num, c_num


# UI Parameters
Nx = 100
x = np.arange(Nx)

tau_fix = 0.9
m_fix = 5
A_fix = 0.0002

for tau in np.arange(0.55, 1.10, 0.05):
    m = m_fix
    A = A_fix

    k = m / Nx
    t, A_t, phi_t = compute_A_and_phi(tau, m, A, Nx)
    Gamma_num, c_num = fit_acoustic_mode(t, A_t, phi_t, k)

    c_ref = 1/np.sqrt(3)
    Gamma_ref = (1/6) * (tau - 0.5) * (k**2)

    plt.figure(1)
    plt.subplot(211)
    plt.plot(tau, Gamma_num / Nx, 'ob')
    plt.plot(tau, Gamma_ref, 'xr')
    plt.ylabel('$\\Gamma$ [1/s]')

    plt.subplot(212)
    plt.plot(tau, c_num / 6, 'ob')
    plt.plot(tau, c_ref, 'xr')
    plt.ylim(0, 1)
    plt.ylabel('$c$ [m/s]')
    plt.xlabel('$\\tau/\\Delta t$ [-]')

plt.tight_layout()
plt.figlegend(['Num','Ref'], loc='upper center', ncol=2, bbox_to_anchor=(0.55, 1.05))
plt.savefig('plots/relaxation_time_study.png', dpi=300, bbox_inches='tight')

for m in np.arange(5, 25, 5):
    tau = tau_fix
    A = A_fix

    k = m / Nx
    t, A_t, phi_t = compute_A_and_phi(tau, m, A, Nx)
    Gamma_num, c_num = fit_acoustic_mode(t, A_t, phi_t, k)

    c_ref = 1/np.sqrt(3)
    Gamma_ref = (1/6) * (tau - 0.5) * (k**2)

    plt.figure(2)
    plt.subplot(211)
    plt.plot(m, Gamma_num / Nx, 'ob')
    plt.plot(m, Gamma_ref, 'xr')
    plt.ylabel('$\\Gamma$ [1/s]')

    plt.subplot(212)
    plt.plot(m, c_num / 6, 'ob')
    plt.plot(m, c_ref, 'xr')
    plt.ylim(0, 1)
    plt.ylabel('$c$ [m/s]')
    plt.xlabel('$m$ [-]')

plt.tight_layout()
plt.figlegend(['Num','Ref'], loc='upper center', ncol=2, bbox_to_anchor=(0.55, 1.05))
plt.savefig('plots/fourier_mode_study.png', dpi=300, bbox_inches='tight')

for A in [0.0002, 0.002, 0.02, 0.2]:
    tau = tau_fix
    m = m_fix

    k = m / Nx
    t, A_t, phi_t = compute_A_and_phi(tau, m, A, Nx)
    Gamma_num, c_num = fit_acoustic_mode(t, A_t, phi_t, k)

    c_ref = 1/np.sqrt(3)
    Gamma_ref = (1/6) * (tau - 0.5) * (k**2)

    plt.figure(3)
    plt.subplot(211)
    plt.semilogx(A, Gamma_num / Nx, 'ob')
    plt.semilogx(A, Gamma_ref, 'xr')
    plt.ylim(0, 0.001)
    plt.ylabel('$\\Gamma$ [1/s]')

    plt.subplot(212)
    plt.semilogx(A, c_num / 6, 'ob')
    plt.semilogx(A, c_ref, 'xr')
    plt.ylim(0, 1)
    plt.ylabel('$c$ [m/s]')
    plt.xlabel('$A$ [kg/m$^3$]')

plt.tight_layout()
plt.figlegend(['Num','Ref'], loc='upper center', ncol=2, bbox_to_anchor=(0.55, 1.05))
plt.savefig('plots/amplitude_study.png', dpi=300, bbox_inches='tight')

plt.show()