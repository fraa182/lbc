#%%
import numpy as np
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

def compute_fft(signal, dt):

    # Compute number of samples
    N = len(signal)

    # Remove mean
    signal = signal - np.mean(signal)

    # FFT
    fft_vals = np.fft.rfft(signal)
    A = np.abs(fft_vals) / N

    # Frequency vector
    f = np.fft.rfftfreq(N, dt)

    return f, A

# UI parameters
R = 1.0
rho_inf = 1.0
U_inf = 1.0

res = 20
Re = 20
tau = 0.7

# Read forces and cut transient
dt_struct = np.dtype([('time', 'f8'), ('L', 'f8'), ('D', 'f8')])
data = np.fromfile(f"sol/forces_res_{res}_Re_{Re:^.2f}_tau_{tau:^.2f}.bin", dtype=dt_struct)

mask = (data['time'] >= 1)

t = data['time'][mask]
L = data['L'][mask]
D = data['D'][mask]

# Compute timestep
dt = t[1] - t[0]

# Compute CL and CD from forces
CL = L / (0.5* rho_inf * U_inf**2 * 2*R)
CD = D / (0.5 * rho_inf * U_inf**2 * 2*R)

# Compute FFT of CL
f, A = compute_fft(CL, dt)

# Find shedding frequency (max of CL FFT)
idx = np.argmax(A)
f_shed = f[idx]

# Convert frequency in Strouhal number (St)
St = f * 2*R / U_inf
St_shed = f_shed * 2*R / U_inf

# Plot CL and CD time history
plt.figure()
plt.subplot(211)
plt.plot(t, CL, 'k')
plt.axhline(y = 0, color='r', linestyle='--')
plt.ylabel("$c_l$ [-]")

plt.subplot(212)
plt.plot(t, CD, 'k')
plt.axhline(y = 2, color='r', linestyle='--')
plt.ylabel("$c_d$ [-]")
plt.xlabel("$t$ [s]")
plt.tight_layout()

# Plot CL FFT vs St
plt.figure()
plt.semilogx(St, A, 'k')
plt.axvline(x=St_shed, color='r', linestyle='--')
plt.xlabel("$St$ [-]")
plt.ylabel("$\\tilde{c_l}$ [-]")
plt.title(f"St shedding = {St_shed:^.2f}")
plt.tight_layout()

plt.show()