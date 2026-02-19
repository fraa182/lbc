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

res_vec = [10, 20, 30, 40, 50]
Re_vec = [20, 100]
tau_vec = [0.85, 0.57]

CD_mean = [0]*len(res_vec)
CD_std = [0]*len(res_vec)

CL_mean = [0]*len(res_vec)
CL_std = [0]*len(res_vec)

St_shed = [0]*len(res_vec)
CD_mean_exp = np.array([2.5, 1.4])

for j in range(len(Re_vec)):
    for i in range(len(res_vec)):

        # Read resolution and Reynolds number
        res = res_vec[i]
        Re = Re_vec[j]

        # Compute tau
        tau = tau_vec[j]#np.min([6*0.2*res/(np.sqrt(3)*Re) + 0.5, 0.9])

        # Read forces and cut transient
        dt_struct = np.dtype([('time', 'f8'), ('L', 'f8'), ('D', 'f8')])
        data = np.fromfile(f"/run/user/1000/gvfs/sftp:host=hpc-legionlogin.polito.it,user=fbellelli/home/fbellelli/lbc/cases/cylinder/sol/forces_res_{res}_Re_{Re:^.2f}_tau_{tau:^.2f}.bin", dtype=dt_struct)

        mask = (data['time'] >= np.mean(data['time']))

        t = data['time'][mask]
        L = data['L'][mask]
        D = data['D'][mask]

        # Compute timestep
        dt = t[1] - t[0]

        # Compute CL and CD from forces
        CL = L / (0.5* rho_inf * U_inf**2 * 2*R)
        CD = D / (0.5 * rho_inf * U_inf**2 * 2*R)

        CL_mean[i] = np.mean(CL)
        CL_std[i] = np.std(CL)

        CD_mean[i] = np.mean(CD)
        CD_std[i] = np.std(CD)

        # Compute FFT of CL
        f, A = compute_fft(CL, dt)

        # Find shedding frequency (max of CL FFT)
        idx = np.argmax(A)
        f_shed = f[idx]

        # Convert frequency in Strouhal number (St)
        St = f * 2*R / U_inf
        St_shed[i] = f_shed * 2*R / U_inf

        # Plot CL and CD time history
        plt.figure()
        plt.subplot(211)
        plt.plot(t, CL, 'k')
        plt.ylabel("$c_l$ [-]")
        plt.suptitle(f"Re = {Re:^.2f} - Res = {res}")

        plt.subplot(212)
        plt.plot(t, CD, 'k')
        plt.ylabel("$c_d$ [-]")
        plt.xlabel("$t$ [s]")
        plt.tight_layout()

        # # Plot CL FFT vs St
        if (Re == 100):
            plt.figure()
            plt.semilogx(St, A, 'k')
            plt.axvline(x=St_shed[i], color='r', linestyle='--')
            plt.xlabel("$St$ [-]")
            plt.ylabel("$\\tilde{c_l}$ [-]")
            plt.title(f"St shedding = {St_shed[i]:^.2f} - Re = {Re:^.2f} - Res = {res}")
            plt.tight_layout()

    plt.figure(101)
    
    plt.subplot(221)
    plt.plot(res_vec,CL_mean,'-o')
    plt.title("$\\bar{C_L}$")

    plt.subplot(222)
    plt.plot(res_vec,CL_std,'-d')
    plt.title("${C_L}_\\sigma$")

    plt.subplot(223)
    plt.plot(res_vec,CD_mean,'-o')
    plt.title("$\\bar{C_D}$")
    plt.xlabel("Res [-]")

    plt.subplot(224)
    plt.plot(res_vec,CD_std,'-d')
    plt.title("${C_D}_\\sigma$")
    plt.xlabel("Res [-]")

    plt.tight_layout()

plt.figure(101)
plt.figlegend([f"Re = {val}" for val in Re_vec], loc='upper center', ncol=len(res_vec), bbox_to_anchor=(0.5, 1.025))

if 100 in Re_vec:
    plt.figure(1001)
    plt.plot(res_vec,St_shed,'-pk')
    plt.axhline(y=0.2, color='r', linestyle='--')
    plt.xlabel("Res [-]")
    plt.ylabel("$St_{shed}$ [-]")
    plt.title(f"Re = {Re}")
    plt.tight_layout()

plt.show()