#%%
import numpy as np
import matplotlib.pyplot as plt

dt_struct = np.dtype([('time', 'f8'), ('L', 'f8'), ('D', 'f8')])
dt_struct_pf = np.dtype([('time', 'f8'), ('D', 'f8'), ('L', 'f8')])

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
rho_inf_pf = 1.1839
U_inf = 1.0
Re_shed = 100
plotTimeHistoryFlag = False

res_vec = np.array([10, 20, 30, 40, 50])
Re_vec = np.array([20, 100])
tau_vec = np.array([0.85, 0.57])
IBB_str = "IBB"

T_conv_proc = 20
T_conv_proc_pf = 5
CD_mean_exp = 1.18 + 6.8/Re_vec**0.89 + 1.96/Re_vec**0.5 - 0.0004*Re_vec/(1 + 3.64e-7*Re_vec**2)

for j in range(len(Re_vec)):
    CD_mean = [np.nan]*len(res_vec)
    CD_std = [np.nan]*len(res_vec)

    CL_mean = [np.nan]*len(res_vec)
    CL_std = [np.nan]*len(res_vec)

    CD_mean_pf = [np.nan]*len(res_vec)
    CD_std_pf = [np.nan]*len(res_vec)

    CL_mean_pf = [np.nan]*len(res_vec)
    CL_std_pf = [np.nan]*len(res_vec)

    St_shed = [np.nan]*len(res_vec)
    St_shed_pf = [np.nan]*len(res_vec)
    for i in range(len(res_vec)):

        # Read resolution and Reynolds number
        res = res_vec[i]
        Re = Re_vec[j]

        # Compute tau
        tau = tau_vec[j]#np.min([6*0.2*res/(np.sqrt(3)*Re) + 0.5, 0.9])

        # Read forces and cut transient
        try:
            data = np.fromfile(f"/run/user/1000/gvfs/sftp:host=hpc-legionlogin.polito.it,user=fbellelli/mnt/beegfs/fbellelli/lbc/cases/cylinder/sol/forces_{IBB_str}_res_{res}_Re_{Re:^.2f}_tau_{tau:^.2f}.bin", dtype=dt_struct)
            data_pf = np.genfromtxt(f"powerflow/forces_Re{Re}_res{res}.txt", dtype=dt_struct_pf, skip_header=22)

            T_conv_last = np.floor(data['time'][-1]*U_inf/(2*R))
            T_conv_last_pf = np.floor(data_pf['time'][-1]*U_inf/(2*R))

            mask = (data['time'] >= (T_conv_last-T_conv_proc)*(2*R)/U_inf) & (data['time'] <= (T_conv_last)*(2*R)/U_inf)
            mask_pf = (data_pf['time'] >= (T_conv_last_pf-T_conv_proc_pf)*(2*R)/U_inf) & (data_pf['time'] <= (T_conv_last_pf)*(2*R)/U_inf)

            t = data['time'][mask]
            L = data['L'][mask]
            D = data['D'][mask]

            t_pf = data_pf['time'][mask_pf]
            L_pf = data_pf['L'][mask_pf]
            D_pf = data_pf['D'][mask_pf]

            t = t - t[0]
            t_pf = t_pf - t_pf[0]

            # Compute timestep
            dt = t[1] - t[0]
            dt_pf = t_pf[1] - t_pf[0]

            # Compute CL and CD from forces
            CL = L / (0.5* rho_inf * U_inf**2 * 2*R)
            CD = D / (0.5 * rho_inf * U_inf**2 * 2*R)

            CL_pf = L_pf / (0.5* rho_inf_pf * U_inf**2 * 2*R)
            CD_pf = D_pf / (0.5 * rho_inf_pf * U_inf**2 * 2*R)

            CL_mean[i] = np.mean(CL)
            CL_std[i] = np.std(CL-CL_mean[i])

            CD_mean[i] = np.mean(CD)
            CD_std[i] = np.std(CD-CD_mean[i])

            CL_mean_pf[i] = np.mean(CL_pf)
            CL_std_pf[i] = np.std(CL_pf)

            CD_mean_pf[i] = np.mean(CD_pf)
            CD_std_pf[i] = np.std(CD_pf)

            if (Re == Re_shed):
                # Compute FFT of CL
                f, A = compute_fft(CL, dt)
                f_pf, A_pf = compute_fft(CL_pf, dt_pf)

                # Find shedding frequency (max of CL FFT)
                idx = np.argmax(A)
                f_shed = f[idx]

                idx_pf = np.argmax(A_pf)
                f_shed_pf = f[idx_pf]

                # Convert frequency in Strouhal number (St)
                St = f * 2*R / U_inf
                St_shed[i] = f_shed * 2*R / U_inf

                St_pf = f_pf * 2*R / U_inf
                St_shed_pf[i] = f_shed_pf * 2*R / U_inf

            if plotTimeHistoryFlag:
                # Plot CL and CD time history
                plt.figure()
                plt.subplot(211)
                plt.plot(t*U_inf/(2*R), CL, 'b')
                plt.axhline(y=0, color='r', linestyle='--')
                #plt.plot(t_pf*U_inf/(2*R), CL_pf, 'r')
                plt.ylabel("$c_l$ [-]")
                plt.suptitle(f"Re = {Re:^.2f} - Res = {res}")

                plt.subplot(212)
                plt.plot(t*U_inf/(2*R), CD, 'b')
                #plt.plot(t_pf*U_inf/(2*R), CD_pf, 'r')
                plt.axhline(y=CD_mean_exp[j], color='r', linestyle='--')
                plt.ylabel("$c_d$ [-]")
                plt.xlabel("$tU_\\infty /D$ [-]")
                #plt.figlegend(['LBC','PowerFLOW'], loc='upper center', ncol=2, bbox_to_anchor=(0.5, 0.95))
                plt.tight_layout()
                plt.savefig(f"plots/time_history_Re{Re}_res{res}.png", dpi=300, bbox_inches='tight')

                # # Plot CL FFT vs St
                if (Re == Re_shed):
                    plt.figure()
                    plt.semilogx(St, A, '-ob')
                    #plt.semilogx(St_pf, A_pf, '--r')
                    plt.xlabel("$St$ [-]")
                    plt.ylabel("$\\tilde{c_l}$ [-]")
                    plt.title(f"Res = {res}")
                    plt.tight_layout()
                    plt.savefig(f"plots/fft_CL_Re{Re}_res{res}.png", dpi=300, bbox_inches='tight')
        except Exception:
            raise Exception

    plt.figure()
    
    plt.subplot(221)
    plt.plot(res_vec,CL_mean,'-ob')
    plt.axhline(y=0, color='r', linestyle='--')
    #plt.plot(res_vec,CL_mean_pf,'--sr')
    plt.title("$\\bar{C_L}$ [-]")

    plt.subplot(222)
    plt.plot(res_vec,CL_std,'-db')
    if (Re == Re_shed):
        plt.axhline(y=0.4, color='r', linestyle='--')
    else:
        plt.axhline(y=0, color='r', linestyle='--')
    #plt.plot(res_vec,CL_std_pf,'--vr')
    plt.title("${C_L}_\\sigma$ [-]")

    plt.subplot(223)
    plt.plot(res_vec,CD_mean,'-ob')
    #plt.plot(res_vec,CD_mean_pf,'--sr')
    plt.axhline(y=CD_mean_exp[j], color='r', linestyle='--')
    plt.title("$\\bar{C_D}$ [-]")
    plt.xlabel("Res [-]")

    plt.subplot(224)
    plt.plot(res_vec,CD_std,'-db')
    #plt.plot(res_vec,CD_std_pf,'--vr')
    plt.axhline(y=0, color='r', linestyle='--')
    plt.title("${C_D}_\\sigma$ [-]")
    plt.xlabel("Res [-]")
    plt.suptitle(f"Re = {Re}")
    #plt.figlegend(['LBC','PowerFLOW'], loc='upper center', ncol=2, bbox_to_anchor=(0.5, 0.95))
    plt.tight_layout()
    plt.savefig(f"plots/grid_convergence_Re{Re}.png", dpi=300, bbox_inches='tight')

    if (Re == Re_shed):
        plt.figure()
        plt.plot(res_vec,St_shed,'-pb')
        #plt.plot(res_vec,St_shed_pf,'--xr')
        plt.axhline(y=0.21, color='r', linestyle='--')
        plt.ylim(0.17,0.23)
        plt.xlabel("Res [-]")
        plt.ylabel("$St_{shed}$ [-]")
        plt.title(f"Re = {Re}")
        plt.tight_layout()
        plt.savefig(f"plots/shedding_Re{Re}.png", dpi=300, bbox_inches='tight')

plt.show()