#%%
import pyvista as pv
import numpy as np
import matplotlib.pyplot as plt
import glob, os, re

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

def get_latest_vtk(folder_path, res=10, tau=0.9):
    # Search for all files matching the pattern
    pattern = os.path.join(folder_path, f"fields_res_{res}_tau_{tau:.2f}_*.vtk")
    files = glob.glob(pattern)
    
    if not files:
        return None

    # Use a regex or string split to extract the number and find the max
    # This lambda finds the integer between '_' and '.'
    latest_file = max(files, key=lambda x: int(re.search(r'fields_res_\d+_tau_[\d\.]+_(\d+)', x).group(1)))
    
    return latest_file

def extract_vertical_line(file, scalar_name='velocity'):

    if not file:
        print("No files found matching the pattern.")
        return []

    # Load the mesh
    mesh = pv.read(file)

    quantity = np.array(mesh.point_data[scalar_name][:, 0])
    quantity = np.reshape(quantity, (mesh.dimensions[1], mesh.dimensions[0]))
    quantity = quantity[:, mesh.dimensions[0]//2]

    return quantity[1:]

def compute_reference(U_top, res):

    H_eff = res - 1
    y_bot = 0.5
    y_top = H_eff

    y = np.arange(1,res)

    uref = U_top*(y - y_bot)/(y_top - y_bot)

    return y, uref, H_eff

# UI parameters (fixed)
U_top = 0.1

velprofile = True
gridstudy = True
taustudy = True

# Relaxation time for the grid resolution study and resolution for the optimal relaxation time study
tau_res = 0.90
res_tau = 11

# Define arrays of resolution and relaxation time
res_vec = np.array([10, 20, 30, 40, 50])
tau_vec = np.concatenate([np.arange(0.51, 1.0, 0.01), np.arange(1, 5.1, 0.1)])

# Initialize L2 error lists
L2err_res = []
L2err_tau = []

res = 10
tau = 0.9

# Comparison against analytical result
if velprofile:
    file = get_latest_vtk("sol", res, tau)

    y, uref, H_eff = compute_reference(U_top, res)

    data = extract_vertical_line(file)

    plt.figure()
    plt.plot(y/H_eff,data/U_top,'-ob',label='LBM')
    plt.plot(y/H_eff,uref/U_top,'--xr',label='Analytic')
    plt.xlabel('$y/y_w$ [-]')
    plt.ylabel('$u/u_w$ [-]')
    plt.tight_layout()
    plt.legend()
    plt.savefig("plots/velocity_profile.png", dpi=300, bbox_inches='tight')

# Grid convergence study
if gridstudy:
    for i in range(len(res_vec)):
        res = res_vec[i]
        file = get_latest_vtk("sol", res, tau_res)

        _, uref, _ = compute_reference(U_top, res)

        data = extract_vertical_line(file)

        L2norm = np.sqrt(np.nansum((data - uref)**2))/np.sqrt(np.sum(uref**2))*100
        L2err_res.append(L2norm)

    plt.figure()
    plt.plot(res_vec,L2err_res,'--db',label="LBM")
    plt.plot(res_vec,(L2err_res[0]) * (res_vec[0])/(res_vec),'-xr',label="$\mathcal{O}(H)$")
    plt.gca().set(xscale="log", yscale="log")
    plt.xlabel("$H$ [-]")
    plt.ylabel("$L_2$ [\\%]")
    plt.title(f"Order of convergence {np.log(L2err_res[-1]/L2err_res[-2])/np.log(res_vec[-2]/res_vec[-1]):^.2f}")
    plt.legend()
    plt.tight_layout()
    plt.savefig("plots/grid_study.png", dpi=300, bbox_inches='tight')

# Optimal relaxation time study
if taustudy:
    for i in range(len(tau_vec)):
        tau = tau_vec[i]
        file = get_latest_vtk("sol", res_tau, tau)

        _, uref, _ = compute_reference(U_top, res_tau)

        data = extract_vertical_line(file)

        L2norm = np.sqrt(np.nansum((data - uref)**2))/np.sqrt(np.sum(uref**2))*100
        L2err_tau.append(L2norm)

    plt.figure()
    plt.plot(tau_vec,L2err_tau,'-k')
    plt.gca().set(yscale="log")
    plt.xlabel("$\\tau/\\Delta t$ [-]")
    plt.ylabel("$L_2$ [\\%]")
    plt.tight_layout()
    plt.savefig("plots/tau_study.png", dpi=300, bbox_inches='tight')

plt.show()