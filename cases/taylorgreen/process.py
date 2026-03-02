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

# UI Parameters
tau = 0.8
U_inf = 0.03
rho_inf = 1.0
Nx = 96
Ny = 72

# Compute derivate quantities
nu = (tau - 0.5) / 3
kx = 2*np.pi / Nx
ky = 2*np.pi / Ny
td = 1 / (nu * (kx**2 + ky**2))
Nt = int(np.floor(td))
t_adim = np.arange(Nt) / td

X, Y = np.meshgrid(np.arange(Nx), np.arange(Ny))

L2_u = np.zeros(Nt)
L2_v = np.zeros(Nt)
L2_p = np.zeros(Nt)

frames = []

for it in range(Nt):
    # Compute analytical solution
    u_analytical = -U_inf * np.sqrt(ky/kx) * np.cos(kx*X) * np.sin(ky*Y) * np.exp(-it/td)
    v_analytical = U_inf * np.sqrt(kx/ky) * np.sin(kx*X) * np.cos(ky*Y) * np.exp(-it/td)
    p_analytical = rho_inf * (287.05 * 298.15) - (rho_inf * U_inf**2 / 4) * (ky/kx * np.cos(2*kx*X) + kx/ky * np.cos(2*ky*Y)) * np.exp(-2*it/td)

    # Read numerical solution
    mesh = pv.read(f"sol/fields_tau_{tau:.6f}_{it:05d}.vtk")

    u = mesh['velocity'][:,0].reshape(mesh.dimensions[1], mesh.dimensions[0])
    v = mesh['velocity'][:,1].reshape(mesh.dimensions[1], mesh.dimensions[0])
    p = mesh['density'].reshape(mesh.dimensions[1], mesh.dimensions[0]) * (287.05 * 298.15)

    L2_u[it] = np.sqrt(np.sum((u - u_analytical)**2))/np.sqrt(np.sum(u_analytical**2))*100
    L2_v[it] = np.sqrt(np.sum((v - v_analytical)**2))/np.sqrt(np.sum(v_analytical**2))*100
    L2_p[it] = np.sqrt(np.sum((p - p_analytical)**2))/np.sqrt(np.sum(p_analytical**2))*100

    if (it % 10 == 0) | (it == Nt-1):
        print(f"Iteration {it} of {Nt-1}")
        plt.figure(1)
        plt.clf()

        plt.subplot(221)
        plt.imshow(np.sqrt(u**2 + v**2)/U_inf, cmap='bwr', vmin=0, vmax=1)
        plt.streamplot(X, Y, u, v, color='k', arrowstyle='-|>', arrowsize=1, density=1, linewidth=1)
        plt.axis('off')
        plt.grid(False)
        plt.title("Velocity (LBM)")

        plt.subplot(222)
        plt.imshow(np.sqrt(u_analytical**2 + v_analytical**2)/U_inf, cmap='bwr', vmin=0, vmax=1)
        plt.streamplot(X, Y, u_analytical, v_analytical, color='k', arrowstyle='-|>', arrowsize=1, density=1, linewidth=1)
        plt.axis('off')
        plt.grid(False)
        plt.title("Velocity (analytical)")

        plt.subplot(223)
        plt.imshow(p/(rho_inf * 287.05 * 298.15), cmap='bwr')
        plt.streamplot(X, Y, u, v, color='k', arrowstyle='-|>', arrowsize=1, density=1, linewidth=1)
        plt.axis('off')
        plt.grid(False)
        plt.title("Pressure (LBM)")

        plt.subplot(224)
        plt.imshow(p_analytical/(rho_inf * 287.05 * 298.15), cmap='bwr')
        plt.streamplot(X, Y, u_analytical, v_analytical, color='k', arrowstyle='-|>', arrowsize=1, density=1, linewidth=1)
        plt.axis('off')
        plt.grid(False)
        plt.title("Pressure (analytical)")

        buf = io.BytesIO()
        plt.savefig(buf, format='png', bbox_inches='tight')
        buf.seek(0)
        frames.append(Image.open(buf))

plt.figure(2)
plt.subplot(311)
plt.plot(t_adim,L2_u,'k')
plt.ylabel("$L_{2,u}$ [\\%]")

plt.subplot(312)
plt.plot(t_adim,L2_v,'k')
plt.ylabel("$L_{2,v}$ [\\%]")

plt.subplot(313)
plt.plot(t_adim,L2_p,'k')
plt.ylabel("$L_{2,p}$ [\\%]")
plt.xlabel("$t/t_d$ [-]")

plt.suptitle(f"$\\tau/\\Delta t$ = {tau:.2f}")
plt.tight_layout()
plt.savefig(f"plots/L2_errors_tau_{tau}.png", dpi=300, bbox_inches='tight')

plt.show()

if frames:
    frames[0].save(
        'plots/solution.gif',
        save_all=True,
        append_images=frames[1:],
        optimize=True,
        duration=100,
        loop=0
    )