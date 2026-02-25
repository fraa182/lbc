#%%
import glob, os
import numpy as np
import pandas as pd
import pyvista as pv
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

def ghia_reader(filename):
    # Read raw file
    with open(filename, "r") as f:
        lines = f.readlines()

    # --- Extract header ---
    header_line = None
    for line in lines:
        if line.strip().startswith("| N"):
            header_line = line
            break

    # Split and clean column names
    raw_cols = header_line.strip().strip("|").split("|")
    cols = [c.strip() for c in raw_cols]

    # Optional: normalize names (recommended)
    cols = [c.replace("=", "").replace(" ", "") for c in cols]

    # --- Extract data ---
    data_lines = []
    for line in lines:
        line = line.strip()
        if line.startswith("|") and not "---" in line and not "N" in line:
            data_lines.append(line)

    # Parse data
    cleaned = [line.strip("|").split("|") for line in data_lines]
    df = pd.DataFrame(cleaned, columns=cols)

    # Convert to numeric
    df = df.map(lambda x: float(x.strip()))
    df["N"] = df["N"].astype(int)

    return df

def lbc_reader(res, Re, tau):

    # Define the pattern with a wildcard (*) for the 'xxxxx' part
    pattern = f"sol/fields_res_{res}_Re_{Re:.6f}_tau_{tau:.6f}_*.vtk"

    # Get a list of all files matching the pattern
    files = glob.glob(pattern)

    if not files:
        print("No matching files found.")
    else:
        # Sort files by modification time (latest last)
        latest_file = max(files, key=os.path.getmtime)

    mesh = pv.read(latest_file)

    u = np.array(mesh.point_data['velocity'][:, 0])
    u = np.reshape(u, (mesh.dimensions[1], mesh.dimensions[0]))
    u = u[:, mesh.dimensions[0]//2]

    v = np.array(mesh.point_data['velocity'][:, 1])
    v = np.reshape(v, (mesh.dimensions[1], mesh.dimensions[0]))
    v = v[mesh.dimensions[1]//2, :]

    x = np.linspace(0, 1, mesh.dimensions[0])
    y = np.linspace(0, 1, mesh.dimensions[1])

    return u, v, x, y

data_u_y = ghia_reader('u_y.txt')
data_v_x = ghia_reader('v_x.txt')

res = 250
tau = 0.9
Re = [100, 400, 1000]#[100, 400, 1000, 3200, 5000, 7500, 10000]

for i in range(len(Re)):

    try:
        u, v, x, y = lbc_reader(res, Re[i], tau)

        plt.figure()
        plt.subplot(211)
        plt.plot(u,y,'-b',label='LBM')
        plt.plot(data_u_y[f'Re{Re[i]}'],data_u_y['y'],'--or')
        plt.ylabel('$y/L$ [-]')
        plt.xlabel('$u/V_\\infty$ [-]')

        plt.subplot(212)
        plt.plot(x,v,'-b')
        plt.plot(data_v_x['x'],data_v_x[f'Re{Re[i]}'],'--or')
        plt.xlabel('$x/L$ [-]')
        plt.ylabel('$v/V_\\infty$ [-]')
        plt.suptitle(f"Re = {Re[i]}")
        plt.tight_layout()
        plt.figlegend(['LBM', 'Ghia et al.'], loc='upper center', ncol=2, bbox_to_anchor=(0.5, 0.95))
        plt.savefig(f"plots/res_{res}_tau_{tau}_Re{Re[i]}.png", dpi=300, bbox_inches='tight')
    except:
        continue

plt.show()