#%%
import pyvista as pv
import numpy as np
import matplotlib.pyplot as plt
import glob
import re

def extract_line(file, x, y, scalar_name='pressure'):

    if not file:
        print("No files found matching the pattern.")
        return [], []

    history = []

    for i in range(len(y)):
        # Load the mesh
        mesh = pv.read(file)
        
        # 2. Find the index of the closest point to your target coordinates
        # Since it's STRUCTURED_POINTS, we can use find_closest_point
        point_idx = mesh.find_closest_point([x,y[i],0])
        
        # 3. Extract the value
        # mesh.point_data returns a dictionary-like object of your SCALARS/VECTORS
        val = mesh.point_data[scalar_name][point_idx]
        
        history.append(val)

    return np.array(history)

Ly = 0.025
Lx = 0.100
res = 50
nu = 1.5e-5
tau = 0.90

dx = Ly / res
dt = (tau - 0.5) * dx**2 / (3 * nu)

dp_dx = 2.0e-5 * (dx / dt**2)

Ny = res
j = np.arange(1, Ny-1)
y = j * dx
H_eff = Ly - dx

x = [Lx/2]

# --- Configuration ---
file = "sol/fields_49900.vtk"  # Adjust to your naming convention
variable = "velocity"          # Can be "pressure", "density", or "velocity"

y_ref = np.linspace(0,Ly,100)
ref = (dp_dx / (2.0 * nu)) * y_ref * (Ly - y_ref)
umax = np.max(ref)

# --- Execution ---
for i in range(len(x)):
    data = extract_line(file, x[i], y, variable)

    plt.plot(y/H_eff,data[:,0]/umax,'o',label=f"x/L = {x[i]/Lx:^.2f}")

plt.plot(y_ref/Ly,ref/umax,'-r',label="Ref")
plt.xlabel("y/H")
plt.ylabel("u/umax")
plt.legend()
plt.show()
# %%
