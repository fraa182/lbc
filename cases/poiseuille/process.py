#%%
import pyvista as pv
import numpy as np
import matplotlib.pyplot as plt
import glob, os, re

def get_latest_vtk(folder_path, res=None):
    # Search for all files matching the pattern
    if res is None:
        pattern = os.path.join(folder_path, "fields_*.vtk")
    else:
        pattern = os.path.join(folder_path, f"fields_res_{res}_*.vtk")
    files = glob.glob(pattern)
    
    if not files:
        return None

    # Use a regex or string split to extract the number and find the max
    # This lambda finds the integer between '_' and '.'
    latest_file = max(files, key=lambda x: int(re.search(r'fields_res_\d+_(\d+)', x).group(1)))
    
    return latest_file

def extract_line(file, x, y, scalar_name='pressure'):

    if not file:
        print("No files found matching the pattern.")
        return [], []

    history = []

    # Load the mesh
    mesh = pv.read(file)

    for i in range(len(y)):
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
nu = 1.5e-5
tau = 0.90
dp_dx_lu = 2.0e-5

x = Lx/2
res_vec = np.array([10, 20, 30, 40, 50])
L2err = []

for i in range(len(res_vec)):
    res = res_vec[i]
    file = get_latest_vtk("sol", res)

    dx = Ly / res
    dt = (tau - 0.5) * dx**2 / (3 * nu)

    dp_dx = dp_dx_lu*(dx / dt**2)

    Ny = res
    j = np.arange(1, Ny-1)
    y = j * dx

    H_eff = Ly - 2*dx
    y0 = 0.5 * dx

    uref = (dp_dx / (2.0 * nu)) * (y - y0) * (H_eff - (y - y0))
    umax = np.max(uref)

    data = extract_line(file, x, y, "velocity")

    L2err.append(np.linalg.norm((data[:,0] - uref)/uref))

    plt.figure()
    plt.plot(y/H_eff,data[:,0]/umax,'--ob',label=f"LBM (H = {res})")
    plt.plot(y/H_eff,uref/umax,'-+r',label="Ref")
    plt.xlabel("y/H")
    plt.ylabel("u/umax")
    plt.legend()

plt.figure()
plt.plot(res_vec,L2err,'--db',label="LBM")
plt.plot(res_vec,(L2err[0]) * (res_vec[0]**2)/(res_vec**2),'-xr',label="O(H^2)")
plt.gca().set(xscale="log", yscale="log")
plt.xlabel("H")
plt.ylabel("L2 error")
plt.legend()

plt.show()
# %%
