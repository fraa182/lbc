import pyvista as pv
import numpy as np
import matplotlib.pyplot as plt
import glob
import re

def extract_time_history(file_pattern, point_coords, scalar_name='pressure'):
    """
    Parses VTK files and extracts the value of a variable at a specific point.
    """
    # 1. Get and sort files numerically (e.g., output_10 before output_100)
    files = glob.glob(file_pattern)
    files.sort(key=lambda f: int(re.sub('\D', '', f)))
    
    if not files:
        print("No files found matching the pattern.")
        return [], []

    history = []
    time_steps = []

    print(f"Processing {len(files)} files...")

    for i, file in enumerate(files):
        # Load the mesh
        mesh = pv.read(file)
        
        # 2. Find the index of the closest point to your target coordinates
        # Since it's STRUCTURED_POINTS, we can use find_closest_point
        point_idx = mesh.find_closest_point(point_coords)
        
        # 3. Extract the value
        # mesh.point_data returns a dictionary-like object of your SCALARS/VECTORS
        val = mesh.point_data[scalar_name][point_idx]
        
        history.append(val)
        time_steps.append(i) # Assuming each file is one time unit

    return np.array(time_steps), np.array(history)

def get_fft_analysis(data, dt):
    """Computes the FFT and frequency bins for a given signal."""
    n = len(data)
    # Remove DC offset (mean) to focus on fluctuations
    data_ac = data - np.mean(data)
    
    # Compute FFT
    fft_values = np.fft.fft(data_ac)
    fft_freq = np.fft.fftfreq(n, d=dt)
    
    # Only take the positive half of the spectrum
    pos_mask = fft_freq > 0

    return fft_freq[pos_mask], np.abs(fft_values[pos_mask])

def plot_probe_location(file_pattern, point_coords, scalar_name='pressure'):
    """
    Plots the last snapshot of the simulation and marks the probe location.
    """
    files = glob.glob(file_pattern)
    files.sort(key=lambda f: int(re.sub('\D', '', f)))
    
    if not files:
        return

    # Load the last file
    mesh = pv.read(files[-1])
    
    # Create the plotter
    plotter = pv.Plotter(window_size=[800, 600])
    
    # Add the scalar field (e.g., pressure or density)
    plotter.add_mesh(mesh, scalars=scalar_name, cmap='viridis', show_edges=False)
    
    # Add a point/marker at the probe location
    # We use a sphere or a cross to make it visible
    probe_marker = pv.Sphere(radius=mesh.length * 0.01, center=point_coords)
    plotter.add_mesh(probe_marker, color='red', label='Probe Location')
    
    # Setup view (assuming XY plane for 2D)
    plotter.view_xy()
    plotter.add_legend()
    plotter.add_axes()
    plotter.show(title=f"Probe Location on {scalar_name} Field")

# --- Configuration ---
FILE_PATTERN = "sol/fields_*.vtk"  # Adjust to your naming convention
TARGET_POINT = [0.015, 0.01, 0]   # (x, y, z) location to probe
VARIABLE = "velocity"          # Can be "pressure", "density", or "velocity"
cx = 0.001/10
cu = 10/0.1
DT = cx/cu

# --- Execution ---
t, data = extract_time_history(FILE_PATTERN, TARGET_POINT, VARIABLE)

# If plotting velocity (which is a vector), take the magnitude or a component
if VARIABLE == "velocity":
    data = data[:,1] #np.linalg.norm(data, axis=1)

plot_probe_location(FILE_PATTERN, TARGET_POINT, VARIABLE)

if True:
    # Cut transient
    t_trans = 100
    data = data[t >= t_trans]
    t = t[t >= t_trans]

    # --- Analysis ---
    freq, magnitude = get_fft_analysis(data, DT)
    freq = freq*0.002/10

    # --- Plotting ---
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(10, 8))
    plt.subplots_adjust(hspace=0.4)

    # Plot 1: Time History
    ax1.plot(t, data, color='tab:blue', linewidth=1.5)
    ax1.set_title(f"Time History of {VARIABLE} at {TARGET_POINT}")
    ax1.set_xlabel("Time")
    ax1.set_ylabel(VARIABLE)
    ax1.grid(True, alpha=0.3)

    # Plot 2: Frequency Domain (FFT)
    ax2.semilogy(freq, magnitude, color='tab:red', linewidth=1.5)
    ax2.set_title(f"Frequency Spectrum (FFT)")
    ax2.set_xlabel("Frequency (St)")
    ax2.set_ylabel("Amplitude (Log Scale)")
    ax2.grid(True, which="both", alpha=0.3)

    # Annotate Peak Frequency
    peak_idx = np.argmax(magnitude)
    ax2.annotate(f'Peak: St = {freq[peak_idx]:.3f}', 
                xy=(freq[peak_idx], magnitude[peak_idx]),
                xytext=(freq[peak_idx]*1.1, magnitude[peak_idx]*1.1),
                arrowprops=dict(arrowstyle='->'))

    plt.show()