#%%
import numpy as np
import pandas as pd
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

data = pd.read_csv('/run/user/1000/gvfs/sftp:host=hpc-legionlogin.polito.it,user=fbellelli/home/fbellelli/lbc/cases/cylinder/benchmark_results.csv')
data_oldsave = pd.read_csv('/run/user/1000/gvfs/sftp:host=hpc-legionlogin.polito.it,user=fbellelli/home/fbellelli/lbc/cases/cylinder/benchmark_results_oldsave.csv')
data_nosave = pd.read_csv('/run/user/1000/gvfs/sftp:host=hpc-legionlogin.polito.it,user=fbellelli/home/fbellelli/lbc/cases/cylinder/benchmark_results_nosave.csv')

plt.figure()
plt.plot(data['Cores'],data['Real_Sec'][0]/data['Real_Sec'],'--dc',label="Code (save every 100 iterations)")
#plt.plot(data_oldsave['Cores'],data_oldsave['Real_Sec'][0]/data_oldsave['Real_Sec'],'--vg',label="Code (old save)")
plt.plot(data_nosave['Cores'],data_nosave['Real_Sec'][0]/data_nosave['Real_Sec'],'--ob',label="Code (no save)")
plt.plot(data['Cores'],np.arange(1,len(data)+1),'-r',label="Linear")
plt.xlabel('Cores')
plt.ylabel('Speed-up')
plt.legend()
#plt.savefig("plots/scalability_cpu_skylake.png", dpi=300, bbox_inches='tight')

plt.show()