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

data_sapphire = pd.read_csv('/run/user/1000/gvfs/sftp:host=hpc-legionlogin.polito.it,user=fbellelli/mnt/beegfs/fbellelli/lbc/cases/cylinder/benchmark_results_sapphire.csv')
data_skylake = pd.read_csv('/run/user/1000/gvfs/sftp:host=hpc-legionlogin.polito.it,user=fbellelli/mnt/beegfs/fbellelli/lbc/cases/cylinder/scalability_cpu_skylake/benchmark_results.csv')
data_oldsave_skylake = pd.read_csv('/run/user/1000/gvfs/sftp:host=hpc-legionlogin.polito.it,user=fbellelli/mnt/beegfs/fbellelli/lbc/cases/cylinder/scalability_cpu_skylake/benchmark_results_oldsave.csv')
data_nosave_skylake = pd.read_csv('/run/user/1000/gvfs/sftp:host=hpc-legionlogin.polito.it,user=fbellelli/mnt/beegfs/fbellelli/lbc/cases/cylinder/scalability_cpu_skylake/benchmark_results_nosave.csv')

plt.figure()
plt.plot(data_sapphire['Cores'],data_sapphire['Real_Sec'][0]/data_sapphire['Real_Sec'],'--dc',label="Sapphire (no save)")
#plt.plot(data_skylake['Cores'],data_skylake['Real_Sec'][0]/data_skylake['Real_Sec'],'--pk',label="Skylake (save)")
#plt.plot(data_oldsave_skylake['Cores'],data_oldsave_skylake['Real_Sec'][0]/data_oldsave_skylake['Real_Sec'],'--vg',label="Skylake (old save)")
plt.plot(data_nosave_skylake['Cores'],data_nosave_skylake['Real_Sec'][0]/data_nosave_skylake['Real_Sec'],'--ob',label="Skylake (no save)")
plt.plot(data_sapphire['Cores'],np.arange(1,len(data_sapphire)+1),'-r',label="Linear")
plt.xlabel('Cores')
plt.ylabel('Speed-up')
plt.legend()
#plt.savefig("plots/scalability_cpu_skylake.png", dpi=300, bbox_inches='tight')

plt.figure()
plt.plot(data_sapphire['Cores'],data_sapphire['Real_Sec'][0]/(data_sapphire['Real_Sec']*data_sapphire['Cores']),'--dc',label="Sapphire (no save)")
#plt.plot(data_skylake['Cores'],data_skylake['Real_Sec'][0]/data_skylake['Real_Sec'],'--pk',label="Skylake (save)")
#plt.plot(data_oldsave_skylake['Cores'],data_oldsave_skylake['Real_Sec'][0]/data_oldsave_skylake['Real_Sec'],'--vg',label="Skylake (old save)")
plt.plot(data_nosave_skylake['Cores'],data_nosave_skylake['Real_Sec'][0]/(data_nosave_skylake['Real_Sec']*data_nosave_skylake['Cores']),'--ob',label="Skylake (no save)")
plt.axhline(y=1, color='r', linestyle='-', label='Linear')
plt.xlabel('Cores')
plt.ylabel('Parallel efficiency')
plt.legend()

plt.figure()
plt.semilogy(data_sapphire['Cores'],data_sapphire['Real_Sec'],'--dc',label="Sapphire (no save)")
plt.semilogy(data_nosave_skylake['Cores'],data_nosave_skylake['Real_Sec'],'--ob',label="Skylake (no save)")
#plt.plot(data_sapphire['Cores'],np.arange(1,len(data_sapphire)+1),'-r',label="Linear")
plt.xlabel('Cores')
plt.ylabel('Wall time')
plt.legend()

plt.show()