# ----------------------------------------------------------------------------------------------------------- #
#                                          IMPORT MODULES                                                     #
# ----------------------------------------------------------------------------------------------------------- #

import numpy as np
import pyvista as pv
from pathlib import Path
from scipy.io import savemat

# ----------------------------------------------------------------------------------------------------------- #
#                                           CALCULATIONS                                                      #
# ----------------------------------------------------------------------------------------------------------- #

sol_folder = "explicit"
for SPL in [130, 145]:
	p_a = 89 if SPL == 130 else 503
	for f_exc in [800, 1000, 1400, 2000]:
		print(f'SPL = {SPL} dB, f = {f_exc} Hz:')
		
		print('Reading...')
		foldername = Path('.')
		file_list = sorted(foldername.glob(f"{sol_folder}/fields_{f_exc}Hz_{p_a}Pa_*.vtk"))
		mesh = pv.read(file_list)
		print('Reading done!')

		print('Processing...')
		v_in = np.zeros([mesh[0].dimensions[1], mesh[0].dimensions[0]])
		v_out = np.zeros([mesh[0].dimensions[1], mesh[0].dimensions[0]])
		for i in range(len(file_list)):
			if i % 2 == 0:
				v_in += mesh[i]['velocity'][:, 0].reshape(mesh[0].dimensions[1], mesh[0].dimensions[0])
			else:
				v_out += mesh[i]['velocity'][:, 0].reshape(mesh[0].dimensions[1], mesh[0].dimensions[0])

		v_in /= len(file_list) / 2
		v_out /= len(file_list) / 2

		data_to_save = {
			'v_in': v_in.astype(np.float64),
			'v_out': v_out.astype(np.float64),
		}
		print('Processing done!')

		print('Saving...')
		savemat(f'v_in_out_{f_exc}Hz_{SPL}dB_{sol_folder}.mat', data_to_save)
		print('Saving done!')
		print(' ')
