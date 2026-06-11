#!/bin/sh

#SBATCH --job-name=nit_orifices
#SBATCH --partition=cpu_skylake_ext
#SBATCH --time=5-00:00:00
#SBATCH --nodes=1
#SBATCH --cpus-per-task=32

clear
export LC_NUMERIC=C

echo "Compiling..."
make clean
make
echo "Done!"
echo ""

export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK

for p_a in 89 503; do
    for f_exc in 800 1000 1400 2000; do
    	echo ""
    	echo "Pressure amplitude: $p_a Pa, frequency: $f_exc Hz"
    	./nit_explicit.o $p_a $f_exc
    	echo "-------------------------------------------------"
    done
done
echo ""

mv sol explicit

echo "Done!"
echo ""
