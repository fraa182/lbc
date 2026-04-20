#!/bin/sh
clear
export LC_NUMERIC=C

echo "Compiling..."
make clean
make
echo "Done!"
echo ""

for p_a in 89 503; do
    for f_exc in 800 1000 1400 2000; do
    echo ""
    echo "Pressure amplitude: $p_a Pa, frequency: $f_exc Hz"
    ./nit.o $p_a $f_exc
    echo "-------------------------------------------------"
    done
done
echo ""