#!/bin/sh
clear
export LC_NUMERIC=C

echo "Compiling..."
make clean
make
echo "Done!"
echo ""

m_fix=5
A_fix=0.0002
tau_fix=0.9

for tau in $(seq 0.55 0.05 1.05); do
    clear
    echo "tau = $tau"
    ./wave.o $tau $m_fix $A_fix
done

for m in $(seq 5 5 20); do
    clear
    echo "m = $m"
    ./wave.o $tau_fix $m $A_fix
done

for A in 0.0002 0.002 0.02 0.2; do
    clear
    echo "A = $A"
    ./wave.o $tau_fix $m_fix $A
done