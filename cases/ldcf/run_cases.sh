#!/bin/sh
clear
export LC_NUMERIC=C

echo "Compiling..."
make clean
make
echo "Done!"
echo ""

T_max=15
res=250
tau=0.9

for Re in 100 400 1000 #3200 5000 7500 10000
do
    echo ""
    ./ldcf.o $res $Re $tau $T_max
    echo "------------------------------"
done
echo ""