#!/bin/sh
clear
export LC_NUMERIC=C

echo "Compiling..."
make clean
make
echo "Done!"
echo ""

decay_time=1
init_iter=1

for tau in $(seq 0.55 0.05 1.50)
do
    ./taylorgreen.o $tau $decay_time $init_iter
    echo "------------------------------"
done
echo ""