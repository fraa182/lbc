#!/bin/sh
clear
export LC_NUMERIC=C

echo "Compiling..."
make clean
make
echo "Done!"
echo ""

gridstudy=$1
taustudy=$2

itermax=500000
tau_res=0.9
res_tau=11

if [ $gridstudy = true ]; then
    echo "Grid resolution study..."
    for res in $(seq 10 10 50)
    do
        ./couette.o $res $tau_res $itermax
        echo "------------------------------"
    done
    echo ""
fi

if [ $taustudy = true ]; then
    echo "Optimal relaxation time study..."
    for tau in $(seq 0.51 0.01 0.99)
    do
        ./couette.o $res_tau $tau $itermax
        echo "------------------------------"
    done

    for tau in $(seq 1 0.1 5)
    do
        ./couette.o $res_tau $tau $itermax
        echo "------------------------------"
    done
    echo ""
fi

echo "Done!"