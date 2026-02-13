#!/bin/sh
clear
export LC_NUMERIC=C

echo "Compiling..."
make clean
make
echo "Done!"
echo ""

itermax=1000000
Re=150

echo "Grid resolution study..."
for res in $(seq 10 10 50)
do
    echo ""
    ./cylinder.o $res $Re $itermax
    echo "------------------------------"
done
echo ""

echo "Done!"