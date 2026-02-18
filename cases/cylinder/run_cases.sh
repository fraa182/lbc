#!/bin/sh
clear
export LC_NUMERIC=C

echo "Compiling..."
make clean
make
echo "Done!"
echo ""

itermax=100000

Re=20
echo "Grid resolution study at Re = $Re..."
for res in $(seq 10 10 50)
do
    echo ""
    tau=$(awk -v r="$res" -v RE="$Re" 'BEGIN { print (6 * 0.2 * r) / (sqrt(3) * RE) + 0.5 }')
    ./cylinder.o $res $Re $tau $itermax
    echo "------------------------------"
done
echo ""

Re=100
echo "Grid resolution study at Re = $Re..."
for res in $(seq 10 10 50)
do
    echo ""
    tau=$(awk -v r="$res" -v RE="$Re" 'BEGIN { print (6 * 0.2 * r) / (sqrt(3) * RE) + 0.5 }')
    ./cylinder.o $res $Re $tau $itermax
    echo "------------------------------"
done
echo ""

echo "Done!"