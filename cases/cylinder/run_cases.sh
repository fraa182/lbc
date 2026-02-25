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
tau=0.85
echo "Grid resolution study at Re = $Re..."
for res in $(seq 10 10 50)
do
    echo ""
    #tau=$(awk -v r="$res" -v RE="$Re" 'BEGIN { val = (6 * 0.2 * r) / (sqrt(3) * RE) + 0.5; print (val > 0.9 ? 0.9 : val)}')
    ./cylinder.o $res $Re $tau $itermax 0
    ./cylinder.o $res $Re $tau $itermax 1
    echo "------------------------------"
done
echo ""

Re=100
tau=0.57
echo "Grid resolution study at Re = $Re..."
for res in $(seq 10 10 50)
do
    echo ""
    #tau=$(awk -v r="$res" -v RE="$Re" 'BEGIN { val = (6 * 0.2 * r) / (sqrt(3) * RE) + 0.5; print (val > 0.9 ? 0.9 : val)}')
    ./cylinder.o $res $Re $tau $itermax 0
    ./cylinder.o $res $Re $tau $itermax 1
    echo "------------------------------"
done
echo ""

echo "Done!"