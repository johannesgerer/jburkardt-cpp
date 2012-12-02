#!/bin/bash
#
#  Compile the program with G++
#
/usr/local/bin/g++ -fopenmp helmholtz.cpp -lm
#
#  Compile the program with ICPC.
#
#icpc -openmp -parallel helmholtz.cpp -lm
#
mv a.out helmholtz
#
#  Run with 1, 2, and 4 threads.
#
echo "Run with 1 thread."
export OMP_NUM_THREADS=1
./helmholtz > helmholtz_local_output.txt
#
echo "Run with 2 threads."
export OMP_NUM_THREADS=2
./helmholtz >> helmholtz_local_output.txt
#
echo "Run with 4 threads."
export OMP_NUM_THREADS=4
./helmholtz >> helmholtz_local_output.txt
#
#  Discard the executable file.
#
rm helmholtz
#
echo "Program output written to helmholtz_local_output.txt"
