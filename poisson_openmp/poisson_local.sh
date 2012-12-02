#!/bin/bash
#
#  Compile the program with G++
#
g++ -fopenmp poisson_openmp.cpp
#
#  Compile the program with ICPC.
#
#icpc -openmp -parallel poisson_openmp.cpp
#
mv a.out poisson
#
#  Run with 1, 2, 4 and 8 threads.
#
echo "Run with 1 thread."
export OMP_NUM_THREADS=1
./poisson > poisson_local_output.txt
#
echo "Run with 2 threads."
export OMP_NUM_THREADS=2
./poisson >> poisson_local_output.txt
#
echo "Run with 4 threads."
export OMP_NUM_THREADS=4
./poisson >> poisson_local_output.txt
#
echo "Run with 8 threads."
export OMP_NUM_THREADS=8
./poisson >> poisson_local_output.txt
#
#  Discard the executable file.
#
rm poisson
#
echo "Program output written to poisson_local_output.txt"
