#!/bin/bash
#
#  Compile the program with G++
#
/usr/local/bin/g++ -fopenmp quad2d_openmp.cpp -lm
#
mv a.out quad2d
#
#  Run with 1, 2, and 4 threads.
#
echo "Run with 1 thread."
export OMP_NUM_THREADS=1
./quad2d > quad2d_local_g++_output.txt
#
echo "Run with 2 threads."
export OMP_NUM_THREADS=2
./quad2d >> quad2d_local_g++_output.txt
#
echo "Run with 4 threads."
export OMP_NUM_THREADS=4
./quad2d >> quad2d_local_g++_output.txt
#
#  Discard the executable file.
#
rm quad2d
#
echo "Program output written to quad2d_local_g++_output.txt"
