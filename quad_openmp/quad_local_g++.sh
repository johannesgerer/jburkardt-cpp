#!/bin/bash
#
#  Compile the program with G++
#
/usr/local/bin/g++ -fopenmp quad_openmp.cpp -lm
#
mv a.out quad
#
#  Run with 1, 2, and 4 threads.
#
echo "Run with 1 thread."
export OMP_NUM_THREADS=1
./quad > quad_local_g++_output.txt
#
echo "Run with 2 threads."
export OMP_NUM_THREADS=2
./quad >> quad_local_g++_output.txt
#
echo "Run with 4 threads."
export OMP_NUM_THREADS=4
./quad >> quad_local_g++_output.txt
#
#  Discard the executable file.
#
rm quad
#
echo "Program output written to quad_local_g++_output.txt"
