#!/bin/bash
#
g++ -fopenmp random_openmp.cpp
mv a.out random_openmp
#
echo "Run with 8 threads."
export OMP_NUM_THREADS=8
./random_openmp > random_openmp_output.txt
#
#  Discard the executable file.
#
rm random_openmp
#
echo "Program output written to random_openmp_output.txt."
