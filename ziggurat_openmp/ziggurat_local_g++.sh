#!/bin/bash
#
#  Compile the program with G++
#
/usr/local/bin/g++ -fopenmp ziggurat_openmp.cpp -lm
#
mv a.out ziggurat
#
#  Run with 1, 2, and 4 threads.
#
echo "Run with 1 thread."
export OMP_NUM_THREADS=1
./ziggurat > ziggurat_local_g++_output.txt
#
echo "Run with 2 threads."
export OMP_NUM_THREADS=2
./ziggurat >> ziggurat_local_g++_output.txt
#
echo "Run with 4 threads."
export OMP_NUM_THREADS=4
./ziggurat >> ziggurat_local_g++_output.txt
#
#  Discard the executable file.
#
rm ziggurat
#
echo "Program output written to ziggurat_local_g++_output.txt"
