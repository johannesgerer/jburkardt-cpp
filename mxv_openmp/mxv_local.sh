#!/bin/bash
#
#  Compile the program with G++.
#
/usr/local/bin/g++ -fopenmp mxv_openmp.cpp
#
#  Compile the program with ICPC.
#
#icpc -openmp -parallel mxv_openmp.cpp
#
mv a.out mxv
#
#  Run with 1, 2, 4 and 8 threads.
#
echo "Run with 1 thread."
export OMP_NUM_THREADS=1
./mxv > mxv_local_output.txt
#
echo "Run with 2 threads."
export OMP_NUM_THREADS=2
./mxv >> mxv_local_output.txt
#
echo "Run with 4 threads."
export OMP_NUM_THREADS=4
./mxv >> mxv_local_output.txt
#
echo "Run with 8 threads."
export OMP_NUM_THREADS=8
./mxv >> mxv_local_output.txt
#
#  Discard the executable file.
#
rm mxv
#
echo "Program output written to mxv_local_output.txt"
