#!/bin/bash
#
#  Compile the program with G++
#
/usr/local/bin/g++ -fopenmp mxm.cpp -lm
#
#  Compile the program with ICPC.
#
#icpc -openmp -parallel mxm.cpp -lm
#
mv a.out mxm
#
#  Run with 1, 2, and 4 threads.
#
echo "Run with 1 thread."
export OMP_NUM_THREADS=1
./mxm > mxm_local_output.txt
#
echo "Run with 2 threads."
export OMP_NUM_THREADS=2
./mxm >> mxm_local_output.txt
#
echo "Run with 4 threads."
export OMP_NUM_THREADS=4
./mxm >> mxm_local_output.txt
#
#  Discard the executable file.
#
rm mxm
#
echo "Program output written to mxm_local_output.txt"
