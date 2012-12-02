#!/bin/bash
#
#  Compile the program with G++
#
/usr/local/bin/g++ -fopenmp compute_pi.cpp -lm
#
#  Compile the program with ICPC.
#
#icpc -openmp -parallel compute_pi.cpp -lm
#
mv a.out compute_pi
#
#  Run with 1, 2, and 4 threads.
#
echo "Run with 1 thread."
export OMP_NUM_THREADS=1
./compute_pi > compute_pi_local_output.txt
#
echo "Run with 2 threads."
export OMP_NUM_THREADS=2
./compute_pi >> compute_pi_local_output.txt
#
echo "Run with 4 threads."
export OMP_NUM_THREADS=4
./compute_pi >> compute_pi_local_output.txt
#
#  Discard the executable file.
#
rm compute_pi
#
echo "Program output written to compute_pi_local_output.txt"
