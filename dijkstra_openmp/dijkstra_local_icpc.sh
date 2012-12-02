#!/bin/bash
#
#  Compile the program with ICPC.
#
icpc -openmp -parallel dijkstra_openmp.cpp -lm
#
mv a.out hello
#
#  Run with 1, 2, and 4 threads.
#
echo "Run with 1 thread."
export OMP_NUM_THREADS=1
./hello > dijkstra_local_icpc_output.txt
#
echo "Run with 2 threads."
export OMP_NUM_THREADS=2
./hello >> dijkstra_local_icpc_output.txt
#
echo "Run with 4 threads."
export OMP_NUM_THREADS=4
./hello >> dijkstra_local_icpc_output.txt
#
#  Discard the executable file.
#
rm hello
#
echo "Program output written to dijkstra_local_icpc_output.txt"
