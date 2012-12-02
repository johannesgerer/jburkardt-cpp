#!/bin/bash
#
#  Compile the program with ICPC.
#
icpc -openmp -parallel schedule_openmp.cpp
#
mv a.out schedule
#
#  Run with 2 threads.
#
echo "Run with 2 threads."
export OMP_NUM_THREADS=2
./schedule > schedule_local_icpc_output.txt
#
#  Discard the executable file.
#
rm schedule
#
echo "Program output written to schedule_local_icpc_output.txt"
