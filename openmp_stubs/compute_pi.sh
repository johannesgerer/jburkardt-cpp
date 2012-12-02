#!/bin/bash
#
#  Copy the STUBS include file temporarily.
#
cp openmp_stubs.hpp omp.h
#
#  Compile the program.
#
g++ ../openmp/compute_pi.cpp -I. -L$HOME/libcpp/$ARCH -lopenmp_stubs
mv a.out compute_pi
rm omp.h
#
#  Run the program.
#
./compute_pi > compute_pi_output.txt
rm compute_pi
#
echo "Normal end of execution."
