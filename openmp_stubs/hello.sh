#!/bin/bash
#
#  Copy the STUBS include file temporarily.
#
cp openmp_stubs.hpp omp.h
#
#  Compile the program.
#
g++ ../openmp/hello.cpp -I. -L$HOME/libcpp/$ARCH -lopenmp_stubs
mv a.out hello
rm omp.h
#
#  Run the program.
#
./hello > hello_output.txt
rm hello
#
echo "Normal end of execution."
