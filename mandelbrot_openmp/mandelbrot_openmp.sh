#!/bin/bash
#
g++ -fopenmp mandelbrot_openmp.cpp
mv a.out mandelbrot_openmp
#
echo "Run with 8 threads."
export OMP_NUM_THREADS=8
./mandelbrot_openmp > mandelbrot_openmp_output.txt
#
#  Discard the executable file.
#
rm mandelbrot_openmp
#
echo "Program output written to mandelbrot_openmp_output.txt."
