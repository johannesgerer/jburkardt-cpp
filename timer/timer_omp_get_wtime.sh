#!/bin/bash
#
/usr/local/bin/g++ -fopenmp timer_omp_get_wtime.cpp >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling timer_omp_get_wtime.cpp"
  exit
fi
rm compiler.txt
#
mv a.out timer_omp_get_wtime
#
#  Run the program.
#
./timer_omp_get_wtime > timer_omp_get_wtime_output.txt
rm timer_omp_get_wtime
#
echo "Program output written to timer_omp_get_wtime_output.txt"
