#!/bin/bash
#
g++ -c monte_carlo_mpi.cpp >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling monte_carlo_mpi.cpp"
  exit
fi
rm compiler.txt
#
g++ monte_carlo_mpi.o -lmpi
if [ $? -ne 0 ]; then
  echo "Errors loading monte_carlo_mpi.o"
  exit
fi
rm monte_carlo_mpi.o
#
mv a.out monte_carlo
mpirun -v -np 4 ./monte_carlo > monte_carlo_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running monte_carlo"
  exit
fi
rm monte_carlo
#
echo "The monte_carlo test problem has been executed."
