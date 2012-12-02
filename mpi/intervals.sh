#!/bin/bash
#
g++ -c intervals_mpi.cpp >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling intervals_mpi.cpp"
  exit
fi
rm compiler.txt
#
g++ intervals_mpi.o -lmpi
if [ $? -ne 0 ]; then
  echo "Errors loading intervals_mpi.o"
  exit
fi
rm intervals_mpi.o
#
mv a.out intervals
mpirun -v -np 4 ./intervals > intervals_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running intervals"
  exit
fi
rm intervals
#
echo "The intervals test problem has been executed."
