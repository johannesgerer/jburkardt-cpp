#!/bin/bash
#
g++ -c search_mpi.cpp >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling search_mpi.cpp"
  exit
fi
rm compiler.txt
#
g++ search_mpi.o -lmpi
if [ $? -ne 0 ]; then
  echo "Errors loading search_mpi.o"
  exit
fi
rm search_mpi.o
#
mv a.out search
mpirun -v -np 4 ./search > search_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running search"
  exit
fi
rm search
#
echo "The search test problem has been executed."
