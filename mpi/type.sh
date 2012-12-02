#!/bin/bash
#
g++ -c type_mpi.cpp >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling type_mpi.cpp"
  exit
fi
rm compiler.txt
#
g++ type_mpi.o -lmpi
if [ $? -ne 0 ]; then
  echo "Errors loading type_mpi.o"
  exit
fi
rm type_mpi.o
#
mv a.out type
mpirun -v -np 4 ./type > type_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running type"
  exit
fi
rm type
#
echo "The type test problem has been executed."
