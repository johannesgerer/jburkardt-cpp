#!/bin/bash
#
mpiCC -c -g communicator_mpi.cpp >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling communicator_mpi.cpp"
  exit
fi
rm compiler.txt
#
mpiCC communicator_mpi.o -lm
if [ $? -ne 0 ]; then
  echo "Errors linking and loading communicator_mpi.o."
  exit
fi
#
rm communicator_mpi.o
#
mv a.out communicator
mpirun -np 4 ./communicator > communicator_local_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running communicator"
  exit
fi
rm communicator
#
echo "Program output written to communicator_local_output.txt"
