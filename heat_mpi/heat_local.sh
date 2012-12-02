#!/bin/bash
#
mpic++ heat_mpi.cpp
#
if [ $? -ne 0 ]; then
  echo "Errors compiling heat_mpi.cpp"
  exit
fi
#
#  Rename the executable.
#
mv a.out heat
#
#  Ask MPI to use 8 processes to run your program.
#
mpirun -np 8 ./heat > heat_local_output.txt
#
#  Clean up.
#
rm heat
#
echo "Program output written to heat_local_output.txt"

