#!/bin/bash
#
mpic++ ring_mpi.cpp
#
if [ $? -ne 0 ]; then
  echo "Errors compiling ring_mpi.cpp"
  exit
fi
#
#  Rename the executable.
#
mv a.out ring
#
#  Ask MPI to use 8 processes to run your program.
#
mpirun -np 8 ./ring > ring_local_output.txt
#
#  Clean up.
#
rm ring
#
echo "Program output written to ring_local_output.txt"

