#!/bin/bash
#
g++ -c not_allocated_arrays.cpp >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling not_allocated_arrays.cpp"
  exit
fi
rm compiler.txt
#
g++ not_allocated_arrays.o
if [ $? -ne 0 ]; then
  echo "Errors linking not_allocated_arrays.o."
  exit
fi
#
rm not_allocated_arrays.o
#
mv a.out not_allocated_arrays
./not_allocated_arrays > not_allocated_arrays_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running not_allocated_arrays."
  exit
fi
rm not_allocated_arrays
#
echo "Program output written to not_allocated_arrays_output.txt"
