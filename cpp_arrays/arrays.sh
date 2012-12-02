#!/bin/bash
#
g++ -c arrays.cpp >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling arrays.cpp"
  exit
fi
rm compiler.txt
#
g++ arrays.o
if [ $? -ne 0 ]; then
  echo "Errors linking arrays.o."
  exit
fi
#
rm arrays.o
#
mv a.out arrays
./arrays > arrays_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running arrays."
  exit
fi
rm arrays
#
echo "Program output written to arrays_output.txt"
