#!/bin/bash
#
g++ -c g++_intrinsics.cpp >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling g++_intrinsics.cpp"
  exit
fi
rm compiler.txt
#
g++ g++_intrinsics.o
if [ $? -ne 0 ]; then
  echo "Errors linking g++_intrinsics.o."
  exit
fi
#
rm g++_intrinsics.o
#
mv a.out g++_intrinsics
./g++_intrinsics > g++_intrinsics_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running g++_intrinsics."
  exit
fi
rm g++_intrinsics
#
echo "Program output written to g++_intrinsics_output.txt"
