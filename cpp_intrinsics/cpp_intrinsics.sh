#!/bin/bash
#
g++ -c cpp_intrinsics.cpp >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling cpp_intrinsics.cpp"
  exit
fi
rm compiler.txt
#
g++ cpp_intrinsics.o -lm
if [ $? -ne 0 ]; then
  echo "Errors linking cpp_intrinsics.o."
  exit
fi
#
rm cpp_intrinsics.o
#
mv a.out cpp_intrinsics
./cpp_intrinsics > cpp_intrinsics_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running cpp_intrinsics."
  exit
fi
rm cpp_intrinsics
#
echo "Program output written to cpp_intrinsics_output.txt"
