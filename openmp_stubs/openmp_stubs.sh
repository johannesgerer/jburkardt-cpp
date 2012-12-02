#!/bin/bash
#
cp openmp_stubs.hpp /$HOME/include
#
g++ -c -g openmp_stubs.cpp >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling openmp_stubs.cpp"
  exit
fi
rm compiler.txt
#
ar qc libopenmp_stubs.a openmp_stubs.o
mv libopenmp_stubs.a ~/libcpp/$ARCH/libopenmp_stubs.a
rm openmp_stubs.o
#
echo "A new version of openmp_stubs has been created."
