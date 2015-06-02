#!/bin/bash
#
cp blas0.hpp /$HOME/include
#
g++ -c -O2 blas0.cpp >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling blas0.cpp."
  exit
fi
rm compiler.txt
#
mv blas0.o ~/libcpp/$ARCH/blas0.o
#
echo "Library installed as ~/libcpp/$ARCH/blas0.o"
