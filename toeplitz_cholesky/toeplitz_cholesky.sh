#!/bin/bash
#
cp toeplitz_cholesky.hpp /$HOME/include
#
g++ -c -g -I/$HOME/include toeplitz_cholesky.cpp >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling toeplitz_cholesky.cpp"
  exit
fi
rm compiler.txt
#
mv toeplitz_cholesky.o ~/libcpp/$ARCH/toeplitz_cholesky.o
#
echo "Library installed as ~/libcpp/$ARCH/toeplitz_cholesky.o"
