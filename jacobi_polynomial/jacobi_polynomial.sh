#!/bin/bash
#
cp jacobi_polynomial.hpp /$HOME/include
#
g++ -c -g -I /$HOME/include jacobi_polynomial.cpp >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling jacobi_polynomial.cpp"
  exit
fi
rm compiler.txt
#
mv jacobi_polynomial.o ~/libcpp/$ARCH/jacobi_polynomial.o
#
echo "Library installed as ~/libcpp/$ARCH/jacobi_polynomial.o"
