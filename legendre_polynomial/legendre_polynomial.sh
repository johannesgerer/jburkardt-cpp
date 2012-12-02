#!/bin/bash
#
cp legendre_polynomial.hpp /$HOME/include
#
g++ -c -g -I /$HOME/include legendre_polynomial.cpp >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling legendre_polynomial.cpp"
  exit
fi
rm compiler.txt
#
mv legendre_polynomial.o ~/libcpp/$ARCH/legendre_polynomial.o
#
echo "Library installed as ~/libcpp/$ARCH/legendre_polynomial.o"
