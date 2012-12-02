#!/bin/bash
#
cp hermite_polynomial.hpp /$HOME/include
#
g++ -c -g -I /$HOME/include hermite_polynomial.cpp >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling hermite_polynomial.cpp"
  exit
fi
rm compiler.txt
#
mv hermite_polynomial.o ~/libcpp/$ARCH/hermite_polynomial.o
#
echo "Library installed as ~/libcpp/$ARCH/hermite_polynomial.o"
