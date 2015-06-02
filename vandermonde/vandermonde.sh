#!/bin/bash
#
cp vandermonde.hpp /$HOME/include
#
g++ -c -g -I/$HOME/include vandermonde.cpp >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling vandermonde.cpp"
  exit
fi
rm compiler.txt
#
mv vandermonde.o ~/libcpp/$ARCH/vandermonde.o
#
echo "Library installed as ~/libcpp/$ARCH/vandermonde.o"
