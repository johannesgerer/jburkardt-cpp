#!/bin/bash
#
cp vandermonde_approx_2d.hpp /$HOME/include
#
g++ -c -g -I/$HOME/include vandermonde_approx_2d.cpp >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling vandermonde_approx_2d.cpp"
  exit
fi
rm compiler.txt
#
mv vandermonde_approx_2d.o ~/libcpp/$ARCH/vandermonde_approx_2d.o
#
echo "Library installed as ~/libcpp/$ARCH/vandermonde_approx_2d.o"
