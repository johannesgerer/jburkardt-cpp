#!/bin/bash
#
cp vandermonde_interp_2d.hpp /$HOME/include
#
g++ -c -g -I /$HOME/include vandermonde_interp_2d.cpp >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling vandermonde_interp_2d.cpp"
  exit
fi
rm compiler.txt
#
mv vandermonde_interp_2d.o ~/libcpp/$ARCH/vandermonde_interp_2d.o
#
echo "Library installed as ~/libcpp/$ARCH/vandermonde_interp_2d.o"
