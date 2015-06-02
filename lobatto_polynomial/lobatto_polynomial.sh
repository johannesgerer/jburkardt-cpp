#!/bin/bash
#
cp lobatto_polynomial.hpp /$HOME/include
#
g++ -c -I/$HOME/include lobatto_polynomial.cpp
if [ $? -ne 0 ]; then
  echo "Errors compiling lobatto_polynomial.cpp"
  exit
fi
#
mv lobatto_polynomial.o ~/libcpp/$ARCH/lobatto_polynomial.o
#
echo "Library installed as ~/libcpp/$ARCH/lobatto_polynomial.o"
