#!/bin/bash
#
cp hermite_polynomial.hpp /$HOME/include
#
g++ -c -I /$HOME/include hermite_polynomial.cpp
if [ $? -ne 0 ]; then
  echo "Errors compiling hermite_polynomial.cpp"
  exit
fi
#
mv hermite_polynomial.o ~/libcpp/$ARCH/hermite_polynomial.o
#
echo "Library installed as ~/libcpp/$ARCH/hermite_polynomial.o"
