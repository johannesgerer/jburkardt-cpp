#!/bin/bash
#
cp polynomial.hpp /$HOME/include
#
g++ -c -I/$HOME/include polynomial.cpp
if [ $? -ne 0 ]; then
  echo "Errors compiling polynomial.cpp"
  exit
fi
#
mv polynomial.o ~/libcpp/$ARCH/polynomial.o
#
echo "Library installed as ~/libcpp/$ARCH/polynomial.o"
