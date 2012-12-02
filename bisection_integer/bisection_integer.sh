#!/bin/bash
#
cp bisection_integer.hpp /$HOME/include
#
g++ -c -g bisection_integer.cpp >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling bisection_integer.cpp"
  exit
fi
rm compiler.txt
#
mv bisection_integer.o ~/libcpp/$ARCH/bisection_integer.o
#
echo "Library installed as ~/libcpp/$ARCH/bisection_integer.o"
