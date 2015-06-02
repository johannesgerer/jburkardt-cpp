#!/bin/bash
#
cp lebesgue.hpp /$HOME/include
#
g++ -c -g -I/$HOME/include lebesgue.cpp >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling lebesgue.cpp"
  exit
fi
rm compiler.txt
#
mv lebesgue.o ~/libcpp/$ARCH/lebesgue.o
#
echo "Library installed as ~/libcpp/$ARCH/lebesgue.o"
