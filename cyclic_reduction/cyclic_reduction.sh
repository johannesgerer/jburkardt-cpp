#!/bin/bash
#
cp cyclic_reduction.hpp /$HOME/include
#
g++ -c -g -I /$HOME/include cyclic_reduction.cpp >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling cyclic_reduction.cpp"
  exit
fi
rm compiler.txt
#
mv cyclic_reduction.o ~/libcpp/$ARCH/cyclic_reduction.o
#
echo "Library installed as ~/libcpp/$ARCH/cyclic_reduction.o"
