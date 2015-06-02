#!/bin/bash
#
cp ou.hpp /$HOME/include
#
g++ -c -g -I/$HOME/include ou.cpp >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling ou.cpp"
  exit
fi
rm compiler.txt
#
mv ou.o ~/libcpp/$ARCH/ou.o
#
echo "Library installed as ~/libcpp/$ARCH/ou.o"
