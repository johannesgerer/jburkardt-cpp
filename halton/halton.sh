#!/bin/bash
#
cp halton.hpp /$HOME/include
#
g++ -c -g -I /$HOME/include halton.cpp >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling halton.cpp"
  exit
fi
rm compiler.txt
#
mv halton.o ~/libcpp/$ARCH/halton.o
#
echo "Library installed as ~/libcpp/$ARCH/halton.o"
