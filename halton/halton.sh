#!/bin/bash
#
cp halton.hpp /$HOME/include
#
g++ -c -I /$HOME/include halton.cpp
if [ $? -ne 0 ]; then
  echo "Errors compiling halton.cpp"
  exit
fi
#
mv halton.o ~/libcpp/$ARCH/halton.o
#
echo "Library installed as ~/libcpp/$ARCH/halton.o"
