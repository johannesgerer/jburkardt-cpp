#!/bin/bash
#
cp ziggurat.hpp /$HOME/include
#
g++ -c -g -I /$HOME/include ziggurat.cpp >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling ziggurat.cpp."
  exit
fi
rm compiler.txt
#
mv ziggurat.o ~/libcpp/$ARCH/ziggurat.o
#
echo "Library installed as ~/libcpp/$ARCH/ziggurat.o"
