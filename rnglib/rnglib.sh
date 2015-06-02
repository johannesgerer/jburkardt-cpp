#!/bin/bash
#
cp rnglib.hpp /$HOME/include
#
g++ -c -g -I/$HOME/include rnglib.cpp >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling rnglib.cpp"
  exit
fi
rm compiler.txt
#
mv rnglib.o ~/libcpp/$ARCH/rnglib.o
#
echo "Library installed as ~/libcpp/$ARCH/rnglib.o"
