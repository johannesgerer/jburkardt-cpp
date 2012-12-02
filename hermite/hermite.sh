#!/bin/bash
#
cp hermite.hpp /$HOME/include
#
g++ -c -g -I /$HOME/include hermite.cpp >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling hermite.cpp"
  exit
fi
rm compiler.txt
#
mv hermite.o ~/libcpp/$ARCH/hermite.o
#
echo "Library installed as ~/libcpp/$ARCH/hermite.o"
