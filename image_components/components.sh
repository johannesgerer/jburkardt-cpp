#!/bin/bash
#
cp components.hpp /$HOME/include
#
g++ -c -g -I /$HOME/include components.cpp >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling components.cpp"
  exit
fi
rm compiler.txt
#
mv components.o ~/libcpp/$ARCH/components.o
#
echo "Library installed as ~/libcpp/$ARCH/components.o"
