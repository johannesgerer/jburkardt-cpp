#!/bin/bash
#
cp fd1d_heat_explicit.hpp /$HOME/include
#
g++ -c -g -I /$HOME/include fd1d_heat_explicit.cpp >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling fd1d_heat_explicit.cpp"
  exit
fi
rm compiler.txt
#
mv fd1d_heat_explicit.o ~/libcpp/$ARCH/fd1d_heat_explicit.o
#
echo "Library installed as ~/libcpp/$ARCH/fd1d_heat_explicit.o"
