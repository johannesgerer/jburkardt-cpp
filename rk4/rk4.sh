#!/bin/bash
#
cp rk4.hpp /$HOME/include
#
g++ -c -g -I /$HOME/include rk4.cpp >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling rk4.cpp"
  exit
fi
rm compiler.txt
#
mv rk4.o ~/libcpp/$ARCH/rk4.o
#
echo "Library installed as ~/libcpp/$ARCH/rk4.o"
