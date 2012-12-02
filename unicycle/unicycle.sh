#!/bin/bash
#
cp unicycle.hpp /$HOME/include
#
g++ -c -g -I /$HOME/include unicycle.cpp >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling unicycle.cpp"
  exit
fi
rm compiler.txt
#
mv unicycle.o ~/libcpp/$ARCH/unicycle.o
#
echo "Library installed as ~/libcpp/$ARCH/unicycle.o"
