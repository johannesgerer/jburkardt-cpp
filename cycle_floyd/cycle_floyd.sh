#!/bin/bash
#
cp cycle_floyd.hpp /$HOME/include
#
g++ -c -g -I /$HOME/include cycle_floyd.cpp >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling cycle_floyd.cpp"
  exit
fi
rm compiler.txt
#
mv cycle_floyd.o ~/libcpp/$ARCH/cycle_floyd.o
#
echo "Library installed as ~/libcpp/$ARCH/cycle_floyd.o"
