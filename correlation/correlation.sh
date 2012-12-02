#!/bin/bash
#
cp correlation.hpp /$HOME/include
#
g++ -c -g -I/$HOME/include correlation.cpp >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling correlation.cpp"
  exit
fi
rm compiler.txt
#
mv correlation.o ~/libcpp/$ARCH/correlation.o
#
echo "Library installed as ~/libcpp/$ARCH/correlation.o"
