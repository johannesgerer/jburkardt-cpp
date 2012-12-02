#!/bin/bash
#
cp bins.hpp /$HOME/include
#
g++ -c -g -I /$HOME/include bins.cpp >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling bins.cpp"
  exit
fi
rm compiler.txt
#
mv bins.o ~/libcpp/$ARCH/bins.o
#
echo "Library installed as ~/libcpp/$ARCH/bins.o"
