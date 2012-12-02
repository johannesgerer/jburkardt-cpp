#!/bin/bash
#
cp cell.hpp /$HOME/include
#
g++ -c -g -I /$HOME/include cell.cpp >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling cell.cpp"
  exit
fi
rm compiler.txt
#
mv cell.o ~/libcpp/$ARCH/cell.o
#
echo "Library installed as ~/libcpp/$ARCH/cell.o"
