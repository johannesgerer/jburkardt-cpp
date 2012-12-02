#!/bin/bash
#
cp hex_grid.hpp /$HOME/include
#
g++ -c -g -I /$HOME/include hex_grid.cpp >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling hex_grid.cpp"
  exit
fi
rm compiler.txt
#
mv hex_grid.o ~/libcpp/$ARCH/hex_grid.o
#
echo "Library installed as ~/libcpp/$ARCH/hex_grid.o"
