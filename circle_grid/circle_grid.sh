#!/bin/bash
#
cp circle_grid.hpp /$HOME/include
#
g++ -c -g circle_grid.cpp >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling circle_grid.cpp"
  exit
fi
rm compiler.txt
#
mv circle_grid.o ~/libcpp/$ARCH/circle_grid.o
#
echo "Library installed as ~/libcpp/$ARCH/circle_grid.o"
