#!/bin/bash
#
cp triangle_grid.hpp /$HOME/include
#
g++ -c -g -I /$HOME/include triangle_grid.cpp >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling triangle_grid.cpp"
  exit
fi
rm compiler.txt
#
mv triangle_grid.o ~/libcpp/$ARCH/triangle_grid.o
#
echo "Library installed as ~/libcpp/$ARCH/triangle_grid.o"
