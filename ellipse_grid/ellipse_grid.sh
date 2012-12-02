#!/bin/bash
#
cp ellipse_grid.hpp /$HOME/include
#
g++ -c -g ellipse_grid.cpp >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling ellipse_grid.cpp"
  exit
fi
rm compiler.txt
#
mv ellipse_grid.o ~/libcpp/$ARCH/ellipse_grid.o
#
echo "Library installed as ~/libcpp/$ARCH/ellipse_grid.o"
