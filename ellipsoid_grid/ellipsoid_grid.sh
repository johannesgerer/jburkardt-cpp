#!/bin/bash
#
cp ellipsoid_grid.hpp /$HOME/include
#
g++ -c -g ellipsoid_grid.cpp >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling ellipsoid_grid.cpp"
  exit
fi
rm compiler.txt
#
mv ellipsoid_grid.o ~/libcpp/$ARCH/ellipsoid_grid.o
#
echo "Library installed as ~/libcpp/$ARCH/ellipsoid_grid.o"
