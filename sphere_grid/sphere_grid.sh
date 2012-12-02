#!/bin/bash
#
cp sphere_grid.hpp /$HOME/include
#
g++ -c -g -I /$HOME/include sphere_grid.cpp >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling sphere_grid.cpp"
  exit
fi
rm compiler.txt
#
mv sphere_grid.o ~/libcpp/$ARCH/sphere_grid.o
#
echo "Library installed as ~/libcpp/$ARCH/sphere_grid.o"
