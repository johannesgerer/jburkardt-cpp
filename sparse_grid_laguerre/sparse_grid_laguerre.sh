#!/bin/bash
#
cp sparse_grid_laguerre.hpp /$HOME/include
#
g++ -c -g -I /$HOME/include sparse_grid_laguerre.cpp >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling sparse_grid_laguerre.cpp"
  exit
fi
rm compiler.txt
#
mv sparse_grid_laguerre.o ~/libcpp/$ARCH/sparse_grid_laguerre.o
#
echo "Library installed as ~/libcpp/$ARCH/sparse_grid_laguerre.o"
