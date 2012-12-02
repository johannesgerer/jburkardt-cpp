#!/bin/bash
#
cp sparse_grid_gl.hpp /$HOME/include
#
g++ -c -g -I /$HOME/include sparse_grid_gl.cpp >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling sparse_grid_gl.cpp"
  exit
fi
rm compiler.txt
#
mv sparse_grid_gl.o ~/libcpp/$ARCH/sparse_grid_gl.o
#
echo "Library installed as ~/libcpp/$ARCH/sparse_grid_gl.o"
