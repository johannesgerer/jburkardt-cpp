#!/bin/bash
#
cp sparse_grid_cc.hpp /$HOME/include
#
g++ -c -g -I /$HOME/include sparse_grid_cc.cpp >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling sparse_grid_cc.cpp"
  exit
fi
rm compiler.txt
#
mv sparse_grid_cc.o ~/libcpp/$ARCH/sparse_grid_cc.o
#
echo "Library installed as ~/libcpp/$ARCH/sparse_grid_cc.o"
