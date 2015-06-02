#!/bin/bash
#
cp simplex_grid.hpp /$HOME/include
#
g++ -c -I/$HOME/include simplex_grid.cpp
if [ $? -ne 0 ]; then
  echo "Errors compiling simplex_grid.cpp"
  exit
fi
#
mv simplex_grid.o ~/libcpp/$ARCH/simplex_grid.o
#
echo "Library installed as ~/libcpp/$ARCH/simplex_grid.o"
