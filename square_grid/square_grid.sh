#!/bin/bash
#
cp square_grid.hpp /$HOME/include
#
g++ -c -I/$HOME/include square_grid.cpp
if [ $? -ne 0 ]; then
  echo "Errors compiling square_grid.cpp"
  exit
fi
#
mv square_grid.o ~/libcpp/$ARCH/square_grid.o
#
echo "Library installed as ~/libcpp/$ARCH/square_grid.o"
