#!/bin/bash
#
cp wedge_grid.hpp /$HOME/include
#
g++ -c -I/$HOME/include wedge_grid.cpp
if [ $? -ne 0 ]; then
  echo "Errors compiling wedge_grid.cpp"
  exit
fi
#
mv wedge_grid.o ~/libcpp/$ARCH/wedge_grid.o
#
echo "Library installed as ~/libcpp/$ARCH/wedge_grid.o"
