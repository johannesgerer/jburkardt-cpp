#!/bin/bash
#
cp polygon_grid.hpp /$HOME/include
#
g++ -c -I/$HOME/include polygon_grid.cpp
if [ $? -ne 0 ]; then
  echo "Errors compiling polygon_grid.cpp"
  exit
fi
#
mv polygon_grid.o ~/libcpp/$ARCH/polygon_grid.o
#
echo "Library installed as ~/libcpp/$ARCH/polygon_grid.o"
