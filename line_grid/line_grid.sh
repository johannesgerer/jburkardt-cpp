#!/bin/bash
#
cp line_grid.hpp /$HOME/include
#
g++ -c -I/$HOME/include line_grid.cpp
if [ $? -ne 0 ]; then
  echo "Errors compiling line_grid.cpp"
  exit
fi
#
mv line_grid.o ~/libcpp/$ARCH/line_grid.o
#
echo "Library installed as ~/libcpp/$ARCH/line_grid.o"
