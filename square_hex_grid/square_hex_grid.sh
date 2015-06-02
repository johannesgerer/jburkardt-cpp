#!/bin/bash
#
cp square_hex_grid.hpp /$HOME/include
#
g++ -c -I /$HOME/include square_hex_grid.cpp
if [ $? -ne 0 ]; then
  echo "Errors compiling square_hex_grid.cpp"
  exit
fi
#
mv square_hex_grid.o ~/libcpp/$ARCH/square_hex_grid.o
#
echo "Library installed as ~/libcpp/$ARCH/square_hex_grid.o"
