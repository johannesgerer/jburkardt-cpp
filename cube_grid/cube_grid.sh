#!/bin/bash
#
cp cube_grid.hpp /$HOME/include
#
g++ -c -I/$HOME/include cube_grid.cpp
if [ $? -ne 0 ]; then
  echo "Errors compiling cube_grid.cpp"
  exit
fi
#
mv cube_grid.o ~/libcpp/$ARCH/cube_grid.o
#
echo "Library installed as ~/libcpp/$ARCH/cube_grid.o"
