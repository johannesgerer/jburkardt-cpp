#!/bin/bash
#
cp pyramid_grid.hpp /$HOME/include
#
g++ -c -I/$HOME/include pyramid_grid.cpp
if [ $? -ne 0 ]; then
  echo "Errors compiling pyramid_grid.cpp"
  exit
fi
#
mv pyramid_grid.o ~/libcpp/$ARCH/pyramid_grid.o
#
echo "Library installed as ~/libcpp/$ARCH/pyramid_grid.o"
