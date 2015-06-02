#!/bin/bash
#
cp hypercube_grid.hpp /$HOME/include
#
g++ -c -I/$HOME/include hypercube_grid.cpp
if [ $? -ne 0 ]; then
  echo "Errors compiling hypercube_grid.cpp"
  exit
fi
#
mv hypercube_grid.o ~/libcpp/$ARCH/hypercube_grid.o
#
echo "Library installed as ~/libcpp/$ARCH/hypercube_grid.o"
