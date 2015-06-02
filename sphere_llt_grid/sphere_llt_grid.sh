#!/bin/bash
#
cp sphere_llt_grid.hpp /$HOME/include
#
g++ -c -I/$HOME/include sphere_llt_grid.cpp
if [ $? -ne 0 ]; then
  echo "Errors compiling sphere_llt_grid.cpp"
  exit
fi
#
mv sphere_llt_grid.o ~/libcpp/$ARCH/sphere_llt_grid.o
#
echo "Library installed as ~/libcpp/$ARCH/sphere_llt_grid.o"
