#!/bin/bash
#
cp sphere_llq_grid.hpp /$HOME/include
#
g++ -c -I/$HOME/include sphere_llq_grid.cpp
if [ $? -ne 0 ]; then
  echo "Errors compiling sphere_llq_grid.cpp"
  exit
fi
#
mv sphere_llq_grid.o ~/libcpp/$ARCH/sphere_llq_grid.o
#
echo "Library installed as ~/libcpp/$ARCH/sphere_llq_grid.o"
