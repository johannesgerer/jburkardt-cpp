#!/bin/bash
#
cp sphere_fibonacci_grid.hpp /$HOME/include
#
g++ -c -I/$HOME/include sphere_fibonacci_grid.cpp
if [ $? -ne 0 ]; then
  echo "Errors compiling sphere_fibonacci_grid.cpp"
  exit
fi
#
mv sphere_fibonacci_grid.o ~/libcpp/$ARCH/sphere_fibonacci_grid.o
#
echo "Library installed as ~/libcpp/$ARCH/sphere_fibonacci_grid.o"
