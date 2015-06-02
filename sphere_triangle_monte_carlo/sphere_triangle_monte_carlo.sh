#!/bin/bash
#
cp sphere_triangle_monte_carlo.hpp /$HOME/include
#
g++ -c -I/$HOME/include sphere_triangle_monte_carlo.cpp >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling sphere_triangle_monte_carlo.cpp"
  exit
fi
rm compiler.txt
#
mv sphere_triangle_monte_carlo.o ~/libcpp/$ARCH/sphere_triangle_monte_carlo.o
#
echo "Library installed as ~/libcpp/$ARCH/sphere_triangle_monte_carlo.o"
