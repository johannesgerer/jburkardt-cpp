#!/bin/bash
#
cp sphere_triangle_quad.hpp /$HOME/include
#
g++ -c -g -I /$HOME/include sphere_triangle_quad.cpp >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling sphere_triangle_quad.cpp"
  exit
fi
rm compiler.txt
#
mv sphere_triangle_quad.o ~/libcpp/$ARCH/sphere_triangle_quad.o
#
echo "Library installed as ~/libcpp/$ARCH/sphere_triangle_quad.o"
