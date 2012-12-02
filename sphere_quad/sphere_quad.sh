#!/bin/bash
#
cp sphere_quad.hpp /$HOME/include
#
g++ -c -g -I /$HOME/include sphere_quad.cpp >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling sphere_quad.cpp"
  exit
fi
rm compiler.txt
#
mv sphere_quad.o ~/libcpp/$ARCH/sphere_quad.o
#
echo "Library installed as ~/libcpp/$ARCH/sphere_quad.o"
