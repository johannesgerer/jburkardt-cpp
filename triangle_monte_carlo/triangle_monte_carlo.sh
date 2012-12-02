#!/bin/bash
#
cp triangle_monte_carlo.hpp /$HOME/include
#
g++ -c -g -I /$HOME/include triangle_monte_carlo.cpp >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling triangle_monte_carlo.cpp"
  exit
fi
rm compiler.txt
#
mv triangle_monte_carlo.o ~/libcpp/$ARCH/triangle_monte_carlo.o
#
echo "Library installed as ~/libcpp/$ARCH/triangle_monte_carlo.o"
