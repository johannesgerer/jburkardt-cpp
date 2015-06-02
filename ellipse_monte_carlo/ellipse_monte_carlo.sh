#!/bin/bash
#
cp ellipse_monte_carlo.hpp /$HOME/include
#
g++ -c -I /$HOME/include ellipse_monte_carlo.cpp >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling ellipse_monte_carlo.cpp"
  exit
fi
rm compiler.txt
#
mv ellipse_monte_carlo.o ~/libcpp/$ARCH/ellipse_monte_carlo.o
#
echo "Library installed as ~/libcpp/$ARCH/ellipse_monte_carlo.o"
