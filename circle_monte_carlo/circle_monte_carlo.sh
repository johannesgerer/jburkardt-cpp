#!/bin/bash
#
cp circle_monte_carlo.hpp /$HOME/include
#
g++ -c -g -I /$HOME/include circle_monte_carlo.cpp >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling circle_monte_carlo.cpp"
  exit
fi
rm compiler.txt
#
mv circle_monte_carlo.o ~/libcpp/$ARCH/circle_monte_carlo.o
#
echo "Library installed as ~/libcpp/$ARCH/circle_monte_carlo.o"
