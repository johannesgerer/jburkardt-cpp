#!/bin/bash
#
cp line_monte_carlo.hpp /$HOME/include
#
g++ -c -g -I /$HOME/include line_monte_carlo.cpp >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling line_monte_carlo.cpp"
  exit
fi
rm compiler.txt
#
mv line_monte_carlo.o ~/libcpp/$ARCH/line_monte_carlo.o
#
echo "Library installed as ~/libcpp/$ARCH/line_monte_carlo.o"
