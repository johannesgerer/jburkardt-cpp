#!/bin/bash
#
cp cube_monte_carlo.hpp /$HOME/include
#
g++ -c -g -I /$HOME/include cube_monte_carlo.cpp >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling cube_monte_carlo.cpp"
  exit
fi
rm compiler.txt
#
mv cube_monte_carlo.o ~/libcpp/$ARCH/cube_monte_carlo.o
#
echo "Library installed as ~/libcpp/$ARCH/cube_monte_carlo.o"
