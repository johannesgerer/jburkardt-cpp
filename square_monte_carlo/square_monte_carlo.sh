#!/bin/bash
#
cp square_monte_carlo.hpp /$HOME/include
#
g++ -c -g -I /$HOME/include square_monte_carlo.cpp >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling square_monte_carlo.cpp"
  exit
fi
rm compiler.txt
#
mv square_monte_carlo.o ~/libcpp/$ARCH/square_monte_carlo.o
#
echo "Library installed as ~/libcpp/$ARCH/square_monte_carlo.o"
