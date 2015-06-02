#!/bin/bash
#
cp wedge_monte_carlo.hpp /$HOME/include
#
g++ -c -I/$HOME/include wedge_monte_carlo.cpp
if [ $? -ne 0 ]; then
  echo "Errors compiling wedge_monte_carlo.cpp"
  exit
fi
#
mv wedge_monte_carlo.o ~/libcpp/$ARCH/wedge_monte_carlo.o
#
echo "Library installed as ~/libcpp/$ARCH/wedge_monte_carlo.o"
