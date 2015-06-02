#!/bin/bash
#
cp hyperball_monte_carlo.hpp /$HOME/include
#
g++ -c -I/$HOME/include hyperball_monte_carlo.cpp
if [ $? -ne 0 ]; then
  echo "Errors compiling hyperball_monte_carlo.cpp"
  exit
fi
#
mv hyperball_monte_carlo.o ~/libcpp/$ARCH/hyperball_monte_carlo.o
#
echo "Library installed as ~/libcpp/$ARCH/hyperball_monte_carlo.o"
