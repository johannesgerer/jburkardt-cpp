#!/bin/bash
#
cp hypersphere_monte_carlo.hpp /$HOME/include
#
g++ -c -I /$HOME/include hypersphere_monte_carlo.cpp
if [ $? -ne 0 ]; then
  echo "Errors compiling hypersphere_monte_carlo.cpp"
  exit
fi
#
mv hypersphere_monte_carlo.o ~/libcpp/$ARCH/hypersphere_monte_carlo.o
#
echo "Library installed as ~/libcpp/$ARCH/hypersphere_monte_carlo.o"
