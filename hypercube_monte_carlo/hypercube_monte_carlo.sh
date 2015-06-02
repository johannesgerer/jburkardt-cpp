#!/bin/bash
#
cp hypercube_monte_carlo.hpp /$HOME/include
#
g++ -c -I /$HOME/include hypercube_monte_carlo.cpp
if [ $? -ne 0 ]; then
  echo "Errors compiling hypercube_monte_carlo.cpp"
  exit
fi
#
mv hypercube_monte_carlo.o ~/libcpp/$ARCH/hypercube_monte_carlo.o
#
echo "Library installed as ~/libcpp/$ARCH/hypercube_monte_carlo.o"
