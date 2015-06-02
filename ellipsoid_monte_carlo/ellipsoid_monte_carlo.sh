#!/bin/bash
#
cp ellipsoid_monte_carlo.hpp /$HOME/include
#
g++ -c -I/$HOME/include ellipsoid_monte_carlo.cpp
if [ $? -ne 0 ]; then
  echo "Errors compiling ellipsoid_monte_carlo.cpp"
  exit
fi
#
mv ellipsoid_monte_carlo.o ~/libcpp/$ARCH/ellipsoid_monte_carlo.o
#
echo "Library installed as ~/libcpp/$ARCH/ellipsoid_monte_carlo.o"
