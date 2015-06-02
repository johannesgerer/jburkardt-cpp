#!/bin/bash
#
g++ -c -I/$HOME/include ellipsoid_monte_carlo_prb.cpp
if [ $? -ne 0 ]; then
  echo "Errors compiling ellipsoid_monte_carlo_prb.cpp"
  exit
fi
#
g++ -o ellipsoid_monte_carlo_prb ellipsoid_monte_carlo_prb.o /$HOME/libcpp/$ARCH/ellipsoid_monte_carlo.o -lm
if [ $? -ne 0 ]; then
  echo "Errors linking and loading ellipsoid_monte_carlo_prb.o."
  exit
fi
#
rm ellipsoid_monte_carlo_prb.o
#
./ellipsoid_monte_carlo_prb > ellipsoid_monte_carlo_prb_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running ellipsoid_monte_carlo_prb."
  exit
fi
rm ellipsoid_monte_carlo_prb
#
echo "Program output written to ellipsoid_monte_carlo_prb_output.txt"
