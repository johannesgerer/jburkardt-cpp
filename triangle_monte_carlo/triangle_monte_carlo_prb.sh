#!/bin/bash
#
g++ -c -g -I/$HOME/include triangle_monte_carlo_prb.cpp >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling triangle_monte_carlo_prb.cpp"
  exit
fi
rm compiler.txt
#
g++ triangle_monte_carlo_prb.o /$HOME/libcpp/$ARCH/triangle_monte_carlo.o -lm
if [ $? -ne 0 ]; then
  echo "Errors linking and loading triangle_monte_carlo_prb.o."
  exit
fi
#
rm triangle_monte_carlo_prb.o
#
mv a.out triangle_monte_carlo_prb
./triangle_monte_carlo_prb > triangle_monte_carlo_prb_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running triangle_monte_carlo_prb."
  exit
fi
rm triangle_monte_carlo_prb
#
echo "Program output written to triangle_monte_carlo_prb_output.txt"
