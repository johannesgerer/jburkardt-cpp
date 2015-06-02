#!/bin/bash
#
g++ -c -g -I/$HOME/include circle_monte_carlo_prb.cpp >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling circle_monte_carlo_prb.cpp"
  exit
fi
rm compiler.txt
#
g++ circle_monte_carlo_prb.o /$HOME/libcpp/$ARCH/circle_monte_carlo.o -lm
if [ $? -ne 0 ]; then
  echo "Errors linking and loading circle_monte_carlo_prb.o."
  exit
fi
#
rm circle_monte_carlo_prb.o
#
mv a.out circle_monte_carlo_prb
./circle_monte_carlo_prb > circle_monte_carlo_prb_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running circle_monte_carlo_prb."
  exit
fi
rm circle_monte_carlo_prb
#
echo "Program output written to circle_monte_carlo_prb_output.txt"
