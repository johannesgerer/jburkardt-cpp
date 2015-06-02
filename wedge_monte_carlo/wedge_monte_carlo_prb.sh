#!/bin/bash
#
g++ -c -I/$HOME/include wedge_monte_carlo_prb.cpp
if [ $? -ne 0 ]; then
  echo "Errors compiling wedge_monte_carlo_prb.cpp"
  exit
fi
#
g++ -o wedge_monte_carlo_prb wedge_monte_carlo_prb.o /$HOME/libcpp/$ARCH/wedge_monte_carlo.o -lm
if [ $? -ne 0 ]; then
  echo "Errors linking and loading wedge_monte_carlo_prb.o."
  exit
fi
#
rm wedge_monte_carlo_prb.o
#
./wedge_monte_carlo_prb > wedge_monte_carlo_prb_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running wedge_monte_carlo_prb."
  exit
fi
rm wedge_monte_carlo_prb
#
echo "Program output written to wedge_monte_carlo_prb_output.txt"
