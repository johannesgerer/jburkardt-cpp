#!/bin/bash
#
g++ -c -I/$HOME/include ellipse_monte_carlo_prb.cpp >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling ellipse_monte_carlo_prb.cpp"
  exit
fi
rm compiler.txt
#
g++ ellipse_monte_carlo_prb.o /$HOME/libcpp/$ARCH/ellipse_monte_carlo.o -lm
if [ $? -ne 0 ]; then
  echo "Errors linking and loading ellipse_monte_carlo_prb.o"
  exit
fi
#
rm ellipse_monte_carlo_prb.o
#
mv a.out ellipse_monte_carlo_prb
./ellipse_monte_carlo_prb > ellipse_monte_carlo_prb_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running ellipse_monte_carlo_prb."
  exit
fi
rm ellipse_monte_carlo_prb
#
echo "Program output written to ellipse_monte_carlo_prb_output.txt"
