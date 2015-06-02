#!/bin/bash
#
g++ -c -g -I/$HOME/include ball_monte_carlo_prb.cpp >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling ball_monte_carlo_prb.cpp"
  exit
fi
rm compiler.txt
#
g++ ball_monte_carlo_prb.o /$HOME/libcpp/$ARCH/ball_monte_carlo.o -lm
if [ $? -ne 0 ]; then
  echo "Errors linking and loading ball_monte_carlo_prb.o."
  exit
fi
#
rm ball_monte_carlo_prb.o
#
mv a.out ball_monte_carlo_prb
./ball_monte_carlo_prb > ball_monte_carlo_prb_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running ball_monte_carlo_prb."
  exit
fi
rm ball_monte_carlo_prb
#
echo "Program output written to ball_monte_carlo_prb_output.txt"
