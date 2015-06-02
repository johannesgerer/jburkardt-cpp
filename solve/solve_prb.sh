#!/bin/bash
#
g++ -c -I/$HOME/include solve_prb.cpp
if [ $? -ne 0 ]; then
  echo "Errors compiling solve_prb.cpp"
  exit
fi
#
g++ solve_prb.o /$HOME/libcpp/$ARCH/solve.o -lm
if [ $? -ne 0 ]; then
  echo "Errors linking and loading solve_prb.o."
  exit
fi
#
rm solve_prb.o
#
mv a.out solve_prb
./solve_prb > solve_prb_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running solve_prb."
  exit
fi
rm solve_prb
#
echo "Program output written to solve_prb_output.txt"
