#!/bin/bash
#
g++ -c -g -I/$HOME/include prob_prb.cpp >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling prob_prb.cpp"
  exit
fi
rm compiler.txt
#
g++ prob_prb.o /$HOME/libcpp/$ARCH/prob.o -lm
if [ $? -ne 0 ]; then
  echo "Errors linking and loading prob_prb.o."
  exit
fi
#
rm prob_prb.o
#
mv a.out prob_prb
./prob_prb > prob_prb_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running prob_prb."
  exit
fi
rm prob_prb
#
echo "Program output written to prob_prb_output.txt"
