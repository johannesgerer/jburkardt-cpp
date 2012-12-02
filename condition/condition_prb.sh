#!/bin/bash
#
g++ -c -g -I/$HOME/include condition_prb.cpp >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling condition_prb.cpp"
  exit
fi
rm compiler.txt
#
g++ condition_prb.o /$HOME/libcpp/$ARCH/condition.o /$HOME/libcpp/$ARCH/r8lib.o -lm
if [ $? -ne 0 ]; then
  echo "Errors linking and loading condition_prb.o."
  exit
fi
#
rm condition_prb.o
#
mv a.out condition_prb
./condition_prb > condition_prb_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running condition_prb."
  exit
fi
rm condition_prb
#
echo "Program output written to condition_prb_output.txt"
