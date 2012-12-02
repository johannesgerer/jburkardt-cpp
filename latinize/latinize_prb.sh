#!/bin/bash
#
g++ -c -g -I/$HOME/include latinize_prb.cpp >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling latinize_prb.cpp"
  exit
fi
rm compiler.txt
#
g++ latinize_prb.o /$HOME/libcpp/$ARCH/latinize.o -lm
if [ $? -ne 0 ]; then
  echo "Errors linking and loading latinize_prb.o."
  exit
fi
#
rm latinize_prb.o
#
mv a.out latinize_prb
./latinize_prb > latinize_prb_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running latinize_prb."
  exit
fi
rm latinize_prb
#
echo "Program output written to latinize_prb_output.txt"
