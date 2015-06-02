#!/bin/bash
#
g++ -c -g -I/$HOME/include interp_prb.cpp >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling interp_prb.cpp"
  exit
fi
rm compiler.txt
#
g++ interp_prb.o /$HOME/libcpp/$ARCH/interp.o -lm
if [ $? -ne 0 ]; then
  echo "Errors linking and loading interp_prb.o."
  exit
fi
#
rm interp_prb.o
#
mv a.out interp_prb
./interp_prb > interp_prb_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running interp_prb."
  exit
fi
rm interp_prb
#
echo "Program output written to interp_prb_output.txt"
