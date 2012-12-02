#!/bin/bash
#
g++ -c -g -I/$HOME/include randlc_prb.cpp >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling randlc_prb.cpp"
  exit
fi
rm compiler.txt
#
g++ randlc_prb.o /$HOME/libcpp/$ARCH/randlc.o -lm
if [ $? -ne 0 ]; then
  echo "Errors linking and loading randlc_prb.o."
  exit
fi
#
rm randlc_prb.o
#
mv a.out randlc_prb
./randlc_prb > randlc_prb_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running randlc_prb."
  exit
fi
rm randlc_prb
#
echo "Program output written to randlc_prb_output.txt"
