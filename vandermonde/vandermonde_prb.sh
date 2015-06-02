#!/bin/bash
#
g++ -c -g -I/$HOME/include vandermonde_prb.cpp >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling vandermonde_prb.cpp"
  exit
fi
rm compiler.txt
#
g++ vandermonde_prb.o /$HOME/libcpp/$ARCH/vandermonde.o -lm
if [ $? -ne 0 ]; then
  echo "Errors linking and loading vandermonde_prb.o."
  exit
fi
#
rm vandermonde_prb.o
#
mv a.out vandermonde_prb
./vandermonde_prb > vandermonde_prb_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running vandermonde_prb."
  exit
fi
rm vandermonde_prb
#
echo "Program output written to vandermonde_prb_output.txt"
