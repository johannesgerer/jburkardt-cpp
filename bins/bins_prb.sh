#!/bin/bash
#
g++ -c -g -I/$HOME/include bins_prb.cpp >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling bins_prb.cpp"
  exit
fi
rm compiler.txt
#
g++ bins_prb.o /$HOME/libcpp/$ARCH/bins.o -lm
if [ $? -ne 0 ]; then
  echo "Errors linking and loading bins_prb.o."
  exit
fi
#
rm bins_prb.o
#
mv a.out bins_prb
./bins_prb > bins_prb_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running bins_prb."
  exit
fi
rm bins_prb
#
echo "Program output written to bins_prb_output.txt"
