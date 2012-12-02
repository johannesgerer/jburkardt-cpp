#!/bin/bash
#
g++ -c -g -I/$HOME/include brent_old_prb.cpp >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling brent_old_prb.cpp"
  exit
fi
rm compiler.txt
#
g++ brent_old_prb.o /$HOME/libcpp/$ARCH/brent_old.o -lm
if [ $? -ne 0 ]; then
  echo "Errors linking and loading brent_old_prb.o."
  exit
fi
#
rm brent_old_prb.o
#
mv a.out brent_old_prb
./brent_old_prb > brent_old_prb_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running brent_old_prb."
  exit
fi
rm brent_old_prb
#
echo "Program output written to brent_old_prb_output.txt"
