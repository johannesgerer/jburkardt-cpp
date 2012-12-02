#!/bin/bash
#
g++ -c -g -I/$HOME/include unicycle_prb.cpp >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling unicycle_prb.cpp"
  exit
fi
rm compiler.txt
#
g++ unicycle_prb.o /$HOME/libcpp/$ARCH/unicycle.o -lm
if [ $? -ne 0 ]; then
  echo "Errors linking and loading unicycle_prb.o."
  exit
fi
#
rm unicycle_prb.o
#
mv a.out unicycle_prb
./unicycle_prb > unicycle_prb_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running unicycle_prb."
  exit
fi
rm unicycle_prb
#
echo "Program output written to unicycle_prb_output.txt"
