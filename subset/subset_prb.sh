#!/bin/bash
#
g++ -c -g -I/$HOME/include subset_prb.cpp >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling subset_prb.cpp"
  exit
fi
rm compiler.txt
#
g++ subset_prb.o /$HOME/libcpp/$ARCH/subset.o -lm
if [ $? -ne 0 ]; then
  echo "Errors linking and loading subset_prb.o."
  exit
fi
#
rm subset_prb.o
#
mv a.out subset_prb
./subset_prb > subset_prb_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running subset_prb."
  exit
fi
rm subset_prb
#
echo "Program output written to subset_prb_output.txt"
