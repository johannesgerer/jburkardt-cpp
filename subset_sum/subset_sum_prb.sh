#!/bin/bash
#
g++ -c -g -I/$HOME/include subset_sum_prb.cpp >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling subset_sum_prb.cpp."
  exit
fi
rm compiler.txt
#
g++ subset_sum_prb.o /$HOME/libcpp/$ARCH/subset_sum.o -lm
if [ $? -ne 0 ]; then
  echo "Errors linking and loading subset_sum_prb.o."
  exit
fi
#
rm subset_sum_prb.o
#
mv a.out subset_sum_prb
./subset_sum_prb > subset_sum_prb_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running subset_sum_prb."
  exit
fi
rm subset_sum_prb
#
echo "Program output written to subset_sum_prb_output.txt"
