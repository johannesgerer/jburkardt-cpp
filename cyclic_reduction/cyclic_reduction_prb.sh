#!/bin/bash
#
g++ -c -g -I/$HOME/include cyclic_reduction_prb.cpp >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling cyclic_reduction_prb.cpp"
  exit
fi
rm compiler.txt
#
g++ cyclic_reduction_prb.o /$HOME/libcpp/$ARCH/cyclic_reduction.o -lm
if [ $? -ne 0 ]; then
  echo "Errors linking and loading cyclic_reduction_prb.o."
  exit
fi
#
rm cyclic_reduction_prb.o
#
mv a.out cyclic_reduction_prb
./cyclic_reduction_prb > cyclic_reduction_prb_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running cyclic_reduction_prb."
  exit
fi
rm cyclic_reduction_prb
#
echo "Program output written to cyclic_reduction_prb_output.txt"
