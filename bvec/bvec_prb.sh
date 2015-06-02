#!/bin/bash
#
g++ -c -I/$HOME/include bvec_prb.cpp
if [ $? -ne 0 ]; then
  echo "Errors compiling bvec_prb.cpp"
  exit
fi
#
g++ bvec_prb.o /$HOME/libcpp/$ARCH/bvec.o -lm
if [ $? -ne 0 ]; then
  echo "Errors linking and loading bvec_prb.o."
  exit
fi
#
rm bvec_prb.o
#
mv a.out bvec_prb
./bvec_prb > bvec_prb_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running bvec_prb."
  exit
fi
rm bvec_prb
#
echo "Program output written to bvec_prb_output.txt"
