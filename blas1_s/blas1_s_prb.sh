#!/bin/bash
#
g++ -c -g -I$HOME/include blas1_s_prb.cpp >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling blas1_s_prb.cpp"
  exit
fi
rm compiler.txt
#
g++ blas1_s_prb.o $HOME/libcpp/$ARCH/blas0.o $HOME/libcpp/$ARCH/blas1_s.o -lm
if [ $? -ne 0 ]; then
  echo "Errors linking and loading blas1_s_prb.o."
  exit
fi
#
rm blas1_s_prb.o
#
mv a.out blas1_s_prb
./blas1_s_prb > blas1_s_prb_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running blas1_s_prb."
  exit
fi
rm blas1_s_prb
#
echo "Program output written to blas1_s_prb_output.txt"
