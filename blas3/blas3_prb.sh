#!/bin/bash
#
g++ -c -g -I$HOME/include blas3_prb.cpp >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling blas3_prb.cpp"
  exit
fi
rm compiler.txt
#
g++ blas3_prb.o $HOME/libcpp/$ARCH/blas0.o $HOME/libcpp/$ARCH/blas3.o -lm
if [ $? -ne 0 ]; then
  echo "Errors linking and loading blas3_prb.o."
  exit
fi
#
rm blas3_prb.o
#
mv a.out blas3_prb
./blas3_prb > blas3_prb_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running blas3_prb."
  exit
fi
rm blas3_prb
#
echo "Program output written to blas3_prb_output.txt"
