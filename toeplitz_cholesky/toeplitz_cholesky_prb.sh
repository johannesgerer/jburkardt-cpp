#!/bin/bash
#
g++ -c -g -I/$HOME/include toeplitz_cholesky_prb.cpp >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling toeplitz_cholesky_prb.cpp"
  exit
fi
rm compiler.txt
#
g++ toeplitz_cholesky_prb.o /$HOME/libcpp/$ARCH/toeplitz_cholesky.o -lm
if [ $? -ne 0 ]; then
  echo "Errors linking and loading toeplitz_cholesky_prb.o."
  exit
fi
#
rm toeplitz_cholesky_prb.o
#
mv a.out toeplitz_cholesky_prb
./toeplitz_cholesky_prb > toeplitz_cholesky_prb_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running toeplitz_cholesky_prb."
  exit
fi
rm toeplitz_cholesky_prb
#
echo "Program output written to toeplitz_cholesky_prb_output.txt"
