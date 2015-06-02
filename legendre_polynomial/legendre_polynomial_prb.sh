#!/bin/bash
#
g++ -c -I/$HOME/include legendre_polynomial_prb.cpp
if [ $? -ne 0 ]; then
  echo "Errors compiling legendre_polynomial_prb.cpp"
  exit
fi
#
g++ legendre_polynomial_prb.o /$HOME/libcpp/$ARCH/legendre_polynomial.o -lm
if [ $? -ne 0 ]; then
  echo "Errors linking and loading legendre_polynomial_prb.o"
  exit
fi
#
rm legendre_polynomial_prb.o
#
mv a.out legendre_polynomial_prb
./legendre_polynomial_prb > legendre_polynomial_prb_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running legendre_polynomial_prb."
  exit
fi
rm legendre_polynomial_prb
#
echo "Program output written to legendre_polynomial_prb_output.txt"
