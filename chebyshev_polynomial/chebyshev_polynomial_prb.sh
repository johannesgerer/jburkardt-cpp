#!/bin/bash
#
g++ -c -g -I/$HOME/include chebyshev_polynomial_prb.cpp >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling chebyshev_polynomial_prb.cpp"
  exit
fi
rm compiler.txt
#
g++ chebyshev_polynomial_prb.o /$HOME/libcpp/$ARCH/chebyshev_polynomial.o -lm
if [ $? -ne 0 ]; then
  echo "Errors linking and loading chebyshev_polynomial_prb.o."
  exit
fi
#
rm chebyshev_polynomial_prb.o
#
mv a.out chebyshev_polynomial_prb
./chebyshev_polynomial_prb > chebyshev_polynomial_prb_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running chebyshev_polynomial_prb."
  exit
fi
rm chebyshev_polynomial_prb
#
echo "Program output written to chebyshev_polynomial_prb_output.txt"
