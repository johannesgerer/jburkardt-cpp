#!/bin/bash
#
cp chebyshev_polynomial.hpp /$HOME/include
#
g++ -c -g -I /$HOME/include chebyshev_polynomial.cpp >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling chebyshev_polynomial.cpp"
  exit
fi
rm compiler.txt
#
mv chebyshev_polynomial.o ~/libcpp/$ARCH/chebyshev_polynomial.o
#
echo "Library installed as ~/libcpp/$ARCH/chebyshev_polynomial.o"
