#!/bin/bash
#
cp bernstein_polynomial.hpp /$HOME/include
#
g++ -c -g -I /$HOME/include bernstein_polynomial.cpp >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling bernstein_polynomial.cpp"
  exit
fi
rm compiler.txt
#
mv bernstein_polynomial.o ~/libcpp/$ARCH/bernstein_polynomial.o
#
echo "Library installed as ~/libcpp/$ARCH/bernstein_polynomial.o"
