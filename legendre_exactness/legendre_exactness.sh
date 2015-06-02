#!/bin/bash
#
g++ -c legendre_exactness.cpp >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling legendre_exactness.cpp"
  exit
fi
rm compiler.txt
#
g++ legendre_exactness.o -lm
if [ $? -ne 0 ]; then
  echo "Errors linking and loading legendre_exactness.o"
  exit
fi
rm legendre_exactness.o
#
chmod ugo+x a.out
mv a.out ~/bincpp/$ARCH/legendre_exactness
#
echo "Executable installed as ~/bincpp/$ARCH/legendre_exactness"
