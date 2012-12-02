#!/bin/bash
#
g++ -c -g int_exactness_legendre.cpp >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling int_exactness_legendre.cpp"
  exit
fi
rm compiler.txt
#
g++ int_exactness_legendre.o -lm
if [ $? -ne 0 ]; then
  echo "Errors linking and loading int_exactness_legendre.o"
  exit
fi
rm int_exactness_legendre.o
#
chmod ugo+x a.out
mv a.out ~/bincpp/$ARCH/int_exactness_legendre
#
echo "Executable installed as ~/bincpp/$ARCH/int_exactness_legendre"
