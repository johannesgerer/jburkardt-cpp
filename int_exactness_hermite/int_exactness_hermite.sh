#!/bin/bash
#
g++ -c -g int_exactness_hermite.cpp >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling int_exactness_hermite.cpp"
  exit
fi
rm compiler.txt
#
g++ int_exactness_hermite.o -lm
if [ $? -ne 0 ]; then
  echo "Errors linking and loading int_exactness_hermite.o"
  exit
fi
rm int_exactness_hermite.o
#
chmod ugo+x a.out
mv a.out ~/bincpp/$ARCH/int_exactness_hermite
#
echo "Executable installed as ~/bincpp/$ARCH/int_exactness_hermite"
