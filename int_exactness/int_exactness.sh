#!/bin/bash
#
g++ -c -g int_exactness.cpp >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling int_exactness.cpp"
  exit
fi
rm compiler.txt
#
g++ int_exactness.o -lm
if [ $? -ne 0 ]; then
  echo "Errors linking and loading int_exactness.o"
  exit
fi
rm int_exactness.o
#
chmod ugo+x a.out
mv a.out ~/bincpp/$ARCH/int_exactness
#
echo "Executable installed as ~/bincpp/$ARCH/int_exactness"
