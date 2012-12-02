#!/bin/bash
#
g++ -c -g -I$HOME/include nint_exactness.cpp >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling nint_exactness.cpp"
  exit
fi
rm compiler.txt
#
g++ nint_exactness.o -lm
if [ $? -ne 0 ]; then
  echo "Errors linking and loading nint_exactness.o."
  exit
fi
#
rm nint_exactness.o
#
chmod ugo+x a.out
mv a.out ~/bincpp/$ARCH/nint_exactness
#
echo "Executable installed as ~/bincpp/$ARCH/nint_exactness"
