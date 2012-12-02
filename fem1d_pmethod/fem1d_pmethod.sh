#!/bin/bash
#
g++ -c -g fem1d_pmethod.cpp >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling fem1d_pmethod.cpp"
  exit
fi
rm compiler.txt
#
g++ fem1d_pmethod.o
if [ $? -ne 0 ]; then
  echo "Errors linking and loading fem1d_pmethod.o."
  exit
fi
#
rm fem1d_pmethod.o
#
chmod ugo+x a.out
mv a.out ~/bincpp/$ARCH/fem1d_pmethod
#
echo "Executable installed as ~/bincpp/$ARCH/fem1d_pmethod"
