#!/bin/bash
#
g++ -c -g fem1d_adaptive.cpp >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling fem1d_adaptive.cpp"
  exit
fi
rm compiler.txt
#
g++ fem1d_adaptive.o
if [ $? -ne 0 ]; then
  echo "Errors linking and loading fem1d_adaptive.o."
  exit
fi
#
rm fem1d_adaptive.o
#
chmod ugo+x a.out
mv a.out ~/bincpp/$ARCH/fem1d_adaptive
#
echo "Executable installed as ~/bincpp/$ARCH/fem1d_adaptive"
