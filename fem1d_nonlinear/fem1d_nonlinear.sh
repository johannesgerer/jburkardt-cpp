#!/bin/bash
#
g++ -c -g fem1d_nonlinear.cpp >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling fem1d_nonlinear.cpp"
  exit
fi
rm compiler.txt
#
g++ fem1d_nonlinear.o
if [ $? -ne 0 ]; then
  echo "Errors linking and loading fem1d_nonlinear.o."
  exit
fi
#
rm fem1d_nonlinear.o
#
chmod ugo+x a.out
mv a.out ~/bincpp/$ARCH/fem1d_nonlinear
#
echo "Executable installed as ~/bincpp/$ARCH/fem1d_nonlinear"
