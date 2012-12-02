#!/bin/bash
#
g++ -c lf2cr.cpp >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling lf2cr.cpp"
  exit
fi
rm compiler.txt
#
g++ lf2cr.o
if [ $? -ne 0 ]; then
  echo "Errors linking and loading lf2cr.o."
  exit
fi
#
rm lf2cr.o
mv a.out ~/bincpp/$ARCH/lf2cr
#
echo "Executable installed as ~/bincpp/$ARCH/lf2cr"
