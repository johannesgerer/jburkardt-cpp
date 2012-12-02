#!/bin/bash
#
g++ -c -g fem1d_sample.cpp >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling fem1d_sample.cpp"
  exit
fi
rm compiler.txt
#
g++ fem1d_sample.o
if [ $? -ne 0 ]; then
  echo "Errors linking and loading fem1d_sample.o."
  exit
fi
#
rm fem1d_sample.o
#
chmod ugo+x a.out
mv a.out ~/bincpp/$ARCH/fem1d_sample
#
echo "Executable installed as ~/bincpp/$ARCH/fem1d_sample"
