#!/bin/bash
#
g++ -c mxv.cpp >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling mxv.cpp"
  exit
fi
rm compiler.txt
#
g++ mxv.o
if [ $? -ne 0 ]; then
  echo "Errors linking and loading mxv.o"
  exit
fi
rm mxv.o
#
mv a.out ~/bincpp/$ARCH/mxv
#
echo "Executable installed as ~/bincpp/$ARCH/mxv"
