#!/bin/bash
#
g++ -c cr2lf.cpp >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling cr2lf.cpp"
  exit
fi
rm compiler.txt
#
g++ cr2lf.o
if [ $? -ne 0 ]; then
  echo "Errors linking and loading cr2lf.o."
  exit
fi
#
rm cr2lf.o
mv a.out ~/bincpp/$ARCH/cr2lf
#
echo "Executable installed as ~/bincpp/$ARCH/cr2lf"
