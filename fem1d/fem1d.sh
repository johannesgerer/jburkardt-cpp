#!/bin/bash
#
g++ -c -g fem1d.cpp >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling fem1d.cpp"
  exit
fi
rm compiler.txt
#
g++ fem1d.o
if [ $? -ne 0 ]; then
  echo "Errors linking and loading fem1d.o."
  exit
fi
#
rm fem1d.o
#
chmod ugo+x a.out
mv a.out ~/bincpp/$ARCH/fem1d
#
echo "Executable installed as ~/bincpp/$ARCH/fem1d"
