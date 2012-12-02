#!/bin/bash
#
g++ -c -g diaphony.cpp >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling diaphony.cpp"
  exit
fi
rm compiler.txt
#
g++ diaphony.o -lm
if [ $? -ne 0 ]; then
  echo "Errors linking and loading diaphony.o"
  exit
fi
rm diaphony.o
#
chmod ugo+x a.out
mv a.out ~/bincpp/$ARCH/diaphony
#
echo "Executable installed as ~/bincpp/$ARCH/diaphony"
