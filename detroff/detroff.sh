#!/bin/bash
#
g++ -c -g -I /$HOME/include detroff.cpp >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling detroff.cpp"
  exit
fi
rm compiler.txt
#
g++ detroff.o -lm
if [ $? -ne 0 ]; then
  echo "Errors linking and loading detroff.o."
  exit
fi
rm detroff.o
#
chmod ugo+x a.out
mv a.out ~/bincpp/$ARCH/detroff
#
echo "Executable installed as ~/bincpp/$ARCH/detroff"
