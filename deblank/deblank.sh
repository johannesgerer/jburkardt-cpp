#!/bin/bash
#
g++ -c deblank.cpp >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling deblank.cpp"
  exit
fi
rm compiler.txt
#
g++ deblank.o
if [ $? -ne 0 ]; then
  echo "Errors linking and loading deblank.o."
  exit
fi
#
rm deblank.o
#
chmod ugo+x a.out
mv a.out ~/bincpp/$ARCH/deblank
#
echo "Program installed as ~/bincpp/$ARCH/deblank"
