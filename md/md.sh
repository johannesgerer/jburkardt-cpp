#!/bin/bash
#
g++ -c -O2 md.cpp >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling md.cpp"
  exit
fi
rm compiler.txt
#
g++ md.o
if [ $? -ne 0 ]; then
  echo "Errors linking and loading md.o."
  exit
fi
#
rm md.o
#
chmod ugo+x a.out
mv a.out ~/bincpp/$ARCH/md
#
echo "Executable installed as ~/bincpp/$ARCH/md"
