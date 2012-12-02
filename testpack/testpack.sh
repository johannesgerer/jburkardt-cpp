#!/bin/bash
#
g++ -c -g testpack.cpp >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling testpack.cpp"
  exit
fi
rm compiler.txt
#
g++ testpack.o
if [ $? -ne 0 ]; then
  echo "Errors linking and loading testpack.o."
  exit
fi
#
rm testpack.o
#
chmod ugo+x a.out
mv a.out ~/bincpp/$ARCH/testpack
#
echo "Executable installed as ~/bincpp/$ARCH/testpack"
