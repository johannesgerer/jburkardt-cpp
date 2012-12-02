#!/bin/bash
#
g++ -c -g wrap.cpp >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling wrap.cpp"
  exit
fi
rm compiler.txt
#
g++ wrap.o
if [ $? -ne 0 ]; then
  echo "Errors linking and loading wrap.o."
  exit
fi
#
rm wrap.o
#
chmod ugo+x a.out
mv a.out ~/bincpp/$ARCH/wrap
#
echo "Executable installed as ~/bincpp/$ARCH/wrap"
