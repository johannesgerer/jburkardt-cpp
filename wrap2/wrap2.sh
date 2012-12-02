#!/bin/bash
#
g++ -c wrap2.cpp >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling wrap2.cpp"
  exit
fi
rm compiler.txt
#
g++ wrap2.o
if [ $? -ne 0 ]; then
  echo "Errors linking and loading wrap2.o."
  exit
fi
#
rm wrap2.o
#
chmod ugo+x a.out
mv a.out ~/bincpp/$ARCH/wrap2
#
echo "Executable installed as ~/bincpp/$ARCH/wrap2"
