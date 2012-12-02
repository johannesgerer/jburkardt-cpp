#!/bin/bash
#
g++ -c -g nas.cpp >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling nas.cpp"
  exit
fi
rm compiler.txt
#
g++ nas.o
if [ $? -ne 0 ]; then
  echo "Errors linking and loading nas.o"
  exit
fi
rm nas.o
#
chmod ugo+x a.out
mv a.out ~/bincpp/$ARCH/nas
#
echo "Executable installed as ~/bincpp/$ARCH/nas"
