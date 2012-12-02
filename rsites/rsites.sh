#!/bin/bash
#
g++ -c rsites.cpp
if [ $? -ne 0 ]; then
  echo "Errors compiling rsites.cpp"
  exit
fi
#
g++ rsites.o -lm
if [ $? -ne 0 ]; then
  echo "Errors linking and loading rsites.o"
  exit
fi
#
rm rsites.o
#
chmod ugo+x a.out
mv a.out ~/bincpp/$ARCH/rsites
#
echo "Executable installed as ~/bincpp/$ARCH/rsites"

