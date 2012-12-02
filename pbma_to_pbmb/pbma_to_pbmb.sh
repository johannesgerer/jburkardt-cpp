#!/bin/bash
#
g++ -c -I$HOME/include pbma_to_pbmb.cpp >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling pbma_to_pbmb.cpp"
  exit
fi
rm compiler.txt
#
g++ pbma_to_pbmb.o -lm
if [ $? -ne 0 ]; then
  echo "Errors linking and loading pbma_to_pbmb.o"
  exit
fi
#
rm pbma_to_pbmb.o
#
chmod ugo+x a.out
mv a.out ~/bincpp/$ARCH/pbma_to_pbmb
#
echo "Executable installed as ~/bincpp/$ARCH/pbma_to_pbmb"
