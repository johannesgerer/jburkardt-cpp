#!/bin/bash
#
g++ -c sum_million.cpp >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling sum_million.cpp"
  exit
fi
rm compiler.txt
#
g++ sum_million.o
if [ $? -ne 0 ]; then
  echo "Errors linking and loading sum_million.o"
  exit
fi
rm sum_million.o
#
chmod ugo+x a.out
mv a.out ~/bincpp/$ARCH/sum_million
#
echo "Program installed as ~/bincpp/$ARCH/sum_million"
