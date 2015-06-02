#!/bin/bash
#
g++ -c tsp_brute.cpp >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling tsp_brute.cpp"
  exit
fi
rm compiler.txt
#
g++ tsp_brute.o -lm
if [ $? -ne 0 ]; then
  echo "Errors linking and loading tsp_brute.o"
  exit
fi
rm tsp_brute.o
#
mv a.out ~/bincpp/$ARCH/tsp_brute
#
echo "Executable installed as ~/bincpp/$ARCH/tsp_brute"
