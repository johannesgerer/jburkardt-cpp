#!/bin/bash
#
g++ -c -g -I/$HOME/include dunavant_prb.cpp >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling dunavant_prb.cpp"
  exit
fi
rm compiler.txt
#
g++ dunavant_prb.o /$HOME/libcpp/$ARCH/dunavant.o -lm
if [ $? -ne 0 ]; then
  echo "Errors linking and loading dunavant_prb.o."
  exit
fi
#
rm dunavant_prb.o
#
mv a.out dunavant_prb
./dunavant_prb > dunavant_prb_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running dunavant_prb."
  exit
fi
rm dunavant_prb
#
echo "Program output written to dunavant_prb_output.txt"
