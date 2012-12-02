#!/bin/bash
#
g++ -c -g -I/$HOME/include owens_prb.cpp >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling owens_prb.cpp"
  exit
fi
rm compiler.txt
#
g++ owens_prb.o /$HOME/libcpp/$ARCH/owens.o -lm
if [ $? -ne 0 ]; then
  echo "Errors linking and loading owens_prb.o."
  exit
fi
#
rm owens_prb.o
#
mv a.out owens_prb
./owens_prb > owens_prb_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running owens_prb."
  exit
fi
rm owens_prb
#
echo "Program output written to owens_prb_output.txt"
