#!/bin/bash
#
g++ -c -g -I/$HOME/include filum_prb.cpp >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling filum_prb.cpp"
  exit
fi
rm compiler.txt
#
g++ filum_prb.o /$HOME/libcpp/$ARCH/filum.o -lm
if [ $? -ne 0 ]; then
  echo "Errors linking and loading filum_prb.o."
  exit
fi
#
rm filum_prb.o
#
mv a.out filum_prb
./filum_prb > filum_prb_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running filum_prb."
  exit
fi
rm filum_prb
#
echo "Program output written to filum_prb_output.txt"
