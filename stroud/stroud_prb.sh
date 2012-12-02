#!/bin/bash
#
g++ -c -g -I/$HOME/include stroud_prb.cpp >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling stroud_prb.cpp"
  exit
fi
rm compiler.txt
#
g++ stroud_prb.o /$HOME/libcpp/$ARCH/stroud.o -lm
if [ $? -ne 0 ]; then
  echo "Errors linking and loading stroud_prb.o."
  exit
fi
#
rm stroud_prb.o
#
mv a.out stroud_prb
./stroud_prb > stroud_prb_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running stroud_prb."
  exit
fi
rm stroud_prb
#
echo "Program output written to stroud_prb_output.txt"
