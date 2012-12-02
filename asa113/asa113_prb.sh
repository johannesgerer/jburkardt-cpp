#!/bin/bash
#
g++ -c -g -I/$HOME/include asa113_prb.cpp >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling asa113_prb.cpp"
  exit
fi
rm compiler.txt
#
g++ asa113_prb.o /$HOME/libcpp/$ARCH/asa113.o -lm
if [ $? -ne 0 ]; then
  echo "Errors linking and loading asa113_prb.o."
  exit
fi
#
rm asa113_prb.o
#
mv a.out asa113_prb
./asa113_prb > asa113_prb_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running asa113_prb."
  exit
fi
rm asa113_prb
#
echo "Program output written to asa113_prb_output.txt"
