#!/bin/bash
#
g++ -c -g -I/$HOME/include r8lib_prb.cpp >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling r8lib_prb.cpp"
  exit
fi
rm compiler.txt
#
g++ r8lib_prb.o /$HOME/libcpp/$ARCH/r8lib.o -lm
if [ $? -ne 0 ]; then
  echo "Errors linking and loading r8lib_prb.o"
  exit
fi
#
rm r8lib_prb.o
#
mv a.out r8lib_prb
./r8lib_prb > r8lib_prb_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running r8lib_prb."
  exit
fi
rm r8lib_prb
#
echo "Program output written to r8lib_prb_output.txt"
