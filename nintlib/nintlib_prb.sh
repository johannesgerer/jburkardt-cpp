#!/bin/bash
#
g++ -c -g -I/$HOME/include nintlib_prb.cpp >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling nintlib_prb.cpp"
  exit
fi
rm compiler.txt
#
g++ nintlib_prb.o /$HOME/libcpp/$ARCH/nintlib.o -lm
if [ $? -ne 0 ]; then
  echo "Errors linking and loading nintlib_prb.o."
  exit
fi
#
rm nintlib_prb.o
#
mv a.out nintlib_prb
./nintlib_prb > nintlib_prb_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running nintlib_prb."
  exit
fi
rm nintlib_prb
#
echo "Program output written to nintlib_prb_output.txt"
