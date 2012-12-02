#!/bin/bash
#
g++ -c -g -I/$HOME/include rkf45_prb.cpp >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling rkf45_prb.cpp"
  exit
fi
rm compiler.txt
#
g++ rkf45_prb.o /$HOME/libcpp/$ARCH/rkf45.o -lm
if [ $? -ne 0 ]; then
  echo "Errors linking and loading rkf45_prb.o."
  exit
fi
#
rm rkf45_prb.o
#
mv a.out rkf45_prb
./rkf45_prb > rkf45_prb_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running rkf45_prb."
  exit
fi
rm rkf45_prb
#
echo "Program output written to rkf45_prb_output.txt"
