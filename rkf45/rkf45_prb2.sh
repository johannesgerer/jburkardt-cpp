#!/bin/bash
#
g++ -c -g -I/$HOME/include rkf45_prb2.cpp >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling rkf45_prb2.cpp"
  exit
fi
rm compiler.txt
#
g++ rkf45_prb2.o /$HOME/libcpp/$ARCH/rkf45.o -lm
if [ $? -ne 0 ]; then
  echo "Errors linking and loading rkf45_prb2.o."
  exit
fi
#
rm rkf45_prb2.o
#
mv a.out rkf45_prb2
./rkf45_prb2 > rkf45_prb2_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running rkf45_prb2."
  exit
fi
rm rkf45_prb2
#
echo "Program output written to rkf45_prb2_output.txt"
