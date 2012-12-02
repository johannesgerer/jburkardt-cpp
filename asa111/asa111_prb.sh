#!/bin/bash
#
g++ -c -g -I/$HOME/include asa111_prb.cpp >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling asa111_prb.cpp"
  exit
fi
rm compiler.txt
#
g++ asa111_prb.o /$HOME/libcpp/$ARCH/asa111.o -lm
if [ $? -ne 0 ]; then
  echo "Errors linking and loading asa111_prb.o."
  exit
fi
#
rm asa111_prb.o
#
mv a.out asa111_prb
./asa111_prb > asa111_prb_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running asa111_prb."
  exit
fi
rm asa111_prb
#
echo "Program output written to asa111_prb_output.txt"
