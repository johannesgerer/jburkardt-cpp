#!/bin/bash
#
g++ -c -g -I/$HOME/include asa245_prb.cpp >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling asa245_prb.cpp"
  exit
fi
rm compiler.txt
#
g++ asa245_prb.o /$HOME/libcpp/$ARCH/asa245.o -lm
if [ $? -ne 0 ]; then
  echo "Errors linking and loading asa245_prb.o."
  exit
fi
#
rm asa245_prb.o
#
mv a.out asa245_prb
./asa245_prb > asa245_prb_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running asa245_prb."
  exit
fi
rm asa245_prb
#
echo "Program output written to asa245_prb_output.txt"
