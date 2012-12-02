#!/bin/bash
#
g++ -c -g -I/$HOME/include asa241_prb.cpp >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling asa241_prb.cpp"
  exit
fi
rm compiler.txt
#
g++ asa241_prb.o /$HOME/libcpp/$ARCH/asa241.o -lm
if [ $? -ne 0 ]; then
  echo "Errors linking and loading asa241_prb.o."
  exit
fi
#
rm asa241_prb.o
#
mv a.out asa241_prb
./asa241_prb > asa241_prb_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running asa241_prb."
  exit
fi
rm asa241_prb
#
echo "Program output written to asa241_prb_output.txt"
