#!/bin/bash
#
g++ -c -g -I/$HOME/include asa091_prb.cpp >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling asa091_prb.cpp"
  exit
fi
rm compiler.txt
#
g++ asa091_prb.o /$HOME/libcpp/$ARCH/asa091.o -lm
if [ $? -ne 0 ]; then
  echo "Errors linking and loading asa091_prb.o."
  exit
fi
#
rm asa091_prb.o
#
mv a.out asa091_prb
./asa091_prb > asa091_prb_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running asa091_prb."
  exit
fi
rm asa091_prb
#
echo "Program output written to asa091_prb_output.txt"
