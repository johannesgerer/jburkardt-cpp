#!/bin/bash
#
g++ -c -g -I/$HOME/include asa047_prb.cpp >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling asa047_prb.cpp."
  exit
fi
rm compiler.txt
#
g++ asa047_prb.o /$HOME/libcpp/$ARCH/asa047.o -lm
if [ $? -ne 0 ]; then
  echo "Errors linking and loading asa047_prb.o."
  exit
fi
#
rm asa047_prb.o
#
mv a.out asa047_prb
./asa047_prb > asa047_prb_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running asa047_prb."
  exit
fi
rm asa047_prb
#
echo "Program output written to asa047_prb_output.txt"
