#!/bin/bash
#
g++ -c -g -I/$HOME/include asa152_prb.cpp >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling asa152_prb.cpp"
  exit
fi
rm compiler.txt
#
g++ asa152_prb.o /$HOME/libcpp/$ARCH/asa152.o -lm
if [ $? -ne 0 ]; then
  echo "Errors linking and loading asa152_prb.o."
  exit
fi
#
rm asa152_prb.o
#
mv a.out asa152_prb
./asa152_prb > asa152_prb_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running asa152_prb."
  exit
fi
rm asa152_prb
#
echo "Program output written to asa152_prb_output.txt"
