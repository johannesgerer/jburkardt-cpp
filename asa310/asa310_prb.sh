#!/bin/bash
#
g++ -c -g -I/$HOME/include asa310_prb.cpp >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling asa310_prb.cpp"
  exit
fi
rm compiler.txt
#
g++ asa310_prb.o /$HOME/libcpp/$ARCH/asa310.o -lm
if [ $? -ne 0 ]; then
  echo "Errors linking and loading asa310_prb.o."
  exit
fi
#
rm asa310_prb.o
#
mv a.out asa310_prb
./asa310_prb > asa310_prb_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running asa310_prb."
  exit
fi
rm asa310_prb
#
echo "Program output written to asa310_prb_output.txt"
