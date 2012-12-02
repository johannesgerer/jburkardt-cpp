#!/bin/bash
#
g++ -c -g -I/$HOME/include asa136_prb.cpp >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling asa136_prb.cpp"
  exit
fi
rm compiler.txt
#
g++ asa136_prb.o /$HOME/libcpp/$ARCH/asa136.o -lm
if [ $? -ne 0 ]; then
  echo "Errors linking and loading asa136_prb.o."
  exit
fi
#
rm asa136_prb.o
#
mv a.out asa136_prb
./asa136_prb > asa136_prb_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running asa136_prb."
  exit
fi
rm asa136_prb
#
echo "Program output written to asa136_prb_output.txt"
